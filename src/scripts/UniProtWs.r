library(httr)
library(readr)
library(stringr)
library(yaml)

#' Recursively combine a nested list by alternating AND and OR
#'
#' @param x - Data to collapse
#' @param and - Collapse with AND

r_paste <- function(x, and = TRUE) {
   if(is.list(x)) {
      x <- sapply(x, r_paste, !and)
      has_key <- names(x) != ""
      names(x)[has_key] <- paste0(names(x)[has_key], ":")
      x <- paste0("(", names(x), x, ")", collapse = ifelse(and, "+AND+", "+OR+"))
   }
   return(x)
}

#' Build a UniProt API query string
#'
#' @param query - An unnamed string to be used as-is or a named list of query fields,
#'                list elements of length > 1 will be OR'ed together
#' @param fields - A character vector of returned fields
#' @param stream - Use the stream API? Otherwise batched queries will be made to the search API
#' @param size - Number of results queried per API call. Max of 500
#' @param format - Format to query

uniProtQueryStr <- function(query, fields, stream = FALSE, size = 500, format = c("tsv", "list", "json", "xml", "fasta", "gff", "obo", "rdf", "xlsx")) {
   if(!is.null(names(query))){
      query <- r_paste(query)
   }

   args <- c(
      query = query,
      compressed = "true",
      format = format
   )

   if(!missing(fields))
      args["fields"] <- paste(fields, collapse = ",")

   endpoint <- "stream"
   if(!stream) {
      endpoint <- "search"
      args["size"] <- size
   }
   paste0(endpoint, "/?",  paste0(names(args), "=", args, collapse = "&"))
}

#' Parse link header
#'
#' @param rsp - An httr response object
#'
#' @return URL from 'link' header

nextURL <- function(rsp) {
   hdrs <- headers(rsp)
   if(!"link" %in% names(hdrs))
      return(NULL)

   link <- hdrs$link
   str_extract(link, "(?<=<)[^>]*")
}

#' Query the UniProt API
#'
#' @param query - UniProt API query parameters as a named vector
#' @param fields - Vector of UniProt API return fields
#' @param processBatch - Function to call on `size` entries
#' @param database - Uniprot database to query, one of c("uniprotkb", "proteomes", "taxonomy")
#' @param stream - Use the stream API? Otherwise batched queries will be made to the search API
#' @param showProgress - Ignored if stream = TRUE, display batch progress
#' @param size - Ignored if stream = TRUE, query batch size
#' @param format - Format to query, "tsv" and "list" are supported
#'
#' @return The result of rbind on each return value of processRows(rows)

queryUniProt <- function(query, fields, processBatch = function(batch) batch, database = c("uniprotkb", "proteomes", "taxonomy"),
      stream = FALSE, showProgress = FALSE, size = 500, format = c("tsv", "list")) {
   database <- match.arg(database)
   format <- match.arg(format)

   url <- paste0("https://rest.uniprot.org/", database, "/",
      uniProtQueryStr(query, fields, stream = stream, size = size, format = format))

   showProgress <- showProgress
   result <- NULL
   if(showProgress) {
      if(!stream) {
         cat("\n")
         total <- NULL
         current <- 0
         cat("Progress: 0%\033[0K\r")
      }else {
         cat("Querying uniprot.org, please wait...\n")
         showProgress <- FALSE
      }
   }
   repeat {
      rsp <- GET(url)

      if(showProgress) {
         if(is.null(total)) total <- as.integer(headers(rsp)$`x-total-results`)

         current <- min(current + size, total)
         percent <- as.integer(100 * current / total)
         cat(sprintf("Progress: %d%% [%d/%d]\033[0K\r", percent, current, total))
      }

      if(http_error(rsp)) {
         if(stream && status_code(rsp) == 429) {
            cat("Uniprot stream failed, switching to search\n")
            return(queryUniProt(query, fields, processRows, rows, stream = FALSE, showProgress, size))
         }
         rspContent <- content(rsp, as = "text", enc = "UTF-8")
         if(!is.character(rspContent))
            rspContent <- toString(rspContent)
         stop("Caught http error: ", rspContent, call. = FALSE)
      }else {
         rspContent <- memDecompress(content(rsp, as = "raw"), type = "gzip")
         result <- switch(format,
         "tsv" = {
            rows <- read_tsv(rspContent, progress = FALSE, show_col_types = FALSE)
            result <- rbind(result, processBatch(rows))
         },
         "list" = {
            rspContent <- rawToChar(rspContent)
            entries <- as.vector(str_match_all(rspContent, "[^\n]+")[[1]])
            result <- c(result, entries)
         })
         if(stream || is.null(url <- nextURL(rsp)))
            break
      }
   }

   if(showProgress) cat("\n\n")

   return(result)
}

