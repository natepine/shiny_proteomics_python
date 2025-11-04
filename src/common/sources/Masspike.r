#'--------------------------------------------------------
#'
#' Support for loading data from Masspike data sources
#'
#'--------------------------------------------------------

# ========================== Network ==========================

allowSSLConf <- httr::config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE)

#' A constructor for a connection object to the GFY via the REST API
#'
#' @param url URL of the GFY server.
#' @param username Login name
#' @param key Api key
#' @param verifySSL Should self-signed ssl certificates be rejected
#'
#' @return gfy_api object
#' @export
#' @import httr
#' @import stringr
gfy_api <- function(url, username, key, verifySSL = TRUE) {
   hconf <- c(httr::authenticate(username, key, type = "basic"))
   if(!verifySSL)
      hconf <- c(hconf, allowSSLConf)

   structure(list("url.head" = stringr::str_trim(url), httr.config = hconf), class = "gfy_api")
}


#' A function to construct a url from the header the object is built with and a custom footer
#'
#' @param gfy.obj A gfy_api object set up to query
#' @param url.foot The footer of the url
#' @param sep The separator for the url between the header and footer. Defaults to "/"
#' @return Character string of the url to GET
#' @export
#' @import stringr
make.url <- function(gfy.obj, url.foot, sep = "/") {
   if (is.null(url.foot) || is.na(url.foot)) {
      warning("url.foot is null or na")
      return(gfy.obj$url.head)
   }
   if (!is.character(url.foot)) {
      stop("url.foot is not a character string")
   }
   url.foot <- stringr::str_trim(url.foot)
   if (url.foot == "") {
      warning("url.foot is only whitespace")
   }
   paste(gfy.obj$url.head, url.foot, sep = sep)
}

#' Perform a GET on the gfy_api object
#'
#' @param gfy.obj A gfy_api object set up to query
#' @param url The url to query
#' @param build.url Build the url from the url header that's part of gfy.obj
#' @param parse Logical to indicate if the raw response from GET should be parsed with parsed_content. If set to NULL or FALSE then the response object from GET will be returned
#' @return A response, either parsed if the parse option is set or a httr response object
#' @export
#' @import httr
gfy.get <- function(gfy.obj, url, build.url = TRUE, parse = 'parsed') {
   if (build.url) {
      url <- make.url(gfy.obj, url)
   }
   get.ret <- httr::GET(url, config = gfy.obj$httr.config)
   if (is.null(parse) | (is.logical(parse) && !parse)) {
      return(get.ret)
   } else {
      return(httr::content(get.ret, as = parse))
   }
}

#' Attempt to log in with a username and password
#'
#' @param url URL of the GFY server, "/login"" will be appended
#' @param username Login name
#' @param password Password
#' @param verifySSL Should self-signed ssl certificates be rejected
#' @return A gfy_api object on successful authentication, otherwise NULL
#' @export
#' @import httr
gfy.login <- function(url, username, password, verifySSL = TRUE) {
   hconf <- c()
   if(!verifySSL)
      hconf <- allowSSLConf

   login.ret <- httr::POST(paste0(url, "/login"),
      body = list(username = username, password = password),
      config = hconf)

   login.ret <- httr::content(login.ret)
   if (login.ret$status != "success") {
      if (login.ret$status == "error") {
         warning(paste("Unable to login to", url, "\n", login.ret$error))
      }
      return(NULL)
   } else {
      return(gfy_api(url, username, login.ret$key, verifySSL))
   }
}

# ========================== Masspike API ==========================

getAPIData <- function(gfy.obj, query) {
   if(!curl::has_internet())
      stop("No internet connection")

   tryCatch({
      data <- gfy.get(gfy.obj, query)
      if(!("status" %in% names(data) && data$status == "error"))
         return(data)

      warning(paste0(
         "Query ", gfy.obj$url.head, "/", query, " failed with status: error\n",
         data$error))
      },
      error = function(e){
         warning(paste0(
            "Unable to retrieve api data",
            "\nURL: ", gfy.obj$url.head, "/", query,
            "\nError: ", e))
      }
   )

   stop("Server failed to return data for this id")
}

# ========== Normalization ==========

#' Get normalization options from a quant id
#'
#' @param qid - Quant id to query
#' @param gfy.obj - GFY api options to use
#' @param isSiteQuant - Specifies id type, protein quant or site quant
#'
#' @return A list with elements factors, normBy and normRow or NULL if query fails
#'          Note* factors = NULL means that the response indicates that quant id was not normalized

gfySourceNormalization <- function(qid, gfy.obj, isSiteQuant) {
   noNorm <- list(
      factors = NULL,
      normBy = "none",
      normRow = NA
   )

   if(is.null(gfy.obj))
      return(NULL)

   if(isSiteQuant) {
      module <- "site_quant_info"
   }else {
      module <- "protein_quant_info"
   }

   response <- tryCatch({
      query <- paste0(module, "/", qid, "?all_users=1")
      gfy.get(gfy.obj, query)
   },
   error = function(e) NULL)

   if(is.null(response))
      return(NULL)

   status <- names(response)
   if(!is.null(status) && "error" %in% status)
      return(NULL)

   options <- response[[1]]$options
   if(isSiteQuant) {
      if("normalizationFactors" %in% names(options)) {
         names(options)[names(options) == "normalizationFactors"] <- "normalization_factors"
      }else {
         pqid <- options$normalization
         if(pqid == "NONE" || (is.logical(pqid) && !pqid))
            return(noNorm)

         pqNorm <- getNormalization(gfy.obj, pqid, FALSE)
         if(!is.null(pqNorm)) {
            return(pqNorm)
         }else {
            warning(
               paste("Failed to get protein quant info for ID:", pqid,
               "\nDefaulting to no normalization")
            )
            return(noNorm)
         }
      }
      normBy <- "ALL"
      normRow <- ""
   }else {
      normBy <- options$norm_by
      normRow <- extractUniprot(options$norm_protein)
   }

   normFactors <- options$normalization_factors

   if(all(unlist(normFactors) == 1))
      return(noNorm)

   normBy <- switch(normBy,
      "ALL" = "all",
      "protein" = "row",
      "None" = "none")

   if(!nzchar(normRow)) {
      normRow <- NA
   }

   return(list(
      factors = normFactors,
      normBy = normBy,
      normRow = normRow
   ))
}

# ========== Proteins ==========

#' Get CORE ProteinQuant or SiteQuant data
#'
#' @param id - ProteinQuant or SiteQuant ID to query
#' @param gfy.obj - gfy object used to query CORE server
#' @param isSiteQuant - Should id be used as a SiteQuant ID?
gfySourceProteins <- function(id, gfy.obj, isSiteQuant) {
   if(isSiteQuant){
      module <- "combined_site_quant"
      params <- "all_users=1&peptide_parsimony=UR"
   }else{
      module <- "protein_quant"
      params <- "all_users=1"
   }

   query <- paste0(module, "/", id, "?", params)
   return(getAPIData(gfy.obj, query))
}

# ========== Peptides ==========

#' Get CORE ProteinQuant or SiteQuant peptide data
#'
#' @param id - ProteinQuant or SiteQuant ID to query
#' @param gfy.obj - gfy object used to query CORE server
#' @param isSiteQuant - Should id be used as a SiteQuant ID?
gfySourcePeptides <- function(id, gfy.obj, isSiteQuant) {
   if(isSiteQuant) {
      module <- "site_quant_peptide"
      params <- "all_users=1&peptide_parsimony=UR"
   }else {
      module <- "protein_quant_peptide"
      params <- "all_users=1"
   }

   query <- paste0(module, "/", id, "?", params)
   return(getAPIData(gfy.obj, query))
}

# ============= Data Source =============

gfySource <- createDatasetSource(
   proteins      = function(dataset) gfySourceProteins(     dataset$ID, dataset$gfy.obj, dataset$isSiteQuant),
   peptides      = function(dataset) gfySourcePeptides(     dataset$ID, dataset$gfy.obj, dataset$isSiteQuant),
   normalization = function(dataset) gfySourceNormalization(dataset$ID, dataset$gfy.obj, dataset$isSiteQuant),

   proteins_msg  = function(dataset) paste("Pulling", dataset$idLabel, "data from server. This may take a few seconds, please wait."),
   peptides_msg  = function(dataset) paste("Pulling Peptide data from server. This may take a few seconds, please wait.")
)

gfySrcToCache <- function(cache) {
   src <- gfySource
   src$proteins      <- srcToCache(cache, gfySource$proteins,      "proteins")
   src$peptides      <- srcToCache(cache, gfySource$peptides,      "peptides")
   src$normalization <- srcToCache(cache, gfySource$normalization, "normalization")
   return(src)
}

