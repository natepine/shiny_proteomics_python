suppressPackageStartupMessages({
   source("scripts/UniProtWs.r")
   library(DBI)
   library(dplyr)
   library(data.table)
})

#In-memory cache for remote data
if(!exists("netData"))
   netData <- new.env(parent = baseenv())

GO_types <- paste("GO", c("Molecular Function", "Biological Process", "Cellular Component"))
annot_types <- c(GO_types,
   "Reactome",
   "Interpro",
   "KEGG Pathway"
)

buildFilter <- function(filter, connective = c("AND", "OR")) {
   if(is.null(filter))
      return("")

   if(is.null(names(filter)))
      stop("Non-null filter must have names")

   connective <- sprintf(" %s ", match.arg(connective))
   paste('WHERE', paste0(names(filter), '="', filter, '"', collapse = connective))
}

dbCount <- function(db, table, filter = NULL) {
   count <- sprintf('SELECT Count(*) FROM "%s" %s', table, buildFilter(filter))
   dbGetQuery(db, count)[[1]]
}

dbSubset <- function(db, table, n, filter = NULL) {
   subset <- sprintf(
      'SELECT * FROM "%s" WHERE rowid IN (SELECT rowid FROM "%s" %s ORDER BY RANDOM() LIMIT %d)',
      table, table, buildFilter(filter), n)
   dbGetQuery(db, subset)
}

printStats <- function(db, taxid) {
   cat(sprintf("Database Stats [%d]\n\n", taxid))

   if(taxid == 9606) {
      annot_types <- c(annot_types, "KEGG Disease")
   }

   proteinTable <- paste0("proteins_", taxid)
   cat("Proteins:", dbCount(db, proteinTable), "\n\n")

   cat("Annotations\n")
   annotTable <- paste0("annotations_", taxid)
   counts <- sapply(annot_types, function(type) {
      dbCount(db, annotTable, c(Database = type))
   })
   cat(sep = "\n", sprintf("%s: %d", annot_types, counts))
   total <- dbCount(db, annotTable)
   cat("Total:", total, "\n")

   if(total != sum(counts)) {
      uniqueDBs <- paste("SELECT DISTINCT Database FROM", annotTable)
      uniqueDBs <- dbGetQuery(uniqueDBs)$Database
      unknown <- setdiff(uniqueDBs, annot_types)
      if(length(unknown) > 0)
         warning("Unknown database(s): ", paste(unknown, collapse = ", "))
   }

   cat("\nPlease check parsed values\n\n")

   peek <- bind_rows(
      lapply(annot_types, function(type){
         dbSubset(db, annotTable, 3, c(Database = type))
      })
   )
   print(peek)
}

#' Helper for reading tsv from urls
#'
#' @param url - URL to query for data
#' @param ... - Passed to readr::read_tsv

get_tsv <- function(url, ...) {
   passNull <- function(e) {
      cat("get_tsv failed for", url, "\n")
      cat(e$message, "\n")
      NULL
   }

   tryCatch({
      rsp <- content(GET(url), as = "text", encoding = "UTF-8")
      data <- read_tsv(rsp, ..., show_col_types = FALSE, progress = FALSE)
      as.data.frame(data)
   }, warning = passNull, error = passNull)
}

#' Helper function to insert dataframes into the db
#'
#' @param db - SQLite db connection
#' @param table - A table name in the db with columns matching @keys
#' @param data - A dataframe with colnames matching @keys
#' @param keys - A vector containing the column names of data to insert

dbInsert <- function(db, table, data, keys = NULL) {
   if(is.null(keys)) {
      keys <- colnames(data)
   }else {
      data <- as.data.frame(data)[keys]
   }

   insert <- paste0(
      'INSERT INTO "', table, '" ',
      '(', paste0(keys, collapse = ", "), ') VALUES ',
      '(', paste0("@", keys, collapse = ", "), ')'
   )
   dbExecute(db, insert, params = data)
}

#' Separates @col into new rows by @sep
#'
#' @param data - Dataframe containing columns 'UniprotID' and @col
#' @param col - Name of the column to split
#' @param sep - Character string containing each character to split on
#' @param colName - Name of output column, defaults to DatabaseID
#'
#' @return New data.frame with a row for each split value in @col
#'          or NULL if all values are NA

melt <- function(data, col, sep = ";", colName = "DatabaseID") {
   valid <- !is.na(data[[col]])
   if(!any(valid))
      return(NULL)

   data <- data[valid, c("UniprotID", col)]
   regex <- sprintf("[^ %s][^%s]*", sep, sep)
   rows <- str_extract_all(data[[col]], regex)
   n_rows <- sapply(rows, length)

   data <- data.frame(
      rep(data[["UniprotID"]], n_rows),
      unlist(rows)
   )
   colnames(data) <- c("UniprotID", colName)

   return(data)
}

resetOrgCache <- function(taxid) {
   objs <- c(
      "kegg_pathway",
      "kegg_disease"
   )

   lapply(paste(objs, taxid, sep = "_"), function(cacheObj) {
      if(identical(get0(cacheObj, netData), NA))
         remove(cacheObj, envir = netData)
      NULL
   })

   invisible()
}

parseGO <- function(dbName, data, table, db) {
   data <- melt(data, dbName)

   if(!is.null(data)) {
      matches <- str_match(data$DatabaseID, "(?<Val>.*) \\[(?<ID>GO:[0-9]*)\\]")
      data[c("Annotation", "DatabaseID")] <- matches[, c("Val", "ID")]
      data$Database <- dbName

      dbInsert(db, table, data)
   }

   return(NULL)
}

cacheKEGG <- function(data, taxid, keggDB = c("pathway", "disease")) {
   keggObj <- paste("kegg", keggDB, taxid, sep = "_")
   if(!exists(keggObj, netData)) {
      keggDB <- match.arg(keggDB)
      if(keggDB == "disease") {
         kegg_org <- "hsa"
      }else if(!is.null( kegg <- melt(data, "KEGG") )) {
         kegg_orgs <- str_extract(kegg$DatabaseID, "[^:]*")
         kegg_org <- unique(kegg_orgs)
         if(length(kegg_org) != 1) {
            kegg_org <- sort(table(kegg_orgs), decreasing = TRUE)[[1]]
            warning("Parsed multiple values for kegg_org: ",
               paste0(kegg_orgs, collapse = ", "), "\nDefaulting to ", kegg_org)
         }
      }else {
         return(NULL)
      }

      list_org <- ifelse(keggDB == "disease", "", kegg_org)
      mapURL <- paste("https://rest.kegg.jp/link", keggDB, kegg_org, sep = "/")
      url    <- paste("https://rest.kegg.jp/list", keggDB, list_org, sep = "/")
      keggMap   <- get_tsv(mapURL, c("OrgID", "DatabaseID"))
      keggList <- get_tsv(url, c("DatabaseID", "Annotation"))

      if(!is.null(keggMap) && !is.null(keggList)) {
         keggMap <- as.data.table(keggMap, key = "DatabaseID")
         keggMap <- keggMap[, DatabaseID := str_remove(DatabaseID, "^.*:")]
         keggCache <- keggMap[keggList, on = "DatabaseID", nomatch = NULL]
      }else {
         keggCache <- NA
         warning("Skipping KEGG ", keggDB, " for taxid: ", taxid,
            "\nQuery failed for KEGG organism: ", kegg_org)
      }

      assign(keggObj, keggCache, netData)
   }else {
      keggCache <- get(keggObj, netData)
   }

   return(keggCache)
}

parseKEGG <- function(data, taxid, kegg, keggDB = c("Pathway", "Disease"), table, db) {
   keggDB <- match.arg(keggDB)
   keggCache <- cacheKEGG(data, taxid, tolower(keggDB))

   if(!identical(keggCache, NA)) {
      kegg <- keggCache[kegg, on = "OrgID", nomatch = NULL][, -"OrgID"]
      kegg$Database <- paste("KEGG", keggDB)

      if(keggDB == "Pathway")
         kegg$Annotation <- str_remove(kegg$Annotation, " - .*$")

      dbInsert(db, table, kegg)
   }
}

parseReactome <- function(data, table, db) {
   reactome <- get0("reactome", netData)
   if(!is.null(reactome)) {
      reactome <- reactome[data$UniprotID, nomatch = NULL]
      reactome$Database <- "Reactome"
      dbInsert(db, table, reactome)
   }
}

parseInterpro <- function(data, table, db) {
   interpro <- get0("interpro", netData)
   if(!is.null(interpro)) {
      uniprotMap <- melt(data, "InterPro")
      if(!is.null(uniprotMap)) {
         interpro <- interpro[uniprotMap, on = "DatabaseID", nomatch = NULL]
         interpro$Database <- "Interpro"
         dbInsert(db, table, interpro)
      }
   }
}

initNetData <- function() {
   if(is.null(get0("reactome", netData))) {
      reactomeCols <- c(
         "UniprotID",
         "DatabaseID",
         "link",
         "Annotation",
         "Evidence Code",
         "Species"
      )
      reactomeURL <- "https://reactome.org/download/current/UniProt2Reactome.txt"
      reactome <- get_tsv(reactomeURL, reactomeCols, lazy = TRUE)
      if(!is.null(reactome)) {
         reactome <- as.data.table(reactome[c("UniprotID", "DatabaseID", "Annotation")], key = "UniprotID")
         netData$reactome <- reactome
      }else {
         warning("Failed to pull Reactome data from\n", reactomeURL, "\nContinuing without Reactome annotations")
      }
   }

   if(is.null(get0("interpro", netData))) {
      interproCols <- c(
         "DatabaseID",
         "Entry Type",
         "Annotation"
      )
      interproURL <- "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list"
      interpro <- get_tsv(interproURL, interproCols, lazy = TRUE)
      if(!is.null(interpro)) {
         interpro <- as.data.table(interpro[c("DatabaseID", "Annotation")], key = "DatabaseID")
         netData$interpro <- interpro
      }else {
         warning("Failed to pull Interpro data from\n", interproURL, "\nContinuing without Interpro annotations")
      }
   }
}

#' Populate SQLite table for one organism
#' Overwrites any previous data
#'
#' @param db - SQLite connection
#' @param taxid - Uniprot organism_id to query
#' @param append - Keep current data and replace on collisions, otherwise drop current data
#' @param stream - Use the Uniprot stream endpoint

buildSpeciesDB <- function(db, taxid, append, stream) {
   if(is.na(taxid))
      return(NULL)

   resetOrgCache(taxid)

   toProteome <- list(organism_id = taxid, proteome_type = 1) # Type 1 indicates refererence proteome
   refProteome <- tryCatch(
      queryUniProt(toProteome, database = "proteomes", stream = stream, format = "list"),
      error = function(e) NULL, warning = function(e) NULL
   )

   cat("Updating organism:", taxid)
   query <- list(list(is_isoform = "true", is_isoform = "false"))
   if(!is.null(refProteome) && length(refProteome) == 1) {
      query$proteome <- refProteome
      cat(" [", refProteome, "]", sep = "")
   }else {
      query$organism_id <- taxid
   }
   cat("\n")

   table <- paste0("proteins_", taxid)
   if(!append && dbExistsTable(db, table)) {
      dbExecute(db, sprintf('DROP TABLE "%s"', table))
   }

   proteinSchema <- sprintf(
      'CREATE TABLE IF NOT EXISTS "%s" (
      UniprotID  TEXT PRIMARY KEY ON CONFLICT REPLACE,
      Sequence   TEXT NOT NULL,
      Function   TEXT,
      PDB        TEXT,
      Location   TEXT
      );', table)
   dbExecute(db, proteinSchema)

   annotTable <- paste0("annotations_", taxid)
   if(!append && dbExistsTable(db, annotTable)) {
      dbExecute(db, sprintf('DROP TABLE "%s"', annotTable))
   }

   annotSchema <- sprintf(
      'CREATE TABLE IF NOT EXISTS "%s" (
      DatabaseID TEXT NOT NULL,
      Database   TEXT NOT NULL,
      Annotation TEXT NOT NULL,
      UniprotID  TEXT NOT NULL,
      FOREIGN KEY (UniprotID) REFERENCES "%s" (UniprotID),
      PRIMARY KEY (Database, DatabaseID, UniprotID) ON CONFLICT REPLACE
      );', annotTable, table)
   dbExecute(db, annotSchema)

   storeRows <- function(rows) {
      rows <- rows %>% rename(
         UniprotID               = Entry,
         `GO Molecular Function` = "Gene Ontology (molecular function)",
         `GO Biological Process` = "Gene Ontology (biological process)",
         `GO Cellular Component` = "Gene Ontology (cellular component)",
         Location                = "Subcellular location [CC]",
         Function                = "Function [CC]"
      )

      protCols <- c("UniprotID", "Sequence", "PDB", "Function", "Location")
      dbInsert(db, table, rows, protCols)

      lapply(GO_types, parseGO, rows, annotTable, db)

      kegg <- melt(rows, "KEGG", colName = "OrgID")
      if(!is.null(kegg)) {
         parseKEGG(rows, taxid, kegg, "Pathway", annotTable, db)

         if(taxid == 9606)
            parseKEGG(rows, taxid, kegg, "Disease", annotTable, db)
      }

      parseReactome(rows, annotTable, db)

      parseInterpro(rows, annotTable, db)

      return(NULL)
   }

   fields <- c("accession", "go_f", "go_c", "go_p", "xref_kegg", "xref_interpro", "sequence", "xref_pdb", "cc_subcellular_location", "cc_function")

   updateDB <- tryCatch({
      queryUniProt(query, fields, storeRows, stream = stream, showProgress = TRUE)
      printStats(db, taxid)
      TRUE
   }, error = function(e) {
      cat(sep = "\n",
         "Error, skipping organism",
         e$message,
         "",
         paste("Try again with the command: ./annotations build", taxid),
         ""
      )
      FALSE
   })

   return(updateDB)
}

annotations <- function(append, orgs = NULL, batch = TRUE, dbPath = "data/annotations.db", confPath = "data/conf.yml") {
   if(!curl::has_internet()) {
      cat("No internet connection\n")
      return(invisible(NULL))
   }
   taxids <- orgs
   if(is.null(orgs) || !is.numeric(orgs)) {
      if(file.exists(confPath)) {
         source("common/Configs.r")
         config <- loadConfig(confPath)
         if(!is.null(orgs) && length(orgs) > 0) {
            species_choices <- names(config$ANNOTATED_SPECIES)
            orgs_ <- tolower(orgs)
            choices_ <- tolower(species_choices)
            found <- orgs_ %in% choices_
            if(!all(found)) {
               stop("Unable to map species name(s) to taxid: ", paste0(orgs[!found], collapse = ","), "\n",
                  "Available choices are: ", paste0(species_choices, collapse = ","))
            }
            taxids <- unlist(config$ANNOTATED_SPECIES[choices_ %in% orgs_])
         }else {
            taxids <- unlist(config$ANNOTATED_SPECIES)
            if(all(is.na(taxids))) {
               cat("No configured species to annotate\n")
               return(invisible(NULL))
            }
         }
      }else {
         stop("Numeric taxids not specified and config not found at\n", normalizePath(confPath))
      }
   }

   db <- dbConnect(RSQLite::SQLite(), dbPath)
   on.exit(dbDisconnect(db), add = TRUE)

   taxids <- taxids[!is.na(taxids)]
   if(append && is.null(orgs)) {
      tables <- c("proteins", "annotations")
      tablesExist <- sapply(taxids, function(taxid) {
         all(sapply(paste(tables, taxid, sep = "_"), dbExistsTable, conn = db))
      })
      taxids <- taxids[!tablesExist]
   }

   if(length(taxids) != 0) {
      initNetData()

      metadataSchema <-
         'CREATE TABLE IF NOT EXISTS metadata (
         TaxID       TEXT PRIMARY KEY ON CONFLICT REPLACE,
         LastUpdated DATE NOT NULL
         );'
      dbExecute(db, metadataSchema)
      setMetadata <- 'INSERT INTO metadata (TaxID, LastUpdated) VALUES (?, julianday())'

      cat("Updating database for taxid(s):", paste(taxids, collapse = ", "), "\n")

      dbBegin(db)
      on.exit(dbRollback(db), add = TRUE, after = FALSE)

      lapply(taxids, function(taxid) {
         commit <- buildSpeciesDB(db, taxid, append, stream = !batch)
         if(commit) {
            dbExecute(db, setMetadata, params = taxid)
            dbCommit(db)
            cat("Committed changes to database\n")
         }else {
            dbRollback(db)
            cat("Rolled back changes to database\n")
         }
         dbBegin(db)
      })

      on.exit()
      dbDisconnect(db)
      cat("Updates finished\n")
   }else {
      print("No annotations to update")
   }
}

if(!interactive()) {
   args <- commandArgs(trailingOnly = TRUE)
   action <- ifelse(length(args) == 0, "help", args[[1]])
   if(!any(action == c("append", "build", "update")) || any(c("-h", "help", "-help", "--help") %in% args)) {
      cat(
"
Annotation command line utility
Create and update the IsoParser webapp's protein annotation database for all selected organisms

Usage: ./annotations action [options] [organisms...]

Action
append\tKeeping current data, replacing on ID conflicts
update\tDiscard current data before populating database
build \tSynonym for --no-batch update
help  \tShow this message

Options
+b, --batch\tSlow with low resource demand. Recommended for updating in the background (with webapp running) [Default]
-b, --no-batch\tFast with high resource demand. Recommended for initial setup

Organisms
Taxon ids or organism names in conf.yml
If omitted, organisms will be supplied by conf.yml
   When append is specified, only organisms missing a table will be selected
"
      )
   }else {
      options <- c("+b", "-b", "--batch", "--no-batch")
      batch <- !(any(c("-b", "--no-batch") %in% args) || action == "build")
      args <- setdiff(args, options)

      orgs <- args[-1]
      if(length(orgs) == 0)
         orgs <- NULL

      append <- action == "append"
      annotations(append, orgs = orgs, batch = batch)
   }
}
