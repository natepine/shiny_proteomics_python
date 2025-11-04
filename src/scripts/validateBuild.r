#  Database/config correlation

docker_cmd <- function(cmd) {
   paste("docker compose run --rm -it app", cmd)
}
service_name <- "app"

confPath <- "data/conf.yml"
if(!file.exists(confPath)) {
   cat("App config not found, exiting")
   quit(status = 1)
}

source("common/Configs.r")
config <- loadConfig(confPath)

dbPath <- "data/database.db"
if(file.exists(dbPath)) {
   library(DBI)
   db <- DBI::dbConnect(RSQLite::SQLite(), dbPath)

   warn_unknown <- function(vals, config, label) {
      missing <- setdiff(vals, config)
      if(length(missing) > 1) {
         cat(sep = "\n",
            paste0("WARNING: Database contains unknown ", label, "(s)"),
            paste0(missing, collapse = ", "),
            paste("Delete entries containing these values or configure a", label, "for each")
         )
         exit(1)
      }
   }

   servers <- dbGetQuery(db, "SELECT DISTINCT Server From Metadata;")
   species <- dbGetQuery(db, "SELECT DISTINCT Species From Metadata;")
   dbDisconnect(db)

   warn_unknown(servers, names(config$SERVERS), "server")
   warn_unknown(species, names(config$ANNOTATED_SPECIES), "species")
}


#  BioPlex
bioplexPath <- "data/bioplex_3.0.rds"
if("9606" %in% config$ANNOTATED_SPECIES) {
   cat("\nBioplex data not found\n")

   if(!file.exists(bioplexPath)) {
      cat("Downloading Bioplex data\n")

      noBioPlex <- function(e) {
         cat(sep = "\n",
            "WARNING: Failed to download missing BioPlex data",
            e$message,
            "Interactions tab disabled"
         )
      }
      tryCatch({
         bioplex <- read.delim('https://bioplex.hms.harvard.edu/data/BioPlex_293T_Network_10K_Dec_2019.tsv')
         colnames(bioplex) <- c('GeneA', 'GeneB', 'UniprotA', 'UniprotB', 'SymbolA', 'SymbolB', 'pW', 'pNI', 'pInt')
         saveRDS(bioplex, bioplexPath)
      }, warning = noBioPlex, error = noBioPlex)
   }

   #Unix system calls used (file.symlink and Sys.readlink)
   is_linked <- function(path) Sys.readlink(path) == basename(bioplexPath)

   bioplexLink <- "data/bioplex.rds"
   if(file.exists(bioplexPath)) {
      cat("Linking Bioplex data\n")

      link <- Sys.readlink(bioplexLink)
      has_link <- if(is.na(link)) {
         FALSE
      }else if(link != basename(bioplexPath)) {
         unlink(bioplexLink)
         FALSE
      }else {
         TRUE
      }

      if(!has_link) {
         invisible(capture.output(
            file.symlink(basename(bioplexPath), bioplexLink)))
      }
   }
}


#  Annotations
species <- config$ANNOTATED_SPECIES
if(!identical(species, NA)){
   library(DBI)
   annotPath <- "data/annotations.db"
   if(file.exists(annotPath)) {
      db <- dbConnect(RSQLite::SQLite(), annotPath)

      hasTable <- sapply(species, function(taxid) {
         if(is.na(taxid))
            return(TRUE)
         dbExistsTable(db, paste0("annotations_", taxid))
      })
      if(!all(hasTable)) {
         cat(sep = "\n",
            paste("WARNING: Annotation tables missing for",
               paste0(names(species)[!hasTable], collapse = ", ")),
            "",
            "Generate tables with this command:",
            docker_cmd("annotations build")
         )
      }

      dbDisconnect(db)
   }else {
      cat(sep = "\n",
         "WARNING: Annotation database not found at",
         normalizePath(annotPath),
         "Build the database with the command:",
         docker_cmd("annotations build"),
         "Continuing without enrichment analysis"
      )
   }
}

invisible()
