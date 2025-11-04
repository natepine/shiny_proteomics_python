library(DBI)
library(stringr)

source("scripts/cliPrompts.r")

fields <- c("Server", "Species")

conf_vals <- NULL
confPath <- "data/conf.yml"
if(file.exists(confPath)) {
   source("common/Configs.r")
   library(yaml)

   config <- tryCatch({
      config <- loadConfig(confPath)
      list(
         Server = names(config$SERVERS),
         Species = names(config$ANNOTATED_SPECIES)
      )
   }, error = function(e) {
      cat(sep = "\n", "Error reading conf.yml", e)
      NULL
   })

   setConf <- function(Field, val, newVal) {
      conf <- read_yaml(confPath)
      conf <- switch(Field,
         "Server" = within(conf, names(SERVERS)[names(SERVERS) == val] <- newVal),
         "Species" = within(conf, names(ANNOTATED_SPECIES)[names(ANNOTATED_SPECIES) == val] <- newVal),
         stop("Called setConf with invalid field [", str(Field), "]")
      )
      write_yaml(conf, confPath)
   }
   replaceConf <- function(Field, val, newVal) {
      conf <- read_yaml(confPath)
      conf <- switch(Field,
         "Server" = within(conf, {
            SERVERS[[newVal]] <- SERVERS[[val]]
            SERVERS[[val]] <- NULL
         }),
         "Species" = within(conf, {
            ANNOTATED_SPECIES[[newVal]] <- ANNOTATED_SPECIES[[val]]
            ANNOTATED_SPECIES[[val]] <- NULL
         }),
         stop("Called replaceConf with invalid field [", str(Field), "]")
      )
      write_yaml(conf, confPath)
   }
   printConf <- function(Field, val) {
      conf <- read_yaml(confPath)
      cat(val)
      cat("\n")
      write_yaml(
         switch(Field,
            "Server" = with(conf, SERVERS[[val]]),
            "Species" = with(conf, ANNOTATED_SPECIES[[val]]),
            stop("Called printConf with invalid field [", str(Field), "]")
         ),
         stdout()
      )
   }
}else {
   cat("WARNING: No conf.yml detected\n")
}
if(is.null(config)){
   cat("Continuing without config checks\n\n")
}

check_ws_dups <- function(vals, override = FALSE) {
   function(val) {
      out <- list(
         val = val,
         passed = FALSE
      )
      if(str_detect(val, "\\s")) {
         out$msg <- "Whitespace detected, try again"
      }else if(val %in% vals) {
         cat("Input would cause duplicates\n")
         if(override)
            out$passed <- boolPrompt("Value", val, "will be permanently lost. Continue anyway")
      }else {
         out$passed <- TRUE
      }
      return(out)
   }
}

rename <- function() {
   dbPath <- "data/database.db"
   db <- dbConnect(RSQLite::SQLite(), dbPath)
   on.exit(dbDisconnect(db))

   header("Renaming Utility")
   repeat {
      cat("Fields to rename\n")
      Field <- selectPrompt(c(fields, "Quit"))
      if(Field == "Quit")
         break
      field <- tolower(Field)

      query <- paste("SELECT DISTINCT", Field, "FROM metadata")
      vals <- dbGetQuery(db, query)[[Field]]
      if(length(vals) == 0) {
         cat("Database is empty, exiting\n")
         quit()
      }

      if(!is.null(config)) {
         repeat {
            conf_vals <- config[[Field]]
            unknown <- setdiff(vals, conf_vals)
            if(length(unknown) > 0) {
               cat(sep = "\n",
                  paste("Found", field, "names in database with no matching config"),
                  paste("-", unknown))
               cat("\n")

               unused <- setdiff(conf_vals, vals)
               if(length(unused) > 0) {
                  cat(sep = "\n",
                     paste("Unused config", field, "names"),
                     paste("-", unused))

                  if(boolPrompt("Keep database names and rename a config", field)) {
                     conf_val <- selectPrompt(unused)

                     cat("Rename", conf_val, "to\n")
                     new_conf <- selectPrompt(c(unknown, "Enter manually"))
                     if(new_conf == "Enter manually")
                        new_conf <- prompt("Enter new name or leave blank to cancel",
                           check = check_ws_dups(setdiff(unused, conf_val)))

                     if(nzchar(new_conf)) {
                        setConf(Field, conf_val, new_conf)
                        conf_vals[conf_vals == conf_val] <- new_conf
                     }
                  }else {
                     break
                  }
               }else {
                  break
               }
            }else {
               break
            }
         }
      }

      repeat {
         if(length(vals) > 1) {
            cat("Database", field, "names")
            val <- selectPrompt(vals)
            cat("Leave blank to select a different", field)
            cat("\n")
         }else {
            val <- vals
         }
         newVal <- prompt(paste("Rename", val, "to"),
            check = check_ws_dups(setdiff(vals, val), override = TRUE))
         if(nzchar(newVal))
            break
      }

      query <- paste("UPDATE metadata SET", Field, "= @newVal WHERE", Field, "= @val")
      numUpdated <- dbExecute(db, query, params = list(val = val, newVal = newVal))
      cat("\nSuccessfully updated", numUpdated, "database entries\n\n")

      if(!is.null(config) && val != newVal) {
         conf_vals <- config[[Field]]
         if(val %in% conf_vals) {
            if(newVal %in% conf_vals) {
               cat("Merge conflict: New", field, "already exists in config\n")
               printConf(Field, val)
               cat("\n")
               printConf(Field, newVal)
               if(boolPrompt("Replace existing entry [", newVal, "] with ", val, sep = "")) {
                  replaceConf(Field, val, newVal)
                  conf_vals <- conf_vals[conf_vals != val]
                  cat("Entry replaced\n\n")
               }
            }else {
               setConf(Field, val, newVal)
               conf_vals[conf_vals == val] <- newVal
               cat("Successfully updated config\n\n")
            }
         }else if(!newVal %in% conf_vals) {
            cat("No", field, " '", val, "' found in config\nManual review of conf.yml required\n\n")
         }
      }

      if(Field == "Server") {
         if(dir.exists(val)) {
            if(file.rename(val, newVal)) {
               cat("Successfully renamed cache\n\n")
            }else {
               cat(sep = "",
                  "Unable to move cache data at ", normalizePath(val), "\n",
                  "Please move the appropriate cache to ", normalizePath(newVal), "\n")
            }
         }
      }

      config[[Field]] <- conf_vals
   }
}

if(!interactive()) {
   f <- file("stdin")
   open(f, blocking=TRUE)
   rename()
}
