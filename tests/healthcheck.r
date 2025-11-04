log_levels <- c(
   "silent" = 0,
   "crash"  = 1,
   "error"  = 2,
   "warn"   = 3,
   "info"   = 4
)

usage <- function() {
   cat(sprintf(
"
Healthcheck scans each viewer in the database and attempts to load it with log output for warnings, errors and crashes

Usage: healthcheck [options...]

Options
all         Run the check on all viewers
-h, --help  Show this message

Otherwise, options are filters defined by key-value pairs and include:
before=<YYYY-MM-DD>
after=<YYYY-MM-DD>
n=<max_viewers>
random=<on|off>
verbosity=<%s>

",
   paste(names(log_levels), collapse = "|")))
}

if(!interactive()) {
   args <- commandArgs(trailingOnly=TRUE)
   if(length(args) == 0 || any(args %in% c("--help", "-h"))) {
      usage()
      quit()
   }
}

appPath <- Sys.getenv("APP_DIR")

suppressPackageStartupMessages({
   library(DBI)
   library(jsonlite)

   source(paste0(appPath, "/tests/common.r"))

   source(paste0(appPath, "/IsoParser/global.r"))
   source(paste0(appPath, "/IsoParser/server.r"))
})

dataPath <- paste0(appPath, "/data/")
dbPath   <- paste0(appPath, "/data/database.db")

test_health <- function(verbosity = "error", after = "2000-01-01", before = "3000-01-01", random = FALSE, n = NULL) {
   if(!verbosity %in% names(log_levels))
      stop("Invalid verbosity")

   db <- dbConnect(RSQLite::SQLite(), dbPath)
   on.exit(dbDisconnect(db))

   data <- dbGetQuery(db, "SELECT * FROM metadata WHERE date BETWEEN ? AND ?", c(after, before))
   if(nrow(data) == 0) {
      cat("No viewers selected")
      quit(1)
   }

   if(is.null(n)) {
      n <- nrow(data)
   }else {
      n <- min(n, nrow(data))
   }

   rows <- seq(nrow(data))
   if(random) {
      rows <- sample(rows, n)
   }else {
      rows <- head(rows, n)
   }

   cache <- createCache(dataPath)
   srcs <- cacheSrc(cache)

   results <- list()
   results$total <- n
   results$cached <- 0
   results$peps_cached <- 0

   results$warnings <- createLog()
   results$errors   <- createLog()
   results$crashes  <- createLog()
   for(row in rows) {
      metadata <- parseMetadata(data[row, ])
      streamInfo <- function(msg) {
         if(log_levels[[verbosity]] >= log_levels[["info"]])
            cat(sprintf("[INFO] <%s>\t%s\n", metadata$key, msg))
      }

      proteins <- cache$read(dataset_cache_key(metadata, "proteins"))
      if(is.null(proteins)) {
         streamInfo("Skipping: no proteins")
         next
      }
      srcs$proteins <- function(...) proteins

      streamInfo(sprintf("ID: %d\tServer: %s", metadata$ID, metadata$server))

      nop <- function(...) NULL
      streamFirstLog <- function(v_thresh, log, title, msg) {
         if(log_levels[[verbosity]] >= log_levels[[v_thresh]] &&
            !title %in% names(log)
         ) {
            catln(sprintf("[%s] <%s> %s - %s", toupper(v_thresh), metadata$key, title, msg))
         }
      }
      recordWarning <- function(title, msg, duration) {
         streamFirstLog("warn", results$warnings, title, msg)
         results$warnings <<- logEvent(results$warnings, metadata$key, title, msg)
      }
      recordError <- function(title, msg) {
         streamFirstLog("error", results$errors, title, msg)
         results$errors   <<- logEvent(results$errors, metadata$key, title, msg)
      }

      results$cached <- results$cached + 1
      view <- createLoadDatasetView(nop, recordWarning, recordError, nop, clearMsg = NA)
      peps <- tryCatch({
         loadDataset(metadata, dataPath, NULL, view, srcs = srcs)$peptides
      },
      error = function (e) {
         streamFirstLog("crash", results$crashes, "Test failed", e$message)
         results$crashes <<- logEvent(results$crashes, metadata$key, e$message)
         return(NULL)
      })

      if(!is.null(peps)) {
         results$peps_cached <- results$peps_cached + 1
      }
   }

   return(results)
}

report_health <- function(results) {
   catln("")
   catln("* DETAILS *")
   printLog("Warnings", results$warnings)
   printLog("Errors",   results$errors)
   printLog("Crashes",  results$crashes)

   num_warnings <- summarizeLog(results$warnings)
   num_errors   <- summarizeLog(results$errors)
   num_crashes  <- summarizeLog(results$crashes)

   catln("")
   catln("* SUMMARY *")
   printLog("Warnings", num_warnings, cat_val)
   catln("")
   printLog("Errors",   num_errors,   cat_val)
   catln("")
   printLog("Crashes",  num_crashes,  cat_val)
   catln("")

   warn_cnt  <- countLog(num_warnings)
   err_cnt   <- countLog(num_errors)
   crash_cnt <- countLog(num_crashes)

   catln("")
   catln("* STATS *")
   catln("")
   cat_val("Total   ", results$total)
   cat_val("Cached  ", results$cached)
   cat_val("No Cache", results$total - results$cached)
   catln("")
   cat_val("Peptides Cached   ", results$peps_cached)
   cat_val("No Peptides Cached", results$cached - results$peps_cached)
   catln("")
   cat_val("Warnings", warn_cnt)
   cat_val("Errors  ",   err_cnt)
   cat_val("Crashes ",  crash_cnt)

   return(err_cnt == 0 && crash_cnt == 0)
}

usage <- function() {
   cat(sprintf(
"
Healthcheck scans each viewer in the database and attempts to load it with log output for warnings, errors and crashes

Usage: healthcheck [options...] [-h|--help]

Options are key-value pairs and include:
before=<YYYY-MM-DD>
after=<YYYY-MM-DD>
n=<max_viewers>
random=<on|off>
verbosity=<%s>

",
   paste(names(log_levels), collapse = "|")))
}

if(!interactive()) {
   args <- commandArgs(trailingOnly=TRUE)

   results <- NULL
   if(any(c("-h", "--help") %in% args)) {
      usage()
   }else if("all" %in% args){
      results <- test_health()
   }else {
      args <- as.list(sapply(args, function(arg) {
         l <- strsplit(arg, "=")[[1]]
         if(length(l) != 2) {
            usage()
            stop("arguments must be key-value pairs seperated by '='")
         }
         setNames(l[[2]], l[[1]])
      }, USE.NAMES = FALSE))

      if(!is.null(args$random)) {
         args$random <- args$random == "on"
      }
      if(!is.null(args$n)) {
         args$n <- as.numeric(args$n)
      }

      valid_opts <- c("before", "after", "n", "random", "verbosity")
      valid_args <- names(args) %in% valid_opts
      if(any(!valid_args)) {
         usage()
         stop("unrecognized option(s):", paste(names(args)[!valid_args], collapse = ","))
      }
      cat("Health check...\n")
      results <- do.call(test_health, as.list(args))
   }

   if(!is.null(results)) {
      if(!report_health(results))
         quit(status = 1)
   }
}

