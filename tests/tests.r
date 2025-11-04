app_dir <- Sys.getenv("APP_DIR")

suppressPackageStartupMessages({
   source(paste0(app_dir, "/tests/common.r"))

   source(paste0(app_dir, "/IsoParser/global.r"))
   source(paste0(app_dir, "/IsoParser/lib/InitiateData.r"))
})

testsPath <- paste0(app_dir, "/tests/store/test_graph.csv")

# Check equivalence of `cols` in `a` and `b` based on an ID column, `key`, ignoring duplicate rows and row ordering
df_equiv <- function(a, b, key, cols = NULL, show_n = 3, eps = 0.0001) {
   trunc <- function(vals, n = 5) {
      if(length(vals) > n) {
         paste(paste(head(vals, n), collapse = ", "), "...")
      }else {
         paste(vals, collapse = ", ")
      }
   }

   inequiv <- c()
   if(is.null(cols)) {
      #require identical col names in both
      cols <- union(colnames(a), colnames(b))
   }
   cols <- c(cols, key)

   checkCols <- function(df, desc) {
      missing <- setdiff(cols, colnames(df))
      if(length(missing) != 0) {
         inequiv <<- c(inequiv, paste(desc, "df is missing compared column(s):", trunc(missing)))
      }
   }

   checkCols(a, "Left")
   checkCols(b, "Right")

   a <- unique(a)
   b <- unique(b)

   checkId <- function(df, desc) {
      dups <- duplicated(df[[key]])
      if(any(dups)) {
         dup_vals <- unique(df[[key]][dups])

         inequiv <<- c(inequiv, paste(desc, "df had", length(dup_vals), "duplicate id(s):", trunc(dup_vals)))
      }
   }

   checkId(a, "Left")
   checkId(b, "Right")

   a_row <- nrow(a)
   b_row <- nrow(b)
   if(a_row != b_row) {
      inequiv <- c(inequiv, sprintf("Differing number of rows: %d | %d", a_row, b_row))
   }

   if(length(inequiv) == 0) {
      crosscheckId <- function(x, y, desc_x, desc_y) {
         x_y_ids <- setdiff(x[[key]], y[[key]])
         if(length(x_y_ids) != 0) {
            msg <- sprintf("%s df is missing %d id(s) in %s df: %s", desc_y, length(x_y_ids), desc_x, trunc(x_y_ids))
            inequiv <<- c(inequiv, msg)
         }
      }

      crosscheckId(a, b, "Left", "Right")
      crosscheckId(b, a, "Right", "Left")

      ids <- as.character(intersect(a[[key]], b[[key]]))

      rownames(a) <- as.character(a[[key]])
      rownames(b) <- as.character(b[[key]])
      a <- a[ids, cols, drop = FALSE]
      b <- b[ids, cols, drop = FALSE]

      equiv <- a == b
      cols_num_equiv <- sapply(a, is.numeric) &
                        sapply(b, is.numeric)
      equiv[, cols_num_equiv] <- abs(a[cols_num_equiv] - b[cols_num_equiv]) <= eps

      row_inequiv <- as.numeric(rowSums(!equiv))
      if(sum(row_inequiv) != 0) {
         a[!equiv] <- paste0(a[!equiv], "*")
         b[!equiv] <- paste0(b[!equiv], "*")

         a <- a[row_inequiv != 0, , drop = FALSE]
         b <- b[row_inequiv != 0, , drop = FALSE]
         n_neq <- nrow(a)

         a$`_df_` <- "<-"
         b$`_df_` <- "->"
         a_neq_b <- rbind(a, b)

         col_order <- c("_df_", key, setdiff(cols, key))
         a_neq_b <- a_neq_b[order(a_neq_b[[key]]), col_order, drop = FALSE]
         rownames(a_neq_b) <- NULL

         more <- ""
         if(n_neq > show_n) {
            more <- "\n..."
            a_neq_b <- head(a_neq_b, show_n * 2)
         }
         msg <- paste0("Left df and Right df have ", n_neq, " differing row(s)\n",
                       paste(capture.output(print(a_neq_b)), collapse = "\n"),
                       more)
         inequiv <- c(inequiv, msg)
      }
   }
   return(inequiv)
}

runTests <- function(path = testsPath, types = c(), test_regex = NULL) {
   if(!file.exists(path)) {
      cat_hdr("No Tests")
      cat_val("File not found", path)
      return(NULL)
   }

   tests <- read.csv(path)

   log <- createLog()

   total  <- 0
   passed <- 0
   for(i in seq(nrow(tests))) {
      test <- tests[i, ]

      # Filter by type
      if(length(types) != 0 && !test$type %in% types) {
         next
      }
      # Filter by regex
      if(!is.null(test_regex) && regexpr(test_regex, test$name) == -1) {
         next
      }

      # Run test
      f <- getTest(test$type)
      if(is.function(f)) {
         total <- total + 1
         event <- f(test)
         if(!is.null(event)) {
            log <- logEvent(log, paste(test$name, "-", test$desc), event$severity, event$details)
         }else {
            passed <- passed + 1
         }
      }else {
         cat_val("Skipping unknown test type", test$type)
      }
   }

   #Print stuff
   printLog("DETAILS", log)

   num_events <- summarizeLog(log)
   printLog("SUMMARY",  num_events,  cat)

   cat_hdr("Results")
   cat_val("Tests run", total)
   cat_val("Tests passed", passed)
}

getTest <- function(type) {
   c(
      "load" = loadTest
      #"transform" = transformTest
   )[[type]]
}

loadTest <- function(test) {
   req_file <- function(desc) list(severity = "Invalid Test", details = paste("Missing", desc, "file"))

   store <- createStore(paste0(app_dir, "/tests/store/"))

   metadata <- store$read(c("metadata", paste0(test$name, ".db")))
   if(is.null(metadata))
      return(req_file("metadata"))

   proteins <- store$read(c(test$from, "proteins"))
   if(is.null(proteins))
      return(req_file("proteins"))

   norm <- store$read(c(test$from, "normalization"))
   if(is.null(norm))
      return(req_file("normalization"))

   peptides <- store$read(c(test$from, "peptides"))
   if(is.null(peptides))
      return(req_file("peptides"))

   mosaic_prot <- store$read(c(test$to, "proteins"))
   if(is.null(mosaic_prot))
      return(req_file("expected proteins"))

   mosaic_pep <- store$read(c(test$to, "peptides"))
   if(is.null(mosaic_pep))
      return(req_file("expected peptides"))

   metadata <- parseMetadata(metadata)
   norm <- norm_from_cache(norm)
   result <- NULL
   nop <- function(...) NULL
   view <- createLoadDatasetView(
      nop,
      function(title, msg, duration) result <<- list(severity = "WARN", details = paste(title, "-", msg)),
      function(title, msg, duration) result <<- list(severity = "ERROR", details = paste(title, "-", msg)),
      nop,
      NA
   )
   src <- createDatasetSource(
      proteins = function(...) proteins,
      peptides = function(...) peptides,
      normalization = function(...) norm
   )

   mosaic_actual <- tryCatch({
      loadDataset(metadata, "", NULL, view, src, gfyFormat)
   }, error = function(e){
      result <<- list(severity = "CRASH", details = e$message)
   })

   if(!is.null(result))
      return(result)

   ineq_prot <- df_equiv(mosaic_prot, mosaic_to_cache(mosaic_actual, "proteins"), "MosaicID", cols = colnames(mosaic_prot))
   if(length(ineq_prot) != 0)
      return(list(severity = "ERROR", details = paste0("\nProteins (Left = expected):\n", paste(ineq_prot, collapse = "\n"))))

   ineq_pep <- df_equiv(mosaic_pep, mosaic_to_cache(mosaic_actual, "peptides"), "PeptideId", cols = colnames(mosaic_pep))
   if(length(ineq_pep) != 0)
      return(list(severity = "ERROR", details = paste0("\nPeptides (Left = expected):\n", paste(ineq_pep, collapse = "\n"))))

   return(NULL)
}

if(!interactive()) {
   runTests()
}
