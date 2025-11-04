#' Create a view to pass to `loadDataset`
#'
#' General Notes:
#'  - `duration` parameters may be NULL to signify no timeout
#'  - If clearMsg is NA, then sticky messages will have a duration of 10 instead of NULL
#'
#' @param notify   - A function that takes a message and duration, then returns a notification handle which may be passed to `clearMsg`
#' @param warn     - A function that takes a title, message, and duration and returns a notification handle which may be passed to `clearMsg`
#' @param error    - A function that takes a title and message
#' @param progress - A function that takes a progress value between 0 and 100
#' @param clearMsg - A function that takes a title and message or NA to indicate that the functionality is not available
createLoadDatasetView <- function(notify = NULL, warn = NULL, error = NULL, progress = NULL, clearMsg = NULL) {
   if(is.null(notify)) {
      notify <- function(msg, duration) showNotification(msg, duration = duration)
   }
   if(is.null(warn)) {
      warn <- function(title, msg, duration) showNotification(paste0(title, ": ", msg), duration = duration)
   }
   if(is.null(clearMsg)) {
      clearMsg <- function(msg_id) removeNotification(msg_id)
   }

   if(is.null(error)) {
      error <- function(title, msg) sendSweetAlert(title = title, text = msg, type = "error")
   }

   if(is.null(progress)) {
      progress <- function(val) updateProgressBar(id = "progressBar1", value = val)
   }

   list(
      notify = notify,
      warn = warn,
      error = error,
      progress = progress,
      clearMsg = clearMsg
   )
}

createNotifyProgress <- function(view) {
   function(msg = "", timeStart = NULL, progress = NULL){
      now <- proc.time()
      if(!is.null(timeStart)) {
         delta <- now - timeStart
         secs <- formatC(delta["elapsed"], digits = 2, format = "f")
         msg <- paste(msg, "in", secs, "seconds.")
      }

      if(nzchar(msg)) {
         view$notify(msg, duration = 4)
      }

      if(!is.null(progress)) {
         view$progress(progress)
      }

      return(now)
   }
}

createStickyNotify <- function(view) {
   function(msg, do){
      if(is.null(msg))
         return(do())

      msg_time <- if(is.function(view$clearMsg)) NULL else 10
      msg_id <- view$notify(msg, duration = msg_time)

      x <- do()

      if(is.function(view$clearMsg)) {
         view$clearMsg(msg_id)
      }

      return(x)
   }
}


dumpVars <- function(vars, labels = NULL, msg = NULL) {
   if(!is.null(msg)) cat(msg, strrep("-", nchar(msg)), sep = "\n")
   if(is.null(labels)) {
      labels <- names(vars)
   }else if(!is.null(names(vars))) {
      labels <- paste(names(vars), labels)
   }
   for(i in seq_along(vars)) {
      cat("\n")
      if(!is.null(labels)) cat(labels[[i]], ":\n", sep = "")
      cat(vars[[i]], sep = "\n")
   }
}

# References dataset${order, numSamples}
transformProteins <- function(dataset, proteins, view) {
   # Reorder Columns
   proteins$columns <- proteins$columns[, dataset$order, drop = FALSE]
   protNames <- colnames(proteins$columns)

   # Reorder Classes and Peptide counts
   proteins$numClasses <- length(proteins$classes)
   if(proteins$numClasses > 1) {
      classes <- proteins$classes

      # Assuming each column name begins with a class name. (validated in MosaicDF.r)
      # See MassPike's www/modules:
      #  siteQuant/lib/site_quant_iterator.php @ generateSchema()
      #  protein_quant/lib/protein_quant_iterator.php @ __constructor()

      classRegex <- paste0("^", str_vec_regex(classes))
      colClasses <- str_extract(protNames, classRegex)

      classKeep <- classes %in% colClasses
      pepCountCols <- str_subset(names(proteins$info), "Peptides$")
      if(!all(classKeep)) {
         classes <- classes[classKeep]
         proteins$numClasses <- sum(classKeep)

         proteins$info <- proteins$info[colnames(proteins$info) != pepCountCols[!classKeep]]
         pepCountCols <- pepCountCols[classKeep]
      }

      proteins$classNames <- colClasses
      proteins$classIDs <- match(colClasses, classes)
      pepCountOrder <- unique(proteins$classIDs)

      #Reorder
      proteins$info[pepCountCols] <- proteins$info[pepCountCols][pepCountOrder]
      pepCountInd <- names(proteins$info) %in% pepCountCols
      names(proteins$info)[pepCountInd]  <- pepCountCols[pepCountOrder]
      proteins$classes <- classes[pepCountOrder]
   }else {
      proteins$classNames <- rep(proteins$classes, dataset$numSamples)
      proteins$classIDs <- rep(1, dataset$numSamples)
   }

   return(proteins)
}

# Reference dataset${proteins$columns, proteins$classes, proteins$classNames, classesMap, names}
transformPeptides <- function(dataset, peptides, view, classMap, classedFactors) {
   if(is.null(peptides))
      return(NULL)

   protNames <- colnames(dataset$proteins$columns)
   protClasses <- dataset$proteins$classes

   peptides$info$Class <- as.character(peptides$info$Class)
   pepClasses <- unique(peptides$info$Class)

   # proteins$classes is a subset of pepClasses (dropping ignored columns)
   dumpClasses <- function(msg) dumpVars(list(Proteins = protClasses, Peptides = pepClasses), "classes", msg)
   if(!all(pepClasses %in% protClasses)) {
      dumpClasses("Unknown peptide class(es)")
      view$warn("Unknown peptide classes", paste(pepClasses, ", "), 4)
   }

   if(!any(protClasses %in% pepClasses)) {
      dumpClasses("No protein class found in peptide classes")
      view$error("Peptide class mismatch", "Peptide classes do not intersect protein classes")
      return(NULL)
   }

   protClassNames <- dataset$proteins$classNames

   pepNames <- colnames(peptides$columns)
   anyPepCol <- str_vec_regex(pepNames)
   protCols <- str_remove(protNames, sprintf("^%s[~_]", protClassNames)) # Don't match class names later
   protPepNames <- str_extract(protCols, anyPepCol)
   if(any(is.na(protPepNames))) {
      msg <- "Failed to map peptide column names"
      debugNames <- protNames
      debugNames[is.na(protPepNames)] <- paste(debugNames[is.na(protPepNames)], "***")
      dumpVars(list(Peptide = pepNames, Protein = debugNames), "names", msg)
      view$error(msg, sprintf("[%s]", paste(protNames[is.na(protPepNames)], collapse = ", ")))
      return(NULL)
   }

   pepColumns <- list()
   pepColIDs <- list()
   for(class in pepClasses) {
      classRows <- class == peptides$info$Class
      classCols <- class == protClassNames

      if(!any(classCols)) { # Drop unattested classes
         peptides$info     <- peptides$info[   !classRows, ]
         peptides$columns  <- peptides$columns[!classRows, ]
         next
      }

      columns <- peptides$columns[classRows, , drop = FALSE]

      # Reorder and rename
      classLabel <- class
      pepLabels <- NULL
      if(dataset$classesMap) {
         classOrder <- match(protPepNames[classCols], pepNames)

         # *Should* be impossible since dataset$classesMap is TRUE
         if(any(duplicated(classOrder))) {
            view$error("Invalid state", "Duplicates in peptide column selection")
            return(NULL)
         }

         columns <- columns[, classOrder, drop = FALSE]
         pepLabels <- dataset$names[classCols]

         classLabel <- classMap[[class]]
         peptides$info[classRows, "Class"] <- classLabel
      }

      # Normalize peptides to match proteins
      if(!is.null(classedFactors)) {
         classFactors <- unlist(classedFactors[[class]][colnames(columns)])
         if(length(classFactors) != ncol(columns)) {
            view$error("Dimension mismatch", paste("Unable to normalize peptides of class", classLabel))
            return(NULL)
         }
         columns <- t(t(columns) * classFactors)
      }

      # Finalize names (after indexing classedFactors)
      if(!is.null(pepLabels)) {
         colnames(columns) <- pepLabels
      }

      pepColIDs[[classLabel]] <- which(classCols)

      rowIDs <- as.character(seq(sum(classRows)))
      rownames(columns) <- rowIDs
      pepColumns[[classLabel]] <- columns

      peptides$info[classRows, "SumSN"] <- rowSums(columns, na.rm = TRUE)
      peptides$info[classRows, "data_row"] <- rowIDs
   }
   peptides$columnIDs <- pepColIDs
   peptides$columns <- pepColumns

   return(peptides)
}

dfltSrc <- do.call(createDatasetSource, gfySource)
dfltFmt <- gfyFormat

# Load and setup viewer data
loadDataset <- function(dataset, dataPath, idMaps, view = createLoadDatasetView(), srcs = NULL, format = NULL) {
   notifyProgress <- createNotifyProgress(view)
   stickyNotify <- createStickyNotify(view)

   if(is.null(srcs))
      srcs <- dfltSrc

   if(inherits(srcs, "tmtmosaic_src"))
      srcs <- list(srcs)

   if(is.null(format))
      format <- dfltFmt

   trySrcs <- function(dataset, type, err_handler) {
      data <- NULL
      for(src in srcs) {
         src_get_msg <- src[[paste(type, "msg", sep = "_")]]
         if(is.null(src_get_msg)) {
            src_get_msg <- function(...) NULL
         }

         data <- stickyNotify(src_get_msg(dataset), function() {
            tryCatch({
               src[[type]](dataset)
            },
               error = err_handler
            )
         })

         if(!is.null(data))
            break
      }

      dataset[[type]] <- data

      if(!is.null(data)) {
         # Run formatter
         data <- format[[type]](dataset)

         if(inherits(data, "tmtmosaic_src_state")) {
            if(is.null(dataset$state))
               dataset$state <- list()

            dataset$state[names(data$state)] <- data$state
            data <- data$result
         }

         dataset[[type]] <- data
      }

      return(dataset)
   }

   # Make ID map available
   dataset$state <- list(idMap = idMaps[[dataset$species]])

# Begin
   view$progress(0)
   timeStart <- proc.time()

# Get Protein Information ------------
   caught_err <- FALSE
   dataset <- trySrcs(dataset, "proteins", function(e) {
      cat(e$message, sep = "\n")
      view$error("Failed to load viewer data", e$message)
      caught_err <<- TRUE
      return(NULL)
   })

   if(is.null(dataset$proteins)) {
      if(!caught_err) {
         view$error("Failed to load viewer data", "Data missing or invalid format")
      }
      return(NULL)
   }

   timeProt <- notifyProgress(paste(dataset$idLabel, "Data Loaded"), timeStart, 30)

# Get Normalization Factors
   dataset <- trySrcs(dataset, "normalization", function(e) {
      view$warn("Failed to load normalization", e$message, 4)
      return(NULL)
   })

   hasFactors <- !is.null(dataset$normalization)

   timeNorm <- notifyProgress(ifelse(hasFactors, "Loaded normalization", "No normalization"),
                              if(hasFactors) timeProt else NULL,
                              40)

# Get Peptide Information
   dataset <- trySrcs(dataset, "peptides", function(e) {
      view$warn("Failed to load peptides", e$message, 4)
      return(NULL)
   })
   loadedPeptides <- !is.null(dataset$peptides)

   timePep <- notifyProgress(
      ifelse(loadedPeptides,
         "Peptide Data Loaded", 
         "Skipped Peptides"),
      if(loadedPeptides) timeNorm else NULL,
      70)

# Validate Quant
   ## Proteins
   quant <- as.mosaic_df(dataset$proteins, dataset$peptides, dataset$isSiteQuant, idMap = idMaps[[dataset$species]])
   dataset$proteins <- quant$proteins
   dataset$peptides <- quant$peptides
   rm(quant)

   if(is.null(dataset$proteins))
      return(NULL)

   dataset$proteins <- transformProteins(dataset, dataset$proteins, view)

   # Check plex boundaries, map user input classes to gfy classes
   classMap <- data.frame(userClass = dataset$classes, gfyClassID = dataset$proteins$classIDs) %>% unique()
   dataset$classesMap <- !any(apply(classMap, 2, duplicated))
   classMap$gfyClass <- dataset$proteins$classes[classMap$gfyClassID]
   if(dataset$classesMap) {
      classMap <- with(classMap, setNames(userClass, gfyClass))
   }else {
      cat(sep = "\n", "Class mismatch:", capture.output(print(classMap)))
   }

   # Set default ID
   if(is.na(dataset$initialID) || !(dataset$initialID %in% dataset$proteins$info$MosaicID)){
      view$warn("Default ID not found", sprintf("[%s] missing, showing random ID instead", dataset$initialID), 4)
      row <- sample(1:nrow(dataset$proteins$info), 1)
      dataset$initialID  <- dataset$proteins$info[[row, "MosaicID"]]
   }

   ## Normalization
   # Order factors and handle missing data
   normOpts <- dataset$normalization
   classedFactors <- NULL
   if(!is.null(normOpts)) {
      if(is.null(normOpts$factors)) {
         normOpts$factors <- rep(1, dataset$numSamples)
      }else {
         classedFactors <- normOpts$factors
         factors <- unlist(classedFactors[dataset$proteins$classes])
         normOpts$factors <- factors[dataset$order]
      }
      names(normOpts$factors) <- colnames(dataset$proteins$columns)

      # Translate row ID to MosaicID
      if(!is.na(normOpts$normRow)) {
         index <- dataset$proteins$info$UniprotID == normOpts$normRow
         normOpts$normRow <- ifelse(any(index), dataset$proteins$info[index, "MosaicID"][[1]], NA)
      }
   }else {
      # Set default factors for no normalization
      normOpts <- list(factors = NULL, normRow = NA, normBy = "none")
   }
   dataset$normalization <- normOpts

   ## Peptides
   if(loadedPeptides){
      dataset$peptides <- transformPeptides(dataset, dataset$peptides, view, classMap, classedFactors)

      if(!dataset$classesMap) {
         # No visual link between protein and peptide data
         n <- ncol(dataset$peptides$columns[[1]])
         pepColors <- brewer.pal(pmin(pmax(3, n), 12), "Set3")
         dataset$pepColors <- rep(pepColors, length.out = n)
         view$warn("Custom order breaks plex boundaries", "Peptides will use default colors, class names and structure.", duration = 10)
      }
   }

   timeQuant <- notifyProgress(paste("Validated quant data"), timePep, 80)

# GO Categories ----------------------
   allGO <- NULL
   if(!is.na(dataset$taxid)){
      allGO <- getGOSet(dataPath, dataset$proteins$info$UniprotID, dataset$taxid)

      loadedGO <- FALSE
      if(is.null(allGO)){
         msg <- paste0("Annotations For Species '", dataset$species, "' Unavailable")
      }else if (nrow(allGO) == 0) {
         warning <- sprintf("The species you selected for this dataset is %s. None of the UniprotID's match the %s GO dataset. Are you sure you picked the right species?",
                            dataset$species, dataset$species)
         view$warn("Annotation mismatch", warning, 4)
         allGO <- NULL
         msg <- paste0("No ID Matches In '", dataset$species, "' Annotation Dataset")
      }else {
         loadedGO <- TRUE
         msg <- sprintf("%s annotations loaded", dataset$species)
      }
      notifyProgress(msg, if(loadedGO) timeQuant else NULL, 100)
   }
   dataset$allGO <- allGO

# Finished -----------------------------
   notifyProgress("Finished Loading Data", timeStart, 100)

   # Clear state
   dataset$state <- NULL

   return(dataset)
}

getGOSet <- function(dataPath, UniprotIDs, taxid){
   go_path <- paste0(dataPath, "annotations.db")
   if(!file.exists(go_path))
      return(NULL)

   db <- dbConnect(SQLite(), go_path, flags = SQLITE_RO)
   on.exit(dbDisconnect(db))

   table <- paste("annotations", taxid, sep = "_")
   if(!dbExistsTable(db, table))
      return(NULL)

   # Escape ' literals
   ids <- unique(UniprotIDs)
   ind <- str_detect(ids, "'")
   ids[ind] <- str_replace_all(ids[ind], "'", "''")

   query <- sprintf("SELECT * FROM %s WHERE UniprotID IN (%s)", table, paste0("'", ids, "'", collapse = ","))
   GO <- dbGetQuery(db, query)

   attr(GO, "LastUpdated") <- dbGetQuery(db, "SELECT date(LastUpdated) FROM metadata WHERE TaxID=?", params = taxid)[[1]]

   return(GO)
}

### Database Fnxns
getMetadata <- function(key) parseMetadata(readMetadata(key))

readMetadata <- function(key){
   db <- dbConnect(SQLite(), sqlitePath, flags = SQLITE_RO)
   on.exit(dbDisconnect(db))

   query <- paste0("SELECT * FROM metadata WHERE Key = ?;")
   return(dbGetQuery(db, query, params = key))
}

parseMetadata <- function(metadata){
   if(length(metadata) == 0 || nrow(metadata) != 1) return(NULL)

   metadata <- as.list(metadata)
   for(target in names(JSON_Vars)){
      tryCatch({
         if(target %in% names(metadata) && !is.null(metadata[[target]]) && !is.na(metadata[[target]])) {
            as.type <- JSON_Vars[[target]]
            metadata[[target]] <- metadata[[target]] %>% fromJSON %>% as.type
         } else {
            # Initialize with safe defaults
            if(target == "ColumnIDs") metadata[[target]] <- numeric(0)
            else if(target == "ColumnClasses") metadata[[target]] <- character(0)
            else if(target == "ColumnNames") metadata[[target]] <- character(0)
            else if(target == "GroupIDs") metadata[[target]] <- numeric(0)
            else if(target == "GroupColors") metadata[[target]] <- character(0)
            else if(target == "GroupNames") metadata[[target]] <- character(0)
         }
      }, error = function(e) {
         cat("[ERROR] Failed to parse JSON for", target, ":", e$message, "\n")
         # Initialize with safe defaults as above
         if(target == "ColumnIDs") metadata[[target]] <- numeric(0)
         else if(target == "ColumnClasses") metadata[[target]] <- character(0)
         else if(target == "ColumnNames") metadata[[target]] <- character(0)
         else if(target == "GroupIDs") metadata[[target]] <- numeric(0)
         else if(target == "GroupColors") metadata[[target]] <- character(0)
         else if(target == "GroupNames") metadata[[target]] <- character(0)
      })
   }
   
   # Check critical variables before using them
   if(length(metadata$ColumnIDs) == 0) {
      cat("[ERROR] ColumnIDs is empty or not properly parsed\n")
      return(NULL)
   }
   
   if(length(metadata$ColumnClasses) == 0) {
      cat("[ERROR] ColumnClasses is empty or not properly parsed\n")
      return(NULL)
   }
   
   # Ensure ColumnClasses has enough elements for the indices in ColumnIDs
   if(max(metadata$ColumnIDs) > length(metadata$ColumnClasses)) {
      cat("[ERROR] ColumnIDs contains indices beyond the length of ColumnClasses\n")
      cat("[DEBUG] ColumnIDs:", paste(metadata$ColumnIDs, collapse=", "), "\n")
      cat("[DEBUG] ColumnClasses length:", length(metadata$ColumnClasses), "\n")
      return(NULL)
   }
   
   # Ensure ColumnNames has enough elements for the indices in ColumnIDs
   if(max(metadata$ColumnIDs) > length(metadata$ColumnNames)) {
      cat("[ERROR] ColumnIDs contains indices beyond the length of ColumnNames\n")
      cat("[DEBUG] ColumnIDs:", paste(metadata$ColumnIDs, collapse=", "), "\n")
      cat("[DEBUG] ColumnNames length:", length(metadata$ColumnNames), "\n")
      return(NULL)
   }
   
   # The same for GroupIDs
   if(max(metadata$ColumnIDs) > length(metadata$GroupIDs)) {
      cat("[ERROR] ColumnIDs contains indices beyond the length of GroupIDs\n")
      cat("[DEBUG] ColumnIDs:", paste(metadata$ColumnIDs, collapse=", "), "\n")
      cat("[DEBUG] GroupIDs length:", length(metadata$GroupIDs), "\n")
      return(NULL)
   }
   
   # Safe version of metadata list creation - with checks before indexing
   result <- tryCatch({
      list(
         key          = as.character(metadata$Key),
         ID           = as.numeric(metadata$QID),
         username     = as.character(metadata$Username),
         name         = as.character(metadata$Dataset),
         notes        = as.character(metadata$Notes),
         date         = as.character(metadata$Date),
         
         # Use length of ColumnIDs instead of NumSamples
         numSamples   = length(metadata$ColumnIDs),
         numGroups    = as.numeric(metadata$NumGroups),
         numPlexes    = as.numeric(metadata$NumPlexes),
         plexLevel    = as.numeric(metadata$PlexLevel), 
         server       = as.character(metadata$Server),
         species      = as.character(metadata$Species),
         initialID    = as.character(metadata$InitialID),
         areReps      = as.logical(metadata$AreReps),
         isSiteQuant  = as.logical(metadata$IsSiteQuant),
         
         # Safe indexing with error handling
         classes      = if(length(metadata$ColumnClasses) >= length(metadata$ColumnIDs)) 
                           metadata$ColumnClasses[metadata$ColumnIDs] 
                        else character(length(metadata$ColumnIDs)),
         
         names        = if(length(metadata$ColumnNames) >= length(metadata$ColumnIDs)) 
                           metadata$ColumnNames[metadata$ColumnIDs] 
                        else character(length(metadata$ColumnIDs)),
         
         groups       = if(length(metadata$GroupIDs) >= length(metadata$ColumnIDs)) 
                           metadata$GroupIDs[metadata$ColumnIDs] 
                        else numeric(length(metadata$ColumnIDs)),
         
         groupColors  = metadata$GroupColors,
         groupNames   = metadata$GroupNames,
         
         order        = metadata$ColumnIDs
      )
   }, error = function(e) {
      cat("[ERROR] Failed to create metadata list:", e$message, "\n")
      return(NULL)
   })
   
   if(is.null(result)) {
      return(NULL)
   }
   
   # Process groups safely
   result <- tryCatch({
      within(result, {
         unique_groups <- unique(groups)
         if(length(unique_groups) > 0) {
            groups <- match(groups, unique_groups)
            if(length(groups) > 0 && length(groupColors) >= max(groups)) {
               columnColors <- groupColors[groups]
            } else {
               # Default colors if indices are out of bounds
               columnColors <- rep("#CCCCCC", length(groups))
            }
         } else {
            columnColors <- rep("#CCCCCC", length(groups))
         }
      })
   }, error = function(e) {
      cat("[ERROR] Failed to process groups:", e$message, "\n")
      return(NULL)
   })
   
   # Add taxonomy ID if species exists in the configuration
   if(!is.null(result) && !is.null(result$species) && result$species %in% names(config$ANNOTATED_SPECIES)) {
      result$taxid <- config$ANNOTATED_SPECIES[[result$species]]
   } else {
      result$taxid <- NA
   }
   
   # Set up server connection if server is in configuration
   if(!is.null(result) && !is.null(result$server) && result$server %in% names(config$SERVERS)) {
      server_auth <- config$SERVERS[[result$server]]
      if(!is.null(server_auth)) {
         url <- paste0("https://", result$server, "/gfy/www/modules/api/v1")
         result$gfy.obj <- gfy_api(url, server_auth$User, server_auth$Key, server_auth$VerifySSL)
      }
   }
   
   # Set ID label based on isSiteQuant
   if(!is.null(result)) {
      result$idLabel <- ifelse(result$isSiteQuant, "Site", "Protein")
   }
   
   return(result)
}
