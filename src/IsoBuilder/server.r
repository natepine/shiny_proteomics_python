source("lib/Utils.r")

# Constants
name_length <- 15

options(shiny.maxRequestSize = 250 * 1024^2)

# Define server logic
function(input, output, session) {
   # Session variables
   values <- reactiveValues(
      #Dataset values
      usePreset = TRUE,
      limitClasses = TRUE,
      editClasses = FALSE,

      #Table values
      isAdmin = is.na(config$SUPER_USER_PASSWORD),
      queryAdminTable = 0,
      delete_key = NULL
   )
   metadata <- reactiveValues()
   lastQID <- NULL

   # Authenticate
   user_info <- secure_server(
      check_credentials = check_gfy_credentials(input$server)
   )

   # Dynamic UIs
   output$username <- renderUI({HTML(paste("<b>", req(user_info$name), "</b>"))})

   output$csvFormat <- renderTable({
      req(input$data_type == "CSV/TSV")
      is_site_quant <- input$is_site_quant

      ex <- gfyProteinExample(is_site_quant)
      return(cbind(list(Pattern = "Example Data"), ex))
   })

   output$csvPepFormat <- renderTable({
      req(input$data_type == "CSV/TSV")
      is_site_quant <- input$is_site_quant

      ex <- gfyPeptideExample(is_site_quant)
      return(cbind(list(Pattern = "Example Data"), ex))
   })

   output$ui_show_as_replicates <- renderUI({
      req(metadata$AreReps) # !NULL && TRUE
      checkboxInput("summarize_reps", "Summarize Reps", value = FALSE)
   })

   output$protChoiceUI <- renderUI({
      proteins <- req(values$proteins)
      selected <- metadata$InitialID
      if(is.null(selected)) {
         selected <- proteins$info$MosaicID[[1]]
      }
      selectizeInput("initialID",
         paste("Search for", values$idLabel),
         choices = list(Symbol = proteins$info$MosaicID),
         selected = selected,
         multiple = FALSE)
   })

   output$groupColorUI <- renderUI({
      groupNames <- req(metadata$GroupNames)
      presets <- isolate(req(values$groupColors))
      lapply(seq(length(groupNames)), function(i) {
         id <- paste0("groupColor", i)
         observeEvent(input[[id]], {
            values$groupColors[i] <- req(input[[id]])
         })

         colourInput(inputId = id, label = groupNames[i], presets[i])
      })
   })

   observe({
      updateTextInput(session, "num_groups", value = metadata$NumGroups)
   })
   observe({
      numSamples <- req(metadata$NumSamples)
      numGroups <- req(input$num_groups)
      iNumGroups <- as.integer(numGroups)
      if(is.na(iNumGroups)) {
         updateTextInput(session, "num_groups", value = NA)
      }else {
         validNum <- iNumGroups %>% pmax(1) %>% pmin(numSamples)
         if(validNum != numGroups) {
            updateTextInput(session, "num_groups", value = validNum)
         }
      }
   })

   observe({
      n <- req(metadata$NumGroups)
      preset <- req(input$colorPalPreset)

      if(!values$usePreset) {
         # Reloaded from database
         colors <- isolate(metadata$GroupColors)
      }else {
         colors <- switch(preset,
            Basic     = brewer.pal(n, "Set1"),
            Steve     = c("#F7FFA1", "#756027", "#45FF42", "#FF9100", "#1244D4", "#FF0000",
                              "#1F1F1F", "#FFA1FC", "#9903A1", "#BBFAF5", "#0E661C"),
            Pastel    = brewer.pal(n, "Pastel1"),
            Rainbow   = rainbow(n),
            Warm      = heat.colors(n),
            Earth     = terrain.colors(n),
            Vibrant   = topo.colors(n),
            Light     = cm.colors(n),
            Neon      = c("#E6FB04", "#00FF00", "#FF00CC", "#9D00FF", "#00FFFF", "#FF0000"),
            Dark      = c("#8A2BE2", "#000000", "#008B8B", "#006400", "#483D8B", "#B8C709"),
            Jade      = c("#94E8B4", "#72BDA3", "#5E8C61", "#4E6151", "#3B322C", "#CCAD99"),
            `Beach ball` = c("#25CED1", "#FCEADE", "#FF8A5B", "#EA526F", "#FFFFFF", "#FFC300"),
                        brewer.pal(n, "Set1")) # unknown preset
      }
      # Ensure n colors
      colors <- rep(colors, length.out = n)
      # Update barchart just once
      values$groupColors <- colors
      for(i in seq(n)){
         # Update inputs individually
         updateColourInput(session, paste0("groupColor", i), value = colors[i])
      }
   })

   asIDs <- function(vals) {
      match(vals, unique(vals))
   }

   lastTab <- function(tab) {
      id <- which(navTabs == tab)
      updateTabItems(session, "sidebarMenu", navTabs[[id - 1]])
   }

   nextTab <- function() {
      tab <- input$sidebarMenu
      id <- which(navTabs == tab)
      remainingSubmits <- paste0("submit_", navTabs[-(1:id)])

      shinyjs::enable(remainingSubmits[[1]])
      lapply(remainingSubmits[-1], shinyjs::disable) # Require (re)submission of following tabs
      
      updateProgressBar(session, "progress", value = unname(id), total = length(navTabs))
      updateTabItems(session, "sidebarMenu", navTabs[[id + 1]])
   }

   restoreInputs <- function(savedViewer) {
      restore <- setdiff(names(savedViewer), c("Date"))
      for(col in restore){
         metadata[[col]] <- savedViewer[[col]]
      }

      updateSelectizeInput(session, "initialID",selected = savedViewer$InitialID)
      updateNumericInput(  session, "num_groups",  value = savedViewer$NumGroups)
      updateTextInput(     session, "datasetName", value = savedViewer$Dataset)
      updateTextAreaInput( session, "notes",       value = savedViewer$Notes)
      updateCheckboxInput( session, "are_reps",    value = savedViewer$AreReps)

      # Show saved colors initially, then respond to future input
      values$usePreset <- FALSE
      observeEvent(input$colorPalPreset, values$usePreset <- TRUE, ignoreInit = TRUE, once = TRUE)

      showNotification("Loaded Saved Viewer")
   }

   # Process User Inputs
   lapply(navTabs[-1], function(tab) {
      observeEvent(
         input[[paste0("back_", tab)]],
         lastTab(tab)
      )
   })

   observeEvent(input$submit_source, {
      type <- req(input$data_type)

      isSiteQuant <- as.logical(input$is_site_quant)
      ignoreCache <- as.logical(input$ignore_cache)
      species     <- input$species
      if(is.null(species)) species <- NA

      idMap <- NULL
      if(!is.na(species)){
         idMap <- idMaps[[species]]
      }

      cache <- createCache(dataPath)

      # Configure data source
      QID <- NA
      if(type == "Server") {
         req(validateQID(input$QID))
         QID <- as.numeric(input$QID)
         server <- user_info$gfy_server
         cache_key <- mosaic_cache_key(QID, isSiteQuant, "proteins", server)

         src <- function() {
            if(!ignoreCache) {
               data <- cache$read(cache_key)

               if(!is.null(data)) {
                  return(data)
               }
            }

            data <- tryCatch({
               gfySourceProteins(QID, user_info$gfy.obj, isSiteQuant)
            },
            error = function(c) {
               showModal(modalDialog(
                  title = "Error | QID Not Found",
                  HTML(paste0(
                     "Quant ID could not be loaded, it may be archived or mistyped. Please check and try submitting again.",
                     "<h5>Error msg: ", c$message, "</h5>"
                  ))
               ))
               return(NULL)
            })
         }

         makeCacheData <- function(data) { # delay cache write until viewer creation
            function() {
               cache$write(data, cache_key)
               return(NULL)
            }
         }
      }else if(type == "CSV/TSV") {
         quantCSV <- req(input$quantCSV)

         validType <- function(path, desc) {
            if (endsWith(tolower(path), ".tsv")) {
               return("tsv")
            } else if (endsWith(tolower(path), ".csv")) {
               return("csv")
            } else {
               showModal(modalDialog(
                  title = "Error | Invalid File Type",
                  HTML(sprintf(
                     "%s could not be imported: type '%s' is not supported. Please check and try submitting again.",
                     desc, type
                  ))
               ))
               return(NULL)
            }
         }

         quantType <- req(validType(quantCSV$datapath, "Quant file"))
         genError <- function(desc, type) {
            function(e) {
               showModal(modalDialog(
                  title = sprintf("Error | %s Import Failed", type),
                  HTML(sprintf(
                     "%s %s could not be imported. Please check and try submitting again.<h5>Error msg: %s</h5>",
                     desc, type, e$message
                  ))
               ))
               return(NULL)
            }
         }

         quantPath <- quantCSV$datapath
         src <- function() {
            tryCatch({
               csvSourceProteins(quantPath, quantType)
            },
            error = genError("Quant", toupper(quantType)))
         }

         makeCacheData <- function(data) { # delay cache write until viewer creation
            function() {
               QID <- autoIncrID(sqlitePath)

               cache_key <- mosaic_cache_key(QID, isSiteQuant, "proteins", "Upload")
               cache$write(data, cache_key)
               unlink(quantPath, force = TRUE)

               return(QID)
            }
         }

         server <- "Upload"
      }

      formatError <- function(desc) function(e) {
         cat(e$message, sep = "\n")
         showModal(modalDialog(
            title = "Error | Invalid Format",
            HTML(paste0(
               desc, " could not be parsed",
               "<h5>", e$message, "</h5>"
            ))
         ))
      }

      # Get expression data
      proteins <- NULL
      data <- src()
      if(!is.null(data)) {
         proteins <- tryCatch({
            gfyFormatProteins(data, isSiteQuant)
         },
         parse_failed = function(w) {
            cat(w$message, sep = "\n")
            showModal(modalDialog(
               title = "Error | Parsing Failed",
               HTML(paste0(
                  "Data could not be parsed",
                  "<h5>Unable to locate scaled or sum columns</h5>"
               ))
            ))
            return(NULL)
         },
         error = formatError(paste(values$idLabel, "Quant data")))
      }
      req(proteins)

      cacheData <- makeCacheData(data)

      peptides <- NULL
      if(type == "CSV/TSV" && !is.null(input$pepQuantCSV)) {
         pepPath <- input$pepQuantCSV$datapath

         pepQuantType <- req(validType(pepPath, "Peptide Quant"))

         rawPeptides <- tryCatch({
            csvSourcePeptides(pepPath, pepQuantType)
         }, error = genError("Peptide Quant", toupper(pepQuantType)))
         req(rawPeptides)

         saveAll <- cacheData
         cacheData <- function() {
            QID <- saveAll()
            cache_key <- mosaic_cache_key(QID, isSiteQuant, "peptides", "Upload")
            cache$write(rawPeptides, cache_key)
            unlink(pepPath, force = TRUE)
            return(QID)
         }

         peptides <- tryCatch({
            gfyFormatPeptides(rawPeptides, isSiteQuant)
         },
         error = formatError("Peptide data"))
         req(peptides)
      }

      proteins <- tryCatch({
         as.mosaic_df(proteins, peptides, isSiteQuant, idMap = idMap)$proteins
      },
      error = function(e) {
         cat(e$message, sep = "\n")
         showModal(modalDialog(
            title = "Error | Invalid Format",
            HTML(paste0(
               "Data could not be parsed",
               "<h5>", e$message, "</h5>"
            ))
         ))
         return(NULL)
      })
      req(proteins)
      proteins$cacheData <- cacheData

      values$proteins      <- proteins
      values$editClasses   <- FALSE
      values$limitClasses  <- TRUE
      values$usePreset     <- TRUE

      if(!is.null(metadata$QID) && (is.na(QID) || QID != metadata$QID)){
         # Clear old metadata
         lapply(names(metadata), function(key) metadata[[key]] <- NULL)
      }

      restored <- FALSE
      if(type == "Server" && file.exists(sqlitePath)) {
         server <- user_info$gfy_server
         db <- dbConnect(SQLite(), sqlitePath, flags = SQLITE_RO)
         on.exit(dbDisconnect(db))
         if(!ignoreCache){
            savedViewer <- getViewerEntry(db, QID, isSiteQuant, server)
            if(!is.null(savedViewer)) {
               restoreInputs(savedViewer)
               restored <- TRUE
            }
         }else {
            key <- findKey(db, QID, isSiteQuant, server)
            if(!is.null(key)) {
               # Keep current key
               metadata$Key <- key
            }
         }
      }


      metadata$QID         <- QID
      metadata$Server      <- server
      metadata$IsSiteQuant <- isSiteQuant
      metadata$Species     <- species

      sep <- ifelse(isSiteQuant, "[~_]", "~")
      protNames <- colnames(proteins$columns)

      setColumns <- !restored
      classes <- rep("default", length(protNames))
      columns <- protNames

      classColRegex <- sprintf("^(%s)%s(.*)$", str_vec_regex(proteins$classes), sep)
      validColnames <- all(str_detect(protNames, classColRegex))
      if(!validColnames) {
         values$editClasses <- TRUE # Failed to parse classes, allow user to correct
         warning("Malformed protein colnames")
      }else {
         class_columns <- str_match(protNames, classColRegex)[, -1, drop = FALSE]
         classes <- class_columns[, 1]
         columns <- class_columns[, 2]
         if(restored) {
            if(length(classes) != length(metadata$ColumnClasses)) {
               setColumns <- TRUE
            }else {
               classMap <- unique(data.frame(classes, metadata$ColumnClasses))
               setColumns <- any(apply(classMap, 2, duplicated))
            }
         }
      }
      if(setColumns){
         metadata$ColumnNames    <- columns
         metadata$ColumnClasses  <- classes
         metadata$NumSamples     <- length(columns)
         metadata$NumPlexes      <- length(unique(classes))
         metadata$PlexLevel      <- length(unique(columns))
      }

      classLimits <- NA
      if(validColnames) classLimits <- asIDs(classes)
      values$classLimits <- classLimits

      if(metadata$NumPlexes * metadata$PlexLevel != metadata$NumSamples){
         values$editClasses <- TRUE # Failed to parse classes, allow user to correct
         showNotification("Failed to parse plex structure, please correct manually", duration = 0)
      }

      Use <- TRUE
      if(!is.null(metadata$ColumnIDs)) {
         Use <- seq(metadata$NumSamples) %in% metadata$ColumnIDs
      }

      if(is.null(metadata$GroupIDs)) {
         if(is.null(metadata$GroupNames)) {
            numGroups <- floor(sqrt(metadata$NumSamples)) #Guess Grouping
         }else {
            numGroups <- length(metadata$GroupNames)
         }
         metadata$GroupIDs <- sort(rep(seq(numGroups), length.out = metadata$NumSamples))
         metadata$NumGroups <- numGroups
      }

      column_df <- data.frame(
         Order = seq(metadata$NumSamples),
         Class = metadata$ColumnClasses,
         Name = metadata$ColumnNames,
         Group = metadata$GroupIDs,
         Use = Use,
         stringsAsFactors = FALSE)

      order <- metadata$ColumnIDs
      if(!is.null(order)){
         order <- c(order, setdiff(seq(metadata$NumSamples), order))
         column_df <- column_df[order, ]
      }

      group_df <- data.frame(Group = unique(column_df$Group))
      groupNames <- metadata$GroupNames
      if(is.null(groupNames) || length(groupNames) != nrow(group_df)){
         groupNames <- as.character(seq(nrow(group_df)))
      }
      group_df$`Group Name` <- groupNames

      class_df <- data.frame(Class = sort(unique(column_df$Class)))
      classNames <- sort(unique(metadata$ColumnClasses))
      if(is.null(classNames) || length(classNames) != nrow(class_df)){
         classNames <- symbolIDs(LETTERS, nrow(class_df))
      }
      class_df$`Class Name` <- classNames

      values$column_df <- column_df
      values$group_df <- group_df
      values$class_df <- class_df

      nextTab()
   })

   observeEvent(input$submit_name, {
      column_df <- hot_to_r(req(input$column_table))
      group_df <- hot_to_r(req(input$group_table))
      class_df <- hot_to_r(req(input$class_table))
      column_df$Group <- req(groupSeq())
      column_df$Use <- usedCols()

      passed <- tableStatus() == "Passed"
      if(!all(passed)){
         showModal(modalDialog(easyClose = TRUE,
            title = "Some requirements have not been met",
            HTML(paste(tableStatus()[!passed], collapse = "<br><br>"))
         ))
         return(NULL)
      }

      columnIDs <- column_df$Order[usedCols()]

      if(values$limitClasses) {
         classes <- column_df$Class[usedCols()]
         parsedClassStruct <- values$classLimits[columnIDs]

         u_classes <- unique(classes)
         userClassStruct <- match(classes, u_classes)
         if(!any(is.na(values$classLimits))) {
            classMap <- unique(data.frame(userClassStruct, parsedClassStruct))
            if(any(duplicated(classMap$parsedClassStruct))) {
               userWidth <- nchar(max(userClassStruct))
               parsedWdith <- nchar(max(parsedClassStruct))
               classWidth <- max(nchar(u_classes))
               remap <- paste(collapse = "<br>",
                  paste(
                     str_pad(userClassStruct, userWidth), "->",
                     str_pad(parsedClassStruct, parsedWdith, "right"), ":",
                     str_pad(classes, classWidth), "->", u_classes[parsedClassStruct]
                  )
               )
               showModal(modalDialog(
                  title = "Class structure mismatch",
                  HTML(paste(sep = "<br>",
                     "The supplied class structure seems to cross plex boundaries. <b>Ignore</b> this message if you believe it is an error and <b>Save Names</b> again.", 
                     "*Note* Peptides will not be renamed or reordered in the viewer.",
                     paste0("<pre>", remap, "</pre>")
                  )),
                  footer = list(
                     actionButton("ignoreClassLimit", "Ignore", `data-dismiss`="modal", `data-bs-dismiss`="modal"),
                     actionButton("limitClasses", "Apply Fix",  `data-dismiss`="modal", `data-bs-dismiss`="modal")
                  )
               ))
               return(NULL)
            }
         }
      }else {
         values$limitClasses <- TRUE
      }

      metadata$ColumnIDs <- columnIDs

      column_df <- column_df %>% arrange(Order)

      metadata$GroupIDs    <- column_df$Group
      metadata$ColumnNames <- column_df$Name
      metadata$GroupNames  <- group_df$`Group Name`
      metadata$PlexLevel   <- metadata$NumSamples / metadata$NumPlexes

      classes <- unique(column_df$Class)
      classMap <- setNames(classes, classes)
      classMap[class_df$Class] <- class_df$`Class Name`
      metadata$ColumnClasses <- unname(classMap[column_df$Class])

      nextTab()
   })

   observeEvent(input$submit_color, {
      metadata$GroupColors <- req(values$groupColors)
      metadata$InitialID <- req(input$initialID)

      nextTab()
   })

   observeEvent(input$submit_notes, {
      name <- input$datasetName
      notes <- input$notes
      req(validateNotes(name, notes))

      metadata$Dataset <- str_remove_all(name, "\'")
      metadata$Notes <- str_remove(notes, "\'")

      nextTab()
   })

   output$review_information_1 <- renderUI({
      req(metadata$Dataset, metadata$Notes)

      qidLabel <- paste0(values$idLabel, "Quant ID")
      output <- c(
         "Dataset Name" = metadata$Dataset,
         setNames(nm = qidLabel, object = metadata$QID), #Hidden when NA
         "User" = user_info$name,
         "Species" = metadata$Species, #Hidden when NA
         "Server" = metadata$Server,
         "User" = user_info$name,
         "Species" = metadata$Species, #Hidden when NA
         "Notes" = metadata$Notes)

      if(is.na(output[qidLabel])) {
         output <- output[!names(output) %in% c("Server", qidLabel)]
      }

      HTML(keyValueStr(output, boldKeys = TRUE))
   })

   output$review_information_2 <- renderUI({
      req(metadata$NumSamples, metadata$PlexLevel, metadata$NumPlexes, metadata$NumGroups)
      proteins <- req(values$proteins)

      peps <- str_subset(colnames(proteins$info), "Peptides")
      output <- c(
         setNames(nm = paste0("Num. ", values$idLabel, "s"), nrow(proteins$info)),
         "Num. Peptides" = paste(colSums(proteins$info[peps]), collapse = ", "),
         "Num. Samples"  = metadata$NumSamples,
         "Num. Groups"   = metadata$NumGroups,
         "Num. Plexes"   = metadata$NumPlexes,
         "Plexing Level" = metadata$PlexLevel)
      HTML(keyValueStr(output, boldKeys = TRUE))
   })

   observeEvent(input$submit_review, {
      finalData <- reactiveValuesToList(metadata)
      finalData$Username <- user_info$name

      req(validateFinal(finalData))

      if(is.na(finalData$QID) && !is.null(values$proteins$cacheData)) {
         # Caching data before updating DB to get a QID for DB
         finalData$QID <- values$proteins$cacheData()
         values$proteins$cacheData <- NULL
      }

      stateKey <- saveToDatabase(finalData, sqlitePath)
      if(!is.null(stateKey)){
         # Successfully saved metadata
         values$stateKey <- stateKey
         updateProgressBar(session, "progress", value = 6, total = 6)
         shinyjs::show("finished_box")
         showNotification("Saved Viewer!")

         # Cache viewer data
         if(!is.null(values$proteins$cacheData)) {
            values$proteins$cacheData()
            values$proteins$cacheData <- NULL
         }

         # Update adminTable
         values$queryAdminTable <- values$queryAdminTable + 1
      }
   })

   output$viewer_link <- renderUI({
      key <- req(values$stateKey)
      viewerLink <- paste0("/", config$VIEWER_SUBDIR, "/?sessid=", key)

      HTML(
         "Your viewer was saved as ", key,
         "<br>",
         '<a href="', viewerLink, '" target="_blank">Go To Viewer</a>'
      )
   })

   observe({
      isSiteQuant <- as.logical(input$is_site_quant)
      req(!is.null(input$is_site_quant))

      idLabel <- "Protein"
      status <- "danger"
      if(isSiteQuant) {
         idLabel <- "Site"
         status <- "success"
      }
      session$sendCustomMessage(type = "setLabelText",
         message = c("quantCSV", paste(idLabel, "Quant")))

      updateTextInput(session, "QID",
         label = paste0(idLabel, "Quant ID")
      )
      values$idLabel <- idLabel

      req(!identical(config$ANNOTATED_SPECIES, NA))
      selected <- input$species
      updateRadioGroupButtons(inputId = "species",
         label = "Species",
         choices = names(config$ANNOTATED_SPECIES),
         selected = selected,
         status = status
      )
   })

   usedCols <- reactive({
      # Clean garbage input
      column_df <- hot_to_r(req(input$column_table))
      keep <- as.logical(column_df$Use)
      replace(keep, is.na(keep), FALSE)
   })

   groupSeq <- reactive({
      #Condense groupIDs range from 1:numSamples (sparse) to indices (sequential)
      column_df <- hot_to_r(req(input$column_table))
      groups <- unique(column_df$Group)
      groupMap <- setNames(seq(length(groups)), sort(groups))
      groupMap[as.character(column_df$Group)]
   })

   tableStatus <- reactive({
      col_df <- req(hot_to_r(input$column_table))
      group_df <- req(hot_to_r(input$group_table))
      class_df <- req(hot_to_r(input$class_table))

      col_df$Use <- usedCols()

      validateTables(col_df, group_df, class_df, name_length)
   })

   observe({
      metadata$AreReps <- ifelse(is.null(input$are_reps), FALSE, as.logical(input$are_reps))
   })

   getColumns <- function() {
      column_df <- hot_to_r(req(input$column_table))
      column_df$Use <- used <- usedCols()
      values$column_df <- column_df

      #Drop unused columns
      return(column_df[used, ])
   }

   setColumn <- function(column, data) {
      #Drop unused columns
      values$column_df$Use <- used <- usedCols()
      values$column_df[used, column] <- data
   }

   observeEvent(input$ignoreClassLimit, values$limitClasses <- FALSE)

   observeEvent(input$limitClasses, {
      column_df <- getColumns()

      parsedClassStruct <- values$classLimits[column_df$Order]
      u_classes <- unique(column_df$Class)

      setColumn("Class", u_classes[parsedClassStruct])
   })

   observeEvent(input$setGroups, {
      numGroups <- as.integer(req(input$num_groups))
      column_df <- getColumns()

      groups <- rep(seq(numGroups), length.out = nrow(column_df))
      if(input$seqGroups){
         groups <- sort(groups)
      }

      metadata$NumGroups <- numGroups
      setColumn("Group", groups)
   })

   observeEvent(input$sortTable, {
      sort_by <- req(input$sortBy)

      column_df <- hot_to_r(req(input$column_table))
      column_df$Use <- used <- usedCols()
      values$column_df <- column_df[order(used, decreasing = TRUE), ]
      column_df <- column_df[used, ]

      if(sort_by %in% c("Class", "Order") || input$sortRows) {
         index <- order(column_df[[sort_by]])
         values$column_df$Use <- used <- usedCols()
         values$column_df[used, ] <- column_df[index, ]
      }else {
         setColumn(sort_by, sort(column_df[[sort_by]]))
      }
   })

   forceRowSort <- reactive(!is.null(input$sortBy) && input$sortBy %in% c("Order", "Class"))

   observeEvent(forceRowSort(), {
      if(forceRowSort()) {
         updateCheckboxInput(session, "sortRows", value = TRUE)
         shinyjs::disable("sortRows")
      }else {
         shinyjs::enable("sortRows")
      }
   })

   uniqueNames <- function(names) {
      if(any(duplicates <- duplicated(names))) {
         for(dup in unique(names[duplicates])) {
            matches <- names == dup
            prefixes <- paste0("_", seq(sum(matches)))
            names[matches] <- str_c(str_sub(names[matches], 1, name_length - nchar(prefixes)), prefixes)
         }
      }
      return(names)
   }

   observeEvent(input$makeUnique, {
      group_df <- hot_to_r(req(input$group_table))
      class_df <- hot_to_r(req(input$class_table))
      column_df <- getColumns()

      name <- uniqueNames(column_df$Name)

      group_df$`Group Name` <- uniqueNames(group_df$`Group Name`)
      class_df$`Class Name` <- uniqueNames(class_df$`Class Name`)

      all_names <- list(name, group_df$`Group Name`, class_df$`Class Name`)

      duplicates <- sapply(all_names, duplicated) %>% unlist
      all_names <- unlist(all_names)
      if(any(duplicates)) {
         showNotification(duration = 0,
            paste0("Some names are still duplicated: ",
               paste0(unique(all_names[duplicates]), collapse = ", ")
         ))
      }else{
         showNotification("Successfully made names unique.", duration = 4)
      }

      setColumn("Name", name)
      values$group_df <- group_df
      values$class_df <- class_df
   })

   trimNames <- function(names) {
      too_long <- nchar(names) > name_length
      if(any(too_long)) {
         names[too_long] <- str_sub(names[too_long], 1, name_length)
      }
      return(names)
   }

   observeEvent(input$trimNames, {
      group_df <- hot_to_r(req(input$group_table))
      class_df <- hot_to_r(req(input$class_table))
      column_df <- getColumns()

      name                  <- trimNames(column_df$Name)
      group_df$`Group Name` <- trimNames(group_df$`Group Name`)
      class_df$`Class Name` <- trimNames(class_df$`Class Name`)

      setColumn("Name", name)
      values$group_df <- group_df
      values$class_df <- class_df
   })

   observeEvent(input$autoName, {
      group_df  <- hot_to_r(req(input$group_table))
      column_df <- getColumns()
      groupIDs <- req(groupSeq())[usedCols()]

      autoNames <- group_df$`Group Name`[groupIDs]
      if(metadata$NumPlexes > 1){
         class_df <- hot_to_r(req(input$class_table))
         classMap <- setNames(class_df$`Class Name`, class_df$Class)
         autoNames <- paste(classMap[column_df$Class], autoNames, sep = "~")
      }

      groupIDs <- ave(rep(1, nrow(column_df)), groupIDs, FUN=cumsum)
      setColumn("Name", paste(autoNames, groupIDs, sep = "_"))
   })

   output$tableValidation <- renderUI({
      status <- tableStatus()
      textBold <- c(
         "unique" = "Names must be unique",
         "length" = "Limit names to 15 characters"
      )
      out <- sapply(names(status), function(test) {
         testStatus <- status[[test]]

         color <- "text-success"
         icon <- "ok"
         hover <- ""
         if(testStatus != "Passed") {
            color <- "text-danger"
            icon <- "remove"
            hover <- "<sup class='h5'>Hover for details</sup>"
         }

         paste0(
            "<span data-toggle=\"tooltip\" data-html=\"true\" data-placement=\"right\" title=\"", testStatus,"\">",
               "<b class='", color, "'>", textBold[[test]], "</b>",
               " ",
               "<span class='glyphicon glyphicon-", icon, " ", color, "'></span>",
            "</span>",
            hover
         )
      })
      session$sendCustomMessage(type = "clearTooltips", message = " ")
      tags$div(HTML(paste(out, collapse = "</br>")))
   })

   # observeEvent(input$sendemail, {
   #    flag <- FALSE
   #    if (flag) {
   #       showNotification("ERROR: Email isn't formatted correctly.", duration = 4)
   #    } else {
   #       body <- list(
   #             metadata$QID,
   #             "User"         = user_info$name,
   #             "Dataset Name" = metadata$Dataset,
   #             "Species"      = metadata$Species,
   #             "Server"       = metadata$Server,
   #             "Notes"        = metadata$Notes,
   #             "Link" = paste0("/", config$VIEWER_SUBDIR, "/?sessid=", values$stateKey)
   #       )
   #       names(body)[1] <- ifelse(metadata$IsSiteQuant, "SiteQuant ID", "Protein Quant ID")
   #       body <- keyValueStr(body, collapse = "\n")
   #       to <- paste0("<", input$useremail, ">")
   #       from <- "<noreply.tmtphospho@gmail.com>"
   #       subject <- "TMTPhospho Link"
   #       tryCatch({
   #          sendmail(from, to, subject, body, control = list(smtpServer = "localhost"))
   #          showNotification(paste0("Successfully sent email with link to ", input$useremail, ". If not seen, check your spam folder and click 'allow'."), duration = 8)
   #       },
   #       error = function(e){
   #          print(e)
   #       })
   #    }
   # })

   # Outputs:
   output$column_table <- renderRHandsontable({
      column_df <- req(values$column_df)

      order <- column_df$Order - 1 #range ajusted for JS
      column_df <- column_df %>% arrange(Order)

      hot <- rhandsontable(column_df, width = 600, manualRowMove = order,
            allowEmpty = FALSE, allowInsertRow = FALSE, allowInsertColumn = FALSE) %>%
         hot_col("Order", readOnly = TRUE) %>%
         hot_col("Name", width = 140, halign = "htCenter") %>%
         hot_col("Class", width = 140, readOnly = !values$editClasses, halign = "htCenter") %>%
         hot_col("Use", halign = "htCenter") %>%
         hot_table(highlightCol = TRUE, highlightRow = TRUE) %>%
         hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE) %>%
         hot_col("Group", type = "numeric") %>%
         hot_validate_numeric("Group", min = 1, max = metadata$NumSamples)

      return(hot)
   })

   observe({
      column_df <- hot_to_r(req(input$column_table))
      groups <- sort(unique(column_df$Group[usedCols()]))

      old_df <- hot_to_r(isolate(input$group_table))
      group_df <- data.frame(Group = groups) %>%
         mutate(`Group Name` = Group)

      if(!is.null(old_df) && 
         (nrow(old_df) != length(groups) || any(old_df$Group != groups))) {
         map <- setNames(old_df$`Group Name`, old_df$Group)
         matched <- group_df$Group %in% old_df$Group
         group_df[matched, "Group Name"] <- map[as.character(group_df[matched, "Group"])]

         metadata$NumGroups <- length(groups)
         values$group_df <- group_df
      }
   })

   output$group_table <- renderRHandsontable({
      if(is.null(values$group_df))
         return(NULL)

      rhandsontable(values$group_df, allowEmpty = FALSE, allowInvalid = FALSE,
            allowInsertRow = FALSE, allowInsertColumn = FALSE) %>%
         hot_col("Group", width = 140, halign = "htCenter", readOnly = TRUE) %>%
         hot_col("Group Name", width = 140, halign = "htCenter") %>%
         hot_table(highlightCol = TRUE, highlightRow = TRUE) %>%
         hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
   })

   observe({
      column_df <- hot_to_r(req(input$column_table))
      classes <- sort(unique(column_df$Class[usedCols()]))
      
      old_df <- hot_to_r(isolate(input$class_table))
      class_df <- data.frame(Class = classes) %>%
         mutate(`Class Name` = Class)

      if(!is.null(old_df)) {
         map <- setNames(old_df$`Class Name`, old_df$Class)
         matched <- class_df$Class %in% old_df$Class
         class_df[matched, "Class Name"] <- map[class_df[matched, "Class"]]
      }

      metadata$NumPlexes <- nrow(class_df)
      values$class_df <- class_df
   })

   output$class_table <- renderRHandsontable({
      if(is.null(values$class_df))
         return(NULL)

      rhandsontable(values$class_df, allowEmpty = FALSE, allowInvalid = FALSE,
            allowInsertRow = FALSE, allowInsertColumn = FALSE) %>%
         hot_col("Class", width = 140, halign = "htCenter", readOnly = TRUE) %>%
         hot_col("Class Name", width = 140, halign = "htCenter") %>%
         hot_table(highlightCol = TRUE, highlightRow = TRUE) %>%
         hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
   })

   plotBarchart <- function(columns, groupIDs, groupNames, colNames, colClasses, numPlexes, areReps, protein) {
      proteins <- req(values$proteins)
      if(is.null(protein) ||
         !(protein %in% rownames(proteins$columns) ||
           protein %in% seq(nrow(proteins$columns)))){
            protein <- 1
      }

      groupIDs <- req(groupIDs)[req(columns)] %>% asIDs()
      groupNames <- req(groupNames)[groupIDs]

      colColors <- req(values$groupColors)[groupIDs]
      colNames <- req(colNames)[columns]
      colClasses <- req(colClasses)[columns]

      req(!is.null(areReps))
      summarize_reps <- FALSE
      if(areReps){
         summarize_reps <- ifelse(is.null(input$summarize_reps), FALSE, input$summarize_reps)
      }

      heights <- proteins$columns[protein, columns]

      grouping <- groupColumns(colClasses, numPlexes, groupNames, areReps)
      g <- tmtBarchart(heights, colColors, colNames, grouping, labelGroups = areReps,
         summarize = summarize_reps, title = proteins$info[protein, "GeneSymbol"])

      ggplotly(g, tooltip = "text") %>%
         plotlyDefaults()
   }

   barchart <- reactive({
      plotBarchart(
         columns     = metadata$ColumnIDs,
         groupIDs    = metadata$GroupIDs,
         groupNames  = metadata$GroupNames,
         colNames    = metadata$ColumnNames,
         colClasses  = metadata$ColumnClasses,
         numPlexes   = metadata$NumPlexes,
         areReps     = metadata$AreReps,
         protein     = input$initialID
      )
   })

   output$protein_plot <- renderPlotly(barchart())
   output$plot_review  <- renderPlotly(barchart())

   output$plot_preview_names <- renderPlotly({
      proteins <- req(values$proteins)
      group_df  <- hot_to_r(req(input$group_table))
      column_df <- hot_to_r(req(input$column_table))
      column_df$Group <- req(groupSeq())
      shiny::validate(need(all(tableStatus() == "Passed"), "Input does not meet requirements"))

      req(ncol(proteins$columns) == nrow(column_df))

      columns <- column_df$Order[usedCols()]
      column_df <- column_df %>% arrange(Order)

      plotBarchart(
         columns     = columns,
         groupIDs    = column_df$Group,
         groupNames  = group_df$`Group Name`,
         colNames    = column_df$Name,
         colClasses  = column_df$Class,
         numPlexes   = metadata$NumPlexes,
         areReps     = metadata$AreReps,
         protein     = metadata$InitialID
      )
   })

   adminTable <- reactive({
      #Update when database changes
      values$queryAdminTable
      getAdminTableData(sqlitePath)
   })

   visibleAdminTable <- reactive({
      #Check Admin permissions
      adminTable <- adminTable()
      if(!is.null(adminTable) && !values$isAdmin){
         adminTable <- adminTable %>% filter(Username == user_info$name)
      }
      adminTable
   })

   observe({
      values$isAdmin <- is.na(config$SUPER_USER_PASSWORD) ||
         input$adminpassword == config$SUPER_USER_PASSWORD
   })

   observeEvent(input$delete_viewer, {
      session$sendCustomMessage(type = "resetValue", message = "delete_viewer")
      row <- as.numeric(input$delete_viewer)
      adminTable <- req(visibleAdminTable())
      if(row > 0 && row <= nrow(adminTable)) {
         showModal(modalDialog(easyClose = TRUE,
            title = "Confirm Viewer Deletion",
            tags$div(style = "font-size:14pt",
               adminTable[row, c("Dataset", "QID", "NumSamples", "Notes", "Species", "Server")] %>%
                  as.list %>%
                  keyValueStr(boldKeys = TRUE) %>%
                  HTML("<hr>")
            ),
            fluidRow(align = "center",
               column(6, modalButton("Cancel")),
               column(6,
                  HTML("<button type='button' class='btn btn-danger' data-dismiss=\"modal\" data-bs-dismiss=\"modal\" ",
                     "onclick='Shiny.setInputValue(\"confirm_delete_viewer\", \"1\");'>Delete</button>")
               )
            ),
            footer = NULL
         ))
         values$delete_key <- adminTable[row, "Key"]
      }
   })

   observeEvent(input$confirm_delete_viewer, {
      session$sendCustomMessage(type = "resetValue", message = "confirm_delete_viewer")
      key <- req(values$delete_key)
      deleteViewer(key, sqlitePath)
      values$queryAdminTable <- values$queryAdminTable + 1
      values$delete_key <- NULL
   })

   output$admintable <- renderDataTable({
      adminTable <- visibleAdminTable()
      shiny::validate(need(!is.null(adminTable), "The viewer database is empty, go through the editor to populate it!"))
      shiny::validate(need(nrow(adminTable) > 0, "You don't have any viewers yet, go through the editor to add one!"))
      isAdmin <- values$isAdmin

      adminTable$AreReps <- as.logical(adminTable$AreReps)
      adminTable$IsSiteQuant <- as.logical(adminTable$IsSiteQuant)

      adminTable <- adminTable %>% mutate(
         Link = paste0("<b><a href=/", config$VIEWER_SUBDIR,
         "/?sessid=", Key, " target=\"_blank\" style=\"color:",
         ifelse(IsSiteQuant,
            "#76AE83\">TMT Mosaic PTMs",
            "#214E98\">TMT Mosaic"),
         "</a></b>"),
         Delete = paste0("<button id='a' type='button' class='btn btn-danger action-button' ",
            "onclick='Shiny.setInputValue(\"delete_viewer\", \"", 1:nrow(adminTable), "\");'></button>"))

      cols <- c("Date", "QID", "Username", "Dataset", "Link", "Key", "Notes", "Server", "Species",
         "NumSamples", "NumGroups", "PlexLevel", "AreReps", "IsSiteQuant", "InitialID", "Delete")
      
      # Hide Species column if no species are defined
      if(all(is.na(adminTable$Species))){
         cols <- setdiff(cols, "Species")
      }

      # add colors:
      colors <- brewer.pal(9, "Set1")
      styledCols <- list(
         "QID" = 1,
         "Username" = 2,
         "Dataset" = 3,
         "NumSamples" = 4,
         "NumGroups" = 5,
         "PlexLevel" = 7)

      if(!isAdmin){
         cols <- setdiff(cols, "Username")
         styledCols$Username <- NULL
      }

      dt <- datatable(
         adminTable[cols],
         escape = FALSE, 
         options = list(
            initComplete = datatableTheme(),
            order = list(which(cols == "Date") - 1, "desc"),
            scrollX = TRUE
         ),
         selection = 'none',
         rownames = FALSE) %>%
         formatStyle("NumSamples", backgroundColor = "#EBEED6")

      for(col in names(styledCols)){
         dt <- formatStyle(dt, col, fontWeight = "bold", color = colors[styledCols[[col]]])
      }

      if(isAdmin){
         dt <- dt %>%
         formatStyle("Username", backgroundColor = "#EBEED6")
      }

      return(dt)
   })
}