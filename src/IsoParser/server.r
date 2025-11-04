#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

## Load Sources
src <- function(local_path, ...) source(paste0(app_dir, "/IsoParser/", local_path), ...)
src("lib/GetData.r")
src("lib/InitiateData.r")
src("lib/Outputs.r")

# Define server logic
function(input, output, session) {
   metadata <- reactive({
      queryStr <- parseQueryString(session$clientData$url_search)
      
      # Check if sessid exists in the query string
      if (!("sessid" %in% names(queryStr))) {
         sendSweetAlert(
               title = "Error: No key supplied",
               text = "The url does not contain a viewer key. Please try again.",
               type = "error"
            )
         return(NULL)
      }
      
      key <- queryStr[["sessid"]]
      
      # Additional safety check
      if (is.null(key) || key == "") {
         sendSweetAlert(
               title = "Error: No key supplied",
               text = "The url does not contain a viewer key. Please try again.",
               type = "error"
            )
         return(NULL)
      }

      metadata <- getMetadata(key)
      if (is.null(metadata)) {
         sendSweetAlert(
               title = "Error: Viewer not found",
               text = "No viewer was found for the key supplied in the url.",
               type = "error"
            )
         return(NULL)
      }

      #Log current key for troubleshooting
      warning(paste(sep = "\n\n", "",
         "*** VIEWER KEY ***",
         key,
         "******************", ""))

      return(metadata)
   })

   observeEvent(metadata(), {
      # Add a validation check to ensure metadata is not NULL
      req(metadata())
      
      # Store query string in a local variable
      queryStr <- parseQueryString(session$clientData$url_search)
      
      session$sendCustomMessage(type = "setDatasetName", message = metadata()$name)

      if(metadata()$isSiteQuant){
         updateRadioGroupButtons(session, "numNearest", "Number of Closest Sites")
         updatePrettyCheckbox(session, "show_volcano_tables", "Show Site Tables")
         session$sendCustomMessage("setPTM", list())
      }

      # Check if species exists in idMaps before accessing
      if(!metadata()$species %in% names(idMaps)){
         # Check for unlisted species, for backwards compatibility if ANNOTATED_SPECIES changes
         idMaps[[metadata()$species]] <- loadIdMap(metadata()$taxid)
      }

      cache <- createCache(dataPath)
      srcs <- list(cacheSrc(cache))
      if(!is.null(metadata()$gfy.obj))
         srcs[[2]] <- gfySrcToCache(cache)

      ## Session Constants
      dataset <- tryCatch({
         loadDataset(metadata(), dataPath, idMaps, srcs = srcs)
      }, error = function(e) {
         cat(e$message, sep = "\n")
         sendSweetAlert(
            title = "Error loading viewer",
            text = e$message,
            type = "error"
         )
         NULL
      },
      parse_failed = function(w) {
         cat(w$message, sep = "\n")
         sendSweetAlert(
            title = "Error | Parsing Failed",
            text = paste0(
               "Data could not be parsed: Unable to locate scaled or sum columns"
            ),
            type = "error"
         )
         NULL  # Return NULL on parse failure
      })
      req(dataset)

      ## Set UI state
      shinyjs::hide(selector = ".content > .progress-group")
      shinyjs::show(selector = ".tab-content")

      if(!is.null(dataset$peptides$info)){
         shinyjs::show(selector = "li:has(a[href$='#shiny-tab-peptides'])")
         shinyjs::show("dtable_pep_download")
         #source("tab/Peptides.r")
      }

      if(file.exists(paste0(dataPath, "bioplex.rds"))) {
         if(!is.na(dataset$taxid) && dataset$taxid == 9606) {
            shinyjs::show(selector = "li:has(a[href$='#shiny-tab-bioplex'])")
         }
      }

      if(!is.na(config$ADMIN_EMAIL)){
         shinyjs::show(selector = "li:has(a[href$='#shiny-tab-contact'])")
      }

      ## Reactive Values
      values <- reactiveValues(
         activeMosaicID = dataset$initialID,
         normProtCols = dataset$proteins$columns,
         normPepCols = dataset$peptides$columns
      )

      # Check if "selected" exists in queryStr before accessing it
      if("selected" %in% names(queryStr)) {
         activeMosaicID <- findMosaicID(dataset$proteins$info, queryStr[["selected"]])
         if(!is.null(activeMosaicID))
            values$activeMosaicID <- activeMosaicID
         rm(activeMosaicID)
      }

      proteinTMT <- reactive({
         # Add validation to ensure normProtCols has the activeMosaicID
         req(values$activeMosaicID %in% rownames(values$normProtCols))
         values$normProtCols[values$activeMosaicID, ]
      })

      GO <- goEnrichment(dataset$allGO,
         dataset$proteins$info %>%
            transmute(
               UniprotID,
               uniqueID = MosaicID,
               displayNames = GeneSymbol),
         urls = c(
            "GO" = "https://amigo.geneontology.org/amigo/term/%s",
            "KEGG" = "https://www.kegg.jp/entry/%s",
            "Reactome" = "https://reactome.org/PathwayBrowser/#/%s",
            "Interpro" = "https://www.ebi.ac.uk/interpro/entry/InterPro/%s/"
         ))

      if(!is.null(GO)){
         shinyjs::show(selector = "li:has(a[href$='#shiny-tab-go'])")
      }

      output$idUI <- renderUI({
         label <- ifelse(dataset$isSiteQuant, "SQID:", "PQID:")
         tagList(
            div(paste("Current Session", label, dataset$ID)),
            div(paste("TMT Mosaic", version))
         )
      })

      # Summary and Bioplex tabs share protein selection and history
      proteinSelect <- sharedSelectize(
         choices = dataset$proteins$info$MosaicID,
         selected = reactive({values$activeMosaicID}),
         msg = "The ID you selected was not found. Make sure that it is formatted correctly.",
         callback = updateProteinChoice)

      # Tabs
      source("tab/Explore.r",       local = TRUE)
      source("tab/Peptides.r",      local = TRUE)
      source("tab/Interactions.r",  local = TRUE)
      source("tab/PCA.r",           local = TRUE)
      source("tab/Volcano.r",       local = TRUE)
      source("tab/Correlation.r",   local = TRUE)
      source("tab/GO.r",            local = TRUE)
      source("tab/Edit.r",          local = TRUE)
      source("tab/Advanced.r",      local = TRUE)
      source("tab/Dataset.r",       local = TRUE)

      # Set up plots
      updateProteinChoice(values$activeMosaicID)
   }, once = TRUE)
}
