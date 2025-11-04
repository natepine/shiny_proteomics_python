# Set up UI
if(!is.na(dataset$server) && dataset$server != "Upload")
   shinyjs::show(selector = ".shinyjs-hide:has(#server_info_box)")
if(!is.null(GO))
   shinyjs::show(selector = ".shinyjs-hide:has(#annot_info_box)")

output$dataset_ui <- renderUI({
   info <- with(dataset, c(
      "Dataset"         = name,
      "Species"         = species, #Hidden when NA
      "Plexing Level"   = plexLevel,
      "Num Plexes"      = numPlexes,
      "Default ID"      = initialID,
      "Created"         = date,
      "Notes"           = notes))

   row_count_name <- ifelse(dataset$isSiteQuant, "Site Count", "Protein Count")
   info[row_count_name] <-  nrow(dataset$proteins$info)

   peps <- grep("Peptides", colnames(dataset$proteins$info))
   if(any(peps)) {
      peps_pulled <- colSums(dataset$proteins$info[peps])
      info["Peptide Count"] <- paste(peps_pulled, collapse = ", ")
   }

   info <- info[intersect(
      c("Dataset",
        "Species",
        "Plexing Level",
        "Plex Count",
        "Protein Count", "Site Count",
        "Peptide Count",
        "Default ID",
        "Created",
        "Notes"), names(info))]

   HTML(keyValueStr(info, boldKeys = TRUE))
})

output$server_info_ui <- renderUI ({
   info <- with(dataset, c(
      "Server"    = server,
      "Username"  = username))
   id_name <- ifelse(dataset$isSiteQuant, "SQID", "PQID")
   info <- append(setNames(dataset$ID, id_name), info)

   pepInfo <- dataset$peptides$info
   if(!is.null(pepInfo) && "RunLoadPath" %in% colnames(pepInfo)) {
      files <- paste0(unique(pepInfo$RunLoadPath), collapse = "<br>")
      info[["Raw Files"]] <- sprintf("<br><div style='white-space:nowrap;overflow-x:scroll'>%s</div>", files)
   }

   HTML(keyValueStr(info, boldKeys = TRUE))
})

output$annot_info_ui <- renderUI ({
   info <- with(dataset, c(
      "Taxon ID" = taxid,
      "Last Updated" = attr(allGO, "LastUpdated"),
      "Annotation Types" = paste0("<div class='centered'>", sort(unique(allGO$Database)), "</div>", collapse = "")
   ))

   HTML(keyValueStr(info, boldKeys = TRUE))
})

dtable_filename <- function(tag = "TMTMosaic") {
   function(){
      filename <- input$dtable_filename
      format <- req(input$dtable_format)
      if(filename == ""){
         type <- ifelse(dataset$isSiteQuant, "sq", "pq")
         filename <- paste(sep = "_", format, type, dataset$ID, tag)
      }
      ext <- ifelse(format == "Summary", ".xlsx", ".csv")
      return(paste0(filename, ext))
   }
}

output$dtable_download <- downloadHandler(
   filename = dtable_filename(),
   content = function(file) {
      format <- req(input$dtable_format)

      waiting <- showNotification("Preparing Download...", duration = NULL)
      on.exit(removeNotification(waiting), add = TRUE)

      if(format == "Summary") {
         writeXLSX(dataset, idMaps, file, session$clientData$url_hostname)
      }else {
         data <- data.frame(error = "Unknown error")
         if(format == "Raw") {
            data$error <- "File not found"
            rdsFile <- paste0(ifelse(dataset$isSiteQuant, "sq", "pq"), "_", dataset$ID, ".rds")
            rdsPath <- paste0(dataPath, dataset$server, "/proteins/", rdsFile)
            if(file.exists(rdsPath)) {
               data <- readRDS(rdsPath)
            }
         }else if(format == "Live") {
            columns <- values$normProtCols
            if(!is.null(columns)) {
               #Add Normalization row
               columns <- rbind(advanced$factors, columns)

               info <- dataset$proteins$info
               info <- rbind(c("**NORMALIZATION_FACTORS**", rep(NA, ncol(info)-1)), info)

               colnames(columns) <- theme$names
               data <- cbind(info, columns)
            }
         }
         write.table(data, file, sep = ",", row.names = FALSE, col.names = TRUE)
      }
   }
)

if(!is.null(dataset$peptides$info)){
   observeEvent(input$dtable_format, {
      if(input$dtable_format ==  "Summary") {
         shinyjs::disable("dtable_pep_download")
      }else {
         shinyjs::enable("dtable_pep_download")
      }
   })

   output$dtable_pep_download <- downloadHandler(
      filename = dtable_filename("peptides_TMTMosaic"),
      content = function(file) {
         format <- req(input$dtable_format)
         sep <- ","
         waiting <- showNotification("Preparing Download...", duration = NULL)
         if(format == "Raw") {
            rdsFile <- paste0(ifelse(dataset$isSiteQuant, "sq", "pq"), "_", dataset$ID, ".rds")
            rdsPath <- paste0(dataPath, dataset$server, "/peptides/", rdsFile)
            data <- data.frame(error = "Download Failed: File not found")
            if(file.exists(rdsPath)) {
                data <- readRDS(rdsPath)
            }
         }else if(format == "Live") {
            info <- dataset$peptides$info
            info <- info[order(info$Class), ]

            columns <- values$normPepCols

            colIndex <- 0
            mat <- matrix(ncol = dataset$numSamples, nrow = nrow(info))
            colnames(mat) <- character(ncol(mat))

            factors <- numeric(dataset$numSamples)
            for(class in unique(info$Class)) {
               rows <- class == info$Class
               if(!any(rows)) next

               classData <- columns[[class]]
               classCols <- seq(ncol(classData)) + colIndex
               mat[rows, classCols] <- classData[info[rows, "data_row"], ]

               colIDs <- dataset$peptides$columnIDs[[class]]
               factors[classCols] <- advanced$factors[colIDs]

               colNames <- colnames(classData)
               if(!dataset$classesMap) {
                  colNames <- paste(class, colNames, sep = "~")
               }else {
                  colNames <- theme$names[colIDs]
               }
               colnames(mat)[classCols] <- colNames

               colIndex <- max(classCols)
            }
            info <- info[colnames(info) != "data_row"]
            info <- rbind(c("**NORMALIZATION_FACTORS**", rep(NA, ncol(info)-1)), info)
            mat <- rbind(factors, mat)

            data <- cbind(info, mat)
         }
         write.table(data, file, sep = sep, row.names = FALSE, col.names = TRUE)
         removeNotification(waiting)
      }
   )
}
