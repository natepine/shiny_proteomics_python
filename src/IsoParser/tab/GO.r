observe({
   req(GO)
   sortOrder <- req(input$go_heatmap_sort)
   counts <- GO$goCounts

   # Sort choices
   if(sortOrder == "Name"){
      counts <- counts[order(counts$Annotation), ]
   }else{ # Up or Down
      decr <- sortOrder == "Down"
      counts <- counts[order(counts$AnnotCount, decreasing = decr), ]
   }

   goParams <- setNames(
      object = counts$DatabaseID,
      nm = sprintf("%s|%s > %s [n=%d]", counts$Database, counts$DatabaseID, counts$Annotation, counts$AnnotCount))

   # Display nearest first
   if(!is.null(summary$sigGO)){
      isNear <- goParams %in% summary$sigGO$DatabaseID
      choices <- list(
         `Nearest GO` = goParams[isNear],
         All          = goParams[!isNear]
      )
   }else{
      choices <- list(
         All = goParams
      )
   }

   updateSelectizeInput(session, "goPlotter_annotation", choices = choices, server = TRUE)
})

goHeatmapIDs <- reactive({
   # Update only with 'Generate Heatmap' button
   req(input$button_plot_go)

   isolate({
      selected <- req(input$goPlotter_annotation)
      uniprotIDs <- filter(dataset$allGO, DatabaseID == selected)$UniprotID
      filter(dataset$proteins$info, UniprotID %in% uniprotIDs)$MosaicID
   })
})

goLogFC <- reactive({
   mosaicIDs <- req(goHeatmapIDs())
   calculateLogFC(dataset, values$normProtCols, mosaicIDs)
})

goPlotterHeight <- function(){
   height <- tryCatch({
      goLogFC <- req(goLogFC())
      numProteins <- nrow(goLogFC)
      size <- req(input$go_heatmap_height)
      size <- as.numeric(size)

      header_px <- 61
      axis_label_px <- 10 * max(nchar(theme$names))
      static_px <- header_px + axis_label_px

      width <- req(session$clientData$output_plot_go_heatmap_width)
      row_px <- (width - 340) / dataset$numSamples

      aspect_ratio <- 2 ** input$go_heatmap_height

      max_height <- 50000000 / width
      best_guess <- static_px + row_px * numProteins * aspect_ratio

      output <- best_guess
      if(best_guess > max_height) {
         showNotification("Plot was too long, clamping height")
         output <- max_height
      }

      output
   },
   error = function(e) {
      400
   })

   return(height)
}

output$plot_go_heatmap <- renderPlot({
   goLogFC <- req(goLogFC())

   colors <- c(req(input$go_heatmap_lower_col), req(input$go_heatmap_upper_col))
   fcCutoff <- as.numeric(req(input$go_heatmap_fc))

   highlightedGS <- NULL
   if(!is.null(values$nearestIDs)){
      show_nearest <- input$go_heatmap_show_nearest
      if(show_nearest) {
         highlightedGS <- dataset$proteins$info[values$nearestIDs, "GeneSymbol"]
      }
   }

   dendroHeatmap(goLogFC, colors, fcCutoff, isDotPlot = FALSE, colNames = theme$names, highlightedGS) +
      theme(plot.margin = unit(c(0,0,0,0), "cm"))
}, height = goPlotterHeight)

output$go_plotter_download <- downloadHandler(
   filename = function() {
      file <- tryCatch({
         goLogFC()
         databaseID <- req(input$goPlotter_annotation)
         paste0(dataset$name, "_", databaseID, ".csv")
      }, error = function(e) {
         "downloadError.csv"
      })
      return(file)
   },
   content = function(file) {
      data <- tryCatch({
         goLogFC <- goLogFC()
         colnames(goLogFC) <- dataset$names
         ids <- data.frame("MosaicID" = goHeatmapIDs(), "Gene Symbol" = rownames(goLogFC))
         cbind(ids, goLogFC)
      }, error = function(e) {
         data.frame(error = "No data plotted")
      })

      write.table(data, file, sep = ",", row.names = FALSE, col.names = TRUE)
   }
)
