# Tab Constants
groupRows <- ave(rep(1, dataset$numSamples), dataset$groups, FUN=cumsum)

updateSelectInput(session, "pca_x", choices = 1:dataset$numSamples, selected = 1)
updateSelectInput(session, "pca_y", choices = 1:dataset$numSamples, selected = 2)

pcaData <- reactive({
   columns <- req(values$normProtCols)
   pcaData <- na.omit(columns)
   pcaData[apply(pcaData, 1, var) != 0, ]
})

pca <- reactive({
   pcaData <- req(pcaData())
   prcomp(t(pcaData), scale. = TRUE)
})

output$pca_header <- renderText({
   pcaData <- req(pcaData())
   paste("PCA: Using",
      nrow(pcaData),
      "of",
      nrow(values$normProtCols),
      paste0(tolower(dataset$idLabel), "s"))
})

output$plot_pca <- renderPlotly({
   pca <- req(pca())
   pcaX <- req(as.numeric(input$pca_x))
   pcaY <- req(as.numeric(input$pca_y))
   showNotification("Plotting PCA", duration = 4)

   x_range   <- as.numeric(c(input$pca_x_min, input$pca_x_max))
   y_range   <- as.numeric(c(input$pca_y_min, input$pca_y_max))
   size        <- as.numeric(input$pca_point_size)
   plot_width  <- as.numeric(input$pca_width)
   plot_height <- as.numeric(input$pca_height)
   labels <- c(input$pca_x_lab, input$pca_y_lab)

   scores <- as.data.frame(pca$x)
   variance <- pca$sdev / sum(pca$sdev)

   data <- data.frame(
      x = scores[, pcaX],
      y = scores[, pcaY],
      name = theme$names,
      group = theme$groupNames[dataset$groups])

   dflts <- function(x, dflt, test = is.na) ifelse(test(x), dflt, x)
   not <- function(f) function(...) !f(...)

   pad <- c(-10, 10)
   x_range <- dflts(x_range, range(data$x) + pad)
   y_range <- dflts(y_range, range(data$y) + pad)
   labels  <- dflts(labels, sprintf("Principal Component %d [%.2f%% Var]",
      c(pcaX, pcaY), variance[c(pcaX, pcaY)] * 100),
      test = not(nzchar))

   colorlist <- setNames(theme$groupColors, theme$groupNames)

   g <- ggplot(data) +
      geom_point(
         aes(
            x = x,
            y = y,
            text = name,
            fill = group),
         alpha = .7,
         shape = 21,
         size = size) +
      scale_fill_manual(values = colorlist) +
      labs(
         x = labels[[1]],
         y = labels[[2]]) +
      expand_limits(
         x = x_range,
         y = y_range) +
      guides(fill = "none") +
      theme_bw()

   ggplotly(g, tooltip = "text") %>%
      config(
         modeBarButtons = list(list("toImage")),
         toImageButtonOptions = list(
            filename = paste(dataset$name, "PCA", pcaX, pcaY,  sep = "_"),
            width  = plot_width,
            height = plot_height
         ),
         displaylogo = FALSE)
})

output$pca_loadings <- renderDataTable({
   pca <- req(pca())
   pcaX <- req(as.numeric(input$pca_x))
   pcaY <- req(as.numeric(input$pca_y))

   dt <- as.data.frame(pca$rotation[, c(pcaX, pcaY), drop = FALSE])
   ids <- rownames(dt)
   dt <- cbind(MosaicID = ids, dt, SET = ids)

   js_col_ind <- c(pcaX = 1, pcaY = 2)

   colDefs <- c(
      addShinyButton("SET", "switch_to_mosaicID"),
      addExpFormat(unname(js_col_ind))
   )

   datatable(dt,
      options = list(
         initComplete = datatableTheme(),
         searchHighlight = TRUE,
         columnDefs = colDefs,
         pageLength = 5,
         info = FALSE,
         order = list(js_col_ind[["pcaX"]], "desc")
      ),
      class = "compact cell-border",
      rownames = FALSE,
      selection = 'none',
      escape = FALSE)
})

output$pca_download <- downloadHandler(
   filename = function() {
      pca <- pca()
      if(is.null(pca)) return("PCA_error.txt")

      paste0(dataset$name, "_PCA.zip")
   },
   content = function(file) {
      pca <- pca()
      if(is.null(pca)) {
         write("No PCA data", file)
      }else {
         repeat {
            #Ensure unique path
            d <- paste0(dataPath, "tmp/", stringi::stri_rand_strings(1, 8))
            if(!dir.exists(d))
               break
         }
         dir.create(d, recursive = TRUE)
         lastDir <- setwd(d)
         zipFolder <- "PCA"
         dir.create(zipFolder)

         sep <- ","
         files <- c("scores.csv", "loadings.csv")
         scores <- cbind(Column = theme$names, pca$x)
         write.table(scores, "scores.csv", sep = sep, row.names = FALSE, col.names = TRUE)

         loadings <- pca$rotation
         loadings <- cbind(MosaicID = rownames(loadings), loadings)
         write.table(loadings, "loadings.csv", sep = sep, row.names = FALSE, col.names = TRUE)

         zip(file, files)
         setwd(lastDir)
         unlink(d, recursive = TRUE)
      }
   }
)

output$full_dendro <- renderPlot({
   dend_data <- req(values$normProtCols)
   colnames(dend_data) <- theme$names

   if(!is.null(advanced$bridge)){
      bridges <- advanced$bridge[advanced$bridge != "None"] %>%
         as.integer
      if(length(bridges) != 0) {
         #remove bridge columns
         validRows <- rowSums(is.na(dend_data[, bridges, drop = FALSE])) == 0
         dend_data <- dend_data[validRows, -bridges]
      }
   }

   ggdendrogram(hclust(dist(t(dend_data))), rotate = TRUE) +
      theme(
         panel.grid.major.y   = element_blank(),
         axis.text.x          = element_blank(),
         axis.text.y          = element_text(size = 12)
      )
}, height = function(){dataset$numSamples * 14})

if(dataset$numPlexes > 1){
   updateSelectizeInput(session, "histogrid_plex", choices = unique(dataset$classes), server = TRUE)
   shinyjs::show("histogrid_plex_column")
}

output$histogrid <- renderUI({
   if(input$overlay_histogrid)
      plotlyOutput("histogrid_summary")
   else
      plotOutput("histogrid_detailed")
})

plexCols <- reactive({
   plexCols <- rep(TRUE, ncol(req(values$normProtCols)))

   if(length(unique(dataset$classes)) > 1){
      plex <- req(input$histogrid_plex)
      plexCols <- dataset$classes == plex
   }

   return(plexCols)
})

histogridData <- reactive({
   plexData <- req(values$normProtCols)
   plexCols <- which(plexCols())

   plexData <- plexData[, plexCols]
   plexData <- log2(plexData / rowMeans(plexData, na.rm = TRUE))

   long_data <- NULL
   for(plexCol in 1:length(plexCols)){
      col <- plexCols[[plexCol]]
      long_data <- rbind(long_data,
         data.frame(
            log2fc = plexData[, plexCol],
            name = theme$names[[col]],
            group = theme$groupNames[dataset$groups][[col]],
            row = groupRows[[col]],
            color = theme$colors[[col]]
         )
      )
   }

   long_data <- long_data %>% filter(!is.na(log2fc) & log2fc != log2(0))

   return(long_data)
})

histogridLimit <- reactive({
   data <- req(histogridData())

   cutoff <- 0.03 # Discard most extreme 2% & center data
   xlim <- sort(data$log2fc)[nrow(data) * c(cutoff, 1 - cutoff)] %>%
      abs %>% max

   return(xlim)
})

output$histogrid_detailed <- renderPlot({
   req(!input$overlay_histogrid)
   long_data <- req(histogridData())

   plexCols <- plexCols()
   xlim <- histogridLimit()

   annot_df <- long_data %>%
      group_by(name, row, group) %>%
      summarize(mean = mean(log2fc), var = var(log2fc), .groups = "keep") %>%
      mutate(lab = paste(
         c("Mean", "Var"),
         c(mean, var) %>% formatC(3),
         sep = ": ", collapse = '\n'))

   numGroups <- length(unique(dataset$groups[plexCols]))
   wrapPlots <- numGroups == 1 || numGroups == sum(plexCols)
   if(!wrapPlots){
      annot_df$lab <- with(annot_df, paste(name, lab, sep = "\n"))
   }

   g <- ggplot(long_data) +
      geom_histogram(
         aes(
            x = log2fc,
            fill = I(color)
         ),
         bins = 80,
         color = "black",
         na.rm = TRUE) +
      geom_label(data = annot_df,
         aes(
            x = -Inf,
            y = Inf,
            label = lab),
         vjust = 1,
         hjust = 0) +
      geom_vline(data = annot_df,
         aes(
            xintercept = 0
         )) +
      geom_vline(data = annot_df,
         aes(
            xintercept = mean
         ),
         size = 1,
         linetype = "dashed",
         color = "red")

   if(wrapPlots){
      g <- g + facet_wrap(vars(factor(name, theme$names)), ncol = 4)
   }else{
      g <- g  +
         facet_grid(
            cols = vars(factor(group, theme$groupNames)),
            rows = vars(row))
   }
   # Decided against using plotly here because layout parameters are only applied to the first axes in faceted plots

   g +
      xlim(c(-xlim, xlim)) +
      theme_bw() +
      theme(strip.text.y = element_blank())# Remove row labels after applying the theme
})

output$histogrid_summary <- renderPlotly({
   req(input$overlay_histogrid)
   long_data <- req(histogridData())
   xlim <- histogridLimit()

   g <- ggplot(long_data,
         aes(
            x = log2fc,
            group = name
         )) +
      geom_density(# Outline for visibility of lighter colors
         color = "grey",
         size = 1.2,
         na.rm = TRUE,
      ) +
      geom_density(
         aes(,
            color = I(color),
            text = name
         ),
         na.rm = TRUE
      ) +
      geom_vline(data = data.frame(x = 0),
         aes(
            xintercept = x
         )) +
      xlim(c(-xlim, xlim)) +
      theme_bw()

   ggplotly(g, tooltip = "text") %>%
      plotlyDefaults(showlegend = FALSE)
})
