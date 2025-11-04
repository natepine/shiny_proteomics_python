#Reactives
volcano <- reactiveValues()

#UI
observe({
   choices <- setNames(seq(dataset$numGroups), theme$groupNames)
   selected <- choices[table(dataset$groups) > 1]
   if(length(selected) < 2) {
      selected <- c(selected, setdiff(choices, selected))[1:2] %>% sort
   }
   updateSelectInput(session, "volcano_grp_a", choices = choices, selected = selected[[1]])
   updateSelectInput(session, "volcano_grp_b", choices = choices, selected = selected[[2]])
})

#Calculations
buildVolcano <- observe({
   columns <- req(values$normProtCols)

   input$button_update_volcano
   req(input$volcano_padjust)
   if(dataset$numGroups < 2) {
      showNotification("Sorry, volcano plot cannot be used with less than 2 groups.",
         duration = NULL, type = "warning")
      return()
   }
   a <- req(isolate(input$volcano_grp_a))
   b <- req(isolate(input$volcano_grp_b))

   if (a == b) {
      sendSweetAlert(
         title = "Same Groups...",
         text = "You selected the same group for both Group A and Group B. They cannot be the same.",
         type = "warning"
      )
      return(NULL)
   }

   sub_index <- dataset$groups %in% c(a,b)
   groups <- dataset$groups[sub_index]
   a_index <- groups == a
   b_index <- groups == b
   if(sum(a_index) < 2 || sum(b_index) < 2) {
      showNotification("At least two columns in each group are required",
         duration = NULL, type = "warning")
      return()
   }
   columns <- columns[, sub_index]

   # Remove rows with 0 variance or less than 2 values:
   row_valid <- apply(columns, 1, var, na.rm = TRUE)
   row_valid[is.na(row_valid)] <- 0
   columns <- columns[row_valid != 0, , drop = FALSE]

   a_df <- columns[, a_index] %>% t() %>% as.data.frame() # Transposed for mapply
   b_df <- columns[, b_index] %>% t() %>% as.data.frame()
   showNotification("(1/3) Loaded groups...", duration = 4)

   # Get fold changes:
   fc_df <- dataset$proteins$info[row_valid != 0, c("MosaicID", "GeneSymbol"), drop = FALSE] %>%
      mutate(
         A_MEAN = colMeans(a_df, na.rm = TRUE),
         B_MEAN = colMeans(b_df, na.rm = TRUE))

   fc_df$A_MEAN[fc_df$A_MEAN == 0] <- .001
   fc_df$B_MEAN[fc_df$B_MEAN == 0] <- .001
   showNotification("(2/3) Made Table...", duration = 4)

   # Only compute P-values with sufficient data
   hasData <- colSums(!is.na(a_df)) > 1 &
               colSums(!is.na(b_df)) > 1
   shiny::validate(need(any(hasData), "Too many missing values"))

   # Convert t.test errors to NaN
   t.test_p.val <- function(A, B) {
      tryCatch(
         t.test(A, B, var.equal = TRUE)$p.value,
      error = function(cond){NaN})
   }

   # !hasData is automatically filled with NA
   fc_df$negLog10Pval[hasData] <- mapply(a_df[hasData], b_df[hasData], FUN = t.test_p.val) %>%
      p.adjust(method = input$volcano_padjust) %>%
      log10() * -1

   fc_df <- fc_df %>%
      mutate(log2FC = log2(B_MEAN / A_MEAN)) %>%
      arrange(MosaicID, GeneSymbol, log2FC, negLog10Pval, A_MEAN, B_MEAN)

   showNotification("(3/3) Calculated PValues...", duration = 4)

   # Filter out NaN and fc of 0:
   volcano$data <- fc_df %>% filter(!(is.nan(negLog10Pval) | log2FC == log2(0)))
}, suspended = TRUE)

startVolcano <- observe({
   # Listen for volcano init conditions
   req(input$volcano_grp_a)
   req(input$volcano_grp_b)
   req(values$normProtCols)
   req(input$sidebarMenu == "volcano")

   buildVolcano$resume()
   startVolcano$destroy()
})

volcano_binned <- reactive({
   fc_df <- req(volcano$data)

   x_cutoff <- req(input$volcano_x_cutoff)
   y_cutoff <- req(input$volcano_y_cutoff)

   log2fc_threshold        <-  log2(as.numeric(x_cutoff))
   negLog10Pval_threshold  <- -log10(as.numeric(y_cutoff))

   hasPval <- !is.na(fc_df$negLog10Pval)
   fc_df <- fc_df %>%
    mutate(DownReg = hasPval & log2FC <= -log2fc_threshold & negLog10Pval >= negLog10Pval_threshold,
           UpReg   = hasPval & log2FC >=  log2fc_threshold & negLog10Pval >= negLog10Pval_threshold)

   fc_df$COLOR                <- "NotSig"
   fc_df$COLOR[fc_df$UpReg]   <- "UpReg"
   fc_df$COLOR[fc_df$DownReg] <- "DownReg"
   
   fc_df$pointID <- 1:nrow(fc_df)

   return(fc_df)
})

#Annotations
output$volcano_annot_ui <- renderUI({
   auto_annotate <- isolate(input$volcano_annotation)

   if(is.null(auto_annotate)){
      fc_df <- req(volcano_binned())
      numSig <- sum(fc_df[c("UpReg", "DownReg")])

      if(numSig == 0){
         auto_annotate <- "none"
      }else if(numSig < 20){
         auto_annotate <- "all"
      }else if(numSig < 200){
         auto_annotate <- "randAll"
      }else {
         auto_annotate <- "distAll"
      }
   }

   # Add GO annotation option only if GO data is available
   go_choices <- if(!is.null(GO)) c("goAnnotation") else c()
   go_names <- if(!is.null(GO)) c("GO Annotation") else c()
   
   choices <- c("none", "currentProtein", "nearestProteins", "randUp", "randDown", "randAll", "distUp", "distDown", "distAll", "up", "down", "all", go_choices, "custom")
   geneSymbol <- dataset$proteins$info[values$activeMosaicID, "GeneSymbol"]
   idLabels <- paste0(dataset$idLabel, "s")
   names(choices) <- c(
      "None",
      paste0("Current ", dataset$idLabel, " (", geneSymbol, ")"),
      paste(length(values$nearestIDs), "Closest", idLabels, "to", geneSymbol),
      "Random 10% Upreg.",
      "Random 10% Downreg.",
      "Random 10% Sig.",
      paste("Top Upreg.",   idLabels),
      paste("Top Downreg.", idLabels),
      paste("Top Sig.",     idLabels),
      paste("All Upreg.",   idLabels),
      paste("All Downreg.", idLabels),
      paste("All Sig.",     idLabels),
      go_names,
      paste("Custom",       idLabels))

   selectInput("volcano_annotation", "Annotate Points", selected = auto_annotate, choices = choices)
})

output$volcano_custom_annot_ui <- renderUI({
   req(input$volcano_annotation == "custom")

   updateSelectizeInput(session, "volcano_annot_proteins", choices = dataset$protein$info$MosaicID, selected = NULL, server = TRUE)
   fluidRow(
      column(10, selectInput("volcano_annot_proteins", NULL, choices = NULL, multiple = TRUE)),
      column(2, align = "center", actionButton("volcano_annotate_custom", "Annotate"))
   )
})

output$volcano_go_annot_ui <- renderUI({
   req(input$volcano_annotation == "goAnnotation")
   req(!is.null(GO))
   
   # Get available GO terms from current enrichment results if available
   go_choices <- NULL
   go_names <- NULL
   
   if(!is.null(input$volcano_go) && !is.null(volcano$goSigGO)) {
      # Use current enrichment results
      sigGO <- volcano$goSigGO
      go_choices <- sigGO$DatabaseID
      go_names <- paste0(sigGO$Database, ": ", sigGO$Annotation, " (", sigGO$SubsetHits, " hits)")
      names(go_choices) <- go_names
   } else {
      # Fallback to all available GO terms from the GO calculator
      if(!is.null(volcanoGO)) {
         # Access GO data through the calculator's internal allGO
         # This requires adding a method to expose available GO terms
         tryCatch({
            allTerms <- GO$goCounts %>%
               arrange(desc(AnnotCount)) %>%
               head(100)  # Limit to top 100 to avoid overwhelming UI
            go_choices <- allTerms$DatabaseID
            go_names <- paste0(allTerms$Database, ": ", allTerms$Annotation, " (", allTerms$AnnotCount, " proteins)")
            names(go_choices) <- go_names
         }, error = function(e) {
            go_choices <- c()
            go_names <- c()
         })
      }
   }
   
   if(is.null(go_choices) || length(go_choices) == 0) {
      fluidRow(
         column(12, tags$p("No GO terms available. Please run GO enrichment analysis first.", style = "color: orange;"))
      )
   } else {
      fluidRow(
         column(10, selectInput("volcano_go_term", "Select GO Term:", choices = go_choices, selected = NULL)),
         column(2, align = "center", actionButton("volcano_annotate_go", "Annotate"))
      )
   }
})

observe({
   labels <- isolate(volcano$labels)
   fc_df <- volcano_binned()
   req(fc_df)
   req(input$volcano_annotation != "custom")
   req(input$volcano_annotation != "goAnnotation")

   showNotification("Updating labels...", duration = 1)

   if(is.null(labels) || length(labels) != nrow(fc_df)){
      labels <- setNames(rep(FALSE, nrow(fc_df)), fc_df$MosaicID)
   }

   #Returns the indices of a random 10% of TRUE values
   Ten_Percent <- function(bool){
      num_in <- sum(bool)
      if(num_in == 0) return(NULL)

      total_ind <- which(bool)
      if(num_in == 1) return(total_ind)

      return(sample(total_ind, ceiling(num_in * .10)))
   }

   #Returns the 'numHits' most distant TRUE values
   Edge <- function(bool, numHits = 15){
      if(sum(bool) > numHits){
         set <- fc_df[bool, c("log2FC", "negLog10Pval")] %>% t %>% abs

         x_cutoff <- input$volcano_x_cutoff
         y_cutoff <- input$volcano_y_cutoff
         origin <- c(log2(as.numeric(x_cutoff)), -log10(as.numeric(y_cutoff)))
         set <- set - origin # Center data on cutoff
         
         dist_squared <- colSums(set * set)
         near <- order(dist_squared, decreasing = TRUE)[-(1:numHits)]
         bool[bool][near] <- FALSE
      }
      return(bool)
   }

   labels[labels] <- FALSE
   active <- switch(
      input$volcano_annotation,

      "currentProtein" = fc_df$MosaicID == values$activeMosaicID,
      "nearestProteins" = match(values$nearestIDs, fc_df$MosaicID),

      "up" = fc_df$UpReg,
      "down" = fc_df$DownReg,
      "all" = fc_df$UpReg | fc_df$DownReg,

      "randUp" = Ten_Percent(fc_df$UpReg),
      "randDown" = Ten_Percent(fc_df$DownReg),
      "randAll" = c(Ten_Percent(fc_df$UpReg), Ten_Percent(fc_df$DownReg)),

      "distUp" = Edge(fc_df$UpReg),
      "distDown" = Edge(fc_df$DownReg),
      "distAll" = Edge(fc_df$UpReg) | Edge(fc_df$DownReg),

      NULL
   )

   if(!is.null(active) && any(active)){
      labels[active] <- TRUE
   }

   volcano$labels <- labels
})

observeEvent(input$volcano_annotate_custom, {
   labeled <- req(input$volcano_annot_proteins)

   volcano$labels[volcano$labels] <- FALSE
   volcano$labels[labeled] <- TRUE
})

observeEvent(input$volcano_annotate_go, {
   req(!is.null(GO))
   req(exists("volcanoGO") && !is.null(volcanoGO))
   selected_go <- req(input$volcano_go_term)
   
   tryCatch({
      # Get proteins associated with the selected GO term
      go_uniprots <- volcanoGO$getProteinsForGOTerm(selected_go)
      
      if(length(go_uniprots) > 0) {
         # Map UniprotIDs to MosaicIDs
         go_proteins <- dataset$proteins$info %>%
            filter(UniprotID %in% go_uniprots) %>%
            pull(MosaicID)
         
         if(length(go_proteins) > 0) {
            volcano$labels[volcano$labels] <- FALSE
            volcano$labels[go_proteins] <- TRUE
            
            # Set the annotation color from GO color mapping if available
            if(!is.null(volcano$goColorMapping) && selected_go %in% names(volcano$goColorMapping)) {
               # Update the annotation color input to match the selected GO term color
               updateColourInput(session, "volcano_annotation_color", 
                               value = volcano$goColorMapping[[selected_go]])
            }
            
            showNotification(paste("Annotated", length(go_proteins), "proteins for selected GO term"), duration = 3)
         } else {
            showNotification("No proteins found for selected GO term in current dataset", duration = 3, type = "warning")
         }
      } else {
         showNotification("No proteins found for selected GO term", duration = 3, type = "warning")
      }
   }, error = function(e) {
      showNotification(paste("Error accessing GO data:", e$message), duration = 5, type = "error")
   })
})

observeEvent(input$volcano_go_annotate_term, {
   req(!is.null(GO))
   req(exists("volcanoGO") && !is.null(volcanoGO))
   selected_go <- req(input$volcano_go_annotate_term)
   
   # Validate selected_go
   if(is.na(selected_go) || is.null(selected_go) || selected_go == "") {
      showNotification("Invalid GO term selected", duration = 3, type = "warning")
      return()
   }
   
   tryCatch({
      # Get proteins associated with the selected GO term
      go_uniprots <- volcanoGO$getProteinsForGOTerm(selected_go)
      
      if(length(go_uniprots) > 0) {
         # Map UniprotIDs to MosaicIDs
         go_proteins <- dataset$proteins$info %>%
            filter(UniprotID %in% go_uniprots) %>%
            pull(MosaicID)
         
         if(length(go_proteins) > 0) {
            # Store GO proteins for color-only annotation (no text labels)
            volcano$goProteins <- go_proteins
            # Use the single GO color selector
            volcano$goColor <- if(is.null(input$volcano_go_color)) "#FF6B6B" else input$volcano_go_color
            
            # Clear regular labels to avoid text annotation
            volcano$labels[volcano$labels] <- FALSE
            
            showNotification(paste("Colored", length(go_proteins), "proteins for GO term (color only)"), duration = 3)
         } else {
            showNotification("No proteins found for selected GO term in current dataset", duration = 3, type = "warning")
         }
      } else {
         showNotification("No proteins found for selected GO term", duration = 3, type = "warning")
      }
   }, error = function(e) {
      showNotification(paste("Error accessing GO data:", e$message), duration = 5, type = "error")
   })
})

observeEvent(input$volcano_click, {
   if(!is.null(volcano$labels)){
      index <- as.numeric(input$volcano_click)
      volcano$labels[index] <- !volcano$labels[index]
   }
   session$sendCustomMessage(type = "resetValue", message = "volcano_click")
})

observeEvent(input$volcano_doubleclick, {
   row <- as.numeric(input$volcano_doubleclick)

   mosaicID <- volcano_binned()[[row, "MosaicID"]]
   proteinSelect$addToHistory(mosaicID)
   updateTabItems(session, "sidebarMenu", "summary")

   session$sendCustomMessage(type = "resetValue", message = "volcano_doubleclick")
})

#Plot
output$plot_volcano <- renderPlotly({
   shiny::validate(need(dataset$numGroups > 1, "Volcano plot can't be used with only one group."))
   req(!is.null(volcano$labels))

   down <- req(input$volcano_up_color)
   up <- req(input$volcano_down_color)

   plot_format <- input$volcano_plot_format

   colorlist <- c(
      "NotSig" = "#D3DDDC",
      "DownReg" = req(input$volcano_down_color),
      "UpReg" = req(input$volcano_up_color)
   )

   # Override background colors if GO background color is set
   if(!is.null(input$volcano_go_bg_color)) {
      colorlist["NotSig"] <- input$volcano_go_bg_color
   }

   fc_df <- volcano_binned()

   # Get plot limits:
   xmax <- max(abs(fc_df$log2FC), na.rm = TRUE)
   ymax <- max(fc_df$negLog10Pval, na.rm = TRUE)

   # Make rectangles for cutoffs:
   x_cutoff <- log2(as.numeric(input$volcano_x_cutoff))
   y_cutoff <- -log10(as.numeric(input$volcano_y_cutoff))
   df_rect <- lapply(c(1, 2), function(rect){
      list(
         type = "rect",
         x0 = c(-100, x_cutoff)[rect],
         x1 = c(-x_cutoff, 100)[rect],
         y0 = y_cutoff,
         y1 = 100,
         fillcolor = "#C6B19D",
         line = list(width = 0), #Hide border
         opacity = 0.3
      )
   })

   annotation_color <- input$volcano_annotation_color
   margin_factor <- 1.05 #5% padding

   g <- plot_ly(fc_df, source = "volcano") %>%
      add_trace(
         type = "scattergl",
         mode = "markers",
         x = ~ log2FC,
         y = ~ negLog10Pval,
         hoverinfo = "text",
         text = ~ GeneSymbol,
         colors = colorlist,
         color = ~ COLOR,
         marker = list(
            size = 5,
            line = list(
               color = "black",
               width = 1
            )
         ),
         customdata = ~ pointID
      )

   # Add GO protein coloring (color only, no text labels)
   if(!is.null(volcano$goProteins) && length(volcano$goProteins) > 0) {
      go_data <- fc_df %>% filter(MosaicID %in% volcano$goProteins)
      if(nrow(go_data) > 0) {
         g <- g %>%
            add_trace(
               data = go_data,
               type = "scattergl",
               mode = "markers",
               x = ~ log2FC,
               y = ~ negLog10Pval,
               hoverinfo = "text",
               text = ~ GeneSymbol,
               showlegend = FALSE,
               marker = list(
                  size = 6,
                  color = if(is.null(volcano$goColor)) "#FF6B6B" else volcano$goColor,
                  line = list(
                     color = "black",
                     width = 1
                  )
               )
            )
         
         # Add text labels if requested
         if(!is.null(input$volcano_go_show_labels) && input$volcano_go_show_labels) {
            g <- g %>%
               add_annotations(
                  data = go_data,
                  font = list(color = if(is.null(volcano$goColor)) "#FF6B6B" else volcano$goColor),
                  x = ~ log2FC,
                  y = ~ negLog10Pval,
                  xanchor = ~ ifelse(UpReg, "left", ifelse(DownReg, "right", "center")),
                  yanchor = "center",
                  text = ~ GeneSymbol,
                  bgcolor = "white",
                  bordercolor = "grey",
                  showarrow = TRUE,
                  arrowhead = 0,
                  ax = ~ ifelse(UpReg, 7, ifelse(DownReg, -7, 0)),
                  ay = ~ ifelse(UpReg | DownReg, -5, -15),
                  standoff = 3
               )
         }
      }
   }

   # Add regular annotations with text labels
   if(any(volcano$labels)){
      annotation_data <- fc_df %>% filter(volcano$labels)
      g <- g %>%
         add_trace(
            data = annotation_data,
            type = "scattergl",
            mode = "markers",
            x = ~ log2FC,
            y= ~ negLog10Pval,
            hoverinfo = "skip",
            showlegend = FALSE,
            marker = list(size = 5, color = annotation_color)
         ) %>%
         add_annotations(
            data = annotation_data,
            font = list(color = annotation_color),
            x = ~ log2FC,
            y = ~ negLog10Pval,
            xanchor = ~ ifelse(UpReg, "left", ifelse(DownReg, "right", "center")),
            yanchor = "center",
            text = ~ GeneSymbol,
            bgcolor = "white",
            bordercolor = "grey",
            # opacity = 0.8,
            showarrow = TRUE,
            arrowhead = 0,
            ax = ~ ifelse(UpReg, 7, ifelse(DownReg, -7, 0)),
            ay = ~ ifelse(UpReg | DownReg, -5, -15),
            standoff = 3
         )
   }
   
   title <- "P-value"
   if(input$volcano_padjust != "none") {
      title <- paste(title, "Adj.")
   }
   title <- paste0("-log10(",title, ")")

   g <- g %>%
      layout(
         shapes = df_rect,
         xaxis = list(
            title = "log2(FC)",
#            fixedrange = TRUE,
            zeroline = FALSE,
            range = xmax * c(-1, 1) * margin_factor
         ),
         yaxis = list(
            title = title,
            fixedrange = TRUE,
            zeroline = FALSE,
            range = ymax * c((1 - margin_factor), margin_factor)
         ),
         legend = list(title = FALSE, yanchor = "center", y = 0.5),
         font = list(family = "serif")
      ) %>%
      config(
         displaylogo = FALSE,
         modeBarButtonsToRemove = c(
            "select2d",
            "lasso2d",
            "hoverClosestCartesian",
            "hoverCompareCartesian"),
         toImageButtonOptions = list(
            format = plot_format, # "png",
            filename = paste0(dataset$ID, "_volcano"),
            width = NULL,
            height = NULL
         ),
         edits = list(
            annotationTail = TRUE
         )) %>%
      onRender('function(el) {
         el.on(\'plotly_click\', singleClickHandler);
         el.on(\'plotly_doubleclick\', doubleClickHandler);
      }')

   return(g)
})

output$plot_volcano_groups <- renderPlotly({
   #Visualize the groups used in the volcano plot
   g1 <- as.numeric(req(input$volcano_grp_a))
   g2 <- as.numeric(req(input$volcano_grp_b))

   colors <- theme$colors
   bordervec <- rep("black", length(colors))
   names <- theme$names

   # Update colors to only show relevant groups:
   fade <- !dataset$groups %in% c(g1, g2)
   colors[fade] <- adjustcolor(colors[fade], alpha.f = .1)
   bordervec[fade] <- adjustcolor(bordervec[fade], alpha.f = .2)

   # Save to object
   g <- ggplot() +
      geom_bar(aes(x = names, y = 1), fill = colors, size = .2, color = bordervec, stat = "identity") +
      scale_x_discrete(limits = names) +
      scale_y_continuous(breaks = NULL, labels = NULL) +
      ggtitle(paste("Comparing", theme$groupNames[[g1]], "to", theme$groupNames[[g2]])) +
      theme(
         axis.title   = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y  = element_blank()
      )

   ggplotly(g, tooltip = NULL) %>%
      config(displayModeBar = FALSE) %>%
      layout(xaxis = list(fixedrange = TRUE), yaxis = list(fixedrange = TRUE))
})

### GO Analysis
if(!is.null(GO)){
   shinyjs::show("volcano_GO")
   volcanoGO <- GO$goCalculator()

   observeEvent(volcano_binned(), {
      # Clear selection
      updateRadioGroupButtons(session, "volcano_go", selected = character(0))
   })

   output$volcano_goTable <- renderDataTable({
      subset <- req(input$volcano_go)
      fc_df <- req(volcano_binned())

      IDs <- switch(subset,
         "down" = fc_df[fc_df$DownReg, "MosaicID"],
         "all" = fc_df[fc_df$DownReg | fc_df$UpReg, "MosaicID"],
         "up" = fc_df[fc_df$UpReg, "MosaicID"]
      )
      shiny::validate(need(length(IDs) != 0, "No significant proteins found"))

      volcanoGO$subsetData(IDs)

      minMatches <- as.numeric(input$volcanoMinMatches)
      pvalCutoff <- as.numeric(input$volcanoPvalCutoff)
      p.adj <- input$volcanoPvalAdjust

      sigGO <- volcanoGO$filterGO(minMatches, pvalCutoff, p.adj)
      shiny::validate(need(sigGO, "No significant GO categories found, try widening parameters."))

      # Store significant GO results for use in annotation
      volcano$goSigGO <- sigGO
      
      # Pass tooltips in a reactive
      volcano$goProtTooltips <- volcanoGO$protTooltips(maxTooltips = 40)

      # Get the base GO table and process it safely
      tryCatch({
         # Try to get raw data instead of rendered table
         goTable <- NULL
         
         # First try to get the raw data if available
         if(exists("getSigGO", where = volcanoGO)) {
            goTable <- volcanoGO$getSigGO()
         } else {
            # Fallback to goTable method but handle the result carefully
            goTableResult <- volcanoGO$goTable()
            
            # If it's already a DataTable widget, we can't easily extract the data
            # So we'll return the widget as-is and add our functionality via JavaScript
            if(inherits(goTableResult, c("datatables", "htmlwidget"))) {
               showNotification("GO table returned as widget, using basic display", type = "message")
               return(goTableResult)
            } else {
               goTable <- goTableResult
            }
         }
         
         # Validate and convert goTable
         if(is.null(goTable)) {
            return(datatable(data.frame(Message = "No GO data available"), 
                           options = list(pageLength = 5), rownames = FALSE))
         }
         
         # Try to convert to data frame if it's not already
         if(!is.data.frame(goTable)) {
            if(is.matrix(goTable) || is.list(goTable)) {
               goTable <- as.data.frame(goTable, stringsAsFactors = FALSE)
            } else {
               return(datatable(data.frame(Message = paste("GO data has unsupported format:", class(goTable)[1])), 
                              options = list(pageLength = 5), rownames = FALSE))
            }
         }
         
         if(nrow(goTable) == 0) {
            return(datatable(data.frame(Message = "No GO terms found"), 
                           options = list(pageLength = 5), rownames = FALSE))
         }
         
         # Add a color column for volcano plot annotation
         # Create default colors for each GO term (cycling through a palette)
         default_colors <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", 
                           "#DDA0DD", "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9")
         num_terms <- nrow(goTable)
         
         # Check if we have DatabaseID or similar identifier column
         id_col <- NULL
         possible_id_cols <- c("DatabaseID", "ID", "TermID", "GOTerm", "Term")
         for(col in possible_id_cols) {
            if(col %in% colnames(goTable)) {
               id_col <- col
               break
            }
         }
         
         if(is.null(id_col)) {
            # If no ID column found, create simple row numbers
            goTable$DatabaseID <- paste0("GO_", seq_len(nrow(goTable)))
            id_col <- "DatabaseID"
         }
         
         # Create color input HTML for each row - all grey buttons
         color_inputs <- sapply(1:num_terms, function(i) {
            term_id <- as.character(goTable[[id_col]][i])
            
            # Validate term_id
            if(is.na(term_id) || is.null(term_id) || term_id == "" || term_id == "NA") {
               term_id <- paste0("unknown_term_", i)
            }
            
            # Create a grey button that triggers annotation
            paste0('<div class="volcano-go-button" data-term-id="', term_id, '" ',
                   'style="width: 40px; height: 22px; background-color: #E5E5E5; color: white; ',
                   'border: 1px solid #666; cursor: pointer; border-radius: 3px; margin: 2px auto; ',
                   'display: flex; align-items: center; justify-content: center; font-size: 10px; font-weight: bold;" ',
                   'onclick="Shiny.setInputValue(\'volcano_go_annotate_term\', \'', term_id, '\', {priority: \'event\'});" ',
                   'title="Click to annotate plot">GO</div>')
         })
         
         # Add the button column to the table
         goTable$`Annotate` <- color_inputs
         
         # Reorder columns to put button first for visibility
         col_order <- c("Annotate", setdiff(names(goTable), "Annotate"))
         goTable <- goTable[, col_order]
         
         # Render the datatable with custom options
         datatable(goTable, 
                   options = list(
                      columnDefs = list(
                         list(targets = 0, orderable = FALSE, width = "50px") # Button column
                      ),
                      pageLength = 10,
                      autoWidth = FALSE,
                      dom = 'tip'
                   ),
                   escape = FALSE, # Allow HTML in button column
                   rownames = FALSE)
                   
      }, error = function(e) {
         showNotification(paste("Error processing GO table:", e$message), type = "error")
         # Return a simple error table as fallback
         datatable(data.frame(Error = paste("GO table processing failed:", e$message)), 
                   options = list(pageLength = 5), 
                   rownames = FALSE)
      })
   })
}

# Simple script for GO table interaction
output$volcano_go_color_script <- renderUI({
   tags$script(HTML('
      $(document).ready(function() {
         // Simple click handler for GO terms - no color picker needed
         console.log("GO table interaction ready - using single color selector");
      });
   '))
})

#Tables
buildVolcanoTable <- function(data){
   cols <- c(
      "UniprotID",
      "GeneSymbol",
      "Description",
      "Sequence",
      "NumPeps",
      "neglogPVal",
      "log2FC",
      "Site")

   info <- dataset$proteins$info[data$MosaicID, ]
   tbl <- cbind(
      data %>% rename(neglogPVal = negLog10Pval),
      info %>% select(-c("MosaicID", "GeneSymbol")) #Drop duplicate columns
   ) %>% mutate(NumPeps = getNumPepsFromID(info, MosaicID))

   tbl[intersect(cols, colnames(tbl))]
}

getDT <- function(data, color){
   data$neglogPVal <- round(as.numeric(data$neglogPVal), 3)
   data$log2FC <- round(as.numeric(data$log2FC), 3)
   
   if(!is.atomic(data$NumPeps)){
      data$NumPeps <- apply(data$NumPeps, 1, function(peps){
         paste0("{", paste0(peps, collapse = ","), "}")
      })
   }

   # Add the buttons to change to the proteins
   data$SET <- rownames(data) # Rownames are MosaicIDs

   if(dataset$isSiteQuant) {
      protSeqs <- getProteinSequence(dataset, data$UniprotID)
      data <- data %>% mutate(
         ProtSize = nchar(protSeqs),
      ) %>% mutate(
         seqTitle = paste0("Site Pos: ", Site, "/", ProtSize, "<br>", trunc_seqs(protSeqs, Sequence)),
      )
   }

   # Reorder data
   cols <- c(
      "UniprotID",
      "GeneSymbol",
      "Description",
      "Sequence",
      "ProtSize",
      "NumPeps",
      "neglogPVal",
      "log2FC",
      "GO",
      "SET",
#hidden
      "Site",
      "seqTitle"
   )
   data <- data[intersect(cols, colnames(data))]

   # Render data
   colDefs <- c(
      addLink("https://www.uniprot.org/uniprotkb/", "UniprotID", "/entry", c(color = color)),
      addShinyButton("SET", "switch_to_mosaicID"))
   if(dataset$isSiteQuant)
      colDefs <- c(colDefs, addPosMarker(data, "ProtSize", "Site", "seqTitle", enableHTML = TRUE))

   dt <- data %>%
      datatable(
         options = list(
            sDom = '<"top">tp<"bottom">',
            columnDefs = colDefs,
            initComplete = datatableTheme(color),
            searchHighlight = TRUE,
            lengthMenu = c(10, 25, 50, 100, nrow(data)),
            pageLength = 10,
            autoWidth = FALSE),
         class = "compact cell-border",
         rownames = FALSE,
         selection = 'none',
         escape = FALSE)
   return(dt)
}

observeEvent(volcano_binned(), {
   fc_df <- volcano_binned()

   volcano$tableDown <- NULL
   if (any(fc_df$DownReg)) {
      volcano$tableDown <- buildVolcanoTable(fc_df[fc_df$DownReg, ])
   }

   volcano$tableUp <- NULL
   if (any(fc_df$UpReg)) {
      volcano$tableUp <- buildVolcanoTable(fc_df[fc_df$UpReg, ])
   }
})

output$volcano_down_count <- renderUI({
   fc_df <- req(volcano_binned())
   color <- req(input$volcano_down_color)
   tags$h4(paste(sum(fc_df$DownReg), paste0("Down Regulated ", dataset$idLabel, "s")), style = paste0("color:", color))
})

output$volcano_up_count <- renderUI({
   fc_df <- req(volcano_binned())
   color <- req(input$volcano_up_color)
   tags$h4(paste(sum(fc_df$UpReg), paste0("Up Regulated ", dataset$idLabel, "s")), style = paste0("color:", color))
})

addEnrichment <- function(table, type){
   if(!is.null(input$volcano_go) && !is.null(volcano$goProtTooltips) &&
         input$volcano_go %in% c(type, "all")){
      table$GO <- volcano$goProtTooltips[table$UniprotID]
   }
   return(table)
}

output$table_volcano_down <- renderDataTable({
   shiny::validate(need(volcano$tableDown, paste0("No downregulated ", tolower(dataset$idLabel), "s found with current cutoffs.")))
   color <- req(input$volcano_down_color)

   table <- addEnrichment(volcano$tableDown, "down")

   getDT(table, color)
})

output$table_volcano_up <- renderDataTable({
   shiny::validate(need(volcano$tableUp, paste0("No upregulated ", tolower(dataset$idLabel), "s found with current cutoffs.")))
   color <- req(input$volcano_up_color)

   table <- addEnrichment(volcano$tableUp, "up")

   getDT(table, color)
})

output$volcano_go_download <- downloadHandler(
   filename = function() {
      return(paste(dataset$ID, input$volcano_go, input$volcano_grp_a, "vs", input$volcano_grp_b, "volcano_GO.csv", sep = "_"))
   },
   content = function(file) {
      req(volcanoGO)
      write.table(volcanoGO$getSigGO(), file, sep = ",", row.names = FALSE, col.names = TRUE)
   }
)

output$volcano_download <- downloadHandler(
   filename = function() {
      filename <- input$volcano_filename

      if (filename == "") {
         filename <- paste0(dataset$ID, "_volcano_", input$volcano_choice)
      }
      
      ext <- input$volcano_table_format

      return(paste0(filename, ".", ext))
   },
   content = function(file) {
      if (input$volcano_choice == "upreg") {
         if (is.null(volcano$tableUp)) {
            df <- data.table(ERROR = c(paste("No significant upregulated", tolower(dataset$idLabel))))
         } else {
            df <- volcano$tableUp
         }
      } else if (input$volcano_choice == "downreg") {
         if (is.null(volcano$tableDown)) {
            df <- data.table(ERROR = c(paste("No significant downregulated", tolower(dataset$idLabel))))
         } else {
            df <- volcano$tableDown
         }
      } else if (input$volcano_choice == "upanddown") {
         if (is.null(volcano$tableUp) && is.null(volcano$tableDown)) {
            df <- data.table(ERROR = c(paste("No significant", tolower(dataset$idLabel))))
         } else {
            df <- rbind(volcano$tableUp, volcano$tableDown)
         }
      } else if (input$volcano_choice == "full"){
         df <- buildVolcanoTable(volcano_binned())
      }
      # sep:
      if (input$volcano_table_format == "csv") {
         write.csv(df, file)
      } else if (input$volcano_table_format == "tsv") {
         write.table(df, file, sep = "\t")
      }
   }
)

observeEvent(input$volcano_go_color, {
   # Update the GO color when the color selector changes
   if(!is.null(volcano$goProteins) && length(volcano$goProteins) > 0) {
      volcano$goColor <- input$volcano_go_color
      showNotification("GO annotation color updated", duration = 1)
   }
})

observeEvent(input$volcano_go_bg_color, {
   # Background color change triggers plot update automatically via reactive
   showNotification("Background color updated", duration = 1)
})

observeEvent(input$volcano_go_show_labels, {
   # Label toggle triggers plot update automatically via reactive
   if(input$volcano_go_show_labels) {
      showNotification("GO labels enabled", duration = 1)
   } else {
      showNotification("GO labels disabled", duration = 1)
   }
})
