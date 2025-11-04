# Tab reactives
summary <- reactiveValues()

if(!is.null(GO)){
   goCalc <- GO$goCalculator()
   shinyjs::show("GO_header")
}else{
   goCalc <- NULL
}

updateProteinChoice <- function(mosaicID) {
   if(!is.null(mosaicID)){
      values$activeMosaicID <- mosaicID
      updateNearest(proteinTMT(), isProtein = TRUE)
   }else{ #Set to search results
      validSearch <- updateNearest(getSliderValues(), isProtein = FALSE)

      if(validSearch){
         closest <- names(summary$nearestProteins)[1]
         proteinSelect$addToHistory(closest, useCallback = FALSE) #Avoid callback loop
         values$activeMosaicID <- closest
      }
   }
}

updateNearest <- function(pattern, isProtein){
   searchable <- !is.na(pattern)
   if(sum(searchable) <= 1){
      sendSweetAlert(
         title = "Missing values",
         text = "At least two values are required for distance calculations.",
         type = "error",
         timer = 2000
      )
      return(FALSE)
   }
   if(input$distance_func == "corr" && sd(pattern[searchable]) == 0){
      sendSweetAlert(
         title = "Invalid search",
         text = "Standard deviation cannot be zero for correlation distance calculations.",
         type = "error",
         timer = 2000
      )
      return(FALSE)
   }

   columns <- values$normProtCols
   colnames(columns) <- theme$names
   subset <- subsetColumns(columns, pattern)
   n <- isolate(as.numeric(input$numNearest))

   nearestProteins <- getRowDistance(
      mat = t(subset$columns[, subset$keep, drop = FALSE]),
      vec = subset$key[subset$keep],
      method = input$distance_func)

   summary$nearestProteins <- nearestProteins[1:min(n, length(nearestProteins))]
   summary$isProtein <- isProtein
   values$nearestIDs <- names(summary$nearestProteins)

   subset$columns <- subset$columns[values$nearestIDs, , drop = FALSE]
   summary$nearestSubset <- subset

   # Update GO
   if(!is.null(goCalc)){
      goCalc$subsetData(values$nearestIDs)
      summary$sigGO <- goCalc$filterGO(
         as.numeric(input$goMinMatches),
         as.numeric(input$goPvalCutoff),
         input$goPvalAdjust)
   }

   return(TRUE)
}

output$ui_select_protein <- proteinSelect$addInput("select_protein", label = paste("Search for", dataset$idLabel))

if("Sequence" %in% names(dataset$proteins$info)) {
   output$current_seq <- renderUI({
      id <- values$activeMosaicID
      info <- dataset$proteins$info[id, ]
      seq <- info[["Sequence"]]

      color_residues <- TRUE
      mod_count <- str_count(seq, "#")
      if(mod_count > 0) {
         if("Site" %in% names(info)) {
            sites <- str_split(info[["Site"]], ";")[[1]]
            if(mod_count > length(sites)) {
               #Need to filter modified residues (#) by AA index in Site column
               color_residues <- FALSE

               sites <- as.integer(sites)
               c_vec <- strsplit(seq, "")[[1]]
               mod_ind <- which(c_vec == "#")
               c_vec <- c_vec[c_vec != "."]
               targets <- which(c_vec == "#") - seq(mod_count)
               pep_pattern <- paste0(c_vec[c_vec != "#"], collapse = "")

               prot_seq <- getProteinSequence(dataset, info[["UniprotID"]])
               locations <- str_locate_all(prot_seq, pep_pattern)[[1]]
               for(start in locations[, "start"]) {
                  offsets <- sites - start + 1
                  if(all(offsets %in% targets)) {
                     mod_found <- targets %in% offsets
                     extra_mods <- mod_ind[which(!mod_found)]

                     str_sub(seq, extra_mods, extra_mods) <- ""
                     color_residues <- TRUE
                     break
                  }
               }
            }
         }

         #Center on all modified sites
         site <- str_extract(seq, ".#(.*#)?")
         short_seq <- trunc_seqs(seq, site, width = 23, elipsis = "...", tags = NULL, ignore = NULL)

         if(color_residues) {
            colorRegex <- '<b class=site-color>\\1</b>'
            short_seq <- str_replace_all(short_seq, "(.)#", colorRegex)
         }

         seq <- HTML(sprintf('<div data-toggle="toggle" title="%s">%s</div>', seq, short_seq))
      }

      return(seq)
   })
}

observe({
   mosaicID <- req(input$switch_to_mosaicID)
   cat("SET button clicked with MosaicID:", mosaicID, "\n")
   proteinSelect$addToHistory(mosaicID)
   updateTabItems(session, "sidebarMenu", "summary")
   cat("Switched to summary tab and updated protein selection\n")
})

observeEvent(input$button_search_pattern, {
   updateProteinChoice(NULL)
})

observeEvent({
   input$button_update_parameters
}, {
   if(summary$isProtein){
      updateProteinChoice(values$activeMosaicID)
   }else{
      updateProteinChoice(NULL)
   }
})

observe({
   summary$data_heatmap <- calculateLogFC(dataset, values$normProtCols, values$nearestIDs)
})

buildPatternUI <- function(id, label, setId){
   renderUI({
      req(dataset$numGroups)
      labels <- c(rep(label, dataset$numGroups), "ALL")
      input_list <- lapply(1:(dataset$numGroups + 1), function(i) {
         list(
            actionButton(
               inputId = paste(id, i, sep = "_"),
               label = HTML(labels[i]),
               class = "minibttn",
               onclick = paste0('Shiny.setInputValue(\"', setId, '\", ', i, ')'),
               style = paste0("background-color:", c(theme$groupColors, "white")[i], ";"))
         )
      })
      tagList(input_list)
   })
}

output$pattern_ui_up <- buildPatternUI("up_but", "&uarr;", "button_group_up")
output$pattern_ui_mid <- buildPatternUI("mid_but", "-", "button_group_mid")
output$pattern_ui_down <- buildPatternUI("down_but", "&darr;", "button_group_down")

output$pattern_ui_use <- renderUI({
   req(dataset$numGroups)

   ui <- tagList()
   for (i in 1:(dataset$numGroups)) {
      ui[[i]] <- prettyCheckbox(
         inputId = paste0("checkbox_use_group_", i),
         label = NULL,
         value = TRUE,
         inline = TRUE,
         shape = "round",
         bigger = FALSE,
         width = "16.02px")
   }
   ui[[length(ui) + 1]] <- actionButton("button_toggle_use_up", "ALL", class = "minibttn usebttn")
   ui[[length(ui) + 1]] <- actionButton("button_toggle_use_down", "None", class = "minibttn usebttn", style = "width:40px;")

   ui
})

output$pattern_sliders <- renderUI({
   if(input$allSliders) {
      lapply(1:dataset$numSamples, function(i) {
         id <- paste0("patternSlider", i)
         value <- isolate(input[[id]])
         if(is.null(value)) {
            value <- 1
         }
         noUiSliderInput(
            inputId = id,
            min = 0, max = 1,
            value = value,
            height = 60,
            width = 15,
            orientation = "vertical",
            direction = "rtl",
            color = theme$colors[i],
            tooltips = FALSE,
            update_on = "end",
            behaviour = "snap",
            inline = TRUE
         )
      })
   }else {
      lapply(1:dataset$numGroups, function(i) {
         id <- paste0("patternSliderGroup", i)
         value <- isolate(input[[id]])
         if(is.null(value)) {
            value <- 1
         }
         noUiSliderInput(
            inputId = id,
            min = 0, max = 1,
            value = value,
            height = 60,
            width = 15,
            orientation = "vertical",
            direction = "rtl",
            color = theme$groupColors[i],
            tooltips = FALSE,
            update_on = "end",
            behaviour = "snap",
            inline = TRUE
         )
      })
   }
})

observe({
   lapply(seq(dataset$numGroups), function(group) {
      observeEvent({
         input[[paste0("checkbox_use_group_", group)]]
      }, {
         use <- input[[paste0("checkbox_use_group_", group)]]
         ids <- paste0("patternSliderGroup", group)
         if(isolate(input$allSliders)) {
            ids <- paste0("patternSlider", which(dataset$groups == group))
         }
         req(all(ids %in% names(input)))
         for (id in ids) {
            updateNoUiSliderInput(
                  inputId = id,
                  value = as.numeric(use),
                  disable = !use
               )
         }
      }, ignoreInit = TRUE)
   })
})

getSliderValues <- function() {
   active <- getActiveSliderGroups()
   vals <- rep(NA, dataset$numSamples)
   if(input$allSliders) {
      activeSliders <- which(active[dataset$groups])
      vals[activeSliders] <- sapply(activeSliders, function(slider) {
         input[[paste0("patternSlider", slider)]]
      })
   }else {
      active <- which(active)
      groupVals <- sapply(paste0("patternSliderGroup", 1:dataset$numGroups),
         function(id) {
            input[[id]]
         })
      activeCols <- dataset$groups %in% active
      vals[activeCols] <- groupVals[dataset$groups[activeCols]]
   }
   vals
}

getActiveSliderGroups <- function() {
   sapply(paste0("checkbox_use_group_", seq(dataset$numGroups)),
      function(id) {
         input[[id]]
      })
}

setSliders <- function(inputName, value){
   group <- as.numeric(input[[inputName]])
   activeGroups <- getActiveSliderGroups()

   setGroup <- function(group) {
      if(activeGroups[[group]]){
         if(input$allSliders){
            sliders <- paste0("patternSlider", which(group == dataset$groups))
            lapply(sliders, function(slider) {
               updateNoUiSliderInput(
                  inputId = slider,
                  value = value
               )
            })
         }else {
            updateNoUiSliderInput(
               inputId = paste0("patternSliderGroup", group),
               value = value
            )
         }
      }
   }

   if(group == dataset$numGroups + 1) { #ALL
      lapply(1:dataset$numGroups, setGroup)
   }else {
      setGroup(group)
   }
   session$sendCustomMessage(type = "resetValue", message = inputName)
}

observeEvent(input$button_group_up,   setSliders("button_group_up", 1))
observeEvent(input$button_group_mid,  setSliders("button_group_mid", 0.5))
observeEvent(input$button_group_down, setSliders("button_group_down", 0))

observeEvent(input$button_reset_pattern, {
   if(input$allSliders){
      for (i in 1:dataset$numSamples) {
         updateNoUiSliderInput(
            inputId = paste0("patternSlider", i),
            value = 1)
      }
   }else {
      for (i in 1:dataset$numGroups) {
         updateNoUiSliderInput(
            inputId = paste0("patternSliderGroup", i),
            value = 1)
      }
   }
   for (i in 1:dataset$numGroups) {
      updatePrettyCheckbox(
         inputId = paste0("checkbox_use_group_", i),
         value = TRUE)
   }
   showNotification(duration = 4, "Successfully reset pattern.")
})

observeEvent(input$button_toggle_use_up, {
   for (i in 1:dataset$numGroups) {
      updatePrettyCheckbox(
         inputId = paste0("checkbox_use_group_", i),
         value = TRUE)
   }
   showNotification(duration = 4, "Successfully turned on all groups.")
})

observeEvent(input$button_toggle_use_down, {
   for (i in 1:dataset$numGroups) {
      updatePrettyCheckbox(
         inputId = paste0("checkbox_use_group_", i),
         value = FALSE)
   }
   showNotification(duration = 4, "Successfully turned off all groups.")
})

observeEvent(input$allSliders, {
   sapply(1:dataset$numGroups, function(group){
      use <- input[[paste0("checkbox_use_group_", group)]]
      ids <- paste0("patternSlider", which(dataset$groups == group))
      groupId <- paste0("patternSliderGroup", group)
      if(input$allSliders) {
         value <- input[[groupId]]
      }else {
         value <- mean(sapply(ids, function(id){input[[id]]}))
         ids <- groupId
      }
      for (id in ids) {
         updateNoUiSliderInput(
               inputId = id,
               value = value,
               disable = !use
            )
      }
   })
}, ignoreInit = TRUE)

observeEvent(input$toBioplex,    updateTabItems(session, "sidebarMenu", "bioplex"))
observeEvent(input$toGOplotter,  updateTabItems(session, "sidebarMenu", "go"))
observeEvent(input$toVolcano,    updateTabItems(session, "sidebarMenu", "volcano"))
observeEvent(input$toPCA,        updateTabItems(session, "sidebarMenu", "pca"))
observeEvent(input$toPeptides,   updateTabItems(session, "sidebarMenu", "peptides"))
observeEvent(input$toFig,        updateTabItems(session, "sidebarMenu", "tmtBarChart"))

output$table_nearest_proteins <- renderDataTable({
   req(summary$nearestProteins)

   nearestIDs <- values$nearestIDs
   info <- dataset$proteins$info[nearestIDs, ]
   df <- info %>%
      transmute(
         Rank = 1:length(values$nearestIDs),
         UniprotID,
         GeneSymbol,
         Description,
         Dist = formatC(summary$nearestProteins, digits = 2, format = "e"),
         NumPeps = getNumPepsFromID(info, MosaicID)
      )

   if(dataset$isSiteQuant){
      protSeqs <- getProteinSequence(dataset, df$UniprotID)
      df <- df %>% mutate(
         Sequence = info$Sequence,
         ProtSize = nchar(protSeqs),
         site = info$Site
      ) %>% mutate(
         seqTitle = paste0("Site Pos: ", site, "/", ProtSize, "<br>", trunc_seqs(protSeqs, Sequence)),
      )
   }

   # Check if any GO categories were found. If so, create the list
   if (!is.null(summary$sigGO)) {
      df$GO <- goCalc$protTooltips()
   }

   # Add the buttons to change to the proteins
   df$SET <- nearestIDs

   # Style Columns
   #           Do not modify columns after this point
   colDefs <- c(
      addColumnClasses(intersect(names(df), c("GO", "SET")), 'dt-center'),
      addLink("https://www.uniprot.org/uniprotkb/", "UniprotID", "/entry"),
      addShinyButton("SET", "switch_to_mosaicID"))

   if(dataset$isSiteQuant)
      colDefs <- c(colDefs, addPosMarker(df, "ProtSize", "site", "seqTitle", enableHTML = TRUE))

   dt <- datatable(df, options = list(
      initComplete = datatableTheme(),
      columnDefs = colDefs,
      searchHighlight = TRUE,
      pageLength = 10),
      class = "compact cell-border",
      rownames = FALSE,
      selection = 'none',
      escape = FALSE)

   return(dt)
})

output$table_nearest_go <- renderDataTable({
   req(goCalc)
   shiny::validate(need(summary$sigGO, "No significant GO categories found, try widening parameters."))
   goCalc$goTable()
})

output$summaryPlots <- renderUI({
   req(!is.null(input$extend_plots))
   
   copy_val <- function(id) {
      val <- input[[id]]
      ifelse(is.null(val), FALSE, val)
   }

   controls <- tagList()
   width <- 12
   if(dataset$isSiteQuant) {
      if(dataset$areReps) {
         width <- 6
      }
      controls <- tagList(controls,
         column(width, align = "center",
            checkboxInput("show_protein_sites", "All Protein Sites",
               value = copy_val("show_protein_sites")))
      )
   }
   if(dataset$areReps) {
      controls <- tagList(controls,
         column(width, align = "center",
            checkboxInput("summarize_reps", "Summarize Reps",
               value = copy_val("summarize_reps")))
      )
   }

   # Reformat ui to extend plots
   if(input$extend_plots){
      hm_height <- paste0(max(400, 30 * length(values$nearestIDs)) * input$heatmap_size, "px")

      ui <- span(align = "left",
         fluidRow(controls),
         plotlyOutput("plot_current_protein"),
         plotlyOutput("plot_nearest_proteins"),
         plotOutput("plot_heatmap", height = hm_height)
      )
   }else {
      ui <- fluidRow(
         column(4, align = "left",
            fluidRow(controls),
            plotlyOutput("plot_current_protein")
         ),
         column(4, align = "left", plotlyOutput("plot_nearest_proteins")),
         column(4, align = 'left', plotOutput("plot_heatmap"))
      )
   }

   return(ui)
})

observe({
   input$extend_plots # Resize plot when extended
   output$plot_current_protein <- renderPlotly(current_protein())
})

current_protein <- reactive({
   mosaicID <- req(values$activeMosaicID)

   summarize_reps <- FALSE
   if(dataset$areReps){
      req(!is.null(input$summarize_reps))
      summarize_reps <- input$summarize_reps
   }

   geneSymbol <- dataset$proteins$info[mosaicID, "GeneSymbol"]
   title <- paste0(geneSymbol, " using {",
      paste(getNumPepsFromID(dataset$proteins$info, mosaicID), collapse = ","),
      "} peptides")
   yLab <- ifelse(advanced$isRA, "TMT RA", "TMT S/N")

   groups <- groupColumns(dataset$classes, dataset$numPlexes, theme$groupNames[dataset$groups], dataset$areReps)
   g <- tmtBarchart(proteinTMT(), theme$colors, theme$names, groups, dataset$areReps,
      summarize_reps, title, yLabel = yLab)
   data_ids <- mosaicID

   if(dataset$isSiteQuant) {
      req(!is.null(input$show_protein_sites))
      if(input$show_protein_sites) {
         uniprotID <- dataset$proteins$info[mosaicID, "UniprotID"]
         full_protein <- dataset$proteins$info[dataset$proteins$info$UniprotID == uniprotID, "MosaicID"]
         if(length(full_protein) > 1) {
            data_ids <- full_protein
            protein_data <- data.frame()
            for(site in setdiff(full_protein, mosaicID)) {
               protein_data <- rbind(
                  protein_data,
                  data.frame(
                     x = tmt_scale_x(groups, summarize_reps),
                     RA = unlist(values$normProtCols[site, ]),
                     MosaicID = site,
                     Name = dataset$proteins$info[[site, "GeneSymbol"]],
                     color = theme$colors,
                     group = groups
                  ) %>% arrange(x)
               )
            }

            if(summarize_reps) {
               protein_data <- protein_data %>%
                  group_by(x, group, MosaicID, color, Name) %>%
                  summarize(RA = mean(RA, na.rm = TRUE, Name = Name, color = color))
            }

            g <- g +
               geom_line(data = protein_data,
                  aes(
                     x = x,
                     y = RA,
                     group = MosaicID
                  ),
                  color = "black",
                  na.rm = TRUE,
                  show.legend = FALSE) +
               geom_point(data = protein_data,
                  aes(
                     x = x,
                     y = RA,
                     group = MosaicID,
                     text = paste0("Site: ", Name, "\nClick to switch to site"),
                     fill = I(color),
                     customdata = MosaicID
                  ),
                  shape = 21,
                  color = "black",
                  size = 2,
                  na.rm = TRUE,
                  show.legend = FALSE)
         }
      }
   }

   d <- dataset$proteins$columns[data_ids, , drop = FALSE]
   colnames(d) <- theme$names

   # Save downloadable data:
   summary$tmt_plot <- g
   summary$tmt_data <- cbind(dataset$proteins$info[data_ids, , drop = FALSE], d)

   ggplotly(g,
         tooltip = "text",
         source = "current_protein") %>%
      plotlyDefaults(showlegend = FALSE) %>%
      event_register("plotly_click")
})

observe({
   id <- event_data("plotly_click", priority = "event", source = "current_protein")$customdata
   if(!is.null(id)) {
      proteinSelect$addToHistory(id)
   }
})

observe({
   input$extend_plots # Resize plot when extended
   output$plot_nearest_proteins <- renderPlotly(nearest_proteins())
})

nearest_proteins <- reactive({
   subset <- req(summary$nearestSubset)
   tmt_values <- subset$columns
   numNearest <- nrow(tmt_values)

   if(!summary$isProtein){
      title <- "Pattern"
      highlight <- "green"

      # Treat the pattern as a protein and draw it on top
      tmt_values <- rbind(
         PATTERN = subset$key,
         tmt_values
      )
   }else{
      title <- dataset$proteins$info[values$activeMosaicID, "GeneSymbol"]
      highlight <- "red"
   }

   positions <- seq(ncol(tmt_values))
   long_data <- data.frame()
   for (row in 1:nrow(tmt_values)) {
      long_data <- rbind(
         long_data,
         data.frame(
            means = tmt_values[row, ],
            ProteinId = rownames(tmt_values)[row],
            x = positions
         )
      )
   }

   maxY <- 100
   if(any(tmt_values > maxY, na.rm = TRUE)){
      showNotification(paste("Some values were too big and have been plotted as", maxY), duration = 2)
   }

   g <- ggplot() +
      aes(
         x = x,
         y = pmin(means, maxY), # clamp values to preserve view
         group = factor(ProteinId, levels = rownames(tmt_values)),
         text = paste0(ProteinId, ": ", formatC(means, digits = tooltipDigits), "<br>Click to switch to ", tolower(dataset$idLabel)),
         customdata = ProteinId
      ) %>% suppressWarnings +
      geom_line(
         data = long_data[-positions, ],
         color = "#8819A8",
         linewidth = 0.2,
         alpha = 0.5,
         na.rm = TRUE,
         show.legend = FALSE) +
      geom_line(
         data = long_data[positions, ],
         color = highlight,
         linewidth = 1,
         na.rm = TRUE,
         show.legend = FALSE) +
      expand_limits(y = 0) +
      scale_x_continuous(breaks = seq(dataset$numSamples), labels = theme$names) +
      ggtitle(paste0(numNearest, " Closest ", dataset$idLabel, "s To ", title)) + # Should always be plural
      xlab(NULL) +
      ylab("TMT RA")

   distances <- summary$nearestProteins
   if(!summary$isProtein){
      distances <- c(0, distances)
   }
   d <- cbind(dataset$proteins$info[rownames(tmt_values), , drop = FALSE],
      data.frame(Distance = distances), tmt_values)
   d$MosaicID <- rownames(tmt_values) # Copy 'Pattern' into the info section

   # Save for download option:
   summary$nearest_plot <- g
   summary$nearest_data <- d

   ggplotly(g, tooltip = "text", source = "nearest") %>%
      plotlyDefaults() %>%
      layout(showlegend = FALSE)
})

observe({
   mosaicID <- req(event_data("plotly_click", source = "nearest"))$customdata
   req(mosaicID != "Pattern")
   proteinSelect$addToHistory(mosaicID)
})

output$plot_heatmap <- renderPlot({
   logFC <- req(summary$data_heatmap)
   if(input$heatmap_use_group_colors){
      colors <- theme$colors
   }else{
      colors <- c(input$heatmap_lower_col, input$heatmap_upper_col)
   } 
   fcCutoff <- as.numeric(input$heatmap_fc)
   isDotPlot <- input$summary_heatmap_type == "Dotplot"
   highlightedGS <- dataset$proteins$info[values$activeMosaicID, "GeneSymbol"]

   g <- dendroHeatmap(logFC, colors, fcCutoff, isDotPlot, colNames = theme$names, highlightedGS)

   d <- logFC
   colnames(d) <- theme$names

   # Save for download option:
   summary$heatmap_plot <- g
   summary$heatmap_data <- cbind(dataset$proteins$info[isolate(values$nearestIDs), ], d)

   return(g)
})

output$summary_plot_download <- downloadHandler(
   filename = function() {
      filename <- input$summary_plot_filename
      type <- input$summary_plot_choice

      if (filename == "") {
         filename <- switch(type,
            "tmt" = values$activeMosaicID,
            "nearest" = paste0(length(values$nearestIDs), "_nearest_prots_to_", values$activeMosaicID),
            "heatmap" = paste0(values$activeMosaicID, "_heatmap"),
            "all" = values$activeMosaicID)
         filename <- paste0(dataset$ID, "_", filename)
      }

      ext <- tolower(input$summary_plot_format)
      if(type == "all") {
         ext <- "zip"
      }

      return(paste0(filename, ".", ext))
   },
   content = function(file) {
      type <- input$summary_plot_choice
      format <- tolower(input$summary_plot_format)

      if(type == "all") {
         repeat {
            #Ensure unique path
            d <- paste0("data/tmp/", stringi::stri_rand_strings(1, 8))
            if(!dir.exists(d))
               break
         }
         dir.create(d, recursive = TRUE)
         lastDir <- setwd(d)
         zipFolder <- paste0(dataset$ID, "_", values$activeMosaicID)
         dir.create(zipFolder)

         types <- paste0(downloadables$summary, "_plot")
         files <- paste0(zipFolder, "/", types, ".", format)
         names(files) <- types
         lapply(types, function(type) {
            plot <- summary[[type]]
            if(!is.null(plot)) {
               ggsave(files[[type]], plot, format)
            }
         })

         zip(file, files)
         setwd(lastDir)
         unlink(d, recursive = TRUE)
      }else {
         type <- paste0(type, "_plot")
         plot <- req(summary[[type]])
         ggsave(file, plot, format)
      }
   }
)

output$summary_data_download <- downloadHandler(
   filename = function() {
      filename <- input$summary_data_filename
      type <- input$summary_data_choice

      if (filename == "") {
         filename <- switch(type, 
            "tmt" = values$activeMosaicID,
            "nearest" = paste0(length(values$nearestIDs), "_nearest_prots_to_", values$activeMosaicID),
            "heatmap" = paste0(values$activeMosaicID, "_heatmap")
         )
         filename <- paste0(dataset$ID, "_", filename)
      }

      return(paste0(filename, ".", tolower(input$summary_data_format)))
   },
   content = function(file) {
      type <- paste0(input$summary_data_choice, "_data")
      format <- tolower(input$summary_data_format)

      sep <- switch(format,
         "CSV" = ",",
         "TSV" = "\t",
         ","
      )

      data <- req(summary[[type]])
      write.table(data, file, sep = sep, row.names = FALSE, col.names = TRUE)
   }
)

output$go_download <- downloadHandler(
   filename = function() {
      return(paste(dataset$ID, values$activeMosaicID, "GO.csv", sep = "_"))
   },
   content = function(file) {
      req(goCalc)
      write.table(goCalc$getSigGO(), file, sep = ",", row.names = FALSE, col.names = TRUE)
   }
)

