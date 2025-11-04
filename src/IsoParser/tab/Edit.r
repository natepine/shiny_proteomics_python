theme <- reactiveValues(
   groupColors = dataset$groupColors,
   groupNames = dataset$groupNames,
   colors = dataset$columnColors,
   names = dataset$names
)

# Tab variables
plotDownload <- list()
themeInputIds <- list(
   groupColors = paste0("theme_color_", 1:dataset$numGroups),
   names = paste0("theme_name_", 1:dataset$numSamples),
   groupNames = paste0("theme_group_name_", 1:dataset$numGroups)
)

getFigureParam <- function(input_val, default, escape) {
   if (is.null(input_val)) {
      return(default)
   }
   else if (input_val == escape) {
      return(default)
   }
   else {
      return(input_val)
   }
}

output$ui_tmtBarChart <- renderUI({
   plotOutput("tmtBarChart", width = input$tmtBarChart_width, height = input$tmtBarChart_height)
})

observeEvent(input$apply_theme, {
   setInputs <- function(inputType) {
      idSet <- themeInputIds[[inputType]]
      if(!is.null(input[[idSet[[1]]]])){
         theme[[inputType]] <- sapply(idSet, function(id) {input[[id]]}) %>% unname
         return(TRUE)
      }
      return(FALSE)
   }
   if(setInputs("groupColors")){
      theme$colors <- theme$groupColors[dataset$groups]
   }
   setInputs("names")
   setInputs("groupNames")
})

observeEvent(input$reset_theme_colors, {
   idSet <- themeInputIds$groupColors
   sapply(1:dataset$numGroups, function(i) {
      updateColourInput(session, idSet[[i]], value = dataset$groupColors[[i]])
   })
})

observeEvent(input$reset_theme_names, {
   idSet <- themeInputIds$names
   sapply(1:dataset$numSamples, function(i) {
      updateTextInput(session, idSet[[i]], value = dataset$names[[i]])
   })
})

observeEvent(input$reset_theme_group_names, {
   idSet <- themeInputIds$groupNames
   sapply(1:dataset$numGroups, function(i) {
      updateTextInput(session, idSet[[i]], value = dataset$groupNames[[i]])
   })
})

output$theme_colors_ui <- renderUI({
   names <- theme$groupNames
   colors <- theme$groupColors
   idSet <- themeInputIds$groupColors
   lapply(1:dataset$numGroups, function(i) {
      colourInput(inputId = idSet[[i]], label = names[[i]], colors[[i]])
   })
})

output$theme_names_ui <- renderUI({
   names <- theme$names
   idSet <- themeInputIds$names
   lapply(1:dataset$numSamples, function(i) {
      textInput(inputId = idSet[[i]], label = names[[i]], value = names[[i]])
   })
})
output$theme_group_names_ui <- renderUI({
   names <- theme$groupNames
   idSet <- themeInputIds$groupNames
   lapply(1:dataset$numGroups, function(i) {
      textInput(inputId = idSet[[i]], label = names[[i]], value = names[[i]])
   })
})

output$tmtBarChart <- renderPlot({
   mosaicID <- req(values$activeMosaicID)
   label_size <- req(input$tmtBarChart_label_size)
   title_size <- req(input$tmtBarChart_title_size)

   # Independent radioButtons could be added to this tab instead of input$summarize_reps
   summarize_reps <- FALSE
   if(dataset$areReps){
      req(!is.null(input$summarize_reps))
      summarize_reps <- input$summarize_reps
   }

   # Use custom values for title, names, colors, etc.
   numPeps <- paste(getNumPepsFromID(dataset$proteins$info, mosaicID), collapse = ",")
   geneSymbol <- dataset$protein$info[values$activeMosaicID, "GeneSymbol"]
   title <- getFigureParam(input$tmtBarChart_title, paste0(geneSymbol, " using {", numPeps, "} peptides"), "")
   xlab  <- getFigureParam(input$tmtBarChart_title_x, NULL, "")
   ylab  <- getFigureParam(input$tmtBarChart_title_y, ifelse(advanced$isRA, "TMT RA", "TMT S/N"), "")

   groups <- groupColumns(dataset$classes, dataset$numPlexes, theme$groupNames[dataset$groups], dataset$areReps)
   g <- tmtBarchart(proteinTMT(), theme$colors, theme$names, groups, dataset$areReps,
      summarize_reps, title, xlab, ylab)

   g <- g + theme(
      axis.text = element_text(size = label_size),
      title     = element_text(size = title_size)
   )

   # Save params for downloading:
   plotDownload$tmtBarChart <<- g
   plotDownload$filename <<- mosaicID

   return(g)
})

output$tmtBarChart_download <- downloadHandler(
   filename = function() {
      filename <- input$tmtBarChart_filename
      if (is.null(filename) || filename == "") {
         filename <- paste0(dataset$ID, "_", plotDownload$filename)
      }
      return(paste0(filename, ".", tolower(input$tmtBarChart_format)))
   },
   content = function(file) {
      if(!is.null(plotDownload$tmtBarChart)){
         ggsave(file, plot = plotDownload$tmtBarChart, width = input$tmtBarChart_width / 70,
            height = input$tmtBarChart_height / 70, device = tolower(input$tmtBarChart_format))
      }else{
         showNotification("Error downloading plot")
      }
   }
)
