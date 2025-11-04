# Constants
defaults <- list(
   factors = rep(1, dataset$numSamples),
   norm = "none",
   row = NULL,
   scaledSum = 100 * dataset$proteins$numClasses,
   baseColumns = NULL
)

if(!is.null(dataset$normalization)) {
   defaults$norm <- dataset$normalization$normBy
   defaults$row <- dataset$normalization$normRow
   if(!is.null(dataset$normalization$factors)) {
      defaults$factors <- dataset$normalization$factors
   }
}

hasClasses <- dataset$proteins$numClasses > 1
rowType <- ifelse(dataset$isSiteQuant, "Site", "Protein")
normChoices <- setNames(c("skip", "none", "row", "all", "manual"), c("Default", "None", rowType, "Column Sums", "Manual"))

suffix <- dataset$proteins$suffixes
if(suffix == "sum") {
   defaults$baseColumns <- dataset$proteins$columns
}

hasSum <- "sum" %in% dataset$proteins$suffixes_available

# Tab reactives
advanced <- reactiveValues(
   bridge = rep("None", dataset$numPlexes),
   isRA = TRUE,
   factors = defaults$factors,
   steps = paste("Using", suffix, "columns")
)

# Set up UI
reset <- c(
   #inputID = reset_func()
   norm_by = function() { updateRadioGroupButtons(session, "norm_by", choices = normChoices, selected = "skip") },
   norm_row = function() { updateSelectizeInput(session, "norm_row", rowType, choices = dataset$proteins$info$MosaicID, selected = defaults$row, server = TRUE) },
   rescale = function() { updateSwitchInput(session, "rescale", value = FALSE) }
)
lapply(paste0("bridge_norm_", unique(dataset$classes)), function(id) {
   reset[[id]] <<- function() { updateSelectInput(session, id, selected = "None") }
})

if(hasClasses) {
   reset$norm_plex <- function() { updatePrettyCheckbox(session, "norm_plex", value = FALSE) }
   reset$class_rescale <- function() { updateSwitchInput(session, "class_rescale", value = TRUE) }
   shinyjs::show("norm_plex")
   shinyjs::show("class_rescale_ui")
}

output$bridge_norm_ui <- renderUI({
   lapply(unique(dataset$classes), function(plex) {
      plexIndex <- dataset$classes == plex
      id <- paste0("bridge_norm_", plex)
      choices <- setNames(c("None", which(plexIndex)), c("None", theme$names[plexIndex]))
      selected <- input$id
      if(is.null(selected) || !selected %in% choices) {
         selected <- "None"
      }
      selectInput(inputId = id,
         label = paste0("Bridge for ", plex),
         choices = choices,
         selected = selected
      )
   })
})

output$custom_col_norm_ui <- renderUI({
   if(hasSum) {
      ui <- tagList(
         div(textOutput("default_norm_by"), class = "h4"),
         hr(style = "border-top: 4px solid #497EA5;"),
         h3("Live Normalization Factors"),
         plotlyOutput("plot_norm_factors", height = 230),
         h3("Custom Normalization"),
         fluidRow(
            column(7,
               radioGroupButtons("norm_by", "Normalize each row by", normChoices),
            ),
            column(5,
               shinyjs::disabled(
                  selectInput("norm_row", rowType, choices = defaults$row, selected = defaults$row)
               )
            )
         ),
         uiOutput("manual_norm_ui")
      )

      if(!is.null(dataset$peptides)) {
         ui <- append(ui, list(
            radioGroupButtons("col_norm_peps", "Peptide Normalization",
               c("Apply" = TRUE, "Skip" = FALSE), select = TRUE)
         ))
         reset[["col_norm_peps"]] <- function() {
            updateRadioGroupButtons(session, "col_norm_peps", selected = FALSE)
         }
      }

      observeEvent(input$norm_row, once = TRUE, reset$norm_row())
   }else {
      ui <- h3("*Sum data is unavailable, renormalization disabled")
   }

   return(ui)
})

output$manual_norm_ui <- renderUI({
   if(input$norm_by != "manual")
      return(NULL)

   return(rHandsontableOutput("manual_norm"))
})

output$manual_norm <- renderRHandsontable({
   table <- t(data.frame(Factor = defaults$factors))
   colnames(table) <- theme$names
   rhandsontable(table, allowEmpty = FALSE, allowInvalid = FALSE,
         allowInsertRow = FALSE, allowInsertColumn = FALSE) %>%
      hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
})

# Tab logic ------------------------------------------

toDefaults <- function(ids = names(reset)) {
   for(id in ids) {
      reset[[id]]()
   }
}

observeEvent(input$col_norm_info, {
   content <- list(paste("If sum data is available, new normalization factors can be applied dynamically.",
      "These factors are multiplied with", rowType, "expression data to readjust each column globally."),
      tags$ul(
         tags$li("Default - use data as-is, skips column normalization"),
         tags$li("None - removes any known normalization"),
         tags$li(paste(rowType, "- uses factors that are the inverse of existing data so that the results are ratios to the selected", rowType)),
         tags$li("Column Sums - generates factors so that after this step the sum of each column is equal to the max of the original column sums"),
         tags$li("Manual - allow `arbitrary normalizations to be loaded into the viewer")
      )
   )

   showModal(modalDialog(title = "Column Normalization", h4(content), easyClose = TRUE))
})

observeEvent(input$bridge_info, {
   rowTypes <- paste0(rowType, "s")
   content <- paste("Within each class, every", rowType, "will be scaled relative to its 'bridge' column.",
      "The transformation preserves the original ratios between columns within each class",
      "for individual", rowTypes, "but adjusts the classes one", rowType, " at a time so that the bridge",
      "columns are always the same. The final value of the bridge columns after this step",
      "will be the mean of all the original bridges across all", rowTypes, "to facilitate",
      "direct comparison across multiple plexes.")

   showModal(modalDialog(title = "Bridge Normalization", h4(content), easyClose = TRUE))
})

observeEvent(input$row_scaling_info, {
   content <- list(paste("Every", rowType, "will be rescaled relative to its sum.",
      "After this step, the sum of each", rowType, "will be equal to the number of distinct classes times 100.",
      "In the case of multiple classes, the 'Scale within classes' toggle controls how rescaling is performed."),
      tags$ul(
         tags$li("If enabled, each class is scaled to 100 individually"),
         tags$li("Otherwise, the full", rowType, " is normalization together so that the sum of all classes is fixed, but the sum of any given class is not")
      )
   )

   showModal(modalDialog(title = "Row Scaling", h4(content), easyClose = TRUE))
})

observeEvent(input$reset_advanced, {
   toDefaults()
   advanced$bridge[unique(dataset$classes)] <- "None"

   applyNormalization("skip", NA)
})

observeEvent(input$apply_advanced, {
   norm_by <- input$norm_by
   if(is.null(norm_by)) {
      showNotification("Please select a normalization method")
      return(NULL)
   }

   norm_row <- NA
   if(norm_by == "row") {
      norm_row <- input$norm_row
      if(is.null(norm_row) || is.na(norm_row) || nchar(norm_row) == 0) {
         showNotification(paste("Please select a", tolower(rowType),  "to normalize"))
         return(NULL)
      }
   }

   rescale <- input$rescale
   rescale <- if(is.null(rescale)) FALSE else as.logical(rescale)

   if(rescale) {
      class_rescale <- input$class_rescale
      class_rescale <- if(is.null(class_rescale)) FALSE else as.logical(class_rescale)
   }else {
      class_rescale <- FALSE
   }

   col_norm_peps <- input$col_norm_peps
   col_norm_peps <- if(is.null(col_norm_peps)) FALSE else as.logical(col_norm_peps)

   advanced$bridge <- sapply(unique(dataset$classes), function(id) { input[[paste0("bridge_norm_", id)]] })

   applyNormalization(norm_by, norm_row, rescale, class_rescale, col_norm_peps)
})

observe({
   norm_by <- req(input$norm_by)
   if(norm_by == "row") {
      shinyjs::enable("norm_row")
   }else{
      shinyjs::disable("norm_row")
   }
})

observe({
   rescale <- input$rescale
   disabled <- is.null(rescale) || !as.logical(rescale)
   updateSwitchInput(session, "class_rescale", disabled = disabled)
})

output$norm_steps <- renderUI({
   steps <- advanced$steps
   if(length(steps) == 1) {
      steps <- paste(steps, "(no action taken)")
   }
   step_sequence <- paste(steps, collapse = " > ")
   ui <- sprintf("<strong>Steps to Live Data:</strong> %s", step_sequence)
   HTML(ui)
})

beginNotify <- function() {
   advanced$steps_tmp <- c()
}

notifyNorm <- function(msg, ...) {
   msg <- paste(msg, ...)
   showNotification(msg, duration = 6)
   advanced$steps_tmp <- c(advanced$steps_tmp, msg)
}

endNotify <- function() {
   advanced$steps <- advanced$steps_tmp
}

applyNormalization <- function(norm_by, norm_row, rescale = FALSE, class_rescale = rescale, col_norm_peps = FALSE) {
   tryCatch({
      beginNotify()

      protein_data <- dataset$proteins$columns
      pep_data <- dataset$peptides$columns
      colFactors <- defaults$factors

      live_suffix <- suffix
      doColNorm <- hasSum && norm_by != "skip"
      if(doColNorm) {
         if(is.null(defaults$baseColumns)) {
            showNotification("Loading sum columns")
            baseData <- src$proteins
            showNotification("Loading sum columns 2")
            if(is.null(baseData))
               stop("Unable to column normalize, failed to load sum data")

            defaults$baseColumns <- baseData$columns[, dataset$order]
            showNotification("Loading sum columns 3")
         }
         protein_data <- defaults$baseColumns
         live_suffix <- "sum"
      }
      notifyNorm("Using", live_suffix, "columns")

      if(doColNorm) {
         factors <- columnNormFactors(protein_data, norm_by, norm_row)
         if(!is.null(factors)) {
            protein_data <- t(t(protein_data) * factors)
            colFactors <- defaults$factors * factors

            if(!is.null(dataset$peptides) && col_norm_peps) {
               pep_colIDs <- dataset$peptides$columnIDs
               for(class in names(pep_data)) {
                  class_colIDs <- pep_colIDs[[class]]
                  class_factors <- colFactors[class_colIDs]
                  pep_data[[class]] <- t(t(pep_data[[class]]) * class_factors)
               }

               notifyNorm("Applied column normalization to peptides")
            }
         }
      }

      hasBridge <- advanced$bridge != "None"
      if(any(hasBridge)) {
         bridge_cols <- as.integer(advanced$bridge[hasBridge])
         target <- mean(protein_data[, bridge_cols], na.rm = TRUE)

         u_classes <- unique(dataset$classes)
         for(plex in u_classes[hasBridge]) {
            bridge <- as.integer(advanced$bridge[plex])
            plex_cols <- plex == dataset$classes

            # Normalize plex_cols so that bridge is *target* in every row
            protein_data[, plex_cols] <- target * protein_data[, plex_cols] / protein_data[, bridge]
         }
         notifyNorm("Bridge normalized row-wise to bridge mean per class")
      }

      isRA <- live_suffix == "scaled"
      if(rescale) {
         msg <- "Scaled row sums to"
         if(class_rescale) {
            classes <- dataset$proteins$classIDs
            max_class_cols <- max(table(classes))
            row_total <- 0
            for(class in unique(classes)) {
               class_cols <- class == classes
               row_sum <- 100 * sum(class_cols) / max_class_cols

               class_data <- protein_data[, class_cols]
               protein_data[, class_cols] <- row_sum * class_data / rowSums(class_data, na.rm = TRUE)
               row_total <- row_total + row_sum
            }
            msg <- sprintf("%s %d within classes (100 per %d columns)", msg, row_total, max_class_cols)
         }else {
            protein_data <- defaults$scaledSum * protein_data / rowSums(protein_data, na.rm = TRUE)
            msg <- paste(msg, defaults$scaledSum)
         }
         notifyNorm(msg)
         isRA <- TRUE
      }

      #Update reactives
      advanced$isRA <- isRA
      advanced$factors <- colFactors
      values$normProtCols <- protein_data
      values$normPepCols <- pep_data

      #Update nearest
      lastID <- NULL
      if(summary$isProtein) {
         lastID <- values$activeMosaicID
      }
      updateProteinChoice(lastID)

      endNotify()
   }, error = function(e) {
      showModal(modalDialog(title = "Failed to normalize", e$message, easyClose = TRUE))
   })
}

columnNormFactors <- function(protein_data, norm_by, norm_row) {
   toRatios <- function(factors) max(factors, na.rm = TRUE) / factors
   factors <- NULL
   if(norm_by == "skip") {
      msg <- "Using default column normalization"
   }else if(norm_by == "none") {
      if(defaults$norm != "none") {
         normLabel <- names(normChoices)[normChoices == defaults$norm]
         factors <- 1 / defaults$factors
         msg <- paste("Removed column normalization by", normLabel)
      }else {
         msg <- "Using default (none) column normalization"
      }
   }else if(norm_by == "row") {
      if(!is.na(norm_row)) {
         if(!identical(defaults$row, norm_row)) {
            if(norm_row %in% rownames(protein_data)) {
               factors <- toRatios(protein_data[norm_row, ])
               msg <- paste0("Normalized columns to ", norm_row)
            }else {
               stop("Invalid MosaicID")
            }
         }else {
            msg <- paste("Using default normalization to", normRow)
         }
      }else {
         stop("No row selected")
      }
   }else if(norm_by == "all") {
      msg <- "Normalized to current column sums"
      factors <- toRatios(colSums(protein_data, na.rm = TRUE))
   }else if(norm_by == "manual") {
      msg <- "Normalized to manual normalization factors"
      norm_df <- as.numeric(hot_to_r(req(input$manual_norm)))
      factors <- norm_df / defaults$factors
   }else {
      stop("Unknown normalization option")
   }

   notifyNorm(msg)
   return(factors)
}

output$default_scale <- renderText({
   msg <- "Unknown"
   if(suffix == "sum") {
      msg <- "Unscaled"
   }else if(suffix == "scaled") {
      msg <- "Scaled to 100"
      if(hasClasses) {
         msg <- paste(msg, "within Classes")
      }
   }
   sprintf("Default: %s", msg)
})

output$plot_default_factors <- renderPlotly({
   factors <- dataset$normalization$factors
   shiny::validate(need(!is.null(factors), "*Note: Server did not return normalization factors, assuming data is un-normalized"))
   normPlot <- tmtBarchart(factors, theme$colors, theme$names, groups = NULL, title = NULL, yLabel = "Multiplied by")

   ggplotly(normPlot, tooltip = "text") %>%
      plotlyDefaults()
})

output$default_norm_by <- renderText({
   msg <- switch(dataset$normalization$normBy,
      "all" = "Column Sums",
      "row" = dataset$normalization$normRow,
      "none" = "None",
      "Unknown")
   sprintf("Default: %s", msg)
})

output$plot_norm_factors <- renderPlotly({
   factors <- req(advanced$factors)
   normPlot <- tmtBarchart(factors, theme$colors, theme$names, groups = NULL, title = NULL, yLabel = "Multiplied by")

   ggplotly(normPlot, tooltip = "text") %>%
      plotlyDefaults()
})

