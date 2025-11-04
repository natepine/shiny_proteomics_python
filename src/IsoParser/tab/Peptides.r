# Reactive Values
peptides <- reactiveValues()

activeUniprotID <- reactive({
   dataset$proteins$info[[values$activeMosaicID, "UniprotID"]]
})

# Tab constants
if(!is.null(dataset$peptides)) {
   unique_classes <- unique(dataset$peptides$info$Class)
}

output$sequenceTitle <- renderText({
   dataset$proteins$info[values$activeMosaicID, c("GeneSymbol", "UniprotID")] %>%
      mutate(GeneSymbol = str_remove(GeneSymbol, "__.*$")) %>%
      paste(collapse = " - ")
})

# Set up UI
if(dataset$isSiteQuant) {
   shinyjs::show("peps_by_uniprot")
}

# Filter tables
byUniprot <- reactive(
   dataset$isSiteQuant &&
   !is.null(input$peps_by_uniprot) &&
   input$peps_by_uniprot
)

pepIndex <- reactive({
   info <- req(dataset$peptides$info)
   index <- rep(TRUE, nrow(info))

   if(length(unique_classes) > 1) {
      plex <- input$peps_plex_choice
      req(!is.null(plex))

      if(plex != "all") {
         index <- info$Class == plex
      }
   }

   if(byUniprot()) {
      currentUniprotID <- activeUniprotID()
      index[index] <- info[index, "UniprotID"] == currentUniprotID
   }else {
      currentMosaicID <- values$activeMosaicID
      index[index] <- info[index, "MosaicID"] == currentMosaicID
   }
   return(index)
})

pepInfo <- reactive({
   index <- req(pepIndex())
   dataset$peptides$info[index, ]
})

pepData <- reactive({
   pepInfo <- req(pepInfo())
   pepData <- list()
   for(class in unique_classes) {
      classInfo <- pepInfo[pepInfo$Class == class, ]
      classData <- values$normPepCols[[class]][classInfo$data_row, , drop = FALSE]
      pepData[[class]] <- classData
   }
   return(pepData)
})

normPepData <- reactive({
   norm <- input$normalize_peps
   req(!is.null(norm))

   pepData <- req(pepData())
   if(norm){
      for(class in names(pepData)) {
         classData <- pepData[[class]]
         pepData[[class]] <- classData / rowSums(classData) * 100
      }
   }
   return(pepData)
})

protSequence <- reactive({
   getProteinSequence(dataset, activeUniprotID())
})

# Peptide selection
observe({
   row <- event_data("plotly_click", priority = "event", source = "peptides")$customdata
   if(!is.null(row))
      peptides$row <- row
})

observe({
   row <- input$table_peptides_rows_selected
   if(!is.null(row))
      peptides$row <- row
})

observe({
   pepInfo <- req(pepInfo())
   # Default to highest S/N peptide
   peptides$row <- which.max(pepInfo$SumSN)
})

# Residue selection
observeEvent(input$residueGroupA, {
   A <- input$residueGroupA
   B <- input$residueGroupB
   if(any(A %in% B)){
      updateSelectInput(session, "residueGroupB", selected = setdiff(B, A))
   }
})

observeEvent(input$residueGroupB, {
   A <- input$residueGroupA
   B <- input$residueGroupB
   if(any(B %in% A)){
      updateSelectInput(session, "residueGroupA", selected = setdiff(A, B))
   }
})

output$plot_peptides <- renderPlotly({
   shiny::validate(need(dataset$peptides, "Sorry, no peptide data provided."))

   pepInfo <- pepInfo()
   pepData <- normPepData()

   # We need at least one row to display:
   numPeps <- nrow(pepInfo)
   shiny::validate(need(numPeps > 0, "Sorry, no peptide data available for the requested class."))

   # Clean up sequences
   pep_names <- pepInfo$PeptideSequence %>%
      str_replace("^-", "DASH") %>%
      str_replace_all("[*-]", "_")

   # Labels
   title <- dataset$proteins$info[[values$activeMosaicID, "GeneSymbol"]]
   if(input$peps_by_uniprot) {
      title <- str_remove(title, "__.*$")
   }

   if(input$normalize_peps){
      y_type <- "RA"
   }else{
      y_type <- "S/N"
   }

   title <- sprintf("%s using %d %s", title, numPeps,
      ifelse(numPeps == 1,
         "peptide",
         "peptides"))

   # Reformat data
   x_offset <- 0
   long_data <- data.frame()
   classIDs <- dataset$proteins$classIDs
   classes <- dataset$proteins$classes[classIDs]
   for(class in unique_classes){
      classRows <- pepInfo$Class == class
      if(!any(classRows)) next

      classData <- pepData[[class]]
      if(dataset$classesMap) {
         colIDs <- match(colnames(classData), dataset$names)
         Names <- theme$names[colIDs]
         Colors <- theme$colors[colIDs]
      }else {
         Names <- paste(class, colnames(classData), sep = "~")
         Colors <- dataset$pepColors
      }
      x <- x_offset + seq(ncol(classData))
      x_offset <- x_offset + ncol(classData) + 0.5

      for(row in which(classRows)) {
         dataRow <- pepInfo[[row, "data_row"]]
         class_data <- data.frame(
            sn = unlist(classData[dataRow, ]),
            x = x,
            name = Names,
            peptide = pep_names[[row]],
            peptideId = pepInfo[[row, "PeptideId"]],
            SumSN = round(pepInfo[[row, "SumSN"]], 0),
            color = Colors,
            rowID = row)
         long_data <- rbind(long_data, class_data)
      }
   }

   min_y <- ifelse(input$peptide_fit_y, NA, 0)
   g <- ggplot(
      data = long_data,
      aes(
         x = x,
         y = sn,
         group = rowID)) +
      geom_line(
         color = "black",
         linewidth = .2,
         na.rm = TRUE,
         show.legend = FALSE) +
      geom_point(
         aes(
            fill = I(color),
            text = paste0(peptide, "\n", y_type, ": ", sn, "\nPeptide Id: ", peptideId, "\nSumSN: ", SumSN, "\nClick to select peptide"),
            customdata = rowID
         ),
         color = "black",
         shape = 21,
         size = 2,
         na.rm = TRUE,
         show.legend = FALSE) %>% suppressWarnings() +
      scale_x_discrete(limits = unique(long_data$x), labels = unique(long_data$name)) %>% suppressWarnings() +
      scale_y_continuous(limits = c(min_y, NA)) +
      ggtitle(title) +
      xlab(NULL) +
      ylab(paste("TMT", y_type))

   gp <- ggplotly(g,
         source = "peptides",
         tooltip = "text") %>%
      config(displayModeBar = FALSE) %>%
      layout(
         xaxis = list(fixedrange = TRUE),
         yaxis = list(fixedrange = TRUE),
         showlegend = FALSE) %>%
      event_register("plotly_click")

   return(gp)
})

output$plot_selected_peptide <- renderPlotly({
   row <- peptides$row
   shiny::validate(need(row, "Select a row below to view the individual peptides."))

   rowInfo <- pepInfo()[row, ]
   class <- rowInfo$Class

   heights <- pepData()[[class]][rowInfo$data_row, , drop = TRUE]
   if(dataset$classesMap){
      index <- dataset$classes == class
      names <- theme$names[index]
      colors <- theme$colors[index]
   }else {
      names <- names(heights)
      colors <- dataset$pepColors
   }

   # Plot
   g <- tmtBarchart(heights, colors, names, groups = NULL, title = pepInfo()[[row, "PeptideSequence"]], yLabel = "TMT S/N")

   gp <- ggplotly(g, tooltip = "text") %>%
      config(displayModeBar = FALSE) %>%
      layout(
         xaxis = list(fixedrange = TRUE),
         yaxis = list(fixedrange = TRUE)
      )

   return(gp)
})

colorResidues <- function(seq, resGroups, groupColors, highlight, concat = TRUE){
   numGroups <- length(resGroups)
   if(numGroups != length(groupColors) || numGroups == 0)
      return(seq)

   colorMap <- c()
   for(i in 1:numGroups){
      group <- resGroups[[i]]
      if(!is.null(group)){
         uniqueGroup <- setdiff(group, names(colorMap))# Remove possible duplicate residues
         colorMap[uniqueGroup] <- groupColors[i]
      }
   }
   seq <- strsplit(seq, "")
   class <- ""
   key <- "color"
   if(highlight) {
      class <- "class='seq-highlight'"
      key <- "border-color"
   }
   tagOpen <- sprintf("<b><span %s style=\'%s:", class, key)
   sapply(seq, function(pep){
      index <- pep %in% names(colorMap)
      pep[index] <- paste0(tagOpen, colorMap[pep[index]], "\';>", pep[index], "</span></b>")
      if(concat){
         pep <- paste0(pep, collapse = "")
      }
      return(pep)
   })
}

output$table_peptides <- renderDataTable({
   shiny::validate(need(dataset$peptides, paste0("Sorry, the peptide information for this ", dataset$idLabel, "Quant ID is unavailable.")))

   sequence <- protSequence()

   pepInfo <- pepInfo()
   numPeps <- nrow(pepInfo)
   req(numPeps > 0)

   if (!is.null(sequence)) {
      pepSeq <- pepInfo$PeptideSequence
      if(dataset$isSiteQuant){
         pepSeq <- str_remove(pepSeq, "#")
      }
      locations <- str_locate(sequence, 
                     str_extract(pepSeq, "(?<=\\.).*(?=\\.)"))
      pepInfo$Start <- locations[, "start"]
      pepInfo$End   <- locations[, "end"]
   }

   pepInfo$SumSN <- round(pepInfo$SumSN, 0)

   # Add links and set column names
   classes <- unique(pepInfo$Class)
   colors <- rep_len(suppressWarnings(brewer.pal(length(classes), "Dark2")), length(classes))
   names(colors) <- classes
   pepInfo <- pepInfo %>% mutate(
      Link = paste0('<a href="https://', dataset$server, '/gfy/www/modules//sdig/index.php?search_id=', SearchId,
         '&peptide_id=', PeptideId,
         '&scanf=', FirstScanNumber, '" ',
         'target="_blank" ',
         'style="color:', colors[Class],
         '">Link</a>')) %>%
      rename("Scan#" = FirstScanNumber, "File" = RunLoadPath)

   # Color Residues
   pepInfo$PeptideSequence <- colorResidues(pepInfo$PeptideSequence, list(input$residueGroupA, input$residueGroupB),
      highlight = FALSE, list(input$residueA_col, input$residueB_col))

   # Truncate path
   pepInfo$File <- str_extract(pepInfo$File, "([^/]+$)")

   # Drop extra cols
   cols <- c("PeptideSequence", "Link", "Class", "SumSN", "Start", "End", "File", "Scan#")

   # Show DT:
   d <- datatable(pepInfo[cols[cols %in% names(pepInfo)]],
      options = list(
         initComplete = datatableTheme(),
         order = list(which(cols == "Start") - 1, "asc"),
         searchHighlight = TRUE,
         lengthMenu = c(10, 25, 50, 100, numPeps), pageLength = 50,
         autoWidth = TRUE),
      class = "compact cell-border",
      rownames = FALSE,
      selection = list(mode = 'single'),
      escape = FALSE)
   return(d)
})

output$ui_peps_choice <- renderUI({
   req(dataset$peptides)
   req(length(unique_classes) > 1)
   req(protSequence())
   choices <- unique_classes

   numpeps <- getNumPepsFromID(dataset$proteins$info, values$activeMosaicID)
   numpeps <- numpeps[order(names(numpeps))]

   names(choices) <- paste0(choices, " - (", numpeps, " peptides)")
   choices <- choices[numpeps != 0]

   choices[[paste0("All - (", sum(unlist(numpeps)), " peptides)")]] <- "all"

   selectInput("peps_plex_choice", "Class", choices = choices, selected = "all")
})

getIndexStr <- function(i) {
   if (i < 10) {
      return(paste0(i, "--|"))
   }
   else if (i < 100) {
      return(paste0(i, "-|"))
   }
   else if (i < 1000) {
      return(paste0(i, "|"))
   }
   return(i)
}

getIndexStrTotal <- function(str) {
   if (nchar(str) > 1000 && is.null(peptides$row)) {
      return(paste0(paste(sapply(seq(1, 1000, 50), getIndexStr), collapse = "<br>")))
   } else {
      return(paste0(paste(sapply(seq(1, nchar(str), 50), getIndexStr), collapse = "<br>")))
   }
}

output$current_peptide_index <- renderUI({
   req(dataset$peptides)
   seq <- protSequence()
   req(seq)
   return(HTML(getIndexStrTotal(seq)))
})

output$current_peptide_sequence <- renderUI({
   shiny::validate(need(protSequence(), "Sorry, sequence for this UniprotID not found..."))

   pepSequences <- pepInfo()$PeptideSequence
   req(pepSequences)

   highlight <- input$highlight_residues
   req(!is.null(highlight))

   row <- peptides$row
   prot_seq <- protSequence()
   seq_sep <- colorResidues(prot_seq, list(input$residueGroupA, input$residueGroupB), list(input$residueA_col, input$residueB_col),
      highlight, concat = FALSE)[, 1]

   if(dataset$isSiteQuant){
      pepSequences <- str_remove_all(pepSequences, "#")
   }

   tagSeq <- function(seq_sep, pepseq, tagStart, tagEnd){
      pattern <- str_extract(pepseq, "(?<=\\.).*(?=\\.)")
      location <- str_locate(prot_seq, pattern)
      found <- !is.na(location[, "start"])
      if(any(found)){
         start <- location[found, "start"]
         end   <- location[found, "end"]
         num_seq <- sum(found)
         if(num_seq > 1){
            # Merge overlapping intervals
            acs_order <- order(start)
            start <- start[acs_order]
            end <- end[acs_order]
            i <- 1
            while(i < num_seq){
               if(start[[i + 1]] <= end[[i]] + 1 && end[[i + 1]] >= end[[i]]){
                  end[[i]] <- end[[i + 1]]
                  start <- start[-(i + 1)]
                  end <- end[-(i + 1)]
                  num_seq <- num_seq - 1
               }else{
                  i <- i + 1
               }
            }
         }
         
         seq_sep[start] <- paste0(tagStart, seq_sep[start])
         seq_sep[end] <- paste0(seq_sep[end], tagEnd)
      }
      return(seq_sep)
   }
   
   # Underline
   seq_sep <- seq_sep %>% tagSeq(unique(pepSequences), "<u>", "</u>")

   # Highlight
   if (!is.null(row)) {
      seq_sep <- seq_sep %>% tagSeq(pepSequences[row], "<span style=\"color:red\">", "</span>")
      # Crop:
      #if(length(seq_sep > 1000))
      #{
      #  center <- floor((start + end)/2)
      #  window <- c(center - 500, center + 500)
      #  begin_str <- "..."
      #  end_str <- "..."
      #  if(window[1] < 1)
      #  {
      #    window[1] = 1
      #    begin_str <- NULL
      #  }
      #  if(window[2] > length(seq_sep))
      #  {
      #    window[2] = length(seq_sep)
      #    end_str <- NULL
      #  }
      #  seq_sep <- c(begin_str, seq_sep[window[1]:window[2]], end_str)
      #}
   } else {
      # Crop:
      if (length(seq_sep) > 1000) {
         seq_sep <- c(seq_sep[1:998], "...")
      }
   }
   # Add breaks:
   index <- seq(1, length(seq_sep), by = 10)
   breaks <- rep("  ", length(index))
   breaks[seq(6, length(breaks), by = 5)] <- "<br>"
   seq_sep[index] <- paste0(breaks, seq_sep[index])

   out <- paste(seq_sep, collapse = "")
   return(HTML(out))
})
