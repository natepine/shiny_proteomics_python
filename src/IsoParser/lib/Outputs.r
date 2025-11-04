# Reusable output templates

#' Creates a heatmap with a dendrogram on the y-axis.
#'
#' @param logFC A dataframe with log2 foldchange data. Plot labels will be copied from dimnames(logFC)
#' @param colors Either c(lowColor, highColor) or a color vector mapped to the columns of logFC
#' @param fcCutoff Clamps logFC to c(-1, 1) * log2(fcCutoff)
#' @param isDotPlot Plot points rather than tiles
#' @param colNames Optionally override colnames(logFC) (allows duplicates)
#' @param highlightedRows A subset of rownames(logFC) to be bolded and colored blue
#' @depends ggplot2, ggdendro
dendroHeatmap <- function(logFC, colors, fcCutoff, isDotPlot, colNames = NULL, highlightedRows = NULL){
   halfWidth <- 0.5
   numProteins <- nrow(logFC)
   numColumns <- ncol(logFC)
   colorScaled <- length(colors) == 2

   if(!colorScaled && length(colors) != numColumns){
      warning("Expected 2 colors or a color for each column, unable to plot dendro-heatmap")
      return(NULL)
   }
   if(numColumns != length(colNames)) {
      warning("Column names do not match number of columns, unable to plot dendro-heatmap")
      return(NULL)
   }
   
   if (numProteins > 1) {
      dendroData <- logFC
      dendroData[is.na(dendroData)] <- 0 # CHANGE THIS IF YOU DON'T WANT TO JUST ADD 0's...

      dend <- as.dendrogram(hclust(dist(dendroData)))
      dend_data <- dendro_data(dend)

      # Flip coordinates to plot a vertical dendrogram
      segment_data <- with(segment(dend_data),
         data.frame(
            x = y + numColumns + halfWidth,
            y = x,
            xend = yend + numColumns + halfWidth,
            yend = xend
         ))

      dendroMatch <- match(rownames(logFC), dend_data$labels$label)
      y_centers <- dend_data$labels$x[dendroMatch]
      
   } else {
      y_centers <- 1 # just 1 row
   }

   if(is.null(colNames)) {
      colNames <- colnames(logFC)
   }
   heatmap_data <- data.frame(
      logFC     = as.vector(t(logFC)),
      x_centers = 1:numColumns,
      y_centers = rep(y_centers, each = numColumns)
   )

   if(!colorScaled){
      inverseColors <- getInverseColors(colors)
      heatmap_data$COLOR <- ifelse(heatmap_data$logFC >= 0, colors, inverseColors)
   }

   # Apply FC:
   log2Cutoff <- log2(fcCutoff)
   heatmap_data$logFC[heatmap_data$logFC > log2Cutoff] <- log2Cutoff
   heatmap_data$logFC[heatmap_data$logFC < -log2Cutoff] <- -log2Cutoff

   # Limits for the vertical axes
   padding <- 0.1
   yRange <- range(heatmap_data$y_centers) + (halfWidth + padding) * c(-1, 1)
   # Limits for color and alpha scales
   limits <- c(-log2Cutoff, log2Cutoff)
   breaks <- round(c(-log2Cutoff, 0, log2Cutoff), 2)

   plot <- ggplot(heatmap_data,
      aes(
         x = x_centers,
         y = y_centers
      ))

   # Choose graph type:
   if(colorScaled){
      if(isDotPlot){
         plot <- plot +
            geom_point(
               aes(size = abs(logFC), fill = logFC),
               color = "black",
               shape = 21)
      }else{
         plot <- plot +
            geom_tile(
               aes(height = 1, fill = logFC),
               color = "black",
               size = 1)
      }

      plot <- plot +
         scale_fill_gradient2(
            low = colors[[1]],
            high = colors[[2]],
            mid = "black",
            na.value = "gray",
            limits = limits,
            breaks = breaks,
            guide = guide_colorbar(
               title = "log2 ( value / rowmean )",
               barwidth = 10,
               title.position = "top",
               label.position = "bottom"))
   }else{
      # length(colors) is equal to numColumns
      if(isDotPlot){
         plot <- plot +
            geom_point(
               aes(size = abs(logFC)),
               fill = heatmap_data$COLOR,
               color = heatmap_data$COLOR,
               shape = ifelse(heatmap_data$logFC >= 0, 24, 25))
      }else{
         plot <- plot +
            geom_tile(
               aes(height = 1),
               fill = "black",
               color = "black",
               size = 1) +
            geom_tile(
               aes(height = 1, alpha = logFC),
               fill = heatmap_data$COLOR,
               color = "black",
               size = 1) +
            scale_alpha_continuous(
               limits = limits,
               breaks = breaks,
               guide = guide_legend(title = "log2 ( value / rowmean )", barwidth = 10, title.position = "top", label.position = "bottom")
            )
      }
   }

   # Add dendrogram
   if(numProteins > 1){
      plot <- plot +
         geom_segment(data = segment_data,
            aes(
               x = x,
               y = y,
               xend = xend,
               yend = yend),
            color = "black")
   }
   
   # Highlight rows
   highlight <- rownames(logFC) %in% highlightedRows

   plot +
      coord_equal() +
      guides(size = "none") +
      scale_x_continuous(
         breaks = 1:numColumns,
         labels = colNames,
         expand = c(0, 0),
         limits = c(0, NA)) +
      scale_y_continuous(
         breaks = y_centers,
         labels = rownames(logFC),
         limits = yRange,
         expand = c(0, 0)) +
      labs(x = NULL, y = NULL) +
      theme_void() +
      theme(
         legend.position = "top",
         text = element_text(family = "serif"),
         axis.text.x = element_text(hjust = 1, vjust = 0.19, angle = 90),
         axis.text.y = element_text(hjust = 1, vjust = 0.19,
            color = ifelse(highlight, "blue", "black"),
            face  = ifelse(highlight, "bold", "plain")),
         plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"), # margin: top, right, bottom, and left
         legend.title.align = 0.5)
}

#' Creates a sharedSelectize pool
#' 
#' @param choices All valid inputs. (NULL, "", FALSE, length 0 values, NA) will be removed. See shiny::isTruthy
#' @param selected A reactive value to update the inputs when changed ex: reactiveVal() or reactive({...})
#' @param msg Displayed when an input is invalid
#' @param callback A function given then current input when a new selection is made (!= history[1])
#' @param historyLength Max history length before new entries cause old ones to be discarded
sharedSelectize <- function(choices, selected, msg = NULL, callback = NULL, historyLength = 10){
   validateChoices <- function(choices){
      truthy <- sapply(choices, isTruthy)
      if(any(!truthy)){
         choices <- choices[truthy]
         warning("Some sharedSelectize values are invalid. They have been removed.")
      }
      return(choices)
   }

   choices <- validateChoices(choices)

   # Shared function variables
   sortedChoices <- sort(choices)
   history <- c(isolate(selected()))
   input <- getDefaultReactiveDomain()$input

   getChoices <- function(){
      if(length(history) > 1){
         setNames(
            list(history, setdiff(sortedChoices, history)),
            c("History", "Symbol"))
      }else{
         list(Symbol = sortedChoices)
      }
   }

   #' Creates a selectize input with shared history, listens to all
   #' inputs in the pool and adds the new value to the history list
   #'
   #' @param id The input id for a new shared selectizeInput
   #' @param ... Passed on selectizeInput
   addInput <- function(id, ...){
      observe({
         # Inputs seem to be set to "" momentarily after each update when server = TRUE, ignoring with req
         choice <- req(input[[id]])
         isNew <- choice != history[1]
         added <- addToHistory(choice, useCallback = isNew)
         if(!added && !is.null(msg)){
            sendSweetAlert(
               title = "Not Found...",
               text = msg,
               type = "error")
         }
      })
      observe({
         updateSelectizeInput(inputId = id, choices = getChoices(), selected = selected(), server = TRUE)
      })
      renderUI({
         isolate({ #Initialize once
            updateSelectizeInput(inputId = id, choices = getChoices(), selected = selected(), server = TRUE)
            selectizeInput(id, choices = NULL, selected = NULL, ...)
         })
      })
   }

   #' Add to the input's history
   #'
   #' @param choice The value to be added if valid
   #' @param useCallback Should the callback function be called?
   #' @return Wheather choice was successfully added to the history
   addToHistory <- function(choice, useCallback = TRUE){
      if(choice %in% choices){
         # Remove any duplicates
         history <- history[history != choice]
         history <- c(choice, history)
         # Limit history to 10 entries
         history <<- history[1:min(length(history), historyLength)]
         if(!is.null(callback) && useCallback){
            isolate(callback(choice))
         }
         return(TRUE)
      }
      return(FALSE)
   }
   
   list(
      addInput = addInput,
      addToHistory = addToHistory,
      class = "sharedSelectizePool"
   )
}

#' Go calculator generator
#'
#' @param allGO A GO dataframe of unique rows with columns {"UniprotID", "Annotation", "Database", "DatabaseID"}
#' @param idMap A dataframe of all observations with columns {"uniqueID", "UniprotID", "displayName"}
#'                where uniqueID is required if the UniprotID column contains duplicates
#'                and displayName is an optional column for renaming table tooltips
#' @param urls If not NULL, a named vector of format strings passed to sprintf with one %s to be supplied a DatabaseID
#'                Names will be used as regex to select from unique values of Database in @allGO, first match will be used
#' @return An object of class goGenerator or NULL if data is unusable

goEnrichment <- function(allGO, idMap, urls = NULL) {
   if(is.null(allGO) || is.null(idMap)) {
      warning("Unable to create goGenerator: NULL parameter(s)")
      return(NULL)
   }

   if(is.null(idMap$uniqueID)) {
      idMap$uniqueID <- idMap$UniprotID
   }

   if(any(duplicated(idMap$uniqueID))) {
      warning("Unable to create goGenerator: duplicated IDs detected in idMap")
      return(NULL)
   }

   missingCols <- setdiff(c("UniprotID", "Annotation", "Database"), names(allGO))
   if(length(missingCols) != 0) {
      warning("Unable to create goGenerator: missing required columns ", paste(missingCols, sep =  ", "))
      return(NULL)
   }

   rownames(idMap) <- idMap$uniqueID
   allGO <- allGO %>% filter(UniprotID %in% idMap$UniprotID)

   if(nrow(allGO) == 0){
      warning("Unable to create goGenerator, no matching UniprotIDs.")
      return(NULL)
   }

   if(!"DatabaseID" %in% names(allGO))
      allGO$DatabaseID <- allGO$Annotation # Use annotation as a backup unique id

   #Count observed and annotated ids
   databases <- do.call(rbind, lapply(unique(allGO$Database), function(database) {
      urlFmt <- NA
      if(!is.null(urls)) {
         urlInd <- str_which(database, names(urls))
         if(length(urlInd) != 0)
            urlFmt <- urls[[urlInd[[1]]]]
      }
      data.frame(
         numObs = sum(idMap$UniprotID %in% filter(allGO, Database == database)$UniprotID),
         urlFmt = urlFmt,
         row.names = database
      )
   }))

   #Count observations per Annotation
   idCounts <- table(idMap$UniprotID)
   if(all(idCounts==1)) {
      goCounts <- allGO %>%
         count(DatabaseID, Database, Annotation) %>%
         rename(AnnotCount = n)
   } else {
      goCounts <- allGO %>%
         mutate(n = idCounts[UniprotID]) %>%
         count(DatabaseID, Database, Annotation, wt = n) %>%
         rename(AnnotCount = n)
   }

   #' A GO calculator (1 per UI) to store GO data and build GO ui elements
   #'
   #' @return An objecto of class goCalculator
   goCalculator <- function(){
      # Internal variables
      subsetIDs <- NULL
      goSubset <- NULL
      colors <- brewer.pal(n = 9, "Set1")
      colors <- c(colors[1:6], "#808000", colors[7:9])

      # External variables
      sigGO <- NULL

      #' Filter the GO data by uniqueID for further analysis
      #'
      #' @param uniqueIDs IDs used to filter GO data
      #' @return Whether any GO data was found for @uniprotIDs
      subsetData <- function(uniqueIDs){
         subsetIDs <<- uniqueIDs
         ids <- idMap[subsetIDs, "UniprotID"]
         goSubset <<- allGO %>%
            filter(UniprotID %in% ids) %>%
            mutate(n = table(ids)[UniprotID])

         if(nrow(goSubset) == 0){
            goSubset <<- NULL
            sigGO <<- NULL
            return(FALSE)
         }

         return(TRUE)
      }

      #' Calculate p-values and determine significance of last subsetData()
      #'
      #' @param minMatches Cutoff for number of GO matches
      #' @param pvalCutoff Cutoff for a significant p-value
      #' @param p.adj      Method to be passed to p.adjust
      #' @return Significant annotations
      filterGO <- function(minMatches, pvalCutoff, p.adj){
         if(is.null(goSubset)){
            return(NULL)
         }

         # Calculate probabilities:
         sigGO <<- goSubset %>% count(DatabaseID, Database, wt = n) %>%
            left_join(goCounts, by = c("DatabaseID", "Database")) %>%
            transmute(
               Database,
               DatabaseID,
               Annotation,
               AnnotCount,
               TotalObs = databases[Database, "numObs"],
               SubsetHits = n,
               SubsetSize = length(subsetIDs)) %>%
            mutate(PVal = phyper(SubsetHits, AnnotCount, TotalObs - AnnotCount, SubsetSize, lower.tail = FALSE)) %>%
            mutate(pCorr = p.adjust(PVal, method = p.adj)) %>%
            filter(SubsetHits >= minMatches, pCorr <= pvalCutoff) %>%
            distinct(DatabaseID, .keep_all = TRUE) %>%
            arrange(pCorr)

         if(nrow(sigGO) == 0) sigGO <<- NULL

         return(sigGO)
      }

      #' Generate tooltips based on most recent subsetData() and filterGO() results
      #'
      #' @param maxTooltips Max |s used in protTooltips() before truncating
      #' @return A list of tooltips to be used in a table named by UniprotID
      protTooltips <- function(maxTooltips = 60){
         if(is.null(sigGO) || is.null(goSubset)){
            return(NULL)
         }

         # Limit to maxTooltips params
         overflow <- nrow(sigGO) > maxTooltips
         nTooltips <- min(nrow(sigGO), maxTooltips)
         nAnnots <- sigGO[1:nTooltips, c("DatabaseID", "Annotation")]

         sapply(idMap[subsetIDs, "UniprotID"], function(id) {
            dbIDs <- filter(goSubset, UniprotID == id)$DatabaseID
            exists <- nAnnots$DatabaseID %in% dbIDs
            rowColors <- ifelse(exists, colors, "lightgray")
            tooltips <- rep("", length(exists))
            tooltips[exists] <- sprintf("data-toggle='tooltip' title='%s'", html_escape(nAnnots[exists, "Annotation"], "'"))

            x <- paste0(
               "<span class=\'GO-bar\' style=\'",
                  "color:", rowColors, ";",
                  "background-color:", rowColors, ";",
               "\'",
               tooltips,
               ">|</span>")
            # Indicate output has been truncated
            if(overflow) x <- append(x, "<span>...</span>")
            return(paste(x, collapse = ""))
         })
      }

      annotationList <- function(collapse = "<br>") {
         sapply(sigGO$DatabaseID, function(annotID){
            uniprotIDs <- filter(goSubset, DatabaseID == annotID)$UniprotID
            mappedIDs <- filter(idMap, UniprotID %in% uniprotIDs)$uniqueID
            ids <- intersect(subsetIDs, mappedIDs)

            if(!is.null(idMap$displayNames)) {
               tooltip <- idMap[ids, "displayNames"]
            } else {
               tooltip <- idMap[ids, "UniprotID"]
            }
            paste(tooltip, collapse = collapse)
         })
      }

      #' Build a datatable to display significant annotations
      goTable <- function(){
         if(is.null(sigGO) || is.null(goSubset)){
            return(NULL)
         }

         # Add tooltips
         sigGO$hits <- annotationList()

         # Add Links
         urlFmts <- databases[sigGO$Database, "urlFmt"]
         hasUrl <- !is.na(urlFmts)
         if(any(hasUrl)) sigGO[hasUrl, "url"] <- sprintf(urlFmts[hasUrl], sigGO[hasUrl, "DatabaseID"])

         # Format Table:
         sigGO$PVal <- formatC(sigGO$PVal, format = "e", digits = 2)
         sigGO$pCorr <- formatC(sigGO$pCorr, format = "e", digits = 2)

         # Render data
         colDefs <- addTooltip(sigGO, "SubsetHits", "hits", enableHTML = TRUE)
         if(any(hasUrl)) colDefs <- c(colDefs, addTableLink(sigGO, "DatabaseID"))

         datatable(sigGO,
            options = list(
               initComplete = datatableTheme(),
               searchHighlight = TRUE,
               autoWidth = FALSE,
               columnDefs = colDefs),
            class = "compact cell-border",
            rownames = FALSE,
            selection = 'none',
            escape = FALSE
         ) %>% formatStyle(
            'DatabaseID',
            target = 'row',
            backgroundColor = styleEqual(sigGO$DatabaseID,
               adjustcolor(
                  rep(colors, length.out = nrow(sigGO)),
                  alpha.f = 0.3))
         )
      }

      #' Get proteins that belong to a specific GO term
      #'
      #' @param goTerm The GO DatabaseID to search for
      #' @return Vector of UniprotIDs that belong to the specified GO term
      getProteinsForGOTerm <- function(goTerm){
         if(is.null(allGO)) return(c())
         
         proteins <- allGO %>%
            filter(DatabaseID == goTerm) %>%
            pull(UniprotID) %>%
            unique()
         
         return(proteins)
      }

      return(structure(
         list(
            subsetData = subsetData,
            filterGO = filterGO,
            getSigGO = function(){
               data <- sigGO
               data$Hits <- annotationList("|")
               return(data)
            },
            goTable = goTable,
            protTooltips = protTooltips,
            getProteinsForGOTerm = getProteinsForGOTerm
         ),
         class = "goCalculator"
      ))
   }

   return(structure(
      list(
         goCalculator = goCalculator,
         goCounts = goCounts
      ),
      class = "goGenerator"
   ))
}

# Helper functions to generate columnDefs for DT::datatable

#' Convert a column to hyperlinks by appending table data to a url
#'
#' @param url     Base URL for each link
#' @param targets Which columns to convert to links
#' @param postfix URL string added after table data
#' @return        A list containing a columnDef entry

addLink <- function(url, targets, postfix = "", style) {
   styleAttr <- ifelse(missing(style), "",
      paste0("style=\"", paste(names(style), style, sep = ":", collapse = " "),"\""))
   list(list(
      targets = targets,
      render = JS(sprintf(
         "function(data, type) {
            if(type === 'display') {
               return '<a href=\"%s' + data + '%s\" target=\"_blank\" %s>' + data + '</a>';
            }
            return data;
         }", url, postfix, styleAttr))
   ))
}

#' Add custom button column
#'
#' @param targets    Which columns to add tooltips to
#' @param inputID    Shiny input id to set to table data
#' @return           A list containing columnDef entries

addShinyButton <- function(targets, inputID) {
   list(list(
      targets = targets,
      searchable = FALSE,
      render = JS(sprintf(
         "function(data, type) {
            if(type === 'display') {
               return \"<button id='a' type='button' class='btn btn-default action-button updatebttn' \" +
                  \"onclick='Shiny.setInputValue(&quot;%s&quot;, &quot;\" + data + \"&quot;);'></button>\";
            }
            return data;
         }", inputID))
   ))
}

#' Convert colunm names to 0 indexed column positions
#'
#' @param df   A data.frame with column names
#' @param c    Column names to convert
#' @return     Vector of matched column indexes

jsCol <- function(df, c) {
   which(names(df) %in% c) - 1
}

#' Add a class to columns
#'
#' @param data    Data.frame containing table data
#' @param targets Which columns to add @classes to
#' @param classes Vector of css classes to apply
#' @return        A list containing columnDef entries

addColumnClasses <- function(targets, classes) {
   list(list(
      targets = targets,
      className = paste(classes, collapse = " ")
   ))
}

#' Convert a column to scientific notation
#'
#' @param targets Which columns to display in scientific notation
#' @return        A list containing columnDef entries

addExpFormat <- function(targets, digits = 2) {
   defs <- list(list(
      targets = targets,
      render = JS(sprintf(
         "function(data, type) {
            if(typeof data !== 'number') {
               data = parseFloat(data);
            }
            if(type === 'display') {
               return data.toExponential(%d);
            }
            return data;
         }", digits))
   ))
   return(defs)
}

#' Convert a column to hyperlinks by appending table data to an optionally hidden column
#'
#' @param data    Data.frame containing table data
#' @param targets Which columns to convert to links
#' @param urlCol  Column name or index containing the base URL for each link or NA for no url, default 'url'
#' @param hidden  Should @urlCol be hidden, default TRUE
#' @return        A list containing columnDef entries

addTableLink <- function(data, targets, urlCol = 'url', hidden = TRUE) {
   if(!is.numeric(urlCol))
      urlCol <- jsCol(data, urlCol)

   defs <- list(list(
      targets = targets,
      render = JS(sprintf(
         "function(data, type, row) {
            if(type === 'display') {
               url = row[%d];
               if(url === null) return data;
               return '<a href=\"' + url + '\" target=\"_blank\">' + data + '</a>';
            }
            return data;
         }", urlCol))
   ))
   if(hidden)
      defs <- c(defs, list(list(targets = urlCol, visible = FALSE, searchable = FALSE)))
   return(defs)
}

#' Add tooltips displaying table data to column on cell hover
#'
#' @param data       Data.frame containing table data
#' @param targets    Which columns to add tooltips to
#' @param tooltipCol Column name or index containing the tooltip contents, default 'tooltip'
#' @param enableHTML Render tooltip contents as HTML, default FALSE
#' @param hidden     Should @tooltipCol be hidden, default TRUE
#' @return           A list containing columnDef entries

addTooltip <- function(data, targets, tooltipCol = 'tooltip', enableHTML = FALSE, hidden = TRUE) {
   if(!is.numeric(tooltipCol))
      tooltipCol <- jsCol(data, tooltipCol)

   defs <- list(list(
      targets = targets,
      render = JS(sprintf(
         "function(data, type, row) {
            if(type === 'display') {
               return '<div class=outline-seq data-toggle=tooltip data-html=%s title=\"' + row[%d] + '\">' + data + '</div>';
            }
            return data;
         }", ifelse(enableHTML, "true", "false"), tooltipCol))
   ))
   if(hidden)
      defs <- c(defs, list(list(targets = tooltipCol, visible = FALSE, searchable = FALSE)))
   return(defs)
}

#' Add visual bars at each mark in a cell displaying the total
#'
#' @param data       Data.frame containing table data
#' @param targets    Which columns to add tooltips to
#' @param marks      Column name or index containing the marks
#' @param title      Optional column name or index containing tooltip content
#' @param enableHTML Render tooltip contents as HTML, default FALSE
#' @param sep        String to split @marks column on for multiple marks, default ";"
#' @param default    Value to display for missing data, default "Unknown"
#' @param hidden     Character vector containing column names to be hidden or NULL to show all, default c(title, marks)
#' @return           A list containing columnDef entries

addPosMarker <- function(data, targets, marks, title = NULL, enableHTML = FALSE, sep = ";", default = "Unknown", hidden = c(title, marks)) {
   if(!is.numeric(marks)) {
      marks <- jsCol(data, marks)
   }
   tagAttr <- "''"
   if(!is.null(title)) {
      if(!is.numeric(title)) {
         title <- jsCol(data, title)
      }
      tagAttr <- sprintf("(title = row[%d]) === null ? '' : ' data-toggle=tooltip data-html=%s title=\"' + title + '\"'", title,
         ifelse(enableHTML, "true", "false"))
   }
   defs <- list(list(
      targets = targets,
      defaultContent = default,
      createdCell = JS("(td) => $(td).addClass('seq-cell')"),
      render = JS(sprintf("function(data, type, row) {
         if(data !== null && type === 'display') {
            marks = row[%d];
            mark_divs = marks.split('%s').map(
               (mark) => '<div class=seq-site style=left:' + parseInt(mark) / data * 100 + '%%></div>'
            ).join('');
            attr = %s;
            return mark_divs + '<div class=outline-seq' + attr + '>' + data + '</div>';
         }
         return data;
      }", marks, sep, tagAttr))
   ))

   if(!is.null(hidden)) {
      if(is.character(hidden)) {
         hidden <- intersect(names(data), hidden)
      }
      defs <- c(defs, list(list(targets = hidden, visible = FALSE, searchable = FALSE)))
   }
   return(defs)
}

#' Truncate both sides of sequences around peptides
#'
#' @param seqs    Full sequences
#' @param peps    Sequences to match, '.' and '#' will be removed
#' @param width   Total number of positions to show, default 150
#' @param elipsis String to indicate truncation, default '...'
#' @param tags    Vector of length 2 to paste around pep or NULL, default is underlined: c("<u>", "</u>")
#' @param ignore  Vector of strings to remove from peps while matching to seqs or NULL, default c('.', '#')
#' @return        If match found, truncated sequence with original peptide inserted
#'                otherwise the original peptide

trunc_seqs <- function(seqs, peps, width = 150, elipsis = "...", tags = c("<u>", "</u>"), ignore = c(".", "#")) {
   result <- peps

   pep_regex <- peps
   if(!is.null(ignore)) {
      ignore <- str_vec_regex(ignore)
      pep_regex <- str_remove_all(peps, ignore)
   }
   pep_regex <- str_escape(pep_regex)

   pep_locs <- str_locate(seqs, pep_regex)
   found <- !is.na(pep_locs[, "start"])
   if(any(found)) {
      pep_locs    <- pep_locs[found, , drop = FALSE]
      peps        <- peps[found]
      seqs        <- seqs[found]
      max_chars   <- pmax(width - nchar(peps), 0) / 2

      pep_locs[, "start"] <- pep_locs[, "start"] - 1

      #pep:                           ALC
      #seq:           MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGN
      # |-left_extra-|                       >>>                        |------------|
      # |                            width                              |
      # |           left_max         |   |          right_max           |
      # | delta  |           start ->|   |<- end               |  delta |
      # (elipsis)|left_end           |pep|            right_end|(elipsis)

      left_end    <- rep(0, sum(found))
      right_end   <- nchar(seqs)
      left_max    <- floor(  max_chars)
      right_max   <- ceiling(max_chars)
      left_delta  <- left_max  -              pep_locs[, "start"]
      right_delta <- right_max - (right_end - pep_locs[, "end"])

      left_extra  <- pmax(left_delta, 0)
      right_extra <- pmax(right_delta, 0)

      left_max    <- left_max    + right_extra
      right_delta <- right_delta - right_extra + left_extra
      right_max   <- right_max   + left_extra
      left_delta  <- left_delta  - left_extra + right_extra

      too_long <- left_delta < 0
      left_end[too_long] <- pep_locs[too_long, "start"] - pmax(left_max[too_long] - nchar(elipsis), 0)
      left <- str_sub(seqs, left_end + 1, pep_locs[, "start"])
      left[too_long] <- paste0(elipsis, left[too_long])

      too_long <- right_delta < 0
      right_end[too_long] <- pep_locs[too_long, "end"] + pmax(right_max[too_long] - nchar(elipsis), 0)
      right <- str_sub(seqs, pep_locs[, "end"] + 1, right_end)
      right[too_long] <- paste0(right[too_long], elipsis)

      if(!is.null(tags))
         peps <- paste0(tags[[1]], peps, tags[[2]])

      result[found] <- paste0(left, peps, right)
   }

   return(result)
}
