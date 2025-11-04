### Get Data Fnxns
getProteinSequence <- function(dataset, UniprotIDs){
   seqs <- NULL
   annotDB <- paste0(dataPath, "annotations.db")
   if(file.exists(annotDB)) {
      db <- dbConnect(SQLite(), annotDB, SQLITE_RO)
      on.exit(dbDisconnect(db))

      table <- paste0("proteins_", dataset$taxid)
      if(dbExistsTable(db, table)) {
         seqs <- dbGetQuery(db, paste("SELECT UniprotID, Sequence FROM", table, "WHERE UniprotID IN (@ids);"), params = list(ids = UniprotIDs))
      }
   }

   if(is.null(seqs) && curl::has_internet()) {
      query <- paste0('accession_id:', unique(UniprotIDs), collapse = '+OR+')
      url <- paste0('https://rest.uniprot.org/uniprotkb/stream?compressed=false&',
         'format=fasta&query=', query)

      failed <- function(cond) {
         message(paste("Uniprot query failed:", url))
         return(NULL)
      }

      data <- tryCatch(
         httr::GET(url) %>%
            httr::content(
               as = "text",
               type = "text/x-fasta",
               encoding = "UTF-8"
            ),
         error = failed,
         warning = failed
      )

      if(!is.null(data)) {
         reg <- ">[^\\|]*."
         data <- str_remove(data, reg) %>%
            str_split(reg) %>%
            unlist

         seqs <- data.frame(
            Sequence = str_remove(data, ".*\n") %>% str_remove_all("\n"),
            UniprotID = str_extract(data, "[^\\|]*")
         )
      }
   }

   if(!is.null(seqs)) {
      seqs <- setNames(seqs$Sequence, seqs$UniprotID)[UniprotIDs]
      names(seqs) <- UniprotIDs
   }else {
      seqs <- setNames(NA, UniprotIDs)
   }

   return(seqs)
}

getInverseColors <- function(hex) {
   inverse <- 255 - col2rgb(hex)
   apply(inverse, 2, function(x)
      rgb(x[1], x[2], x[3], maxColorValue = 255)
   )
}

getFigureParam <- function(input_val, default, escape) {
   if (is.null(input_val)) {
      return(default)
   } else if (input_val == escape) {
      return(default)
   } else {
      return(input_val)
   }
}

#' @param columns - Expression data as a matrix, each protein is a row
#' @param key - A vector of expression data for subsetting, could be a row from columns.
#'              NAs will be used to drop unwanted columns from the dataset.
subsetColumns <- function(columns, key) {
   keep <- !is.na(key)

   key[keep] <- 100 * key[keep] / sum(key, na.rm = TRUE)
   columns <- 100 * columns / rowSums(columns[, keep, drop = FALSE], na.rm = TRUE)

   # Remove proteins with all NA values
   columns <- columns[rowSums(!is.na(columns[, keep, drop = FALSE])) != 0, , drop = FALSE]

   return(list(
      keep = keep,
      key = key,
      columns = columns
   ))
}

#' @param mat - Expression data as a matrix, colnames contains the protein identifiers
#' @param vec - A vector of expression data to compute distances from
#' @param method - The distance function to use
#'
#' @return A named vector of distances in descending order
getRowDistance <- function(mat, vec, method = c("euc", "corr", "cos")) {
   method <- match.arg(method)

   colVecLengths <- function(x){
      (x ^ 2) %>% colSums(na.rm = TRUE) %>% sqrt
   }
   vecLength <- function(x){
      (x ^ 2) %>% sum(na.rm = TRUE) %>% sqrt
   }

   #Ignore zero values in mat and rows with less than two values
   zeroes <- mat == 0
   if(any(zeroes, na.rm = TRUE)) {
      mat[zeroes] <- NA
      updatedCols <- which(colSums(zeroes, na.rm = TRUE) > 0)
      invalidCols <- updatedCols[colSums(!is.na(mat[,updatedCols, drop = FALSE])) < 2]
      if(any(invalidCols)){
         mat <- mat[, -invalidCols, drop = FALSE]
      }
   }

   if (method == "euc") {
      # Renormalize to ignore NAs
      wt <- nrow(mat) / colSums(!is.na(mat))
      vec <- vec %*% t(wt)
      distVec <- colVecLengths(mat - vec)
   } else if (method == "corr") {
      distVec <- 1 - cor(mat, vec)
      distVec <- setNames(as.vector(distVec), rownames(distVec))
   } else if (method == "cos") {
      dotProducts <- (mat * vec) %>% colSums(na.rm = TRUE)
      lenProducts <- colVecLengths(mat) * vecLength(vec)
      distVec <- 1 - (dotProducts / lenProducts)
   }

   return(
      sort(distVec)
   )
}

findMosaicID <- function(info, id) {
   id <- tolower(id)

   index <- FALSE
   find <- function(column) {
      index <<- id == tolower(info[[column]])
      return(any(index))
   }

   if(!find("MosaicID") && !find("GeneSymbol") && !find("UniprotID")){
      return(NULL)
   }
   
   #Return first TRUE index
   return(info$MosaicID[index][[1]])
}

getNumPepsFromID <- function(protInfo, mosaicID) {
   peps <- str_subset(colnames(protInfo), "Peptides")
   return(protInfo[mosaicID, peps])
}

checkIfDuplicatedProt <- function(recent_proteins, new_protein) {
   recent_proteins <- recent_proteins[recent_proteins != new_protein]      #remove any duplicates
   recent_proteins <- recent_proteins[1:min(length(recent_proteins), 9)]   #Limit to 10 recent proteins
   recent_proteins <- c(new_protein, recent_proteins)
   return(recent_proteins)
}

calculateLogFC <- function(dataset, columnData, mosaicIDs) {
   columnData <- columnData[mosaicIDs, , drop = FALSE]

   classIDs <- dataset$proteins$classIDs
   for(classID in unique(classIDs)) {
      classCols <- classIDs == classID
      classData <- columnData[, classCols, drop = FALSE]
      columnData[, classCols] <- classData / rowMeans(classData, na.rm = TRUE)
   }

   # Set FC of 0 to 0
   columnData <- log2(columnData) %>% as.data.frame
   columnData[columnData == log2(0)] <- 0

   # Set names:
   colnames(columnData) <- seq(ncol(columnData)) # Map indices back to theme$names (duplicates names possible)
   rownames(columnData) <- dataset$proteins$info[mosaicIDs, "GeneSymbol"]

   return(columnData)
}

writeXLSX <- function(dataset, idMaps, file, hostname) {
   wb <- createWorkbook(
      creator = paste0("TMT Mosaic ", version),
      title = dataset$name
   )
   addWorksheet(wb, "Data Summary")

   boldStyle <- createStyle(textDecoration = "bold")
   if(!is.null(dataset$peptides)){
      addWorksheet(wb, "Raw Files")
      rawFiles <- unique(dataset$peptides$info$RunLoadPath)
      rawFiles <- str_extract(rawFiles, "([^/]+$)")
      writeData(wb, "Raw Files", startRow = 1, "Raw Files:")
      addStyle(wb, "Raw Files", boldStyle, rows = 1, cols = 1)
      writeData(wb, "Raw Files", startRow = 3, rawFiles)
   }

   addWorksheet(wb, "Legend")

   getSumText <- function() {
      sumText <- "Summed S/N"
      if(dataset$normalization$normBy != "none") {
         sumText <- paste(sep = ", ", sumText,
            switch(dataset$normalization$normBy,
               "row" = paste0("Column Normalized to ", dataset$normalization$normRow),
               "all" = "Column Normalized",
               "Unknown Normalization"))
      }
      return(sumText)
   }

   if(dataset$isSiteQuant){
      format <- format_translator(list(
         GeneSymbol        = "gene_symbol",
         Site              = "Site Position",
         UniprotID         = "Protein Id",
         Description       = "prot_description",
         Motif             = "Motif",
         MaxScore          = "Max Score",
         Redundancy        = "redundancy",
         Sequence          = "sequence",
         Peptides          = list(pattern = "num_quant", regex = TRUE)
      ))

      sumText <- getSumText()
      headers <- c(
         scaled = paste(sumText, "Scaled to 100 per Site", sep = ", "),
         sum = sumText
      )

      query <- paste0("site_quant/", dataset$ID, "?all_users=1&peptide_parsimony=UR")
      colTypes <- writeSheet(wb, dataset, idMaps, "Site Quant", headers, query = query, format = format)

      format <- format_translator(list(
         GeneSymbol        = "geneSymbol",
         Site              = "sitePosStr",
         UniprotID         = "proteinID",
         Description       = "protDesc",
         Motif             = "motifPeptideStr",
         MaxScore          = "maxScoreStr",
         Redundancy        = "redundancyStr",
         Sequence          = "sequence",
         Peptides          = list(pattern = "num_quant", regex = TRUE)
      ))
      query <- paste0("composite_site_quant/", dataset$ID, "?all_users=1&peptide_parsimony=UR")
      comp_colTypes <- writeSheet(wb, dataset, idMaps, "Composite", headers, query = query, format = format)

      tabs <- data.frame(
         c("Site Quant", "Composite"),
         c("", "2 or more sites on the same peptide")
      )

      colTypes <- union(colTypes, comp_colTypes)
      if(!is.null(colTypes) && !any(colTypes == "Num Quant") && any(str_detect(colTypes, "Num Quant"))){
         colTypes <- c(colTypes, "Num Quant")
      }
   }else {
      sumText <- getSumText()
      headers <- c(
         scaled = paste(sumText, "Scaled to 100 per Protein", sep = ", "),
         sum = sumText
      )
      colTypes <- writeSheet(wb, dataset, idMaps, "Protein Quant", headers)

      tabs <- data.frame("Protein Quant", "")
   }

# Legend Tab

   rowDelta <- 2
   nextRow <- 1
   if(nrow(tabs) > 1) {
      writeData(wb, "Legend", "Tabs:")
      writeData(wb, "Legend", tabs, startCol = 2, colNames = FALSE)
      nextRow <- nrow(tabs) + rowDelta
   }

   legend <- c(
      ProteinID      = "Protein ID",
      GeneSymbol     = "Gene Symbol",
      Description    = "Description",
      Site           = "Mod site residue position #",
      Motif          = "Mod site flanked by 6 residues on either end",
      MaxScore       = "ModScore: probability of localization (>13 means >95% confidently localized)",
      Redundancy     = "Unique site or seen more than once in proteome (Redundant)",
      Sequence       = "Peptide sequence",
      `Num Quant`    = "Number of spectra quantified for this site",
      Peptides       = "Number of peptides quantified for this protein"
   )

   sumSelect <- str_detect(colTypes, "sum")
   legend_elem <- intersect(names(legend), colTypes[!sumSelect])
   writeData(wb, "Legend", "Columns:", startRow = nextRow)
   writeData(wb, "Legend", data.frame(legend_elem, legend[legend_elem]), startCol = 2, startRow = nextRow, colNames = FALSE)
   nextRow <- nextRow + length(legend_elem) + rowDelta

   addStyle(wb, "Legend", boldStyle, rows = 1:nextRow, cols = 1:2, gridExpand = TRUE)
   setColWidths(wb, "Legend", 1:3, widths = "auto")

# Data Summary Tab

   summary <- c(
      Dataset = dataset$name,
      ID = dataset$ID,
      Date = "--DATE--",
      `Number of Samples` = dataset$numSamples,
      Species = dataset$species
   )
   names(summary)[which(names(summary)=="ID")] <- paste0(dataset$idLabel, "Quant ID")
   dateCell <- c(2, which(names(summary)=="Date"))
   writeData(wb, "Data Summary", data.frame(names(summary), summary), colNames = FALSE)
   writeData(wb, "Data Summary", Sys.Date(), xy = dateCell)
   addStyle(wb, "Data Summary", createStyle(numFmt = "DATE"), cols = dateCell[[1]], rows = dateCell[[2]])
   nextRow <- length(summary) + rowDelta

   class_regex <- sprintf("^(%s)[_~]", str_vec_regex(dataset$proteins$classes))
   tags <- str_remove(colTypes[sumSelect], class_regex)
   tags <- str_remove(tags, "_sn_sum$")[dataset$order]
   writeData(wb, "Data Summary", startRow = nextRow, data.frame(
      `Sample Name` = dataset$names,
      Tag = tags,
      check.names = FALSE
   ))
   addStyle(wb, "Data Summary", boldStyle, rows = nextRow, cols = 1:2, gridExpand = TRUE)
   nextRow <- nextRow + length(dataset$names) + rowDelta

   addStyle(wb, "Data Summary", boldStyle, rows = 1:length(summary), cols = 1, gridExpand = TRUE)
   setColWidths(wb, "Data Summary", 1:2, widths = "auto")

   url <- paste0("https://", hostname, "/", config$VIEWER_SUBDIR, "/?sessid=", dataset$key)
   names(url) <- url
   class(url) <- "hyperlink"
   writeData(wb, "Data Summary", "Interactive Viewer", startRow = nextRow)
   addStyle(wb, "Data Summary", boldStyle, rows = nextRow, cols = 1)
   writeData(wb, "Data Summary", url, startRow = nextRow + 1)
   nextRow <- nextRow + 2 + rowDelta

   saveWorkbook(wb, file)
}

writeSheet <- function(wb, dataset, idMaps, sheetName, headers, query = NULL, format = NULL) {
   addWorksheet(wb, sheetName)

   cache <- createCache(dataPath)
   fmtProteins <- function(data) {
      if(is.null(data))
         return(NULL)

      data <- gfyFormatProteins(data, dataset$isSiteQuant, format = format,
         suffixes = c("scaled", "sum"), suffix_select = "in")

      if(is.null(data))
         return(NULL)

      as.mosaic_df(data, NULL, dataset$isSiteQuant, idMap = idMaps[[dataset$species]])$proteins
   }
   data <- tryCatch({
      data <- NULL

      if(is.null(query)) { #default queries may be cached
         rawData <- cache$read(dataset_cache_key(dataset, "proteins"))
         data <- fmtProteins(rawData)
      }

      if(is.null(data)) {
         if(is.null(query)) {
            rawData <- gfySourceProteins(dataset$ID, dataset$gfy.obj, dataset$isSiteQuant)
         }else {
            rawData <- getAPIData(dataset$gfy.obj, query)
         }
         data <- fmtProteins(rawData)
      }

      data
   },
      error = function(e) {
         msg <- paste("ERROR:", e$message)
         writeData(wb, sheetName, msg, startRow = 2)
         showNotification(msg, duration = NULL)
         NULL
      }
   )

   if(is.null(data))
      return(NULL)

   data$info <- data$info %>%
      rename(ProteinID = UniprotID) %>%
      select(-MosaicID)

   pepCols <- str_detect(colnames(data$info), "Peptides")

   pepOrder <- order(desc(rowSums(data$info[pepCols])))
   data$info <- data$info[pepOrder, ]
   data$columns <- data$columns[pepOrder, ]

   replacement <- ifelse(dataset$isSiteQuant, "Num Quant", "Peptides")
   if(sum(pepCols) > 1) {
      colnames(data$info)[pepCols] <- str_replace(colnames(data$info)[pepCols], "_Peptides", paste0(" ", replacement))
   }else {
      colnames(data$info)[pepCols] <- replacement
   }

   writeData(wb, sheetName, data$info, startRow = 2)
   freezePane(wb, sheetName, firstActiveRow = 3, firstActiveCol = 3)

   filterSuffix <- function(suffix) {
      index <- endsWith(colnames(data$columns), suffix)
      if(any(index)) {
         columns <- data$columns[, index][, dataset$order]
         return(columns)
      }
      return(NULL)
   }

   data$scaled <- filterSuffix("scaled")
   data$sum <- filterSuffix("sum")

   if(!is.null(data$sum)) {
      if(is.null(data$scaled) || ncol(data$sum) != ncol(data$sum)) {
         if(length(data$classes) > 1) {
            sum_cols <- colnames(data$sum)
            sum_dim <- dim(data$sum)
            data$scaled <- matrix(numeric(prod(sum_dim)), sum_dim[[1]])
            for(class in data$classes) {
               class_cols <- startsWith(sum_cols, class)
               if(any(class_cols)) {
                  col_sums <- data$sum[, class_cols]
                  data$scaled[, class_cols] <- 100 * col_sums / rowSums(col_sums, na.rm = TRUE)
               }
            }
         }else {
            data$scaled <- 100 * data$sum / rowSums(data$sum, na.rm = TRUE)
         }
      }
   }

   if(!is.null(data$scaled)) colnames(data$scaled) <- dataset$names
   if(!is.null(data$sum)) colnames(data$sum) <- dataset$names

   nextCol <- ncol(data$info) + 1
   nextCol <- populateSheet(wb, data$scaled, sheetName, nextCol, "#e2f0d9", headers[["scaled"]])
   nextCol <- populateSheet(wb, data$sum, sheetName, nextCol, "#deebf7", headers[["sum"]])

   headerStyle <- createStyle(textDecoration = "bold", halign = "center")
   addStyle(wb, sheetName, headerStyle, rows = 1:2, cols = 1:nextCol, gridExpand = TRUE, stack = TRUE)
   setColWidths(wb, sheetName, 1:nextCol, width = "auto", ignoreMergedCells = TRUE)

   return(c(colnames(data$info), colnames(data$columns)))
}

populateSheet <- function(wb, data, sheetName, startCol, color, header) {
   if(is.null(data))
      return(startCol)

   endCol <- startCol + ncol(data) - 1
   colorStyle <- createStyle(fgFill = color)

   writeData(wb, sheetName, data, startRow = 2, startCol = startCol)
   mergeCells(wb, sheetName, rows = 1, cols = startCol:endCol)
   addStyle(wb, sheetName, colorStyle, rows = 1:2, cols = startCol:endCol, gridExpand = TRUE)
   writeData(wb, sheetName, header, startCol = startCol)
   return(endCol + 1)
}
