# ========== Translation ==========

gfy_schemas <- list(
   combined_site_quant = list(
      UniprotID = "Protein Id",
      GeneSymbol  = "gene_symbol",
      Description = "prot_description",
      Peptides    = list(pattern = "num_quant", regex = TRUE),
      Site        = "Site Position",
      Sequence    = "sequence"
   ),
   protein_quant = list(
      UniprotID   = "Protein Id",
      GeneSymbol  = "Gene Symbol",
      Description = "Description",
      Peptides    = list(pattern = "(Number of p|P)eptides", regex = TRUE)
   ),
   site_quant_peptide = list(
      UniprotID         = "protein_id",
      GeneSymbol        = "GeneSymbol",
      Class             = "class",
      SearchId          = "search_id",
      PeptideId         = "peptide_id",
      PeptideSequence   = "sequence",
      FirstScanNumber   = "scan",
      RunLoadPath       = "file",
      Site              = "site_position"
   ),
   protein_quant_peptide = list(
      UniprotID         = "ProteinId",
      GeneSymbol        = "GeneSymbol",
      Class             = "Class",
      SearchId          = "SearchId",
      PeptideId         = "PeptideId",
      PeptideSequence   = "PeptideSequence",
      FirstScanNumber   = "FirstScanNumber",
      RunLoadPath       = "RunLoadPath"
   )
)

gfy_formats <- lapply(gfy_schemas, format_translator)

valid_suffixes <- c("sum", "scaled", "median", "mean", "std_dev")

columnsRegex <- "1[23][0-9][nceNCE]?[dD]?_sn"

# ========== Validation ==========

parse_warning <- function(msg, call = sys.call(-1)) {
   w <- structure(
      list(message = msg, call = call),
      class = c("parse_failed", "warning", "condition")
   )
   warning(w)
}
# ========== Normalization ==========

gfyFormatNormalization <- identity

# ========== Proteins ==========

#' Translate CORE ProteinQuant or SiteQuant data to a common format
#'
#' @param data - Raw data to validate and translate
#' @param isSiteQuant - Is `data` ProteinQuant or SiteQuant
#' @param format - A format_translator for the columns returned by CORE.
#' @param suffixes - Select columns. Must be a subset of c("sum", "scaled", "median", "mean", "std_dev")
#' @param suffix_select - How to perform column selection.
#'          first: Attempts to select each value, only returns the first match and sets suffixes to that match.
#'          in: Select suffixes from @suffixes or return NULL if no columns match. Allows missing suffixes if at least one column is selected.
#'          all: Select each suffix from @suffixes or NULL if any suffix is missing.
#'          none: Skip column selection (return full dataset) and set suffixes to all available suffixes.
#'
#' @return A data.frame with standardized column names. Attributes include data for suffixes, counts of contaminants, reverses,
#'         and zero-peptide rows that have been filtered out of the dataset
#'
#'         Or NULL to indicate no data available

gfyFormatProteins <- function(data, isSiteQuant, format = NULL,
      suffixes = c("scaled", "sum"), suffix_select = "first") {

   if(!is.null(format)) {
      stopifnot(inherits(format, "format_translator"))
   }else {
      module <- if(isSiteQuant) "combined_site_quant" else "protein_quant"
      format <- gfy_formats[[module]]
   }

   if(!(length(suffix_select) == 1 && suffix_select %in% c("first", "in", "all", "none"))) {
      msg <- paste0(suffix_select, collapse = ",")
      stop(sprintf("Invalid suffix_select (%s)", msg))
   }
   if(!all(suffixes %in% valid_suffixes)) {
      msg <- paste0(suffixes, collapse = ",")
      stop(sprintf("Invalid suffixes (%s)", msg))
   }

   # Format
   data <- translate_format(data, format)

   if(isSiteQuant)
      data$Site <- as.character(data$Site)
   data$GeneSymbol <- as.character(data$GeneSymbol)

   # Select columns
   dataCols <- str_subset(colnames(data), columnsRegex)

   non_suffix_reg <- paste0(".*", columnsRegex, "_")
   suffixes_available <- unique(str_remove(dataCols, non_suffix_reg))

   if(suffix_select != "none") {
      suffix_found <- suffixes %in% suffixes_available
      if(suffix_select == "all") {
         # Select every suffix
         if(!all(suffix_found)) {
            msg <- paste(sep = "\n",
               "Failed to select ALL suffixes:", paste0(suffixes, collapse = "\n"),
               "\nMissing suffixes:",
               paste(suffixes[!suffix_found], collapse = "\n")
            )
            parse_warning(msg)
            return(NULL)
         }

         selected <- c()
         for(s in suffixes) {
            s_ind <- endsWith(dataCols, s)
            selected <- c(selected, dataCols[s_ind])
         }
      }else if(suffix_select == "in") {
         # Select intersection of suffixes and suffixes_available
         if(!any(suffix_found)) {
            msg <- paste(sep = "\n",
               "Failed to select any suffixes IN:", paste0(suffixes, collapse = "\n"),
               "\nAvailable suffixes:",
               paste(suffixes_available, collapse = "\n")
            )
            parse_warning(msg)
            return(NULL)
         }

         suffixes <- suffixes[suffix_found]
         selected <- c()
         for(s in suffixes) {
            s_ind <- endsWith(dataCols, s)
            selected <- c(selected, dataCols[s_ind])
         }
      }else if(suffix_select == "first") {
         # Select first match
         if(!any(suffix_found)) {
            msg <- paste(sep = "\n",
               "Failed to select FIRST suffix from:", paste0(suffixes, collapse = "\n"),
               "\nAvailable suffixes:",
               paste(suffixes_available, collapse = "\n")
            )
            parse_warning(msg)
            return(NULL)
         }

         suffixes <- suffixes[suffix_found][[1]]
         s_ind <- endsWith(dataCols, paste0("_", suffixes))
         dataCols <- dataCols[s_ind]
      }
   }else {
      suffixes <- suffixes_available
   }

   if(!all(suffixes %in% valid_suffixes)) {
      msg <- paste("One or more parsed suffixes are invalid:", suffixes, sep = "\n")
      parse_warning(msg)
      return(NULL)
   }

   infoColsRegex <- paste0(str_vec_regex(names(format)))
   data <- data[, c(str_subset(colnames(data), infoColsRegex), dataCols), drop = FALSE]

   attr(data, "is_quant") <- colnames(data) %in% dataCols
   attr(data, "suffixes") <- suffixes
   attr(data, "suffixes_available") <- suffixes_available

   return(data)
}

# ========== Peptides ==========

#' Translate CORE ProteinQuant or SiteQuant peptide data to a common format
#'
#' @param data - Raw data to validate and translate
#' @param isSiteQuant - Is `data` ProteinQuant or SiteQuant
#'
#' @return A data.frame with standardized column names and an attribute describing which columns 
#'         are quantitative columns
gfyFormatPeptides <- function(data, isSiteQuant) {
   module <- if(isSiteQuant) "site_quant_peptide"
                       else  "protein_quant_peptide"
   format <- gfy_formats[[module]]

   data <- translate_format(data, format)
   allColsRegex <- paste0(str_vec_regex(names(format)), "|", columnsRegex)
   data <- data[, str_subset(colnames(data), allColsRegex), drop = FALSE]

   columnsRegex <- paste0(columnsRegex, "$")
   is_quant <- str_detect(colnames(data), columnsRegex)

   if(!"GeneSymbol" %in% colnames(data))
      data$GeneSymbol <- NA

   data$Class <- as.character(data$Class)
   data$GeneSymbol <- as.character(data$GeneSymbol)
   attr(data, "is_quant") <- is_quant

   return(data)
}

# ============= Data Formats =============

# Add missing GeneSymbols column to site quant peptide output using the mapping in protein output
withGSMap <- function(fmt) {
   gsFmt <- fmt

   gsFmt$proteins <- function(dataset) {
      if(!dataset$isSiteQuant)
         return(fmt$proteins(dataset))

      schema <- gfy_schemas[["combined_site_quant"]]
      gsMap <- setNames(dataset$proteins[[schema$GeneSymbol]], dataset$proteins[[schema$UniprotID]])

      return(createSrcState(
         fmt$proteins(dataset),
         list(gsMap = gsMap)
      ))
   }

   gsFmt$peptides <- function(dataset) {
      if(!dataset$isSiteQuant)
         return(fmt$peptides(dataset))

      schema <- gfy_schemas[["site_quant_peptide"]]
      symbols <- dataset$state$gsMap[dataset$peptides[[schema$UniprotID]]]
      index <- !is.na(symbols)
      if(any(index)) {
         dataset$peptides[[schema$GeneSymbol]][index] <- symbols[index]
      }else {
         dataset$peptides[[schema$GeneSymbol]] <- NA_character_
      }

      return(createSrcState(
         fmt$peptides(dataset),
         list(gsMap = NULL) # delete state
      ))
   }

   return(gsFmt)
}

gfyFormat <- withGSMap(createDatasetFormat(
   proteins      = function(dataset) gfyFormatProteins(dataset$proteins, dataset$isSiteQuant),
   normalization = function(dataset) gfyFormatNormalization(dataset$normalization),
   peptides      = function(dataset) gfyFormatPeptides(dataset$peptides, dataset$isSiteQuant)
))

gfyProteinExample <- function(isSiteQuant) {
   row <- if(isSiteQuant) {
      c(
         "sp|O60343|TBCD4_HUMAN",
         "TBC1D4",
         "TBC1 domain family member 4",
         "3",
         "370",
         "FEINLISPDTKSV",
         "107.08"
      )
   }else {
      c(
         "sp|P55011|S12A2_HUMAN",
         "SLC12A2",
         "Solute carrier family 12 member 2",
         "3",
         "203.94"
      )
   }

   fmt <- gfy_formats[[if(isSiteQuant) "combined_site_quant" else "protein_quant"]]
   fmt$Columns <- list(pattern = paste0(columnsRegex, "_(scaled|sum)"), regex = TRUE)

   m <- matrix(c("Example Data", row), dimnames = list(c("Pattern", translator_str(fmt))))
   as.data.frame(t(m))
}

gfyPeptideExample <- function(isSiteQuant) {
   row <- if(isSiteQuant) {
      c(
         "sp|Q8VE97|SRSF4_MOUSE",
         "SRSF4",
         "default",
         "16312",
         "11",
         "K.DTDHSRS#PSR.S",
         "557",
         "path/to/ec03881.raw",
         "401",
         "305.91"
      )
   }else {
      c(
         "sp|P55011|S12A2_HUMAN",
         "SLC12A2",
         "default",
         "23066",
         "11558",
         "K.DWLQADMR.D",
         "19315",
         "path/to/eb10744.raw",
         "94.34"
      )
   }

   fmt <- gfy_formats[[if(isSiteQuant) "site_quant_peptide" else "protein_quant_peptide"]]
   fmt$Columns <- list(pattern = paste0(columnsRegex, "_(scaled|sum)"), regex = TRUE)

   m <- matrix(c("Example Data", row), dimnames = list(c("Pattern", translator_str(fmt))))
   as.data.frame(t(m))
}

