library("readr")

#'------------------------------------------------
#'
#' Support for loading data from delimited sources
#'
#'------------------------------------------------

#' Get data from a delimited file
#'
#' @param file - Data file
#' @param type - File type of `file`. Either "csv" or "tsv"
#'
#' @return The file's data

src_delim <- function(file, type = c("csv", "tsv")) {
   type <- match.arg(type)
   read <- if(type == "csv") {
      readr::read_csv
   }else if(type == "tsv") {
      readr::read_tsv
   }
   read(file)
}

# ========== Normalization ==========

csvSourceNormalization <- src_delim

# ========== Proteins ==========

csvSourceProteins <- src_delim

# ========== Peptides ==========

csvSourcePeptides <- src_delim

# ============= Data Source =============

csvSource <- function(proteins, peptides = NULL, normalization = NULL){
   proteins <- function(dataset) csvSourceProteins(proteins$file, protiens$type)

   if(!is.null(peptides))
      peptides <- function(dataset) csvSourcePeptides(peptides$file, peptides$type)

   if(!is.null(normalization))
      normalization <- function(dataset) csvSourceNormalization(normalization$file, normalization$type)

   msg <- function(desc, type)
      function(dataset) paste("Loading", desc, "data, please wait.")

   createDatasetSource(
      proteins = proteins,
      peptides = peptides,
      normalization = normalization,

      proteins_msg  = msg("protein"),
      peptides_msg  = msg("peptide"),
      normalization_msg  = msg("normalization")
   )
}

csvSrcToCache <- function(cache, ...) {
   src <- csvSource(...)
   src$proteins      <- srcToCache(cache, src$proteins,      "proteins")
   src$peptides      <- srcToCache(cache, src$peptides,      "peptides")
   src$normalization <- srcToCache(cache, src$normalization, "normalization")
   return(src)
}

