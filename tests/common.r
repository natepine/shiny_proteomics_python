createStore <- function(path) createCache(
   path,
   function(...) write.csv(..., row.names = FALSE),
   function(...) read.csv(..., check.names = FALSE),
   ext = ".csv",
   isValidData = is.data.frame
)

# Cache conversions
norm_to_cache <- function(norm) {
   norm$factors <- toJSON(norm$factors, auto_unbox = TRUE)
   return(as.data.frame(norm))
}

norm_from_cache <- function(norm) {
   norm <- as.list(norm)
   norm$factors <- parse_json(norm$factors)
   return(norm)
}

mosaic_to_cache <- function(data, type = c("proteins", "peptides"), class_label) {
   type <- match.arg(type)

   if(type == "proteins") {
      cbind(data$proteins$info, data$proteins$columns)
   }else if(type == "peptides") {
      flat_data <- data$peptides$info

      classes <- unique(data$peptides$info$Class)
      for(class in classes) {
         class_index <-data$peptides$info$Class == class 
         class_info <- data$peptides$info[class_index, , drop = FALSE]
         columns <- data$peptides$columns[[class]][class_info$data_row, , drop = FALSE]
         flat_data[class_index, colnames(columns)] <- columns
      }
      return(flat_data)
   }
}

# Logging Utils
catln <- function(...) cat(..., sep = "\n")
cat_val <- function(label, val) catln(paste0(label, ":\t", val))
cat_hdr <- function(hdr) {
   catln("",
         hdr,
         strrep("=", nchar(hdr)))
}

createLog <- list

logEvent <- function(log, key, type, details = "") {
   if(!type %in% names(log))
      log[[type]] <- list()
   log[[type]] <- append(log[[type]], sprintf("<%s> %s", key, details))
   return(log)
}

printLog <- function(hdr, log, writeLog = catln) {
   if(length(log) == 0)
      return(NULL)

   cat_hdr(hdr)
   for(type in names(log)) {
      writeLog(type, unlist(log[[type]]))
   }
}

summarizeLog <- function(log) sapply(log, length)

countLog <- function(log) {
   if(length(log) == 0)
      return(0)
   sum(log)
}

