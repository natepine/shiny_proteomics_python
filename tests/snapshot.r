appPath <- Sys.getenv("APP_DIR")

suppressPackageStartupMessages({
   library(DBI)

   source(paste0(appPath, "/tests/common.r"))

   source(paste0(appPath, "/IsoParser/global.r"))
   source(paste0(appPath, "/IsoParser/server.r"))
})

# Saves results of loadDataset
summarize_proteins <- function(data) {
   return(head(data))
}

summarize_peptides <- function(data, class_label) {
   stopifnot(is.data.frame(data))

   if(!class_label %in% colnames(data)) {
      stop(sprintf("Data is missing column: '%s'", class_label))
   }

   classes <- data[class_label]
   do.call(rbind, lapply(classes, function(class)
      head(data[data[[class_label]] == class, , drop = FALSE])))
}

snap <- function(name, type, input, desc) {
   cache <- createCache(paste0(app_dir, "/data/"))
   store <- createStore(paste0(app_dir, "/tests/store/"))

   req <- function(x, msg) if(is.null(x)) stop(msg) else x

   if(type == "load") {
      key <- input
      raw_metadata <- readMetadata(key)
      if(nrow(raw_metadata) == 0)
         stop("No viewer for key: ", key)

      raw_metadata$Key      <- "N/A"
      raw_metadata$Dataset  <- name
      raw_metadata$Notes    <- desc
      raw_metadata$Date     <- "N/A"
      raw_metadata$Species  <- "Other"

      raw_metadata$Username  <- "test"

      metadata <- parseMetadata(raw_metadata)

      src <- cacheSrc(cache)
      key <- c("gfyRaw", name)

      proteins <- req(src$proteins(metadata), "No cached proteins")
      proteins <- summarize_proteins(src$proteins(metadata))

      norm <- req(src$normalization(metadata), "No cached normalization")
      raw_norm <- norm_to_cache(norm)

      peptides <- req(src$peptides(metadata), "No cached peptides")
      class_label <- if(metadata$isSiteQuant) "class" else "Class"
      peptides <- summarize_peptides(peptides, class_label)

      prot_key <- c(key, "proteins")
      norm_key <- c(key, "normalization")
      peps_key <- c(key, "peptides")
      meta_key <- c("metadata", paste0(name, ".db"))

      raw_metadata$Server    <- "N/A"
      raw_metadata$QID       <- 0

      # Write
      store$write(proteins,     prot_key, overwrite = TRUE)
      store$write(raw_norm,     norm_key, overwrite = TRUE)
      store$write(peptides,     peps_key, overwrite = TRUE)
      store$write(raw_metadata, meta_key, overwrite = TRUE)

      # Read (use rounded values)
      proteins     <- store$read(prot_key)
      raw_norm     <- store$read(norm_key)
      peptides     <- store$read(peps_key)
      raw_metadata <- store$read(meta_key)

      norm <- norm_from_cache(raw_norm)
      metadata <- parseMetadata(raw_metadata)

      nop <- function(...) NULL
      pass_exception <- function(title, msg, ...) stop(msg)
      view <- createLoadDatasetView(nop, pass_exception, pass_exception, nop, NA)
      src <- createDatasetSource(
         proteins = function(...) proteins,
         normalization = function(...) norm,
         peptides = function(...) peptides,
      )
      mosaic <- loadDataset(metadata, "", NULL, view, src, gfyFormat)
      mosaic_key <- c("mosaic", name)
      mosaic_proteins <- mosaic_to_cache(mosaic, "proteins")
      mosaic_peptides <- mosaic_to_cache(mosaic, "peptides")

      entry <- list(
         name = name,
         desc = desc,
         type = type,
         from = paste(key, collapse = "/"),
         to   = paste(mosaic_key, collapse = "/")
      )
      graph <- store$read("test_graph")
      if(is.null(graph)) {
         graph <- data.frame()
      }
      graph <- graph[!(graph$name == name & graph$type == type)]
      graph <- rbind(graph, entry)

      raw_metadata$InitialID <- mosaic$proteins$info$MosaicID[[1]]

      store$write(raw_metadata, meta_key, overwrite = TRUE)
      store$write(mosaic_proteins, c(mosaic_key, "proteins"), overwrite = TRUE)
      store$write(mosaic_peptides, c(mosaic_key, "peptides"), overwrite = TRUE)
      store$write(graph, "test_graph", overwrite = TRUE)
   }
}

if(!interactive()) {
   types <- c("load")#, "source", "format")

   usage <- function() {
      cat(sprintf(
"
Create a new test by saving a summary of a dataset's structure

Usage: snapshot <name> <type> <input> <description>

Types available - %s

load: <input> must be a valid viewer key
",
         paste(types, collapse = ", ")
      ))
      quit(save = "no", status = 1)
   }

   args <- commandArgs(trailingOnly = TRUE)
   if(length(args) <= 3)
      usage()

   name  <- args[1]
   type  <- args[2]
   input <- args[3]
   desc  <- paste(args[-seq(3)], collapse = " ")

   if(!type %in% types)
      usage()

   tryCatch({
      snap(name, type, input, desc)
      cat("Success", sep = "\n")
   }, error = function(e) {
      cat("Snapshot failed: ", e$message, "\n")
   })
}

