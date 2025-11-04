createDatasetSource <- function(
      proteins = NULL,
      proteins_msg = NULL,
      normalization = NULL,
      normalization_msg = NULL,
      peptides = NULL,
      peptides_msg = NULL
) {
   nop <- function(...) NULL
   nop_dflt <- function(f) if(is.null(f)) nop else f

   src <- lapply(list(
      proteins          = proteins,
      proteins_msg      = proteins_msg,
      normalization     = normalization,
      normalization_msg = normalization_msg,
      peptides          = peptides,
      peptides_msg      = peptides_msg),
      nop_dflt
   )

   structure(
      src,
      class = "tmtmosaic_src"
   )
}

createDatasetFormat <- function(proteins = NULL, normalization = NULL, peptides = NULL) {
   id_dflt <- function(f) if(is.null(f)) identity else f

   lapply(list(
      proteins      = proteins,
      normalization = normalization,
      peptides      = peptides
   ), id_dflt)
}

createSrcState <- function(result, state) {
   stopifnot(is.list(state))
   stopifnot(!is.null(names(state)))

   structure(list(
      result = result,
      state = state
   ),
   class = "tmtmosaic_src_state")
}

dataset_cache_key <- function(dataset, type) mosaic_cache_key(dataset$ID, dataset$isSiteQuant, type, dataset$server)

srcToCache <- function(cache, src, type) {
   function(dataset) {
      rawData <- src(dataset)
      write_status <- cache$write(rawData, dataset_cache_key(dataset, type))
      if(startsWith(write_status, "[WARN]")) {
         cat(write_status, sep = "\n")
      }

      return(rawData)
   }
}

cacheReadSrc <- function(cache, type) {
   function(dataset) {
      rawData <- cache$read(dataset_cache_key(dataset, type))
      return(rawData)
   }
}

cacheSrc <- function(cache) {
   createDatasetSource(
      proteins      = cacheReadSrc(cache, "proteins"),
      normalization = cacheReadSrc(cache, "normalization"),
      peptides      = cacheReadSrc(cache, "peptides")
   )
}
