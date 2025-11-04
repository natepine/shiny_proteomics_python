createCache <- function(cacheDir, write = saveRDS, read = readRDS, ext = ".rds", isValidData = function(data) TRUE) {
   if(!startsWith(cacheDir, "/")) {
      cacheDir <- paste(getwd(), cacheDir, sep = "/") # Absolute path
   }

   if(!dir.exists(cacheDir)) {
      dir.create(cacheDir, recursive = TRUE)
   }
   cacheDir <- normalizePath(cacheDir, mustWork = TRUE)

   cachePath <- function(key) {
      if(length(key) == 0) {
         return(NULL)
      }

      path <- cacheDir

      filename <- key[length(key)]
      if(!is.null(ext)) {
         filename <- paste0(filename, ext)
      }

      new_dir <- NULL
      dirs <- key[seq_len(length(key) - 1)]
      # Create dirs as long as they are below the `cacheDir`
      for(dir in dirs) {
         path <- normalizePath(paste(path, dir, sep = "/"), mustWork = FALSE)
         if(startsWith(path, cacheDir)) { # Prevent path traversal outside of sandbox
            if(!file.exists(path)) {
                  new_dir <- c(path, new_dir)

                  # If no file or directory exists, make a new directory in the path
                  dir.create(path)
            }
         }else {
            for(empty_dir in new_dirs) {# Clean up empty directories
               unlink(new_dirs)
            }

            # If we allow keys like c("..", "..", "best_guess_of_parent_dir", "<cacheDir>") leak system details
            # Best to deny anything that crosses `cacheDir`
            return(NULL)
         }
      }

      if(dir.exists(path)) { # Ensure that we end on a directory, not a file
         path <- normalizePath(paste(path, filename, sep = "/"), mustWork = FALSE)
         if(startsWith(path, cacheDir) && !dir.exists(path)) {
            # Sandbox path traversal and disallow references to existing directories, key may be user-supplied
            return(path)
         }
      }
      return(NULL)
   }

   structure(list(
      write = function(data, key, overwrite = FALSE) {
         if(is.null(data))
            return("[INFO] No data to write")

         path <- cachePath(key)
         if(is.null(path))
            return("[WARN] Invalid cache key")

         if(!isValidData(data))
            return("[WARN] Invalid data")

         if(!overwrite && file.exists(path))
            return("[INFO] Overwrite disallowed")

         write(data, path)
         return("")
      },
      read = function(key) {
         data <- NULL
         path <- cachePath(key)
         if(!is.null(path) && file.exists(path)) {
            data <- read(path)
         }
         return(data)
      },
      delete = function(key) {
         path <- cachePath(key)
         if(!is.null(path) && file.exists(path)) {
            unlink(path)
            return(TRUE)
         }
         return(FALSE)
      }
   ), class = "tmtmosaic_cache")
}

mosaic_cache_key <- function(id, isSiteQuant, type = c("proteins", "peptides", "normalization"), origin) {
   type <- match.arg(type)
   prefix <- if(isSiteQuant) "sq" else "pq"
   filename <- paste(sep = "_", prefix, id)
   return(c(origin, type, filename))
}

