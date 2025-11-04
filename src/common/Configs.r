library(yaml)

#keyword_vals <- c("REQUIRED", "OPTIONAL")
keyword_names <- c("WILD_CARD")

as.path <- function(path, ...) {
   paste(path, ..., sep = ">")
}

setRequired <- function(struct) {
   if(is.list(struct)) {
      struct <- lapply(struct, setRequired)
      attr(struct, "REQ") <-
         any(sapply(struct, attr, "REQ"))
   }else {
      if(is.null(struct)) {
         stop("Invalid config template: NULL disallowed")
      }
      attr(struct, "REQ") <- identical(struct, "REQUIRED")
   }
   return(struct)
}

getDefaults <- function(struct, env) {
   if(is.list(struct)) {
      wild <- names(struct) == "WILD_CARD"
      if(any(wild)) {
         recurse <- !wild | sapply(struct, attr, "REQ")
         if(!any(recurse)) {
            return(NA)
         }
         struct <- struct[recurse]
      }
      defaults <- lapply(struct, getDefaults, env)
      names(defaults)[names(defaults) %in% keyword_names] <- "default"
      return(defaults)
   }else {
      default <- struct
      attr(default, "REQ") <- NULL
      if(is.character(default)) {
         if(default == "OPTIONAL") {
            default <- NA
         }else if(substr(default, 1, 1) == "~") {
            obj <- substr(default, 2, nchar(default))
            if(exists(obj, env)) {
               default <- get(obj,env)
            }else {
               default <- NA
            }
         }
      }
      return(default)
   }
}

validate <- function(conf, struct, env = as.environment(-1L), path = "") {
   if(is.list(struct)) {
      if(!is.list(conf)) {
         stop("Expected a list: ", path)
      }

      dup <- duplicated(names(conf))
      if(any(dup)) {
         warning("Ignoring duplicate keys: ",
            as.path(path, unique(names(conf)[dup])))
         conf <- conf[!dup]
      }

      is_keyword <- names(struct) %in% keyword_names
      keywords <- names(struct)[is_keyword]
      unconsumed <- names(conf)
      final_conf <- list()
      build_conf <- function(keys, nextStruct = NULL) {
         new_vals <- lapply(keys, function(key) {
            if(is.null(nextStruct))
               nextStruct <- struct[[key]]
            validate(conf[[key]], nextStruct, env, as.path(path, key))
         })
         names(new_vals) <- keys
         c(final_conf, new_vals)
      }

      if(any(!is_keyword)) {
         struct_names <- names(struct)[!is_keyword]
         valid <- names(conf) %in% struct_names
         if(!"WILD_CARD" %in% keywords && any(!valid)) {
            warning("Ignoring unexpected config: ",
               as.path(path, paste(names(conf)[!valid], collapse = "/")))
            conf <- conf[valid]
         }

         required <- sapply(struct[struct_names], attr, "REQ")
         if(any(required)) {
            req_names <- struct_names[required]
            found <- req_names %in% names(conf)
            if(!all(found)) {
               stop("Missing required config: ",
                  as.path(path, req_names[!found]))
            }
         }

         unconsumed <- names(conf)
         struct_match <- unconsumed %in% struct_names
         if(any(struct_match)) {
            final_conf <- build_conf(unconsumed[struct_match])
            unconsumed <- unconsumed[!struct_match]
         }

         if(any(!required)) {
            opt_names <- struct_names[!required]
            found <- opt_names %in% names(conf)
            if(any(!found)) {
               opt_missing <- opt_names[!found]
               new_vals <- lapply(struct[opt_missing], getDefaults, env)
               names(new_vals) <- opt_missing
               final_conf <- c(final_conf, new_vals)
            }
         }
      }

      if(any(is_keyword)) {
         #if("WILD_CARD" %in% keywords) {

         wild_card <- struct[["WILD_CARD"]]
         if(attr(wild_card, "REQ") && length(unconsumed) == 0) {
            stop("At least one value is required: ", path)
         }
         final_conf <- build_conf(unconsumed, wild_card)
         rm(unconsumed)

         #}
      }
      return(final_conf)
   }else if(length(conf) != 1) {
      stop("Expected scalar value: ",
         path)
   }
   return(conf)
}

#                 Config Parser
# Provide a template with [Name: val] pairs for automated config
# validation and providing missing defaults
#
#                 --- NAMES ---
# WILD_CARD - Allow any name
#     REQUIRED: Must contain at least one value
#     OPTIONAL: May be omitted (parsed as NA)
# non-keyword - Literal match
#
# Duplicates ignored
#
#
#                 --- VALUES ---
# REQUIRED - Stop if missing
# OPTIONAL - Allow missing value (default to NA)
# begins with `~` - Default value supplied by R environment
# non-keyword - Default value (implies optional)
#
# NULL and empty values are not allowed
#

configs_env <- new.env(parent = baseenv())
configs_env$dfl_temp <- yaml.load(
"
SERVERS:
   WILD_CARD:
      User: REQUIRED
      Key: REQUIRED
      Label: OPTIONAL
      VerifySSL: true
ANNOTATED_SPECIES:
   WILD_CARD: OPTIONAL
SUPER_USER_PASSWORD: OPTIONAL
VIEWER_SUBDIR: ~subdir
ADMIN_EMAIL: OPTIONAL
")
configs_env$dfl_env <- (function(){
   appName <- Sys.getenv("APP_NAME")
   if(!nzchar(appName)) appName <- "tmtmosaic"
   subdir <- paste(appName, "IsoParser", sep = "/")

   e <- new.env()
   assign("subdir", subdir, e)
   return(e)
})()

#' Load the YAML file at @configPath and validate it using @template
#'
#' @param configPath - Configuration file to be loaded and validated
#' @param template   - A list containing a configuration template defining the expected schema and validation strategy.
#'                      See "Config Parser" for template syntax details
#' @param env        - An R environment used to inject dynamically defined default values into the template.
#'                      For example, the template value "~x" will be interpreted as the value of the variable x in @env.
#'
#' @return A list with the structure defined by @template and REQUIRED values from the YAML file at @configPath.
#'          All other values may be provided by the template as defaults. Invalid or non-existant config files will throw an error with stop().

loadConfig <- function(configPath, template = configs_env$dfl_temp, env = configs_env$dfl_env) {
   #Read configs
   if(!file.exists(configPath)){
      stop(paste("No config was found at", configPath))
   }
   conf <- read_yaml(configPath)

   conf_struct <- setRequired(template)

   validate(conf, conf_struct, env)
}
