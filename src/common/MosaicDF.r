#----------------------------------------------------
#                     Mosaic Core
#
# Validated, unified internal data format
# Only contains columns with a defined datatype validator in 'core_format'(see below)
# Validated columns will be added to 'info'. Remaining columns will be separated as 'columns'
# and validated with is.numeric
# =================== Translators ===================
invalid_format <- function(msg) stop(paste("Check Format:", msg), call. = FALSE)

gen_with_unique_names <- function(err) {
   function(l, desc, f) {
      if(!is.list(l))
         invalid_format(paste(desc, "must be a list"))

      l_names <- names(l)
      if(is.null(l_names) || !all(nzchar(l_names)))
         invalid_format(paste(desc, "must be a named list"))

      dup_names <- duplicated(l_names)
      if(any(dup_names)) {
         dups <- paste(l_names[dup_names], collapse = ", ")
         invalid_format(sprintf("%s must be a list with unique names - (%s) duplicated", desc, dups))
      }

      Map(f, l_names, l)
   }
}

#' Create a new format translator
#'
#' Convenience class that defines column name transformation for a data.frame.
#'
#' @details
#' The specification is a list where names are the final columns which must be unique and
#' the values are either a string (interpreted as list(match = <string>, regex = FALSE))
#' or a list with elements:
#' - pattern: String to match against original columns
#' - regex: Whether `pattern` is to be considered a regex pattern. Defaults to FALSE if missing
#'
#' @param translator - Format translator specification
format_translator <- function(translator) {
   invalid_translator <- function(msg) stop(paste("Translate Format:", msg), call. = FALSE)
   with_unique_names <- gen_with_unique_names(invalid_translator)

   translator <- with_unique_names(translator, "Format Translator", function(col, val) {
      label <- sprintf("Column [%s]", col)
      invalid_col <- function(msg) invalid_translator(paste(label, msg))

      fields <- names(val)
      if(is.character(val)) {
         if(length(val) != 1)
            invalid_col("expected a single string")
         spec <- list(pattern = val)
      }else if(is.list(val)) {
         if(!"pattern" %in% fields)
            invalid_col("'pattern' must be provided")

         pattern <- val$pattern
         if(!is.character(pattern) && length(pattern) != 1)
            invalid_col("'pattern' must be a single string")
         spec <- list(pattern = pattern)

         if("regex" %in% fields) {
            regex <- val$regex
            if(!is.logical(regex) || length(regex) != 1)
               invalid_col("'regex' is not a single boolean")
            spec$regex <- regex
         }
      }else {
         invalid_col("Value must be a list or a string")
      }

      if(!"regex" %in% fields)
         spec$regex <- FALSE

      return(spec)
   })

   structure(translator, class = "format_translator")
}

#' Apply a format_translator to data
#'
#' Replace column names according to `translator`. If no column matches, throw an error.
#' Each column is translated at most once, in order of translations. For example, if a column
#' matches a pattern and then the translation would match another pattern, the second translation
#' is not applied. Regex patterns replace only the regex match, not the full column name.
#'
#' @param data - A data.frame to translate
#' @param translator - A `format_translator` object which specifies the translation to apply
#'
#' @return The `data` data.frame with column names updated according to `translator`
translate_format <- function(data, translator) {
   stopifnot(inherits(translator, "format_translator"))

   if(is.null(data))
      return(NULL)

   invalid_translate <- function(msg) stop(paste("Translate:", msg), call. = FALSE)

   if(any(duplicated(names(data))))
      invalid_translate("Data has duplicated column names")

   translated <- rep(FALSE, ncol(data))
   for(key in names(translator)) {
      spec <- translator[[key]]

      desc <- ""
      pattern <- spec$pattern
      if(spec$regex) {
         match <- str_detect( names(data)[!translated],        pattern)
         subset <- names(data)[!translated][match]
         replacement <- str_replace(subset, pattern, key)

         desc <- " (regex)" #clearer error messages

         if(any(duplicated(replacement)))
            invalid_translate(paste("Regex substitution produced duplicated column names from:",
               paste(subset, collapse = ", ")))
      }else {
         match <- pattern == names(data)[!translated]
         replacement <- key
      }

      if(!any(match))
         invalid_translate(sprintf("No matching column was found for {%s}%s", pattern, desc))

      names(data)[!translated][match] <- replacement
      translated[!translated][match] <- TRUE

      dups <- duplicated(names(data))
      if(any(dups))
         invalid_translate(paste("Replacement introduced duplicate column names of:",
            paste(unique(names(data)[dups]), collapse = ",")))
   }

   return(data)
}


#' Convert a 'format_translator' to its character representation
#'
#' @param translator - A `format_translator` object which specifies the translation to convert
#'
#' @return The character vector representing the `translator`
translator_str <- function(translator) {
   stopifnot(inherits(translator, "format_translator"))

   sapply(translator, function(spec){
      if(spec$regex) {
         sprintf("/%s/", spec$pattern)
      }else {
         spec$pattern
      }
   })
}

#' Create a new format validator
#'
#' Column Validation Mini-Language
#' Description of the data interface of the 'info' data.frames of proteins and peptides
#' Required = TRUE means that data can be accessed without checking that the column exists
#'    Throw an error if one of these columns is missing, the data is invalid
#' Required = FALSE means that input data may be missing that column and the application must handle it gracefully
#' If the data exists, it's type will be validated by the function supplied in 'valid'
#'
#' @details
#' Each column exists in exactly one grouping. The full validation scheme is a named list of groupings
#' A grouping is a list with the following elements:
#'    match - A predicate for deciding group membership with parameters: is_site_quant, is_peptide
#'    columns - Named list of column specifications where the name is the column name to validate
#' A column specification is either a function, f, (interpreted as list(valid = f)) or a list with the
#' following elements:
#'    valid - A predicate passed the column data for validation of type / invariants. The column is passed as a vector if regex = FALSE and a data.frame otherwise.
#'            Throw your own error, or return FALSE for a generic "<Column> is invalid" error
#'    optional - If included and TRUE, allows the column to be missing. Defaults to FALSE
#'    regex - If included and TRUE, applies the specification to all columns that match the
#'                 column name as a regex. Defaults to FALSE (literal match)
#'
#' **Note: Any columns added after data has been saved must be marked as optional
#'    for compatibility with existing viewers
#' On the other hand, only relax validation or set required columns as optional when it becomes
#'    necessary to accommodate a new format since it can never be required again
#'
#' @param format - Format specification
#'
#' @returns Validated 'format_validator' object
format_validator <- function(format) {
   with_unique_names <- gen_with_unique_names(invalid_format)

   format <- with_unique_names(format, "Format", function(group_name, group) {
      label <- sprintf("Group [%s]", group_name)
      invalid_group <- function(msg) invalid_format(paste(label, msg))

      if(!is.list(group) || !all(c("match", "columns") %in% names(group)))
         invalid_group("must be a list with 'match' and 'columns'")

      match <- group$match
      if(!is.function(match))
         invalid_group("'match' must be a function")

      group$columns <- with_unique_names(group$columns, paste(label, "columns"), function(col_name, col_spec) {
         invalid_col_spec <- function(msg) invalid_group(sprintf("> column [%s] %s", col_name, msg))

         if(is.function(col_spec)) {
            col_spec <- list(valid = col_spec)
         }else {
            if(!is.list(col_spec))
               invalid_col_spec("must be a function or list")

            if(!"valid" %in% names(col_spec))
               invalid_col_spec("must include 'valid'")

            valid <- col_spec$valid
            if(!is.function(valid))
               invalid_col_spec("'valid' must be a function")

            bool_if_exists <- function(field) {
               if(field %in% names(col_spec)) {
                  val <- col_spec[[field]]
                  if(is.na(val) || !is.logical(val) || length(val) != 1)
                     invalid_col_spec(sprintf("'%s' must be a boolean", field))
               }
            }

            bool_if_exists("optional")
            bool_if_exists("regex")
         }
         return(col_spec)
      })
      return(group)
   })

   structure(format, class = "format_validator")
}

#' Check data according to a 'format_validator'
#'
#' @param data Data to be validated
#' @param format Format specification as a 'format_validator' object
#' @param ... Parameters passed to the group match predicates
check_format <- function(data, format, ...) {
   stopifnot(inherits(format, "format_validator"))
   invalid_data <- function(col_name, err = "invalid") stop(sprintf("Column '%s' is %s", col_name, err), call. = FALSE)

   data_cols <- colnames(data)

   for(group in format) {
      match <- group$match
      if(match(...)) {
         columns <- group$columns
         Map(function(col_name, col_spec) {
            valid <- col_spec$valid
            specs <- names(col_spec)

            with_default <- function(spec, dflt) ifelse(spec %in% specs, col_spec[[spec]], dflt)

            optional <- with_default("optional", FALSE)
            regex    <- with_default("regex", FALSE)

            cols <- col_name
            if(regex) {
               cols <- str_subset(data_cols, col_name)
            }else {
               cols <- intersect(data_cols, col_name)
            }

            if(length(cols) == 0) {
               if(!optional) {
                  msg <- "missing"
                  if(regex)
                     msg <- paste(msg, "[regex]")
                  invalid_data(col_name, msg)
               }
            }else {
               col_data <- if(regex) {
                  data[cols]
               }else {
                  data[[cols]]
               }

               if(!valid(col_data))
                  invalid_data(col_name)
            }
         }, names(columns), columns)
      }
   }

   return(data)
}

# Validated, unified internal data format
# Only contains columns with a defined datatype validator
# Validated columns will be added to 'info'. Remaining columns will be separated as 'columns'
# and validated with is.numeric
core_format <- format_validator(list(
   core = list(
      match = function(...) TRUE, # these columns are always guaranteed
      columns = list(
         UniprotID  = is.character,
         GeneSymbol = is.character
      )
   ),
   protein_core = list(
      match = function(is_peptide = TRUE, ...) !is_peptide,
      columns = list(
         Description = is.character,
         Peptides    = list(
            valid = function(cols) all(sapply(cols, is.numeric)),
            regex = TRUE
         )
      )
   ),
   peptide_core = list(
      match = function(is_peptide = FALSE, ...) is_peptide,
      columns = list(
         Class           = is.character,
         SearchId        = list(valid = is.numeric,   optional = FALSE),# Concat w/ url
         PeptideId       = list(valid = is.numeric,   optional = FALSE),# Displayed and concat w/ url
         PeptideSequence = is.character,
         FirstScanNumber = list(valid = is.numeric,   optional = FALSE),# Displayed and concat w/ url
         RunLoadPath     = list(valid = is.character, optional = FALSE) # Only displayed
      )
   ),
   protein_site = list(
      match = function(is_peptide = TRUE, is_site_quant = FALSE, ...) !is_peptide && is_site_quant,
      columns = list(
         Site = is.character, # numeric seperated by ";"
         Sequence = is.character
      )
   ),
   peptide_site = list(
      match = function(is_peptide = TRUE, is_site_quant = FALSE, ...) is_peptide && is_site_quant,
      columns = list(
         Site = is.numeric
      )
   )
))

validate_mosaic_df <- function(proteins, peptides, isSiteQuant) {
   all_numeric <- function(df, desc) {
      is_numeric <- vapply(df, is.numeric, logical(1))
      if(!all(is_numeric)) {
         cols <- paste(names(df)[!is_numeric], ", ")
         stop(paste(desc, "columns must be numeric. Non-numeric columns:", cols))
      }
   }
   check_data <- function(data, is_peptide) {
      check_format(data, core_format, is_peptide = is_peptide, is_site_quant = isSiteQuant)

      is_quant <- attr(data, "is_quant")
      if(is.null(is_quant) || length(is_quant) != ncol(data) || !is.logical(is_quant))
         stop("Expected 'is_quant' attr to be a logical vector along data columns")

      if(sum(is_quant) == 0)
         stop("Expected at least one data column")

      all_numeric(data[, is_quant, drop = FALSE])
   }

   check_data(proteins, FALSE)

   # Parse classes
   pepCounts <- str_subset(names(proteins), "Peptides")
   classes <- str_remove(pepCounts, "[~_ ]?Peptides$")
   if(length(classes) == 1 && !nzchar(classes)) {
      classes <- "default"
   }

   protCols <- colnames(proteins)[attr(proteins, "is_quant")]

   class_regex <- paste0("^(", str_vec_regex(classes), ")")
   if(length(classes) > 1) {
      # Seems like this is probably specific to MassPike format, but who knows where the code assumes it...
      has_class <- str_detect(protCols, class_regex)
      if(any(!has_class)) {
         stop(paste("One or more protein column names do not begin with a class:",
              paste(protCols[!has_class], collapse = ",")))
      }
   }

   if(!is.null(peptides)) {
      check_data(peptides, TRUE)

      pepCols <- colnames(peptides)[attr(peptides, "is_quant")]

      # Test that there is a peptide quant for each protein quant
      for(class in classes) {
         for(pepCol in pepCols) {
            regex <- sprintf("%s[~_]%s", class, pepCol)
            matches <- str_detect(protCols, regex)
            if(sum(matches) > 1) {
               stop(paste0("Expected regex /", regex, "/ to match one at most one protein quant column. Matched ",
                           sum(matches), ": [",
                           paste(protCols[matches], collapse=","), "]"))
            }

            protCols <- protCols[!matches]
         }
      }
      if(length(protCols) != 0) {
            stop(paste0("Expected each protein quant column to match a peptide quant, class combination. ",
                        "Failed to match: [", paste(protCols[matches], collapse=","), "]"))
      }
   }

   return(classes)
}

#-------------------------------------------------------
#                     Mosaic DF
#
# Finalized, validated data format

#' Normalize various IDs
#'
#' @param idMap - Optional map from default ids to custom ids, applied to the UniprotID column
#' @param uniqueGS - Ensure unique gene symbols with make.unique
formatIDs <- function(data, uniqueGS = TRUE, customId = NULL, idMap = NULL) {
   # Filter IDs
   isContaminant <- str_detect(data$UniprotID, "contaminant")
   data <- data[!isContaminant, ]
   attr(data, "contam") <- sum(isContaminant)

   isReverse <- str_detect(data$UniprotID, "##")
   data <- data[!isReverse, ]
   attr(data, "decoys") <- sum(isReverse)

   # Replace IDs
   if(!is.null(idMap)){
      index <- data$UniprotID %in% names(idMap)
      data$UniprotID[index] <- idMap[data$UniprotID[index]]
   }

   # Format Uniprot-style IDs
   data$UniprotID <- extractUniprot(data$UniprotID)
   
   # Ensure each protein has an id
   index <- is.na(data$UniprotID)
   if(any(index)){
      data$UniprotID[index] <- paste("Unknown_", 1:sum(index))
   }

   # Set empty GeneSymbols to their UniprotIDs
   index <- is.na(data$GeneSymbol)
   if(any(index)){
      data$GeneSymbol[index] <- data$UniprotID[index]
   }

   # Filter GeneSymbols
   isReverse <- str_detect(data$GeneSymbol, "##")
   data <- data[!isReverse, ]
   attr(data, "decoys") <- attr(data, "decoys") + sum(isReverse)

   index <- str_detect(data$GeneSymbol, "Rik.{0,1}")
   if (any(index)) {
      data$GeneSymbol[index] <- paste0("Rik_", data$GeneSymbol[index])
   }

   # Add MosaicIDs
   data <- data %>% mutate(
      MosaicID = paste0(GeneSymbol, "__", UniprotID) %>%
         str_replace_all("-", "_")
   )

   # Tag GeneSymbols and MosaicIDs with a custom id
   if(!is.null(customId)){
      ids <- data[[customId]]
      index <- !is.na(ids) & nzchar(ids)

      if(any(!index))
          stop(paste(customId, "contains NA or empty values"))

      data$GeneSymbol[index] <- paste(data$GeneSymbol[index], ids[index], sep = "__")
      data$MosaicID[index]   <- paste(data$MosaicID[index],   ids[index], sep = "__")
   }

   # Make GeneSymbols unique
   if(uniqueGS){
      data$GeneSymbol <- make.unique(data$GeneSymbol, sep = "_")
   }

   return(data)
}

#' Performs common data validation and normalization
as.mosaic_df <- function(proteins, peptides, isSiteQuant, ...) {
   if(is.null(proteins))
      return(NULL)

   classes <- validate_mosaic_df(proteins, peptides, isSiteQuant)

   customId <- if(isSiteQuant) "Site" else NULL
   proteins <- formatIDs(proteins, customId = customId, ...)

   # Remove rows with 0 peptides
   pepCols <- str_subset(names(proteins), "Peptides")
   zeroPeps <- rowSums(proteins[pepCols]) == 0
   proteins <- proteins[!zeroPeps, , drop = FALSE]

   # Clean descriptions
   desc <- proteins$Description
   desc[is.na(desc)] <- "<NA>"
   index <- str_detect(desc, " OS=")
   if(any(index)){
      desc[index] <- str_extract(desc[index], ".*(?= OS=)")
      proteins$Description <- desc
   }


   # Split data
   is_quant <- c(attr(proteins, "is_quant"), FALSE) # MosaicID column was added
   proteins <- list(
      columns = as.matrix(proteins[, is_quant, drop = FALSE]),
      info = as.data.frame(proteins[, !is_quant, drop = FALSE]),
      classes = classes,
      suffixes = attr(proteins, "suffixes"),
      suffixes_available = attr(proteins, "suffixes_available"),
      contaminants = attr(proteins, "contam"),
      reverses = attr(proteins, "decoys"),
      noPeps = sum(zeroPeps)
   )

   if(length(pepCols) > 1) {
      # Set plexes with 0 peptides to NA
      for(pepCol in pepCols) {
         class <- str_remove(pepCol, ".Peptides")
         classCols <- startsWith(colnames(proteins$columns), class)
         naRows <- proteins$info[pepCol] == 0
         proteins$columns[naRows, classCols] <- NA
      }
   }

   # Index by MosaicID
   rownames(proteins$columns) <- proteins$info$MosaicID
   rownames(proteins$info)    <- proteins$info$MosaicID

   if(!is.null(peptides)) {
      peptides <- formatIDs(peptides, uniqueGS = FALSE, customId = customId, ...)

      is_quant <- c(attr(peptides, "is_quant"), FALSE)
      peptides <- list(
         columns = as.matrix(peptides[, is_quant, drop = FALSE]),
         info = as.data.frame(peptides[, !is_quant, drop = FALSE]),
         contaminants = attr(peptides, "contam"),
         decoys = attr(peptides, "decoys")
      )
   }

   return(list(proteins = proteins, peptides = peptides))
}

# ========== Utilities ==========

# For use with specifics & error messages?
parse_warning <- function(msg, call = sys.call(-1)) {
   w <- structure(
      list(message = msg, call = call),
      class = c("parse_failed", "warning", "condition")
   )
   warning(w)
}

extractUniprot <- function(IDs) {
   # Ex:               (sp|)Q80UM7(|MOGS_MOUSE)
   uniprotRegex <- "(?<=\\|).*(?=\\|)"
   uniprotIDs <- str_extract(IDs, uniprotRegex)
   # Only replace on successful extraction
   isUniprotFormat <- !is.na(uniprotIDs)
   IDs[isUniprotFormat] <- uniprotIDs[isUniprotFormat]
   return(IDs)
}

#' Create a regex for matching the longest value in @patterns first
#'
#' @param patterns - possible values to match
#' @return Fully escaped and sorted regex

str_vec_regex <- function(patterns) {
   byLen <- order(nchar(patterns), decreasing = TRUE)
   patterns <- str_escape(patterns[byLen])
   paste0(patterns, collapse = "|")
}

#Auto load formats on demand?
# Returns only columns that match each pattern in order
# Non-exact matches are computed as regex patterns
nameSelect <- function(data, patterns){
   index <- NULL
   unmatched <- names(data)
   for(pattern in patterns){
      if(pattern %in% unmatched){
         add <- pattern
      }else{
         add <- str_subset(unmatched, pattern)
      }
      unmatched <- setdiff(unmatched, add)
      index <- c(index, add)
   }
   return(data[index])
}

# Replaces any occurances of pattern with names(pattern) in each name
# Each name is only modified once
nameReplace <- function(data, patterns){
   names <- names(data)
   replaced <- rep(NA, ncol(data))
   exact <- patterns %in% names(data)
   if(any(exact)){
      replaced[match(patterns[exact], names)] <- names(patterns)[exact]
   }

   if(any(!exact)){
      regex <- patterns[!exact]
      for(i in 1:length(regex)){
         pattern <- regex[i]
         subset <- is.na(replaced)
         subset[subset]   <- str_detect(names[subset], pattern)
         replaced[subset] <- str_replace(names[subset], pattern, names(pattern))
      }
   }

   new <- !is.na(replaced)
   names(data)[new] <- replaced[new]
   return(data)
}

