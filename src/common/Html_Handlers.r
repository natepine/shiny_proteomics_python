## Add functionality for handling HTML outputs
bold <- function(string) {
   return(paste0("<b>", string, "</b>"))
}

color <- function(string, color) {
   return(paste0("<div style=\'color:", color, "\'>", string, "</div>"))
}

defaultBox <- function(..., status = "primary", width = 12) {
   box(width = width, status = status, solidHeader = TRUE,
      ...
   )
}

centeredCol <- function(..., width = 12) {
   column(width, offset = (12 - width) / 2, align = "center",
      ...
   )
}

keyValueStr <- function(data, sep = ": ", collapse = "<br>", boldKeys = FALSE, boldVals = FALSE){
   data <- data[!is.na(data)]

   keys <- names(data)
   vals <- data

   if(boldKeys) keys <- bold(keys)
   if(boldVals) vals <- bold(vals)

   paste(keys, sep = sep, vals, collapse = collapse)
}

html_escape <- function(str, quote = c(NULL, "'", '"')) {
   quote <- match.arg(quote)

   str <- str_replace_all(str, ">", ";gt")
   str <- str_replace_all(str, "<", ";lt")
   str <- str_replace_all(str, "&", ";amp")

   if(!is.null(quote)) {
      esc <- switch(quote,
         "'" = "&#39",
         '"' = "&quot",
         "")

      str <- str_replace_all(str, quote, esc)
   }
   return(str)
}
