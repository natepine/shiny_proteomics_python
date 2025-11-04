f <- stdin()

header <- function(msg) {
   cat("\n\n---", msg, "---\n\n")
}

intCheck <- function(str) {
   val <- NA
   passed <- str_detect(str, "^[0-9]+$")
   if(passed) val <- as.integer(str)
   
   list(passed = passed, val = val, msg = "Please enter an integer")
}

boolPrompt <- function(..., sep = " ") {
   cat(..., "(y/N) ", sep = sep)
   input <- tolower(readLines(f, n = 1))
   cat("\n")

   input %in% c("y", "yes")
}

prompt <- function(msg, check = NULL, failed_msg = "Invalid input", default = NULL) {
   repeat {
      if(nzchar(msg)) {
         cat(msg)
         if(!is.null(default)) {
            cat(" [Default:", default)
            cat("]")
         }
         cat("\n")
      }
      cat("> ")
      line <- readLines(f, 1)
      cat("\n")

      if(is.null(check)) {
         if(is.null(default) || nzchar(line)) {
            return(line)
         }else {
            return(NA)
         }
      }else {
         out <- check(line)
         if(out$passed) {
            return(out$val)
         }else if(msg %in% names(out)){
            cat(msg, sep = "\n")
         }else {
            cat(failed_msg, "\n")
         }
      }
   }
}

selectPrompt <- function(choices) {
   if(length(choices) == 1)
      return(choices)

   indices <- seq(length(choices))
   labels <- choices
   if(!is.null(names(choices)))
      labels <- names(choices)

   cat("\n")
   cat("Select one:", sep = "\n",
      paste0("[", indices, "] ", labels, collapse = "\n"))

   check <- function(input) {
      index <- as.integer(input)
      passed <- index %in% indices
      val <- NA
      if(passed) {
         val <- choices[[index]]
      }

      list(
         passed = passed,
         val = val
      )
   }

   prompt("", check)
}
