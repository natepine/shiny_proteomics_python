source("scripts/cliPrompts.r")

if(file.exists("data/conf.yml")) {
   cat("\nNote conf.yml already exists\n\n")
}

ping_gfy <- function(url, user, key) {
   handler <- function(e) {
      message("\nNetwork error:\n", e)
      return(NA)
   }

   allowSSLConf <- httr::config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE)
   tryCatch({
      rsp <- GET(url, authenticate(user, key), config = allowSSLConf)
      status_code(rsp)
   }, warning = handler, error = handler)
}

getServer <- function() {
   ask <- list(
      URL   = function(){prompt("Enter the URL of your MassPike instance")},
      User  = function(){prompt("MassPike admin username")},
      Key   = function(){prompt("User's API key")},
      Label = function(){prompt("Server display name", default = "URL")},
      AllowSSL = function(){boolPrompt("Allow self-signed SSL certificates")}
   )
   url <- ask$URL()
   data <- list(
      User = ask$User(),
      Key = ask$Key()
   )
   label <- ask$Label()
   if(!is.na(label))
      data$Label <- label

   if(curl::has_internet()) {
      repeat{
         query <- paste0("https://", url, "/gfy/www/modules/api/v1/")
         status <- ping_gfy(query, data$User, data$Key)
         issues <- c()
         if(is.na(status)) {
            cat("Failed to contact", query)
            issues <- c("URL")
         }else if(status != 200) {
            if(status == 401) {
               cat("User is not authorized")
               issues <- c("User", "Key")
            }else if(status == 404) {
               cat("HTTP GET to", query , "returned 404, not found")
               issues <- c("URL")
            }else {
               cat("HTTP GET to", query , "returned error code", status)
               issues <- c("URL", "User", "Key")
            }
         }else {
            cat("MassPike authentication successful\n")
            break
         }
         cat("\n")

         names(issues) <- paste("Edit", issues)
         issues <- c(issues, "Continue anyway" = "cont")
         pick <- selectPrompt(issues)
         if(pick == "cont") {
            break
         }else {
            val <- ask[[pick]]()
            if(pick == "URL")  url          <- val
            else               data[[pick]] <- val
         }
      }
   }else {
      cat("No internet connection, skipping authentication check\n")
   }

   data$VerifySSL <- ! ask$AllowSSL()

   server <- list(data)
   names(server) <- url
   return(server)
}

getSpecies <- function() {
   tax_id <- prompt("Enter an NCBI taxon id", check = intCheck, default = "No annotations")
   label <- prompt("Enter a display name for this organism")
   data <- list(tax_id)
   names(data) <- label
   return(data)
}

f <- file("stdin")
open(f, blocking=TRUE)

suppressPackageStartupMessages({
   library(yaml)
   library(stringr)
   library(httr)
})

cat("Generating conf.yml\n")

conf <- list()

header("MassPike Servers")
repeat {
   conf$SERVERS <- c(conf$SERVERS, getServer())
   if(!boolPrompt("Add another server"))
      break
}

header("Superuser")
if(!boolPrompt("Grant all users access to other users' data")) {
   su_pass <- prompt("Create superuser password")
   if(nzchar(su_pass))
      conf$SUPER_USER_PASSWORD <- su_pass
}

header("Species Annotations")
if(boolPrompt("Enable annotation enrichment analysis")) {
   cat("An annotation database will be built for each taxon id, enabling the GO enrichment UI\n")
   repeat {
      conf$ANNOTATED_SPECIES <- c(
         conf$ANNOTATED_SPECIES, getSpecies()
      )
      if(!boolPrompt("Add another species"))
         break
   }
   if(!any(is.na(conf$ANNOTATED_SPECIES))) {
      cat("Each species has a taxon id\n")
      if(boolPrompt("Add generic 'Other' species"))
         conf$ANNOTATED_SPECIES$Other <- NA
   }
}

header("Reverse Proxy")
if(boolPrompt("Will this webapp be deployed behind a reverse proxy")) {
   cat("Generated URLs will need to point to the configured address\n")
   subdir <- prompt("Enter the viewer's URL path")
   if(nzchar(subdir)) {
      conf$VIEWER_SUBDIR <- subdir
   }
}

header("Contact Info")
admin_email <- prompt("TMT Mosaic admin email", default = "No contact info")
if(!is.na(admin_email))
   conf$ADMIN_EMAIL <- admin_email

header("Current config")
cat("----------------------\n")
write_yaml(conf, stdout())
cat("----------------------\n")

path <- paste0(getwd(), "/data")
repeat {
   if(dir.exists(path)) {
      path <- paste0(normalizePath(path), "/conf.yml")
      if(!file.exists(path) || boolPrompt("Overwrite ", path)) {
         show_err <- function(e) {
            cat("Failed to save config:",
                e$message,
                "\n", sep = "\n")
            FALSE
         }

         saved <- tryCatch({
            write_yaml(conf, path)
            TRUE
         }, error = show_err, warning = show_err)

         if(saved)
            break
      }
   }else {
      cat(path, "does not exist\n")
   }
   path <- prompt("Path to 'data' dir")
}

cat("\nConfig saved!\n")
flush.console()

