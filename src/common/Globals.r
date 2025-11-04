app_dir <- Sys.getenv("APP_DIR")

src_common <- function(file) source(paste(app_dir, "common", file, sep = "/"))

# Load shared libraries
libraries <- function(pkgs) lapply(pkgs, library, character.only = TRUE)
libraries(c(
   "shiny",
#   "shinyjs", namespace collision with colourpicker -> use shinyjs::
   "shinyWidgets",
   "shinydashboard",
   "DT",
   "dplyr",
   "RColorBrewer",
   "rhandsontable",
   "colourpicker",
   "ggplot2",
   "plotly",
   "stringr",
   "jsonlite",
   "RSQLite"
))

#Source shared helper functions
src_common("Html_Handlers.r")
src_common("Outputs.r")
src_common("MosaicDF.r")
src_common("Formats.r")
src_common("formats/Masspike.r")
src_common("sources/Masspike.r")
src_common("sources/Csv.r")
src_common("Configs.r")
src_common("Cache.r")
rm(src_common)

### Server Globals
dataPath <- paste0(app_dir, "/data/")
sqlitePath <- paste0(dataPath, "database.db")

# Load server config
config <- loadConfig(paste0(dataPath, "conf.yml"))
config$server_choices <- with(config, {
   choices <- sapply(names(SERVERS), function(server) {
      label <- SERVERS[[server]]$Label
      if(is.na(label))
         label <- server
      return(label)
   })
   setNames(names(choices), choices)
})

# Load shared data
loadIdMap <- function(taxid){
   if(!is.na(taxid)) {
      path <- paste0(dataPath, "idMaps/", taxid, ".rds")
      if(file.exists(path)) {
         return(readRDS(path))
      }
   }
   return(NULL)
}

idMaps <- list()
if(!identical(config$ANNOTATED_SPECIES, NA)) {
   idMaps <- sapply(config$ANNOTATED_SPECIES, loadIdMap)
}

# Database definitions
JSON_Vars <- c(
   # Column_name   = as.type
   "ColumnClasses" = as.character,
   "GroupNames"    = as.character,
   "GroupColors"   = as.character,
   "GroupIDs"      = as.integer,
   "ColumnNames"   = as.character,
   "ColumnIDs"   = as.integer
)

#Set default themes
theme_set(
   theme_classic() +
   theme(
      axis.line            = element_blank(),
      panel.grid.major.y   = element_line(colour = "grey"),
      text                 = element_text(family = "serif"),
      axis.text.x          = element_text(angle = 90, hjust = 1),
      plot.title           = element_text(hjust = 0.5),
      axis.title           = element_text(face = "bold"),
      plot.margin          = unit(c(0.2, 0.2, 0.2, 0.2), "cm")# c(Top, Right, Bottom, Left)
))

plotlyDefaults <- function(plot, ...){
   config(plot, displayModeBar = FALSE) %>%
   layout(
      xaxis = list(fixedrange = TRUE),
      yaxis = list(fixedrange = TRUE),
      ...)
}

datatableTheme <- function(bgColor = "#3B3B3B"){
   JS(paste0(
      "function(settings, json) {",
         "$(this.api().table().header()).css({'background-color': '", bgColor, "', 'color': '#fff'});",
      "}"
   ))
}

tooltipDigits <- 4

version <- "v4.0"
contact <- NA
if(!is.na(config$ADMIN_EMAIL))
   contact <- sprintf("Contact administrator at [%s] with questions/issues", config$ADMIN_EMAIL)

