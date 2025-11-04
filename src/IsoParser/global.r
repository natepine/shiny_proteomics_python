app_dir <- Sys.getenv("APP_DIR")
source(paste0(app_dir, "/common/Globals.r"))

libraries(c(
    "visNetwork",
    "ggdendro",
    "htmlwidgets",
    "openxlsx"
))

options("openxlsx.maxWidth" = 50)

downloadables <- list(
    summary = list(
        `TMT Barchart` = "tmt",
        Nearest = "nearest",
        Heatmap = "heatmap"
    )
)
