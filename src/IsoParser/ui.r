# TMT Mosaic
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

## Load Sources
source("lib/RowBuilder.r")

dashboardPage(skin = "blue",
   dashboardHeader(title = "TMT Mosaic",
      tags$li(uiOutput("idUI"), class = "dropdown")),
      dashboardSidebar(
         sidebarMenu(id = "sidebarMenu",
         menuItem("Protein Summary", tabName = "summary", icon = icon("search", lib = "glyphicon")),
         shinyjs::hidden(menuItem("Peptides", tabName = "peptides", icon = icon("link", lib = "glyphicon"))),
         shinyjs::hidden(menuItem("Interactions", tabName = "bioplex", icon = icon("project-diagram"))),
         menuItem("PCA", tabName = "pca", icon = icon("record", lib = "glyphicon")),
         menuItem("Volcano", tabName = "volcano", icon = icon("fire", lib = "glyphicon")),
         menuItem("Correlation", tabName = "correlation", icon = icon("stats", lib = "glyphicon")),
         shinyjs::hidden(menuItem("GO Plotter", tabName = "go", icon = icon("table"))),
         menuItem("Edit", tabName = "tmtBarChart", icon = icon("paint-brush")),
         menuItem("Advanced", tabName = "advanced", icon = icon("lock")),
         menuItem("Dataset", tabName = "dataset", icon = icon("table")),
         shinyjs::hidden(menuItem("Contact", tabName = "contact", icon = icon("address-book")))
      )
   ),
   dashboardBody(
      tags$head(
         shinyjs::useShinyjs(),
         includeCSS("../common/www/shinyProteomics.css"),
         includeScript("www/utils.js"),
         includeCSS("www/parser.css"),
         uiOutput("navbarCSS")
      ),
      progressBar(id = "progressBar1", value = 0, status = "plum", size = "xs"),
      shinyjs::hidden(
         tabItems(
            Build_PlotTabs_Summary(),
            Build_PlotTabs_Peptides(),
            Build_PlotTabs_Bioplex(),
            Build_PlotTabs_Pca(),
            Build_PlotTabs_Volcano(),
            Build_PlotTabs_Correlation(),
            Build_PlotTabs_GO_Plotter(),
            Build_TmtBarchart_Plot(),
            Build_Advanced_Tab(),
            Build_Dataset_Tab(),
            Build_Contact_Tab()
         )
      )
   )
)
