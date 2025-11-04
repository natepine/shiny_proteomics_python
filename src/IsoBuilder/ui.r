#IsoBuilder
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

colorPresets <- c(
   "Basic",
   "Steve",
   "Pastel",
   "Rainbow",
   "Warm",
   "Earth",
   "Vibrant",
   "Light",
   "Neon",
   "Dark",
   "Jade",
   "Beach ball"
)

tabCount <- 0
tab <- NULL
nextTab <- function() {
   tabCount <<- tabCount + 1
   tab <<- navTabs[tabCount]
   return(tab)
}

navUI <- function(enableBack = TRUE) {
   submitBttn <- actionButton(paste0("submit_", tab), navBttns[[tab]], class = "bttn")
   if(tabCount != 1) {
      submitBttn <- shinyjs::disabled(submitBttn)
   }
   if(enableBack){
      out <- list(
         column(6, align = "left",  actionButton(paste0("back_", tab), "Back", class = "bttn")),
         column(6, align = "right", submitBttn)
      )
   }else{
      out <- centeredCol(submitBttn)
   }
   return(out)
}

navBoxTab <- function(width,  ...) {
   tabItem(tabName = nextTab(),
      centeredCol(width = width,
         defaultBox(title = names(tab),
            ...
         )
      )
   )
}

species_ui <- ""
if(!identical(config$ANNOTATED_SPECIES, NA)) {
   species_ui <- radioGroupButtons(
      inputId = "species",
      label = "Species",
      choices = names(config$ANNOTATED_SPECIES),
      justified = TRUE,
      status = "danger"
   )
}

ui <- dashboardPage(skin = "blue",
   dashboardHeader(title = "TMTEditor",
                  uiOutput("username", class = "username")),
   dashboardSidebar(
      sidebarMenu(id = "sidebarMenu",
         lapply(names(navTabs), function(navTabName) {
            menuItem(navTabName, tabName = navTabs[[navTabName]], icon = icon("circle"))
         }),

         menuItem("My Viewers",  tabName = "myViewers", icon = icon("list"))
      ),
      collapsed = TRUE
   ),
   dashboardBody(
      tags$head(
         shinyjs::useShinyjs(),
         includeScript("www/utils.js"),
         includeCSS("../common/www/shinyProteomics.css"),
         includeCSS("www/builder.css")
      ),
      progressBar("progress", value = 0, size = "xs", total = 6, status = "primary", display_pct = FALSE),
      tabItems(
         navBoxTab(8, # Quant ID
            br(),
            tabsetPanel(id = "data_type",
               #type = 'hidden',
               tabPanel("Server",
                  centeredCol(
                     textInput("QID",
                        label = "ProteinQuant ID",
                        value = NULL,
                        width = "150px",
                        placeholder = "i.e. 16097"
                     ),
                     prettyCheckbox("ignore_cache",
                        "Refresh Cache",
                        value = FALSE,
                        shape = "round",
                        status = "danger")
                  )
               ),
               tabPanel("CSV/TSV",
                  centeredCol(
                     div(tableOutput("csvFormat"),    style = "overflow-x: scroll"),
                     fileInput("quantCSV",    "Protein Quant",            accept = c("text/csv", ".tsv"), width = "250px"),
                     div(tableOutput("csvPepFormat"), style = "overflow-x: scroll"),
                     fileInput("pepQuantCSV", "Peptide Quant (Optional)", accept = c("text/csv", ".tsv"), width = "250px")
                  )
               )
            ),
            radioGroupButtons("is_site_quant",
               label = NULL,
               choices = list("Protein" = FALSE, "Site" = TRUE),
               justified = TRUE,
               status = "primary",
               checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = NULL),
               width = "200px"
            ),
            species_ui,
            br(),
            navUI(enableBack = FALSE)
         ),
         navBoxTab(10, # Names
            fluidRow(style = "text-align: center;",
               column(4,
                  h3("Requirements"),
                  tags$span(style = "font-size: 15px",
                     uiOutput("tableValidation")
                  )
               ),
               column(3,
                  "Number of Groups",
                  fluidRow(
                     column(4, style = "padding:0px", textInput("num_groups", NULL, 2)),
                     column(4, style = "padding:0px", actionButton("setGroups", "Set")),
                     column(1, style = "padding:0px", checkboxInput("seqGroups", "Seq", TRUE))
                  ),
                  fluidRow(
                     actionButton("trimNames", "Trim"),
                     actionButton("makeUnique", "Make Unique"),
                     actionButton("autoName", "Auto Name")
                  ),
                  fluidRow(
                     column(3, style = "padding:0px", actionButton("sortTable", "Sort")),
                     column(5, style = "padding:0px", selectInput("sortBy", NULL, choices = c("Order", "Class", "Name", "Group"), selected = "Name")),
                     column(4, style = "padding:0px", checkboxInput("sortRows", "As rows", TRUE))
                  )
               ),
               column(5,
                  h3("Editing"),
                  tags$span(style = "font-size: 15px",
                     HTML(keyValueStr(sep = " - ", boldKeys = TRUE, c(
                        "Spreadsheet Data" = "Copy and paste from any spreadsheet directly into the cells.",
                        "Reordering samples" = "Click a row to select it, drag and drop to reorder. (Shift-click for multiple rows)"
                     )))
                  )
               )
            ),
            tags$hr(style = "border-color: black;"),
            fluidRow(
               column(7, align = "left",
                  rHandsontableOutput("column_table")
               ),
               column(5, align = "left",
                  rHandsontableOutput("group_table"),
                  tags$hr(style = "border-color: black;"),
                  rHandsontableOutput("class_table"),
                  tags$br(),
                  centeredCol(
                     defaultBox(title = "Plot Preview", collapsible = TRUE,
                        checkboxInput("are_reps", "Perfect Replicates", FALSE),
                        plotlyOutput("plot_preview_names", height = 305)
                     )
                  )
               )
            ),
            column(4, offset = 4,
               navUI()
            )
         ),
         navBoxTab(10, # Colors
            fluidRow(
               column(9, align = "center",
                  br(),
                  plotlyOutput("protein_plot", width = "90%", height = 500),
                  tags$hr(),
                  column(4, offset = 4,
                     navUI()
                  )
               ),
               column(3, align = "center",
                  uiOutput("ui_show_as_replicates"),
                  h4("Default ID:"),
                  uiOutput("protChoiceUI"),
                  tags$hr(),
                  h4("Group Colors:"),
                  selectInput("colorPalPreset", "Preset", multiple = FALSE, choices = colorPresets),
                  uiOutput("groupColorUI")
               )
            )
         ),
         navBoxTab(4, #Notes
            h4("Describe your datset"),
            br(),
            centeredCol(textInput("datasetName", "Dataset Name:", placeholder = "ex. MCT vs WT")),
            centeredCol(
               textAreaInput("notes", "Notes:", placeholder = paste("**REQUIRED** Please provide detailed notes.",
                  "They will provide context for the experiment and help identify the viewer later.\nex. 3 rep x 3 rep x 4 rep, 11 plex..."),
                  height = 100, resize = "both"),
               br()
            ),
            navUI()
         ),
         tabItem(tabName = nextTab(), # Finalize
            defaultBox(title = names(tab), width = 8,
               h4("Please review the following information before continuing. Ensure that all fields look correct."),
               br(),
               fluidRow(
                  column(6, htmlOutput("review_information_1", style = "font-size: 18px")),
                  column(6, htmlOutput("review_information_2", style = "font-size: 18px"))
               ),
               tags$hr(),
               navUI()
            ),
            defaultBox(title = "Plot Preview", width = 4,
               plotlyOutput("plot_review", height = 305)
            ),
            shinyjs::hidden(tags$span(id = "finished_box",
               defaultBox(title = "Viewer Link", status = "success",
                  uiOutput("viewer_link", style = "font-size: 20px;text-align: center;")
                  # ,
                  # column(11, offset = 1,
                  #    fluidRow(
                  #       textInput("useremail", "Email:"),
                  #       actionButton("sendemail", "Send viewer link to email")
                  #    ),
                  #    h6("NOTE: This will send an email from the address 'noreply.tmtmosaic@gmail.com' or 'noreply.tmtphospho@gmail.com', and as such it may end up in Spam (check there and be sure to mark as 'not spam').")
                  # )
               )
            ))
         ),
         if(!is.na(config$SUPER_USER_PASSWORD)){
            tabItem(tabName = "myViewers",
               passwordInput("adminpassword", "Admin:", placeholder = "Show viewers from all users"),
               centeredCol(style = "font-size:80%", dataTableOutput("admintable", width = "100%")))
         }else{
            tabItem(tabName = "myViewers",
               centeredCol(style = "font-size:80%", dataTableOutput("admintable", width = "100%")))
         }
      )
   )
)

footer <- version
if(!is.na(config$ADMIN_EMAIL))
   footer <- paste(sep = " | ", footer, contact)
ui <- div(id = "page", ui, tags$footer(footer))

set_labels(
  language = "en",
  "Please authenticate" = "CORE Login",
  "Username:" = "Username:",
  "Password:" = "Password:",
  "You are not authorized for this application" = "Login Failed"
)

choices <- config$server_choices
nServers <- length(choices)
if(nServers > 1){
   width <- nServers * (max(nchar(names(choices))) + 3) * 9 + 30
   width <- max(width, 200)
   tags_bottom <-
      centeredCol(
         radioGroupButtons(
            inputId = "server", label = "Server:",
            choices = choices,
            justified = TRUE, status = "success",
            width = paste0(width, "px"),
            checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = NULL)
         )
      )
} else {
   server <- names(choices)
   tags_bottom <- 
      column(12, align = "center",
         tags$h4(paste("Server:", server))
      )
}

secure_app(ui, background = "linear-gradient(rgba(56, 125, 140, 0.5),rgba(56, 125, 140, 0.5))", tags_bottom = tags_bottom)

