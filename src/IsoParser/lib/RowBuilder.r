### Build specific rows within Dashboard body
p_adjust_choices <- list("None" = "none", "BH/FDR" = "fdr", "Bonferroni" = "bonferroni")
amino_acids <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
fcCutoffs <- c(1.1, 1.25, 1.5, 1.75, 2, 2.5, 3, 5)
pvalCutoffs <- c(.00001, .0001, .001, .005, .01, .05, .1, 1)
plotFormats <- c("PNG", "JPG", "SVG", "TIFF", "PDF")

moreInfo <- function(id) {
   actionButton(id, label = "",
      icon = icon("question"),
      style="position: absolute;right: 10px;top: 10px; border-radius: 30px;")
}

Build_ProteinHeader_FluidRow <- function(){
   fluidRow(
      column(3, align = "center", 
         uiOutput("ui_select_protein"),
         htmlOutput("current_seq", style = "font-size: 16pt; overflow-x: scroll;")
      ),
      column(6, align = "center",
         fluidRow(
            column(5, align = "center", radioGroupButtons(
               inputId = "numNearest", label = "Number of Closest Proteins", 
               choices = c(10, 15, 20, 25, 30, 50, 100), 
               justified = TRUE, status = "primary",
               selected = 15, size = "xs"
            )),
            column(3, align = "center", radioGroupButtons(
               inputId = "distance_func", label = "Dist. Function",
               choices = c(
                  "Euc" = "euc",
                  "Corr" = "corr",
                  "Cos" = "cos"
               ),
               justified = TRUE, status = "primary", size = "xs"
            )),
            column(4, align = "center", br(), actionButton("button_update_parameters", "Update Parameters", class = "updatebttn"))
         ),
         fluidRow(
            column(4, align = "center", selectInput("goMinMatches", "Min. Num GO Matches:", 1:10, 2)),
            column(4, align = "center", selectInput("goPvalCutoff", "PVal Cutoff:", choices = pvalCutoffs, .01)),
            column(4, align = "center", selectInput("goPvalAdjust", "PVal Correction", choices = p_adjust_choices, selected = "fdr"))
         ),
         style='border-left:1px solid; border-right:1px solid; border-color:#497EA5'
      ),
      column(3, align = "center",
         fluidRow(
            centeredCol( 
               selectInput("heatmap_fc", "Heatmap FC Cutoff:", choices = fcCutoffs, selected = 2),
               # br(), 
               dropdownButton(
                  fluidRow(
                     column(6, colourInput("heatmap_lower_col", "Lower Color:", "#002fffff", showColour = "background")),
                     column(6, colourInput("heatmap_upper_col", "Upper Color:", "#FF00C9", showColour = "background"))
                  ),
                  radioGroupButtons(
                     inputId = "summary_heatmap_type", label = "Graph Choice:", 
                     choices = c("Heatmap", "Dotplot"),
                     justified = TRUE, status = "success",
                     selected = "Heatmap", size = "xs"
                  ),
                  prettyCheckbox("heatmap_use_group_colors", "(Dotplot) Group Colors", FALSE),
                  conditionalPanel("input.extend_plots",
                     sliderInput("heatmap_size", "Heatmap size:", min = .3, max = 1, step = .01, value = 1, ticks = FALSE)
                  ),
                  label = "Advanced Heatmap",
                  width = 300,
                  circle = FALSE
               )
            )
         )
      )
   )
}

Build_PatternHeader_FluidRow <- function(){
   fluidRow(
      column(5,
         fluidRow(
            tags$span(
               style = "font-size:11.5pt;float:left;text-align:right;",
               HTML(paste(
                  "Group Up (&uarr;)",
                  "Group Mid (-)",
                  "Group Down (&darr;)",
                  "Use",
                  sep = "<br>"))),
            column(9, align = "left", class = "scroll-x",
               fluidRow(
                  uiOutput("pattern_ui_up", inline = TRUE),
                  actionButton("button_reset_pattern", "Reset", class = "minibttn updatebttn", style = "background-color:#E5C2C2;width:40px;")
               ),
               fluidRow(
                  uiOutput("pattern_ui_mid", inline = TRUE)
               ),
               fluidRow(
                  uiOutput("pattern_ui_down", inline = TRUE),
                  actionButton("button_search_pattern", "Search", class = "minibttn updatebttn", style = "width:40px;")
               ),
               fluidRow(
                  uiOutput("pattern_ui_use", style = 'padding:2px 5px;height:25;')
               )
            )
         )
      ),
      column(5, align = "center",
         fluidRow(
            tags$span(style = "float:left;text-align:right;",
               HTML(paste(
                  "&uarr;",
                  "-",
                  "&darr;",
                  sep = "<br>"
               ))
            ),
            column(11, align = "left", class = "scroll-x",
               uiOutput("pattern_sliders", inline = TRUE)
            )
         ),
         checkboxInput("allSliders", "All sliders", value = FALSE)
      ),
      column(2,
         dropdownButton(
            radioButtons(
               inputId = "summary_plot_choice",
               label = "Plot:",
               choices = c(downloadables$summary, All = "all")
            ),
            textInput("summary_plot_filename", "Filename:", placeholder = "(leave blank for default)"),
            radioGroupButtons(
               inputId = "summary_plot_format", label = "File Format:",
               choices = plotFormats,
               justified = FALSE,
               status = "info",
               size = "xs"
            ),
            downloadButton("summary_plot_download", "Download"),
            label = "Plots",
            status = "secondary",
            right = TRUE,
            icon = icon("download", lib = "glyphicon"),
            width = 300,
            circle = FALSE
         ),
         dropdownButton(
            radioButtons(
               inputId = "summary_data_choice",
               label = "Data:",
               choices = downloadables$summary
            ),
            textInput("summary_data_filename", "Filename:", placeholder = "(leave blank for default)"),
            radioGroupButtons(
               inputId = "summary_data_format", label = "File Format:",
               choices = c("CSV", "TSV"),
               justified = FALSE,
               status = "info",
               size = "xs"
            ),
            downloadButton("summary_data_download", "Download"),
            label = "Plot Data",
            status = "secondary",
            right = TRUE,
            icon = icon("download", lib = "glyphicon"),
            width = 300,
            circle = FALSE
         ),
         checkboxInput("extend_plots", "Extend Plots", value = FALSE)
      )
   )
}

Build_PlotTabs_Summary <- function() {
   tabItem(tabName = "summary",
      defaultBox(
            Build_ProteinHeader_FluidRow()
      ),
      defaultBox(
            Build_PatternHeader_FluidRow()
      ),
      centeredCol(
         uiOutput("summaryPlots"),
         span(style = "font-size:80%",
            dataTableOutput("table_nearest_proteins"),
            dataTableOutput("table_nearest_go"),
            shinyjs::hidden(fluidRow(id = "GO_header",
               column(6, align = "center",
                  downloadButton("go_download", "GO Table")
               ),
               column(6, align = "center",
                  tags$h4("*Hover SubsetHits to see hits")
               )
            ))
         ),
         br()
      )  )
}

Build_PlotTabs_Peptides <- function(){
   tabItem(tabName = "peptides",
      defaultBox(status = "warning",
         fluidRow(
            column(4, align = "center", plotlyOutput("plot_peptides", width = "100%", height = 300)),
            column(4, align = "center", plotlyOutput("plot_selected_peptide", width = "100%", height = 300)),
            column(4, align = "left",
               tags$div(textOutput("sequenceTitle"), class = "plot-title", style = "margin: 5.6px 0;"),
               fluidRow(
                  uiOutput("ui_peps_choice")
               ),
               fluidRow(
                  column(2, align = "right", htmlOutput(outputId = "current_peptide_index", style="font-size: 12px; font-family: \"Lucida Console\", Monaco, monospace;")),
                  column(10, align = "left", htmlOutput(outputId = "current_peptide_sequence", style="font-size: 12px; font-family: \"Lucida Console\", Monaco, monospace;", class = "scroll-x"))
               )
            )
         )
      ),
      fluidRow(
         column(3, align = "center",
            fluidRow(style = "margin-bottom: -10px",
               column(6, align = "center",
                  checkboxInput("normalize_peps", "Show RA", value = TRUE)               ),
               column(6, align = "center", checkboxInput("peptide_fit_y", "Fit to Values", value = FALSE))
            ),
            shinyjs::hidden(checkboxInput("peps_by_uniprot", "Full Protein", value = FALSE))
         ),
         column(6, align = "center",
            fluidRow(
               column(2, colourInput("residueA_col", "Color A", "#008CFF", showColour = "background")),
               column(4, selectInput("residueGroupA", "Residue Group A", choices = amino_acids, multiple = TRUE, selected = c("R", "K"))),
               column(4, selectInput("residueGroupB", "Residue Group B", choices = amino_acids, multiple = TRUE, selected = "C")),
               column(2, colourInput("residueB_col", "Color B", "#FF9900", showColour = "background"))
            )
         ),
         column(3, align = "center",
            checkboxInput("highlight_residues", "Highlight Residues")
         )
      ),
      centeredCol(style = "font-size:90%", dataTableOutput("table_peptides"))
   )
}

Build_PlotTabs_Bioplex <- function(){
   tabItem(tabName = "bioplex",
      fluidRow(
         defaultBox(width = 3,
            uiOutput("ui_select_protein_bioplex"),
            tags$hr(),
            radioGroupButtons(
               inputId = "bioplex_graph_type", label = "Graph Choice", 
               choices = list(Current = "current", Nearest = "nearest", All = "all", `2` = "2", `4` = "4", `8` = "8"), 
               justified = TRUE, status = "success",
               selected = "current", size = "xs"
            ),
            radioGroupButtons(
               inputId = "bioplex_graph_layout", label = "Graph Layout", 
               choices = list(Circle = "layout_in_circle", Star = "layout_as_star", Frucht = "layout_with_fr", Nice = "layout_nicely"), 
               justified = TRUE, status = "success",
               selected = "layout_nicely", size = "xs"
            ),
            tags$hr(),
            colourInput("bioplex_current_color", "Current Color:", "#D399F0", showColour = "background"),
            colourInput("bioplex_nearest_color", "Nearest Color:", "#FF9100", showColour = "background"),
            colourInput("bioplex_bait_color", "Bait Color (not Nearest):", "#13F73D", showColour = "background"),
            colourInput("bioplex_other_color", "Other Color:", "#F7FFA1", showColour = "background"),
            colourInput("bioplex_na_color", "N/A Color (not in Bioplex):", "#E2E7EB", showColour = "background"),
            tags$hr(),
            radioGroupButtons(
               inputId = "bioplex_format", label = "File Extension (for download):", 
               choices = list("PNG", "JPEG", "PDF"), 
               justified = TRUE, status = "success",
               selected = "png", size = "xs"
            )
         ),
         defaultBox(width = 9,
            visNetworkOutput("plot_bioplex")
         )
      )
   )
}

Build_PlotTabs_Pca <- function(){
   tabItem(tabName = "pca",
      fluidRow(
         defaultBox(width = 6,
            span(textOutput("pca_header"), style = "font-size: 24px;"),
            plotlyOutput("plot_pca"),
            br(),
            fluidRow(
               column(8,
                  div("Loadings", style = "font-size: 24px;", class = "centered"),
                  dataTableOutput("pca_loadings")
               ),
               column(4,
                  selectInput("pca_x", "X PC", choices = 1),
                  selectInput("pca_y", "Y PC", choices = 2),
                  fluidRow(
                     column(6, downloadButton("pca_download", "Download")),
                     column(6, dropdownButton(
                        label = "Plot Settings",
                        circle = FALSE,
                        right = TRUE,
                        width = "230px",
                        textInput("pca_x_lab", "Label x:", NULL),
                        textInput("pca_y_lab", "Label y:", NULL),
                        fluidRow(
                           column(6,
                              numericInput("pca_x_min", "Min x:", NULL),
                              numericInput("pca_y_min", "Min y:", NULL)
                           ),
                           column(6,
                              numericInput("pca_x_max", "Max x:", NULL),
                              numericInput("pca_y_max", "Max y:", NULL)
                           )
                        ),
                        sliderInput("pca_point_size", "Point Size:", min = 1, max = 15, value = 5, step = 1),
                        hr(),
                        h4("Download Parameters"),
                        fluidRow(
                           column(6, numericInput("pca_width", "Width:", min = 50, max = 10000, NULL)),
                           column(6, numericInput("pca_height", "Height:", min = 50, max = 10000, NULL))
                        )
                     )
                  ))
               )
            )
         ),
         column(width = 6,
            defaultBox(
               title = "Column Dendrogram",
               plotOutput("full_dendro", height = NULL)
            )
         )
      ),
      fluidRow(
         defaultBox(title = "Column Histograms", collapsible = TRUE, collapsed = TRUE,
            shinyjs::hidden(column(6, id = "histogrid_plex_column",
               selectInput("histogrid_plex", "Select a Plex", choices = NULL, width = 150)
            )),
            column(6, checkboxInput("overlay_histogrid", "Summarize Histograms")),
            uiOutput("histogrid")
         )
      )
   )
}

Build_PlotTabs_Volcano <- function(){
   tabItem(tabName = "volcano",
      includeScript("www/volcano_click.js"),
      defaultBox(status = "danger",
         fluidRow(
            column(3, align = "center",
               h4("FC = Group B / Group A"),
               fluidRow(
                  column(6, align = "center", selectInput("volcano_grp_a", "Group A:", selected = NULL, choices = 1)),
                  column(6, align = "center", selectInput("volcano_grp_b", "Group B:", selected = NULL, choices = 2))
               ),
               centeredCol(actionButton("button_update_volcano", "Update Plot", style="width:150px;", class = "updatebttn"))
            ),
            column(2, align = "center",
               selectInput(
                  inputId = "volcano_x_cutoff",
                  label = "X Cutoff (FC):", 
                  choices = fcCutoffs, 
                  selected = 1.75
               ),
               selectInput(
                  inputId = "volcano_y_cutoff",
                  label = "Y Cutoff (pval):", 
                  choices = pvalCutoffs,
                  selected = .01
               ),
               style='border-right:1.3px solid; border-left:1.3px solid; border-color:#E74C3C'
            ),
            column(2, align = "left",
               centeredCol(
                  colourInput("volcano_annotation_color", "Select annotation color", "#1E90FF", palette = "limited", showColour = "background"),
                  br()
               ),
               centeredCol(uiOutput("volcano_annot_ui"))
            ),
            column(3, align = "center", plotlyOutput("plot_volcano_groups", width = "100%", height = 200)),
            column(2, align = "center",
               tags$a(href="https://en.wikipedia.org/wiki/Multiple_comparisons_problem", target='_blank',
                  style="position: absolute;right: 20px;", title="Click to learn more", `data-toggle`="tooltip",
                  icon("question")
               ),
               selectInput("volcano_padjust", label = "P Adjust Method", choices = p_adjust_choices, selected = "none"),
               dropdownButton(
                  radioGroupButtons(
                     inputId = "volcano_plot_format", label = "Plot Download Type:",
                     choices = list(PNG = "png", SVG = "svg", JPEG = "jpeg"),
                     justified = TRUE, status = "success",
                     selected = "png", size = "xs"
                  ),
                  radioButtons(
                     inputId = "volcano_choice",
                     label = "Item to Download:",
                     choices = list(`Upregulated Table` = "upreg", `Downregulated Table` = "downreg", `Up + Down Table` = "upanddown", `Full Table` = "full"),
                     selected = "upreg"
                  ),
                  textInput("volcano_filename", "Filename:", placeholder = "(leave blank for default)"),
                  radioGroupButtons(
                     inputId = "volcano_table_format", label = "File Extension (for saving):",
                     choices = list(CSV = "csv", TSV = "tsv"),
                     justified = TRUE, status = "danger",
                     selected = "csv", size = "xs"
                  ),
                  downloadButton("volcano_download", "Download Table"),
                  label = "DOWNLOAD",
                  status = "danger",
                  size = "lg",
                  icon = icon("download", lib = "glyphicon"),
                  width = 300,
                  circle = FALSE
               ),
               tags$br(),
               prettyCheckbox("show_volcano_tables", "Show Protein Tables", value = TRUE, status = "danger")
            ),
            uiOutput("volcano_custom_annot_ui"),
            uiOutput("volcano_go_annot_ui")
         ),
         fluidRow(
            column(2, align = "left", offset = 1, colourInput("volcano_down_color", label = NULL, "#BE3333", showColour = "background")),
            column(4, align = "center", colourInput("volcano_go_bg_color", "Background Points:", value = "#E5E5E5", palette = "limited", showColour = "background")),
            column(2, align = "center", tags$p("*Click to toggle annotations, double click to view details")),
            column(2, align = "right", colourInput("volcano_up_color", label = NULL, "#40A02D", showColour = "background"))
         )
      ),
      centeredCol(plotlyOutput("plot_volcano")),
      fluidRow(
         column(6, align = "center", uiOutput("volcano_down_count")),
         column(6, align = "center", uiOutput("volcano_up_count"))
      ),
      conditionalPanel(condition = "input.show_volcano_tables == 1",
         column(12, align = "center", style = "font-size:80%", dataTableOutput("table_volcano_down")),
         column(12, align = "center", style = "font-size:80%", dataTableOutput("table_volcano_up")),
         shinyjs::hidden(span(id = "volcano_GO",
            defaultBox(status = "danger",
               fluidRow(
                  column(8, align = "center",
                     radioGroupButtons("volcano_go", "GO Enrichment", choices = c(
                        "Down Regulated" = "down",
                        "All Significant" = "all",
                        "Up Regulated" = "up"),
                        selected = character(0)) # No selection
                  ),
                  column(2, align = "center",
                     colourInput("volcano_go_color", "GO Color:", value = "#FF1493", palette = "limited")
                  ),
                  column(2, align = "center",
                     checkboxInput("volcano_go_show_labels", "Show Labels", value = FALSE)
                  )
               ),
               conditionalPanel(condition = "input.volcano_go !== null",
                  column(4, align = "center", selectInput("volcanoMinMatches", "Min. Num GO Matches:", 1:10, 2)),
                  column(4, align = "center", selectInput("volcanoPvalCutoff", "PVal Cutoff:", choices = pvalCutoffs, selected = .01)),
                  column(4, align = "center", selectInput("volcanoPvalAdjust", "PVal Correction", choices = p_adjust_choices, selected = "fdr")),
                  dataTableOutput("volcano_goTable"),
                  uiOutput("volcano_go_color_script"),
                  column(12, align = "center",
                     tags$p("Click any GO term button to annotate volcano plot | Use controls above to customize colors and labels", 
                           style = "font-size: 12px; color: #666; margin: 5px 0;")
                  ),
                  column(6, align = "center",
                     downloadButton("volcano_go_download", "GO Table")
                  ),
                  column(6, align = "center",
                     tags$h4("*Hover SubsetHits to see hits")
                  )
               )
            )
         ))
      )
   )
}

Build_PlotTabs_Correlation <- function(){
   tabItem(tabName = "correlation",
      includeScript("www/correlation_click.js"),
      defaultBox(status = "success",
         fluidRow(
            column(3, align = "center",
               h4("Comparison 1: FC = Group B / Group A"),
               fluidRow(
                  column(6, align = "center", selectInput("corr_grp1_a", "Group A:", selected = NULL, choices = 1)),
                  column(6, align = "center", selectInput("corr_grp1_b", "Group B:", selected = NULL, choices = 2))
               )
            ),
            column(3, align = "center",
               h4("Comparison 2: FC = Group D / Group C"),
               fluidRow(
                  column(6, align = "center", selectInput("corr_grp2_c", "Group C:", selected = NULL, choices = 1)),
                  column(6, align = "center", selectInput("corr_grp2_d", "Group D:", selected = NULL, choices = 3))
               )
            ),
            column(3, align = "center",
               br(),
               centeredCol(actionButton("button_update_correlation", "Update Plot", style="width:150px;", class = "updatebttn"))
            ),
            column(3, align = "center",
               dropdownButton(
                  radioGroupButtons(
                     inputId = "correlation_plot_format", label = "Plot Download Type:",
                     choices = list(PNG = "png", SVG = "svg", JPEG = "jpeg"),
                     justified = TRUE, status = "success",
                     selected = "png", size = "xs"
                  ),
                  textInput("correlation_filename", "Filename:", placeholder = "(leave blank for default)"),
                  downloadButton("correlation_download", "Download Plot"),
                  label = "DOWNLOAD",
                  status = "success",
                  size = "lg",
                  icon = icon("download", lib = "glyphicon"),
                  width = 300,
                  circle = FALSE
               )
            )
         ),
         fluidRow(
            column(12, align = "center", 
               tags$p("*Each point represents a protein. Click a point to view protein details. Regression line shows linear relationship.")
            )
         )
      ),
      centeredCol(plotlyOutput("plot_correlation")),
      fluidRow(
         column(6, align = "center",
            defaultBox(title = "Linear Regression Statistics", status = "success", collapsible = TRUE,
               dataTableOutput("correlation_stats_table")
            )
         ),
         column(6, align = "center",
            defaultBox(title = "Labeling Controls", status = "primary", collapsible = TRUE,
               numericInput("correlation_label_threshold", "Distance Threshold for Labels:", value = 1.5, min = 0, step = 0.1),
               tags$p("Points with distance ≥ threshold will be labeled on the plot.", style = "font-size:12px;"),
               tags$hr(),
               fluidRow(
                  column(6,
                     colourInput("correlation_annotation_color", "Annotation Color", "#1E90FF", palette = "limited", showColour = "background")
                  ),
                  column(6,
                     uiOutput("correlation_annot_ui")
                  )
               ),
               uiOutput("correlation_custom_annot_ui"),
               tags$hr(),
               tags$h5("Interactive Features:", style = "color: #0275d8; margin-bottom: 8px;"),
               tags$p(HTML("<strong>Hover:</strong> Show protein info tooltip<br>
                           <strong>Click:</strong> Select protein and switch to Protein tab"), 
                      style = "font-size:12px; margin-bottom: 5px;"),
               tags$p("💡 All points are interactive - hover for details!", style = "font-size:11px; color: #666; font-style: italic;")
            )
         )
      ),
      fluidRow(
         column(12, align = "center",
            uiOutput("correlation_significant_count")
         )
      ),
      fluidRow(
         column(12, align = "center",
            defaultBox(title = "All Correlation Data", status = "info", collapsible = TRUE,
               dataTableOutput("table_correlation_all")
            )
         )
      ),
      fluidRow(
         column(12, align = "center",
            defaultBox(title = "Export Data", status = "success", collapsible = TRUE,
               fluidRow(
                  column(6,
                     selectInput("correlation_data_format", "Format:",
                               choices = list("CSV" = "csv", "TSV" = "tsv"),
                               selected = "csv")
                  ),
                  column(6,
                     textInput("correlation_data_filename", "Filename:", placeholder = "(auto-generated)")
                  )
               ),
               fluidRow(
                  column(12, align = "center",
                     downloadButton("correlation_data_download", "Download Data", class = "bttn")
                  )
               )
            )
         )
      )
   )
}

Build_PlotTabs_GO_Plotter <- function(){
   tabItem(tabName = "go",
      defaultBox(
         fluidRow(align = "center", 
            column(6,
               selectInput("goPlotter_annotation", "Search for Annotation:", choices = NULL),
               fluidRow(
                  column(4, br(), actionButton("button_plot_go", "Generate Heatmap", style="width:\"100%\"; border-color:#AAAAAA; color:#000000; background-color:#CFE5C2")),
                  column(4, 
                     radioGroupButtons(
                        inputId = "go_heatmap_sort", label = "Sort Categories By:", 
                        choices = list(Name = "Name", `High-Low` = "Down", `Low-High` = "Up"), 
                        justified = TRUE, status = "info",
                        selected = "Name", size = "xs"
                     )
                  ),
                  column(4, selectInput("go_heatmap_fc", "Heatmap FC Cutoff:", choices = fcCutoffs, selected = 2))
               )
            ),
            column(3,
               sliderInput("go_heatmap_height", "Heatmap Aspect Ratio Factor (log2):", min = -1, max = 1, value = 0, step = 0.1, ticks = FALSE),
               checkboxInput("go_heatmap_show_nearest", "Highlight Nearest", TRUE)
            ),
            column(3,
               colourInput("go_heatmap_upper_col", "Upper Color:", "#FF00C9", showColour = "background"),
               colourInput("go_heatmap_lower_col", "Lower Color:", "#002fffff", showColour = "background")
            )
         ),
         fluidRow(
            column(6,
               p("NOTE: You may need to scroll down to see the full heatmap with a large n (n > 100)", style = "font-size:12pt")
            ),
            column(6, align = "center",
               downloadButton("go_plotter_download", "Download Plot Data")
            )
         )
      ),
      plotOutput("plot_go_heatmap", height = NULL)
   )
}

Build_TmtBarchart_Plot <- function(){
   tabItem(tabName = "tmtBarChart",
      defaultBox(
         column(1, 
            h4("Theme:")
         ),
         column(3, align = "center", 
            dropdownButton(
               uiOutput("theme_colors_ui"),
               centeredCol(
                  actionButton("reset_theme_colors", "Reset Colors",
                     style="width:200px;", class = "updatebttn")
               ),
               label = "Change Group Colors",
               width = 300,
               circle = FALSE
            )
         ),
         column(3, align = "center",
            dropdownButton(
               uiOutput("theme_names_ui"),
               centeredCol(
                  actionButton("reset_theme_names", "Reset Labels",
                     style="width:200px;", class = "updatebttn")
               ),
               label = "Change Sample Labels",
               width = 300,
               circle = FALSE
            )
         ),
         column(3, align = "center",
            dropdownButton(
               uiOutput("theme_group_names_ui"),
               centeredCol(
                  actionButton("reset_theme_group_names", "Reset Labels",
                     style="width:200px;", class = "updatebttn")
               ),
               label = "Change Group Labels",
               width = 300,
               circle = FALSE
            )
         ),
         column(2, align = "center",
            actionButton("apply_theme", "Apply Theme", class = "updatebttn")
         )
      ),
      column(3, align = "center", 
         dropdownButton(
            textInput("tmtBarChart_title", "Title:", placeholder = "(leave blank for default)"),
            textInput("tmtBarChart_title_x", "X-Axis Label:", placeholder = "(leave blank for default)"),
            textInput("tmtBarChart_title_y", "Y-Axis Label:", placeholder = "(leave blank for default)"),
            sliderInput("tmtBarChart_title_size", "Title Size:", min = 0, max = 24, step = 1, value = 14, ticks = FALSE),
            sliderInput("tmtBarChart_label_size", "Label Size:", min = 0, max = 24, step = 1, value = 12, ticks = FALSE),
            size = "lg",
            status = "info",
            icon = icon("search"),
            label = "Advanced Options",
            width = 300,
            circle = FALSE
         )
      ),
      column(3, align = "center",
         dropdownButton(
            textInput("tmtBarChart_filename", "Filename:", placeholder = "(leave blank for default)"),
            radioGroupButtons(
               inputId = "tmtBarChart_format", label = "File Format:",
               choices = plotFormats,
               justified = TRUE, status = "danger",
               selected = "png", size = "xs"
            ),
            downloadButton("tmtBarChart_download", "Download Plot"),
            label = "DOWNLOAD",
            status = "danger",
            size = "lg",
            icon = icon("download", lib = "glyphicon"),
            width = 300,
            circle = FALSE
         )
      ),
      column(3, align = "center",
         sliderInput("tmtBarChart_width", "Figure Width:", min = 100, max = 1000, value = 600)
      ),
      column(3, align = "center",
         sliderInput("tmtBarChart_height", "Figure Height:", min = 100, max = 1000, value = 400)
      ),
      uiOutput("ui_tmtBarChart")
   )
}

Build_Advanced_Tab <- function(){
   tabItem(tabName = "advanced",
      div(style = "text-align:center;color:red;",
         h2("WARNING: This tab will change Live Data for ALL plots"),
         h3("Do not use without a clear understanding of its effects")
      ),
      defaultBox(
         htmlOutput("norm_steps", style = "font-size:20px;")
      ),
      fluidRow(
         column(6,
            defaultBox(
               moreInfo("col_norm_info"),
               h2(HTML(bold("Column Normalization"))),
               h3("Default Normalization Factors"),
               plotlyOutput("plot_default_factors", height = 230),
               uiOutput("custom_col_norm_ui")
            )
         ),
         column(6,
            defaultBox(
               moreInfo("bridge_info"),
               h2(HTML(bold("Bridge Normalization"))),
               h4("Compute the ratio to the selected channel for a given class"),
               uiOutput("bridge_norm_ui")
            ),
            defaultBox(
               moreInfo("row_scaling_info"),
               h2(HTML(bold("Row Scaling"))),
               fluidRow(
                  column(6,
                     h4(id = "rescale_label", "Rescale Rows"),
                     switchInput("rescale", value = FALSE),
                  ),
                  shinyjs::hidden(
                     column(6, id = "class_rescale_ui",
                        h4(id = "class_rescale_label", "Scale within classes"),
                        switchInput("class_rescale", value = TRUE, disabled = TRUE),
                     )
                  )
               ),
               hr(),
               div(textOutput("default_scale"), class = "h4")
            )
         )
      ),
      fluidRow(
         column(6, align = "right",
            actionButton("reset_advanced", "Use Defaults")
         ),
         column(6, 
            actionButton("apply_advanced", "Apply")
         )
      )
   )
}

Build_Dataset_Tab <- function() {
   tabItem(tabName = "dataset",
      fluidRow(
         column(6,
            defaultBox(
               h3("Dataset", class = "centered"),
               htmlOutput("dataset_ui", class = "plain-text")
            ),
            shinyjs::hidden(
               defaultBox(id = "annot_info_box",
                  h3("Annotations", class = "centered"),
                  htmlOutput("annot_info_ui", class = "plain-text")
               )
            )
         ),
         column(6,
            defaultBox(
               h3("Data Download", class = "centered"),
               div(class = "plain-text",
                  p(HTML(paste(sep = "<br>",
                     "Raw - Raw input data for TMT Mosaic",
                     "Live - Current data in memory (ordered, named, and normalized)",
                     "Summary - Polished XLSX file describing 'Raw' data returned by the server"
                  )))
               ),
               fluidRow(
                  column(6, textInput("dtable_filename", "Filename:", placeholder = "(leave blank for default)")),
                  column(6,
                     radioGroupButtons(
                        inputId = "dtable_format", label = "File Type:", 
                        choices = c("Raw", "Live", "Summary"), 
                        justified = TRUE, status = "primary",
                        selected = "Live"
                     )
                  )
               ),
               column(12, align = "center",
                  downloadButton("dtable_download", "Download"),
                  shinyjs::hidden(
                     downloadButton("dtable_pep_download", "Download Peptides")
                  )
               )
            ),
            shinyjs::hidden(
               defaultBox(id = "server_info_box",
                  h3("Server Info", class = "centered"),
                  htmlOutput("server_info_ui", class = "plain-text")
               )
            )
         )
      )
   )
}

Build_Contact_Tab <- function(){
   tabItem(tabName = "contact",
      fluidRow(align = "center",
         h3("TMT Mosaic Administrator"),
         h4(contact)
      )
   )
}

