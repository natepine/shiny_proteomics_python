# Load required library for text repel
tryCatch({
   if (!require(ggrepel, quietly = TRUE)) {
      cat("Warning: ggrepel package not found. Text labels will not be available.\n")
      # Define a dummy function if ggrepel is not available
      geom_text_repel <- function(...) {
         cat("geom_text_repel not available - ggrepel package missing\n")
         return(NULL)
      }
   }
}, error = function(e) {
   cat("Error loading ggrepel package:", conditionMessage(e), "\n")
   # Define a dummy function if there's an error
   geom_text_repel <- function(...) {
      cat("geom_text_repel not available due to error:", conditionMessage(e), "\n")
      return(NULL)
   }
})

#Reactives
correlation <- reactiveValues()

#UI
observe({
   choices <- setNames(seq(dataset$numGroups), theme$groupNames)
   selected <- choices[table(dataset$groups) > 1]
   if(length(selected) < 4) {
      # If we don't have 4 groups with replicates, use all available groups
      selected <- choices[1:min(4, length(choices))]
   }
   updateSelectInput(session, "corr_grp1_a", choices = choices, selected = selected[[1]])
   updateSelectInput(session, "corr_grp1_b", choices = choices, selected = selected[[2]])
   updateSelectInput(session, "corr_grp2_c", choices = choices, selected = selected[[1]])
   updateSelectInput(session, "corr_grp2_d", choices = choices, selected = selected[[min(3, length(selected))]])
})

#Calculations
buildCorrelation <- observe({
   columns <- req(values$normProtCols)

   input$button_update_correlation
   if(dataset$numGroups < 2) {
      showNotification("Sorry, correlation plot cannot be used with less than 2 groups.",
         duration = NULL, type = "warning")
      return()
   }
   
   # Get the selected groups
   grp1_a <- req(isolate(input$corr_grp1_a))
   grp1_b <- req(isolate(input$corr_grp1_b))
   grp2_c <- req(isolate(input$corr_grp2_c))
   grp2_d <- req(isolate(input$corr_grp2_d))

   # Validate that we have different groups for meaningful comparisons
   if (grp1_a == grp1_b) {
      sendSweetAlert(
         title = "Same Groups in Comparison 1",
         text = "You selected the same group for both Group A and Group B in Comparison 1. They cannot be the same.",
         type = "warning"
      )
      return(NULL)
   }
   
   if (grp2_c == grp2_d) {
      sendSweetAlert(
         title = "Same Groups in Comparison 2",
         text = "You selected the same group for both Group C and Group D in Comparison 2. They cannot be the same.",
         type = "warning"
      )
      return(NULL)
   }

   showNotification("(1/4) Validating groups...", duration = 4)
   
   # Calculate fold change for comparison 1 (B/A)
   fc1_data <- calculateFoldChange(columns, grp1_a, grp1_b)
   if(is.null(fc1_data)) return()
   
   # Calculate fold change for comparison 2 (D/C)
   fc2_data <- calculateFoldChange(columns, grp2_c, grp2_d)
   if(is.null(fc2_data)) return()
   
   showNotification("(2/4) Calculated fold changes...", duration = 4)
   
   # Merge the data on MosaicID to ensure we're comparing the same proteins
   merged_data <- merge(fc1_data[, c("MosaicID", "GeneSymbol", "log2FC")], 
                       fc2_data[, c("MosaicID", "log2FC")], 
                       by = "MosaicID", 
                       suffixes = c("_comp1", "_comp2"))
   
   # Remove any rows with infinite or NaN values
   merged_data <- merged_data[is.finite(merged_data$log2FC_comp1) & is.finite(merged_data$log2FC_comp2), ]
   
   if(nrow(merged_data) < 3) {
      showNotification("Not enough valid data points for correlation analysis (need at least 3).",
         duration = NULL, type = "warning")
      return()
   }
   
   showNotification("(3/4) Merged datasets...", duration = 4)
   
   # Calculate correlation and linear regression
   correlation_result <- cor(merged_data$log2FC_comp1, merged_data$log2FC_comp2, use = "complete.obs")
   lm_result <- lm(log2FC_comp2 ~ log2FC_comp1, data = merged_data)
   lm_summary <- summary(lm_result)
   
   # Calculate distance and threshold data for table functions
   tryCatch({
      merged_data$distance <- abs(merged_data$log2FC_comp1 - merged_data$log2FC_comp2)
      merged_data$distance <- round(merged_data$distance, 4)
      # Default threshold of 1.5 - will be updated in reactive functions
      merged_data$above_threshold <- merged_data$distance >= 1.5
   }, error = function(e) {
      cat("Error calculating initial distance/threshold:", conditionMessage(e), "\n")
      merged_data$distance <- rep(0, nrow(merged_data))
      merged_data$above_threshold <- rep(FALSE, nrow(merged_data))
   })
   
   # Store results
   correlation$data <- merged_data
   correlation$cor <- correlation_result
   correlation$lm <- lm_result
   correlation$lm_summary <- lm_summary
   
   # Initialize labels for annotation system and clicked points
   correlation$labels <- setNames(rep(FALSE, nrow(merged_data)), merged_data$MosaicID)
   correlation$clicked_points <- setNames(rep(FALSE, nrow(merged_data)), merged_data$MosaicID)
   
   # Safely create comparison names with error handling
   correlation$comparison1_name <- tryCatch({
      paste(theme$groupNames[[as.numeric(grp1_b)]], "vs", theme$groupNames[[as.numeric(grp1_a)]])
   }, error = function(e) {
      cat("Error creating comparison1_name:", conditionMessage(e), "\n")
      paste("Group", grp1_b, "vs Group", grp1_a)
   })
   
   correlation$comparison2_name <- tryCatch({
      paste(theme$groupNames[[as.numeric(grp2_d)]], "vs", theme$groupNames[[as.numeric(grp2_c)]])
   }, error = function(e) {
      cat("Error creating comparison2_name:", conditionMessage(e), "\n")
      paste("Group", grp2_d, "vs Group", grp2_c)
   })
   
   showNotification("(4/4) Completed correlation analysis!", duration = 4)
   
}, suspended = TRUE)

startCorrelation <- observe({
   # Listen for correlation init conditions
   req(input$corr_grp1_a)
   req(input$corr_grp1_b)
   req(input$corr_grp2_c)
   req(input$corr_grp2_d)
   req(values$normProtCols)
   req(input$sidebarMenu == "correlation")

   buildCorrelation$resume()
   startCorrelation$destroy()
})

# Helper function to calculate fold change for a pair of groups
calculateFoldChange <- function(columns, grp_a, grp_b) {
   sub_index <- dataset$groups %in% c(grp_a, grp_b)
   groups <- dataset$groups[sub_index]
   a_index <- groups == grp_a
   b_index <- groups == grp_b
   
   if(sum(a_index) < 1 || sum(b_index) < 1) {
      showNotification("At least one column in each group is required",
         duration = NULL, type = "warning")
      return(NULL)
   }
   
   columns_subset <- columns[, sub_index]
   
   # Remove rows with 0 variance or less than 1 value:
   row_valid <- apply(columns_subset, 1, function(x) sum(!is.na(x))) > 0
   columns_subset <- columns_subset[row_valid, , drop = FALSE]
   
   a_df <- columns_subset[, a_index, drop = FALSE]
   b_df <- columns_subset[, b_index, drop = FALSE]
   
   # Get fold changes:
   fc_df <- dataset$proteins$info[row_valid, c("MosaicID", "GeneSymbol"), drop = FALSE] %>%
      mutate(
         A_MEAN = rowMeans(a_df, na.rm = TRUE),
         B_MEAN = rowMeans(b_df, na.rm = TRUE))

   fc_df$A_MEAN[fc_df$A_MEAN == 0] <- .001
   fc_df$B_MEAN[fc_df$B_MEAN == 0] <- .001
   
   fc_df <- fc_df %>%
      mutate(log2FC = log2(B_MEAN / A_MEAN)) %>%
      arrange(MosaicID, GeneSymbol, log2FC, A_MEAN, B_MEAN)

   # Filter out infinite values
   fc_df <- fc_df %>% filter(is.finite(log2FC))
   
   return(fc_df)
}

# Output: Correlation plot
output$plot_correlation <- renderPlotly({
   tryCatch({
      data <- req(correlation$data)
      lm_result <- req(correlation$lm)
      
      # Ensure labels are reactive (can be NULL initially) - this makes the plot reactive to label changes
      labels <- correlation$labels
      annotation_color <- input$correlation_annotation_color  # Also make reactive to color changes
      
      # Add customdata for click functionality
      data$customdata <- data$MosaicID
      
      # Get threshold value from input (default to 1.5 if not set)
      threshold <- tryCatch({
         if(is.null(input$correlation_label_threshold) || is.na(input$correlation_label_threshold)) {
            1.5
         } else {
            as.numeric(input$correlation_label_threshold)
         }
      }, error = function(e) {
         cat("Error getting threshold value:", conditionMessage(e), "\n")
         1.5  # Default fallback
      })
      
      # Identify points above threshold (distance from origin)
      tryCatch({
         data$distance <- sqrt(data$log2FC_comp1^2 + data$log2FC_comp2^2)
         data$above_threshold <- data$distance >= threshold
      }, error = function(e) {
         cat("Error calculating distance/threshold:", conditionMessage(e), "\n")
         data$distance <- rep(0, nrow(data))
         data$above_threshold <- rep(FALSE, nrow(data))
      })
      
      # Create tooltip text with MosaicID and GeneSymbol
      tryCatch({
         # Create enhanced tooltip with gene symbol prominence
         data$tooltip_text <- paste0(
            "<b>", ifelse(is.na(data$GeneSymbol) | data$GeneSymbol == "", 
                         data$MosaicID, 
                         paste0(data$GeneSymbol, " (", data$MosaicID, ")")), "</b><br>",
            "log2FC Comp1: ", round(data$log2FC_comp1, 3), "<br>",
            "log2FC Comp2: ", round(data$log2FC_comp2, 3), "<br>",
            "Distance: ", round(data$distance, 3)
         )
         
         # Create simplified hover label for dynamic annotation
         data$hover_label <- ifelse(is.na(data$GeneSymbol) | data$GeneSymbol == "", 
                                   data$MosaicID, 
                                   data$GeneSymbol)
      }, error = function(e) {
         cat("Error creating tooltip text:", conditionMessage(e), "\n")
         data$tooltip_text <- paste0("MosaicID: ", data$MosaicID)
         data$hover_label <- data$MosaicID
      })
      
      # Create the scatter plot with tooltip aesthetics and color coding
      tryCatch({
         # Add color coding based on annotation status
         data$point_color <- "grey"  # Default color for all points
         data$point_size <- 1
         
         # Debug: Print labels status
         cat("Labels status - is.null:", is.null(labels), ", any(labels):", if(!is.null(labels)) any(labels) else "N/A", "\n")
         if(!is.null(labels)) {
            cat("Number of TRUE labels:", sum(labels), "\n")
         }
         
         # First apply annotation colors (nearest proteins, etc.)
         if(!is.null(labels) && any(labels)) {
            # Find which points should be annotated based on MosaicID
            annotated_mosaicids <- names(labels)[labels]
            annotated_indices <- which(data$MosaicID %in% annotated_mosaicids)
            
            cat("Annotated MosaicIDs:", paste(annotated_mosaicids, collapse = ", "), "\n")
            cat("Annotated indices:", paste(annotated_indices, collapse = ", "), "\n")
            
            if(length(annotated_indices) > 0) {
               used_annotation_color <- if(!is.null(annotation_color)) {
                     annotation_color
                  } else {
                     "#1E90FF"
                  }
               data$point_color[annotated_indices] <- used_annotation_color
               data$point_size[annotated_indices] <- 1.5  # Make annotated points slightly larger
               
               cat("Applied annotation color", used_annotation_color, "to", length(annotated_indices), "points\n")
            }
         } else {
            cat("No annotation labels to apply\n")
         }
         
         # Then apply clicked point colors (overrides annotation colors)
         clicked_points <- correlation$clicked_points
         if(!is.null(clicked_points) && any(clicked_points)) {
            clicked_mosaicids <- names(clicked_points)[clicked_points]
            clicked_indices <- which(data$MosaicID %in% clicked_mosaicids)
            
            if(length(clicked_indices) > 0) {
               data$point_color[clicked_indices] <- "purple"
               data$point_size[clicked_indices] <- 2
               
               cat("Applied purple color to", length(clicked_indices), "clicked points\n")
            }
         }
      }, error = function(e) {
         cat("Error setting up point colors:", conditionMessage(e), "\n")
         data$point_color <- "#1f77b4"
         data$point_size <- 1
      })
      
      p <- ggplot(
         data, aes(
            x = log2FC_comp1, 
            y = log2FC_comp2
         )
      ) +
         geom_point(aes(
            text = tooltip_text,
            customdata = MosaicID,
            color = I(point_color),
            size = I(point_size)),
            alpha = 0.7) +
         geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linewidth = 1) +
         labs(
            x = paste("log2 FC:", correlation$comparison1_name),
            y = paste("log2 FC:", correlation$comparison2_name),
            title = "Fold Change Correlation"
         ) +
         theme_minimal() +
         theme(
            plot.title = element_text(size = 14),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10)
         )
      
      # Add ggrepel labels for annotated and clicked points
      tryCatch({
         # Make sure to get the latest clicked points state
         clicked_points <- correlation$clicked_points
         
         # Debug: Check which points should be labeled
         annotated_mosaicids <- if(!is.null(labels) && any(labels)) names(labels)[labels] else character(0)
         clicked_mosaicids <- if(!is.null(clicked_points) && any(clicked_points)) names(clicked_points)[clicked_points] else character(0)
         all_labeled_mosaicids <- unique(c(annotated_mosaicids, clicked_mosaicids))
         
         cat("Annotated MosaicIDs for labeling:", paste(annotated_mosaicids, collapse = ", "), "\n")
         cat("Clicked MosaicIDs for labeling:", paste(clicked_mosaicids, collapse = ", "), "\n")
         cat("All labeled MosaicIDs:", paste(all_labeled_mosaicids, collapse = ", "), "\n")
         
         # Create a subset for labeling (both annotated and clicked points)
         # Use MosaicID instead of point_color to be more explicit
         label_data <- data[data$MosaicID %in% all_labeled_mosaicids, ]
         
         if(nrow(label_data) > 0) {
            cat("Adding ggrepel labels for", nrow(label_data), "points\n")
            p <- p + geom_text_repel(
               data = label_data,
               aes(label = hover_label),
               size = 3,
               max.overlaps = Inf,
               box.padding = 0.3,
               point.padding = 0.3,
               segment.color = "grey50",
               segment.size = 0.3
            )
         } else {
            cat("No points to label (no annotated or clicked points)\n")
         }
      }, error = function(e) {
         cat("Error adding ggrepel labels:", conditionMessage(e), "\n")
      })

      
      # Convert to plotly with enhanced tooltip and click functionality
      g <- ggplotly(p, tooltip = "text") %>%
         plotlyDefaults() %>%
         layout(
            annotations = annotated_mosaicids
         ) %>%
         config(
            displayModeBar = TRUE,
            scrollZoom = TRUE,
            displaylogo = FALSE,
            modeBarButtonsToRemove = c(
               "select2d",
               "lasso2d"
            ),
            toImageButtonOptions = list(
               format = "png",
               filename = paste0(dataset$ID, "_correlation"),
               width = NULL,
               height = NULL
            )
         ) %>%
         onRender('function(el) {
            el.on(\'plotly_click\', correlationSingleClickHandler);
            el.on(\'plotly_doubleclick\', correlationDoubleClickHandler);
            
            el.on(\'plotly_hover\', function(eventData) {
               var point = eventData.points[0];
               if (point && point.customdata) {
               }
            });
            
            el.on(\'plotly_unhover\', function(eventData) {
            });
         }')
      
      # Ensure customdata is properly set for all traces
      for(i in 1:length(g$x$data)) {
         if(!is.null(g$x$data[[i]]$mode) && grepl("markers", g$x$data[[i]]$mode)) {
            g$x$data[[i]]$customdata <- data$MosaicID
            # Also add hover label data for enhanced labeling
            g$x$data[[i]]$hoverlabel <- data$hover_label
            # Enhance hover template for better interactivity
            g$x$data[[i]]$hovertemplate <- paste0(
               "<b>%{text}</b>",
               "<br><i>Click to select this protein</i>",
               "<extra></extra>"
            )
            
            # Set annotated points to markers+text mode
            if(!is.null(labels) && any(labels)) {
               annotated_mosaicids <- names(labels)[labels]
               # Find which points in this trace are annotated
               trace_annotated <- g$x$data[[i]]$customdata %in% annotated_mosaicids
               
               if(any(trace_annotated)) {
                  # For annotated points, change mode to markers+text and add text labels
                  g$x$data[[i]]$mode <- "markers+text"
                  g$x$data[[i]]$text <- ifelse(trace_annotated, 
                                               data$hover_label[match(g$x$data[[i]]$customdata, data$MosaicID)], 
                                               "")
                  g$x$data[[i]]$textposition <- "top center"
                  g$x$data[[i]]$textfont <- list(size = 10, color = "black")
               }
            }
         }
      }
      
      return(g)
      
   }, error = function(e) {
      cat("Error in plot_correlation:", conditionMessage(e), "\n")
      print(traceback())
      
      # Return empty plotly object on error
      plotly_empty() %>%
         layout(title = list(text = paste("Error generating plot:", e$message)))
   })
})

# Output: Statistics table
output$correlation_stats_table <- renderDataTable({
   tryCatch({
      lm_summary <- req(correlation$lm_summary)
      cor_value <- req(correlation$cor)
      
      # Extract key statistics with error handling
      slope <- lm_summary$coefficients[2, 1]
      slope_pval <- lm_summary$coefficients[2, 4]
      intercept <- lm_summary$coefficients[1, 1]
      r_squared <- lm_summary$r.squared
      adj_r_squared <- lm_summary$adj.r.squared
      
      # Safely get comparison names
      comp1_name <- if(!is.null(correlation$comparison1_name)) {
         correlation$comparison1_name
      } else {
         "Comparison 1"
      }
      
      comp2_name <- if(!is.null(correlation$comparison2_name)) {
         correlation$comparison2_name
      } else {
         "Comparison 2"
      }
      
      # Create vectors with matching lengths
      statistics <- c("Pearson Correlation (r)", "R-squared", "Adjusted R-squared", 
                     "Slope", "Intercept", "Slope P-value", "Number of Proteins",
                     "Comparison 1", "Comparison 2")
      
      values <- c(
         round(cor_value, 4),
         round(r_squared, 4),
         round(adj_r_squared, 4),
         round(slope, 4),
         round(intercept, 4),
         ifelse(slope_pval < 0.001, "< 0.001", round(slope_pval, 4)),
         nrow(correlation$data),
         comp1_name,
         comp2_name
      )
      
      # Ensure both vectors have the same length
      if(length(statistics) != length(values)) {
         cat("Warning: Statistics and values vectors have different lengths:",
             length(statistics), "vs", length(values), "\n")
         # Trim to the shorter length
         min_length <- min(length(statistics), length(values))
         statistics <- statistics[1:min_length]
         values <- values[1:min_length]
      }
      
      stats_df <- data.frame(
         Statistic = statistics,
         Value = as.character(values),
         stringsAsFactors = FALSE
      )
      
      datatable(stats_df, 
                options = list(
                  dom = 't',
                  pageLength = 20,
                  ordering = FALSE
                ),
                rownames = FALSE) %>%
         formatStyle(columns = 1:2, fontSize = '14px')
         
   }, error = function(e) {
      cat("Error in correlation_stats_table:", conditionMessage(e), "\n")
      print(traceback())
      
      # Return a simple error table
      error_df <- data.frame(
         Statistic = "Error",
         Value = paste("Unable to generate statistics:", conditionMessage(e)),
         stringsAsFactors = FALSE
      )
      datatable(error_df, options = list(dom = 't'), rownames = FALSE)
   })
})

# Download handler for correlation plot
output$correlation_download <- downloadHandler(
   filename = function() {
      tryCatch({
         base_name <- input$correlation_filename
         if(is.null(base_name) || base_name == "") {
            base_name <- paste0("correlation_", correlation$comparison1_name, "_vs_", correlation$comparison2_name)
            base_name <- gsub("[^A-Za-z0-9_-]", "_", base_name)  # Remove special characters
         }
         paste0(base_name, ".", input$correlation_plot_format)
      }, error = function(e) {
         cat("Error generating filename:", conditionMessage(e), "\n")
         paste0("correlation_plot.", ifelse(is.null(input$correlation_plot_format), "png", input$correlation_plot_format))
      })
   },
   content = function(file) {
      tryCatch({
         data <- req(correlation$data)
         
         # Get threshold value for download plot
         threshold <- tryCatch({
            if(is.null(input$correlation_label_threshold) || is.na(input$correlation_label_threshold)) {
               1.5
            } else {
               as.numeric(input$correlation_label_threshold)
            }
         }, error = function(e) {
            cat("Error getting threshold for download:", conditionMessage(e), "\n")
            1.5  # Default fallback
         })
         
         # Calculate distance and identify points above threshold
         tryCatch({
            data$distance <- sqrt(data$log2FC_comp1^2 + data$log2FC_comp2^2)
            data$above_threshold <- data$distance >= threshold
         }, error = function(e) {
            cat("Error calculating distance for download:", conditionMessage(e), "\n")
            data$distance <- rep(0, nrow(data))
            data$above_threshold <- rep(FALSE, nrow(data))
         })
         
         # Add color coding for download plot
         tryCatch({
            data$point_color <- "grey"  # Default color for all points
            data$point_size <- 1
            
            # First apply annotation colors for download
            labels <- correlation$labels
            if(!is.null(labels) && any(labels)) {
               annotated_mosaicids <- names(labels)[labels]
               annotated_indices <- which(data$MosaicID %in% annotated_mosaicids)
               
               if(length(annotated_indices) > 0) {
                  used_annotation_color <- if(!is.null(input$correlation_annotation_color)) {
                     input$correlation_annotation_color
                  } else {
                     "blue"
                  }
                  data$point_color[annotated_indices] <- used_annotation_color
                  data$point_size[annotated_indices] <- 1.5
               }
            }
            
            # Then apply clicked point colors for download (overrides annotation colors)
            clicked_points <- correlation$clicked_points
            if(!is.null(clicked_points) && any(clicked_points)) {
               clicked_mosaicids <- names(clicked_points)[clicked_points]
               clicked_indices <- which(data$MosaicID %in% clicked_mosaicids)
               
               if(length(clicked_indices) > 0) {
                  data$point_color[clicked_indices] <- "purple"
                  data$point_size[clicked_indices] <- 1
               }
            }
         }, error = function(e) {
            cat("Error setting up download colors:", conditionMessage(e), "\n")
            data$point_color <- "grey"
            data$point_size <- 1
         })
         
         # Create hover labels for download plot
         tryCatch({
            data$hover_label <- ifelse(is.na(data$GeneSymbol) | data$GeneSymbol == "", 
                                      data$MosaicID, 
                                      data$GeneSymbol)
         }, error = function(e) {
            cat("Error creating hover labels for download:", conditionMessage(e), "\n")
            data$hover_label <- data$MosaicID
         })
         
         p <- ggplot(data, aes(
            x = log2FC_comp1, 
            y = log2FC_comp2,
         )) +
         geom_point(aes(color = I(point_color), size = I(point_size)), alpha = 0.7) +
         geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
            labs(
               x = paste("log2 FC:", correlation$comparison1_name),
               y = paste("log2 FC:", correlation$comparison2_name),
               title = "Fold Change Correlation"
            ) +
            theme_minimal() +
            theme(
               plot.title = element_text(size = 14),
               axis.title = element_text(size = 12),
               axis.text = element_text(size = 10)
            )
         
         # Add ggrepel labels for annotated and clicked points in download
         tryCatch({
            # Use same logic as interactive plot for consistency
            annotated_mosaicids <- if(!is.null(labels) && any(labels)) names(labels)[labels] else character(0)
            clicked_mosaicids <- if(!is.null(clicked_points) && any(clicked_points)) names(clicked_points)[clicked_points] else character(0)
            all_labeled_mosaicids <- unique(c(annotated_mosaicids, clicked_mosaicids))
            
            # Create a subset for labeling (both annotated and clicked points)
            label_data <- data[data$MosaicID %in% all_labeled_mosaicids, ]
            
            if(nrow(label_data) > 0) {
               cat("Adding ggrepel labels for download:", nrow(label_data), "points\n")
               p <- p + geom_text_repel(
                  data = label_data,
                  aes(label = hover_label),
                  size = 3,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  segment.size = 0.3
               )
            } else {
               cat("No points to label in download (no annotated or clicked points)\n")
            }
         }, error = function(e) {
            cat("Error adding ggrepel labels to download:", conditionMessage(e), "\n")
         })
               
         ggsave(file, plot = p, device = input$correlation_plot_format, 
                width = 8, height = 6, dpi = 300)
                
      }, error = function(e) {
         cat("Error generating correlation download:", conditionMessage(e), "\n")
         print(traceback())
         # Create a simple error plot if download fails
         p <- ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5, label = paste("Error:", e$message)), size = 5) +
            theme_void()
         ggsave(file, plot = p, device = "png", width = 8, height = 6, dpi = 300)
      })
   }
)

# Click event handlers
observeEvent(input$correlation_click, {
   mosaicID <- input$correlation_click
   if(!is.null(mosaicID)) {
      proteinSelect$addToHistory(mosaicID)
   }
})

observeEvent(input$correlation_doubleclick, {
   mosaicID <- input$correlation_doubleclick
   if(!is.null(mosaicID)) {
      proteinSelect$addToHistory(mosaicID)
   }
})

# Custom click handler for plotly
observeEvent(input$correlation_plotly_click, {
   mosaicID <- input$correlation_plotly_click
   if(!is.null(mosaicID)) {
      proteinSelect$addToHistory(mosaicID)
   }
})

# Custom click handler for ggplotly
correlationSingleClickHandler <- function(event) {
   if(!is.null(event$data$customdata)) {
      mosaicID <- event$data$customdata
      if(!is.null(mosaicID) && length(mosaicID) > 0) {
         proteinSelect$addToHistory(mosaicID)
      }
   }
}

# Tables - Build correlation data table
buildCorrelationTable <- function(data) {
   tryCatch({
      cols <- c(
         "MosaicID",
         "UniprotID",
         "GeneSymbol", 
         "Description",
         "Sequence",
         "NumPeps",
         "log2FC_comp1",
         "log2FC_comp2",
         "distance",
         "Site"
      )
      
      info <- dataset$proteins$info[data$MosaicID, ]
      tbl <- cbind(
         data,
         info %>% select(-c("MosaicID", "GeneSymbol")) # Drop duplicate columns
      ) %>% mutate(
         NumPeps = getNumPepsFromID(info, MosaicID),
         distance = round(distance, 4),
         log2FC_comp1 = round(log2FC_comp1, 4),
         log2FC_comp2 = round(log2FC_comp2, 4)
      )
      
      tbl[intersect(cols, colnames(tbl))]
      
   }, error = function(e) {
      cat("Error building correlation table:", conditionMessage(e), "\n")
      data.frame(
         MosaicID = data$MosaicID,
         GeneSymbol = data$GeneSymbol,
         log2FC_comp1 = round(data$log2FC_comp1, 4),
         log2FC_comp2 = round(data$log2FC_comp2, 4),
         distance = round(data$distance, 4)
      )
   })
}

getCorrelationDT <- function(data, color = "#2E86AB") {
   tryCatch({
      data$log2FC_comp1 <- round(as.numeric(data$log2FC_comp1), 3)
      data$log2FC_comp2 <- round(as.numeric(data$log2FC_comp2), 3)
      data$distance <- round(as.numeric(data$distance), 3)
      
      if(!is.atomic(data$NumPeps)) {
         data$NumPeps <- apply(data$NumPeps, 1, function(peps) {
            paste0("{", paste0(peps, collapse = ","), "}")
         })
      }
      
      # Add the buttons to change to the proteins
      data$SET <- data$MosaicID # Use actual MosaicID values
      
      if(dataset$isSiteQuant) {
         protSeqs <- getProteinSequence(dataset, data$UniprotID)
         data <- data %>% mutate(
            ProtSize = nchar(protSeqs),
         ) %>% mutate(
            seqTitle = paste0("Site Pos: ", Site, "/", ProtSize, "<br>", trunc_seqs(protSeqs, Sequence)),
         )
      }
      
      # Reorder data
      cols <- c(
         "UniprotID",
         "GeneSymbol",
         "Description", 
         "Sequence",
         "ProtSize",
         "NumPeps",
         "log2FC_comp1",
         "log2FC_comp2",
         "distance",
         "SET",
         # hidden
         "Site",
         "seqTitle"
      )
      data <- data[intersect(cols, colnames(data))]
      
      # Render data
      colDefs <- c(
         addLink("https://www.uniprot.org/uniprotkb/", "UniprotID", "/entry", c(color = color)),
         addShinyButton("SET", "switch_to_mosaicID")
      )
      if(dataset$isSiteQuant)
         colDefs <- c(colDefs, addPosMarker(data, "ProtSize", "Site", "seqTitle", enableHTML = TRUE))
      
      dt <- data %>%
         datatable(
            options = list(
               sDom = '<"top">tp<"bottom">',
               columnDefs = colDefs,
               initComplete = datatableTheme(color),
               searchHighlight = TRUE,
               lengthMenu = c(10, 25, 50, 100, nrow(data)),
               pageLength = 10,
               autoWidth = FALSE,
               scrollX = TRUE,
               scrollCollapse = TRUE
            ),
            class = "compact cell-border",
            rownames = FALSE,
            selection = 'none',
            escape = FALSE
         )
      return(dt)
      
   }, error = function(e) {
      cat("Error creating correlation DataTable:", conditionMessage(e), "\n")
      # Return simple table on error
      data %>%
         datatable(
            options = list(
               pageLength = 10,
               scrollX = TRUE,
               scrollCollapse = TRUE,
               autoWidth = FALSE
            ),
            class = "compact cell-border",
            rownames = FALSE
         )
   })
}

# Update correlation table when data changes or threshold changes
observeEvent(list(correlation$data, input$correlation_label_threshold), {
   tryCatch({
      cat("observeEvent triggered - data is null:", is.null(correlation$data), "\n")
      if(!is.null(correlation$data)) {
         cat("Data has", nrow(correlation$data), "rows\n")
         
         # Get current threshold value
         threshold <- if(is.null(input$correlation_label_threshold) || is.na(input$correlation_label_threshold)) {
            1.5
         } else {
            as.numeric(input$correlation_label_threshold)
         }
         cat("Using threshold:", threshold, "\n")
         
         # Update above_threshold based on current threshold
         correlation$data$above_threshold <- correlation$data$distance >= threshold
         
         # Build all data table
         cat("Building correlation table...\n")
         correlation$table_all <- buildCorrelationTable(correlation$data)
         cat("Built table with", nrow(correlation$table_all), "rows\n")
         
         cat("Updated correlation table: All =", nrow(correlation$table_all), "rows\n")
      }
   }, error = function(e) {
      cat("Error updating correlation tables:", conditionMessage(e), "\n")
      print(traceback())
   })
})

# Output: Correlation data table for all points
output$table_correlation_all <- renderDataTable({
   cat("=== Rendering table_correlation_all ===\n")
   tryCatch({
      cat("Rendering correlation all table, data available:", !is.null(correlation$table_all), "\n")
      if(!is.null(correlation$table_all)) {
         cat("Table has", nrow(correlation$table_all), "rows and", ncol(correlation$table_all), "columns\n")
         cat("Column names:", paste(colnames(correlation$table_all), collapse = ", "), "\n")
      }
      shiny::validate(need(correlation$table_all, "No correlation data available."))
      result <- getCorrelationDT(correlation$table_all, "#2E86AB")
      cat("Successfully created DataTable\n")
      return(result)
   }, error = function(e) {
      cat("Error rendering correlation table:", conditionMessage(e), "\n")
      print(traceback())
      # Return empty table with error message
      data.frame(Error = paste("Unable to generate table:", e$message)) %>%
         datatable(options = list(dom = 't'), rownames = FALSE)
   })
})

# Output: Count of significant points
output$correlation_significant_count <- renderUI({
   tryCatch({
      if(!is.null(correlation$data)) {
         threshold <- if(is.null(input$correlation_label_threshold) || is.na(input$correlation_label_threshold)) {
            1.5
         } else {
            as.numeric(input$correlation_label_threshold)
         }
         
         count <- sum(correlation$data$above_threshold, na.rm = TRUE)
         tags$h4(paste(count, paste0("Significant Points (Distance ≥ ", threshold, ")")), 
                style = "color: #E63946")
      }
   }, error = function(e) {
      cat("Error rendering correlation count:", conditionMessage(e), "\n")
      tags$h4("Error calculating count", style = "color: red")
   })
})

# Download handler for correlation data
output$correlation_data_download <- downloadHandler(
   filename = function() {
      tryCatch({
         base_name <- input$correlation_data_filename
         if(is.null(base_name) || base_name == "") {
            base_name <- paste0("correlation_", 
                              gsub("[^A-Za-z0-9_-]", "_", correlation$comparison1_name), "_vs_", 
                              gsub("[^A-Za-z0-9_-]", "_", correlation$comparison2_name))
         }
         
         ext <- input$correlation_data_format
         
         paste0(base_name, ".", ext)
      }, error = function(e) {
         cat("Error generating correlation download filename:", conditionMessage(e), "\n")
         paste0("correlation_data.", ifelse(is.null(input$correlation_data_format), "csv", input$correlation_data_format))
      })
   },
   content = function(file) {
      tryCatch({
         # Export all correlation data
         if(is.null(correlation$table_all)) {
            df <- data.frame(ERROR = "No correlation data available")
         } else {
            df <- correlation$table_all
         }
         
         # Write file based on format
         if(input$correlation_data_format == "csv") {
            write.csv(df, file, row.names = FALSE)
         } else if(input$correlation_data_format == "tsv") {
            write.table(df, file, sep = "\t", row.names = FALSE)
         }
         
      }, error = function(e) {
         cat("Error generating correlation data download:", conditionMessage(e), "\n")
         # Write error message to file
         write.csv(data.frame(ERROR = paste("Download failed:", e$message)), file, row.names = FALSE)
      })
   }
)

# Annotations for correlation plot
output$correlation_annot_ui <- renderUI({
   auto_annotate <- isolate(input$correlation_annotation)

   if(is.null(auto_annotate)){
      data <- req(correlation$data)
      
      # Default annotation based on data size
      if(nrow(data) < 20){
         auto_annotate <- "largestDist10"
      } else if(nrow(data) < 200){
         auto_annotate <- "randomLargeDist20"
      } else {
         auto_annotate <- "largestDist10"
      }
   }

   choices <- c("none", "currentProtein", "nearestProteins", "largestDist10", "randomLargeDist20", "custom")
   geneSymbol <- dataset$proteins$info[values$activeMosaicID, "GeneSymbol"]
   idLabels <- paste0(dataset$idLabel, "s")
   names(choices) <- c(
      "None",
      paste0("Current ", dataset$idLabel, " (", geneSymbol, ")"),
      paste("15 Smallest Distance Overall", idLabels, "by TMT quantitative data"),
      "Top 10 Closest to POI (by Euclidean Distance)",
      "Random 20 from High Distance",
      paste("Custom", idLabels)
   )

   selectInput("correlation_annotation", "Annotate Points", selected = auto_annotate, choices = choices)
})

output$correlation_custom_annot_ui <- renderUI({
   req(input$correlation_annotation == "custom")

   updateSelectizeInput(session, "correlation_annot_proteins", choices = dataset$proteins$info$MosaicID, selected = NULL, server = TRUE)
   fluidRow(
      column(10, selectInput("correlation_annot_proteins", NULL, choices = NULL, multiple = TRUE)),
      column(2, align = "center", actionButton("correlation_annotate_custom", "Annotate"))
   )
})

# Annotation logic observer
observe({
   labels <- isolate(correlation$labels)
   data <- req(correlation$data)
   req(input$correlation_annotation != "custom")

   showNotification("Updating correlation labels...", duration = 1)

   if(is.null(labels) || length(labels) != nrow(data)){
      labels <- setNames(rep(FALSE, nrow(data)), data$MosaicID)
   }

   # Helper function to get largest distance proteins (using Euclidean distance from origin)
   getLargestDistance <- function(numHits = 10){
      # Calculate Euclidean distance from origin for annotation purposes
      euclidean_distance <- sqrt(data$log2FC_comp1^2 + data$log2FC_comp2^2)
      
      if(nrow(data) > numHits){
         ordered_indices <- order(euclidean_distance, decreasing = TRUE)[1:numHits]
         return(ordered_indices)
      } else {
         return(1:nrow(data))
      }
   }
   
   # Helper function to get random selection from largest distances
   getRandomLargestDistance <- function(numHits = 20){
      # Calculate Euclidean distance from origin
      euclidean_distance <- sqrt(data$log2FC_comp1^2 + data$log2FC_comp2^2)
      
      # Get top 50% by distance first, then randomly sample
      top_half_cutoff <- quantile(euclidean_distance, 0.5, na.rm = TRUE)
      eligible_indices <- which(euclidean_distance >= top_half_cutoff)
      
      if(length(eligible_indices) > numHits){
         return(sample(eligible_indices, numHits))
      } else {
         return(eligible_indices)
      }
   }

   labels[labels] <- FALSE
   active_indices <- switch(
      input$correlation_annotation,

      "currentProtein" = which(data$MosaicID == values$activeMosaicID),
      "nearestProteins" = {
         if(!is.null(values$nearestIDs) && length(values$nearestIDs) > 0) {
            which(data$MosaicID %in% values$nearestIDs)
         } else {
            integer(0)  # Return empty if nearestIDs not available
         }
      },
      "largestDist10" = getLargestDistance(10),
      "randomLargeDist20" = getRandomLargestDistance(20),

      NULL
   )

   if(!is.null(active_indices) && length(active_indices) > 0){
      labels[active_indices] <- TRUE
      cat("Annotation observer: Set", length(active_indices), "labels to TRUE for annotation type:", input$correlation_annotation, "\n")
   } else {
      cat("Annotation observer: No active indices for annotation type:", input$correlation_annotation, "\n")
   }

   correlation$labels <- labels
   cat("Annotation observer: Updated correlation$labels, total TRUE labels:", sum(labels), "\n")
})

# Custom annotation handler
observeEvent(input$correlation_annotate_custom, {
   tryCatch({
      labeled <- req(input$correlation_annot_proteins)
      
      if(!is.null(correlation$labels)) {
         # Clear existing labels
         correlation$labels[correlation$labels] <- FALSE
         
         # Set new custom labels
         valid_ids <- labeled[labeled %in% names(correlation$labels)]
         if(length(valid_ids) > 0) {
            correlation$labels[valid_ids] <- TRUE
            showNotification(paste("Annotated", length(valid_ids), "proteins"), duration = 2)
         } else {
            showNotification("No valid proteins found for annotation", type = "warning", duration = 3)
         }
      }
   }, error = function(e) {
      cat("Error in custom annotation:", conditionMessage(e), "\n")
      showNotification("Error applying custom annotations", type = "error", duration = 3)
   })
})

# Click handler for toggling clicked points (separate from annotations)
observeEvent(input$correlation_click, {
   if(!is.null(correlation$clicked_points)){
      mosaicID <- input$correlation_click
      if(!is.null(mosaicID) && mosaicID %in% names(correlation$clicked_points)){
         correlation$clicked_points[mosaicID] <- !correlation$clicked_points[mosaicID]
      }
   }
})
