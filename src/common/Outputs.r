# Global Output Templates

#' Set x spacing for tmtBarchart
#'
#' @param groups     Group identity for each column
#' @param summarize  Should groups have the same position?
#'
#' @return Vector of x positions for each column/group in one plex
tmt_scale_x <- function(groups, summarize) {
   u_groups <- unique(groups)
   group_index <- match(groups, u_groups)
   if(summarize) {
      return(group_index)
   }

   return(order(order(group_index)) + (group_index - 1) * 0.5)
}

#' Helper function for tmtBarchart
#'
#' @param plexes     Column plexes
#' @param numPlexes  Number of plexes to output
#' @param groups     Column groups
#' @param areReps    Should data be grouped as replicates.
#'
#' @return Vector of the group asignment for each column
groupColumns <- function(plexes, numPlexes, groups, areReps) {
   if(areReps) {
      return(groups)
   }else {
      if(all(plexes == "default") && numPlexes != 1) {
         plexes <- rep(1:numPlexes, each = length(plexes)/numPlexes)
      }
      return(plexes)
   }
}

#' Build a barchart to display a protein's expression across each column
#'
#' @param columnValues  The protein expression of each column
#' @param columnColors  A vector of colors for each column
#' @param columnLabels  A vector of labels for each column
#' @param groups        A character vector of each column's grouping for spacing and summarizing
#' @param labelGroups   Logical indicating whether x-axis should be labeled by group or column
#' @param summarize     Summarize each group with errorbars
#' @param title         Plot title
#' @param xLabel        X label
#' @param yLabel        Y label
#'
#' @return A ggplot2 barchart
tmtBarchart <- function(columnValues, columnColors, columnLabels, groups, labelGroups = FALSE,
      summarize = FALSE, title = "TMT Barchart", xLabel = NULL, yLabel = "TMT RA") {

   if(is.null(groups)) {
      summarize <- labelGroups <- FALSE
      groups <- rep(1, length(columnValues))
   }

   tmtData <- data.frame(
      RA = columnValues,
      x = tmt_scale_x(groups, summarize),
      Name = columnLabels,
      color = columnColors,
      group = groups
   )

   if(summarize) {
      boxData <- tmtData %>%
         group_by(x, group) %>%
         summarize(sd = sd(RA, na.rm = TRUE), RA = mean(RA, na.rm = TRUE))

      g <- ggplot(
            tmtData,
            aes(
               x = x,
               y = RA,
               text = paste0("Name: ", Name, "\nValue: ", RA)
            )) +
         geom_errorbar(data = boxData,
            aes(
               ymin = RA - sd,
               ymax = RA + sd,
               text = paste0("Group:", group, "\nMean: ", RA, "\nSD: ", sd)
            ),
            width = .5,
            size = .25) +
         geom_jitter(
            aes(fill = I(color)),
            shape = 21,
            size = 3,
            alpha = 0.8,
            show.legend = FALSE) +
         ylim(0, NA)
   }else {
      tmtData$RA[is.na(tmtData$RA)] <- 0

      g <- ggplot(
         tmtData,
         aes(
            x = x,
            y = RA,
            text = paste0("Name: ", Name, "\nValue: ", RA)
         )) +
      geom_bar(
         fill = tmtData$color,
         size = .2,
         color = "black",
         stat = "identity",
         show.legend = FALSE)
   }

   g <- g +
      ggtitle(title) +
      xlab(xLabel) +
      ylab(yLabel)

   limits <- tmtData$x
   labels <- columnLabels
   if(labelGroups) {
      labels <- groups
      if(!summarize) {
         labels <- groups[!duplicated(groups)]
         limits <- sapply(unique(groups), function(group) {
            mean(tmtData[groups == group, "x"])
         })
      }
   }

   # Ignore continuous limits to discrete scale warning
   suppressWarnings(
      g <- g + scale_x_discrete(limits = limits, labels = labels)
   )

   return(g)
}
