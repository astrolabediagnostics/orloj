#' Cell assignments report.
#'
#' Generate plots and CSV files that report on cell assignments to cell subsets.
#' For each level, the report includes the following over all (subset, channel)
#' combinations:
#' - Heatmap of channel median in that subset
#' - Heatmap of channel coefficient of variation in that subset
#' - Joy plots of channel distribution in each subset
#'
#' @param sample An Astrolabe sample.
#' @param channels List of channels to include in report.
#' @param report_levels If true, intermediate level heatmaps will be exported
#' in addition to terminal level heatmaps.
#' @import patchwork
#' @export
reportCellAssignments <- function(sample,
                                  channels = NULL,
                                  report_levels = FALSE) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  exprs <- fcsExprs(sample)
  levels <- getCellSubsetLevels(sample)
  if (is.null(channels)) channels <- sample$cell_assignments$class_channels

  if (!all(channels %in% colnames(exprs))) {
    stop("At least one channel is missing from the exprs data frame")
  }

  if (!report_levels) {
    levels <- dplyr::filter(levels, Level %in% c("Assignment", "Profile"))
  }

  # Generate reports for each labeling level.
  report <- list()
  for (level_row in seq(nrow(levels))) {
    level_col <- levels$Level[level_row]
    parent_col <- levels$Parent[level_row]

    level_indices <- which(!is.na(exprs[[level_col]]))
    level_exprs <- exprs[level_indices, ]
    parent_labels <- unique(exprs[level_indices, parent_col])

    # Iterate over each parent labels. For Assignment and Profiling, there is
    # only one parent label ("Root").
    level_report <- list()
    for (pl_label in parent_labels) {
      # Get data for this parent label and convert to long format.
      pl_indices <- level_exprs[[parent_col]] == pl_label
      pl_exprs <-
        level_exprs[pl_indices, c(channels, level_col)]
      pl_exprs_long <-
        reshape2::melt(pl_exprs,
                       id.vars = level_col,
                       variable.name = "Channel",
                       value.name = "Intensity")
      # Decide on name of value (either parent or assignment/profile).
      name <- pl_label
      if (level_col %in% c("Assignment", "Profile")) {
        name <- level_col
      }

      channel_order <-
        gtools::mixedsort(as.character(unique(pl_exprs_long$Channel)))

      # Figure: Intensity heatmap.
      level_report[[name]] <-
        plotHeatmapAggregate(pl_exprs_long,
                             x = "Channel",
                             y = level_col,
                             value = "Intensity",
                             type = "cluster_labels",
                             title = name,
                             x_axis_order = channel_order)

      # Figure: CV(Intensity) heatmap.
      level_report[[paste0(name, "_cv")]] <-
        plotHeatmapAggregate(pl_exprs_long,
                             x = "Channel",
                             y = level_col,
                             value = "Intensity",
                             func = function(v) sd(v) / mean(v),
                             type = "cluster_labels_cv",
                             title = paste0(name, ": Coefficient of Variation"),
                             x_axis_order = channel_order)

      # Figure: Joy plot of all channels over all subsets.
      if (level_col == "Assignment") {
        level_report[[paste0(name, "_joy")]] <-
          .plotJoy(pl_exprs_long, level_col)
      }
    }

    report <- c(report, level_report)
  }

  report
}

.plotJoy <- function(exprs_long, level_col) {
  # Joy plot object from expression data in long format.
  level_order <- gtools::mixedsort(unique(exprs_long[[level_col]]))
  channel_order <- gtools::mixedsort(as.character(unique(exprs_long$Channel)))

  # Calculate xlim for each channel.
  xlims <- exprs_long %>%
    dplyr::group_by(Channel) %>%
    dplyr::summarize(
      Min = floor(min(Intensity)),
      Max = ceiling(max(Intensity))
    )

  # Generate the figure, iterating over each (level_value, channel) combination.
  plt <- NULL
  for (level_value in level_order) {
    df <- exprs_long[exprs_long[[level_col]] == level_value, ]
  
    for (channel in channel_order) {
      channel_df <- df[df$Channel == channel, ]
      xlim <- with(xlims, c(Min[Channel == channel], Max[Channel == channel]))

      # Generate this panel.
      obj <- 
        ggplot(channel_df, aes(x = Intensity)) +
        geom_density(fill = "grey70") +
        coord_cartesian(xlim = xlim) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x = element_blank(),
              panel.background = element_blank(),
              text = element_text(size = 12))
    
      # Top row: Add the channel as title.
      if (level_value == level_order[1]) {
        obj <- obj +
          labs(title = channel) +
          theme(plot.title = element_text(hjust = 0.5))
      }
    
      # Lefternmost column: Add the level as Y-axis label.
      if (channel == channel_order[1]) {
        level_label <- gsub("_", " ", level_value)
        level_label <- paste(strwrap(level_label, width = 10), collapse = "\n")
        obj <- obj +
          ylab(level_label)
      } else {
        obj <- obj +
          theme(axis.title.y = element_blank())
      }
    
      if (is.null(plt)) {
        plt <- obj
      } else {
        plt <- plt + obj
      }
    }
  }

  # Layout so that each channel is stacked as one column.
  plt <- plt + plot_layout(ncol = length(channel_order))

  # 100x100 tile for each panel.
  width <- length(channel_order) * 100
  height <- length(level_order) * 100
  list(plt = plt, width = width, height = height)
}
