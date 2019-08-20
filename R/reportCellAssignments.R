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
#' @import patchwork
#' @export
reportCellAssignments <- function(sample, channels = NULL) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  exprs <- fcsExprs(sample)
  if (!is.null(sample$cell_subsets)) {
    levels <- colnames(sample$cell_subsets)
  } else {
    # Reverse compatibility: Possible levels prior to cell_subsets update.
    levels <-
      intersect(colnames(exprs), c("Compartment", "Assignment", "Profiling"))
  }
  if (is.null(channels)) channels <- sample$cell_assignments$class_channels
  channels <- gtools::mixedsort(channels)

  if (!all(channels %in% colnames(exprs))) {
    stop("At least one channel is missing from the exprs data frame")
  }

  # Generate reports for each labeling level.
  report <- list()
  for (level in levels) {
    # Get data for this parent label and convert to long format.
    level_exprs <- exprs[, c(channels, level)]
    level_exprs_long <-
      reshape2::melt(level_exprs,
                     id.vars = level,
                     variable.name = "Channel",
                     value.name = "Intensity")

    # Figure: Intensity heatmap.
    report[[level]] <-
      plotHeatmapAggregate(level_exprs_long,
                           x = "Channel",
                           y = level,
                           value = "Intensity",
                           type = "cluster_labels",
                           title = paste0(level, ": Median"),
                           x_axis_order = channels)

    # Figure: CV(Intensity) heatmap.
    report[[paste0(level, "_cv")]] <-
      plotHeatmapAggregate(level_exprs_long,
                           x = "Channel",
                           y = level,
                           value = "Intensity",
                           func = function(v) sd(v) / mean(v),
                           type = "cluster_labels_cv",
                           title = paste0(level, ": Coefficient of Variation"),
                           x_axis_order = channels)

    # Figure: Joy plot of all channels over all subsets.
    # Subsample for purpose of plot generation.
    if (nrow(level_exprs) > 10000) {
      set.seed(12345)
      level_exprs <- level_exprs[sample(seq(nrow(level_exprs)), 10000), ]
      level_exprs_long <-
        reshape2::melt(level_exprs,
                       id.vars = level,
                       variable.name = "Channel",
                       value.name = "Intensity")
    }
    if (level %in% c("Assignment", "Compartment")) {
      report[[paste0(level, "_joy")]] <- .plotJoy(level_exprs_long, level)
    }
  }

  report
}

.plotJoy <- function(exprs_long, level) {
  # Joy plot object from expression data in long format.
  level_order <- gtools::mixedsort(unique(exprs_long[[level]]))
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
    df <- exprs_long[exprs_long[[level]] == level_value, ]
  
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
