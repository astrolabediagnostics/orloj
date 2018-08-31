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

      # Figure: Intensity joy plot.
      pl_exprs_long$Channel <-
        factor(pl_exprs_long$Channel, levels = channel_order)
      joy_plt <- list()
      joy_plt$plt <-
        ggplot(pl_exprs_long, aes_string(x = "Intensity", y = level_col)) +
        ggridges::geom_density_ridges() +
        facet_wrap(~ Channel, nrow = 1) +
        labs(title = name, x = "Cell Subset", y = "Channel") +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.background = element_blank())
      joy_plt$height <- length(unique(pl_exprs_long[[level_col]])) * 20 + 20
      joy_plt$width <- length(channels) * 75 + 20
      level_report[[paste0(name, "_joy")]] <- joy_plt
    }

    report <- c(report, level_report)
  }

  report
}
