#' Cell assignments report.
#'
#' Generate plots and CSV files that report on cell assignments to cell subsets.
#' Report includes a heatmap of all cell subsets in each level, and biaxial
#' plots of each subset separately.
#'
#' @param sample An Astrolabe sample.
#' @param report_levels If true, intermediate level heatmaps will be exported
#' in addition to terminal level heatmaps.
#' @export
reportCellAssignments <- function(sample, report_levels = FALSE) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  exprs <- fcsExprs(sample)
  levels <- getCellSubsetLevels(sample)
  class_channels <- sample$cell_assignments$class_channels

  if (!report_levels) {
    levels <- dplyr::filter(levels, Level %in% c("Assignment", "Profile"))
  }

  # For Profile level, only report terminal subsets.
  if ("Profile" %in% levels$Level) {
    levels$Parent[levels$Level == "Profile"] <- "Level_0"
  }

  # Normalize expression to the 0-1 range across class channels.
  exprs_norm <- exprs
  exprs_norm[, class_channels] <-
    apply(exprs_norm[, class_channels], 2, function(v) {
      (v - min(v)) / (max(v) - min(v))
    })

  # Generate reports for each labeling level.
  report <- list()
  for (level_row in seq(nrow(levels))) {
    level_col <- levels$Level[level_row]
    parent_col <- levels$Parent[level_row]

    level_indices <- which(!is.na(exprs_norm[[level_col]]))
    level_exprs_norm <- exprs_norm[level_indices, ]
    parent_labels <- unique(exprs_norm[level_indices, parent_col])

    # Figure: Heatmap for each parent label.
    level_report <-
      lapply(nameVector(parent_labels), function(pl_label) {
        # Get data for this parent label and convert to long format.
        pl_indices <- level_exprs_norm[[parent_col]] == pl_label
        pl_exprs_norm <-
          level_exprs_norm[pl_indices, c(class_channels, level_col)]
        pl_exprs_norm_long <-
          reshape2::melt(pl_exprs_norm,
                         id.vars = level_col,
                         variable.name = "Channel",
                         value.name = "Intensity")
        # Generate heatmap.
        heatmap_title <- pl_label
        if (level_col %in% c("Assignment", "Profile")) {
          heatmap_title <- level_col
        }
        plotHeatmapAggregate(pl_exprs_norm_long,
                             x = "Channel",
                             y = level_col,
                             value = "Intensity",
                             type = "cluster_labels",
                             title = heatmap_title)
      })

    # Rename level to Assignment/Profile for these levels.
    if (level_col %in% c("Assignment", "Profile")) {
      names(level_report) <- level_col
    }

    report <- c(report, level_report)
  }

  report
}
