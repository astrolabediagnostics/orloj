#' Cell counts report.
#'
#' Generate plots and CSV files that report on cell counts over all cell
#' subsets. Report includes bar plots of cell counts by hierarchy parent.
#'
#' @param sample An Astrolabe samlpe.
#' @param report_levels If true, intermediate level heatmaps will be exported
#' in addition to terminal level heatmaps.
#' @export
reportCellCounts <- function(sample, report_levels = FALSE) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  subset_cell_counts <- sample$aggregate_statistics$subset_cell_counts

  if (!report_levels) {
    subset_cell_counts <-
      dplyr::filter(subset_cell_counts, Parent %in% c("Assignment", "Profile"))
  }

  report <- list()

  # Figure: Bar plot of subset cell counts, by parent.
  parent_labels <- unique(subset_cell_counts$Parent)

  for (parent_label in parent_labels) {
    parent_cell_counts <-
      dplyr::filter(subset_cell_counts, Parent == parent_label)

    # Order cell subsets by name.
    parent_cell_counts$CellSubset <-
      factor(parent_cell_counts$CellSubset,
             levels = gtools::mixedsort(parent_cell_counts$CellSubset))
    plt_list <-
      orloj::plotBarPlot(parent_cell_counts,
                         x = "CellSubset",
                         y = "N",
                         title = parent_label,
                         theme = NULL)
    plt_list$plt <- plt_list$plt + ggplot2::scale_y_continuous(labels = scales::comma)

    report[[parent_label]] <- plt_list
  }

  report
}
