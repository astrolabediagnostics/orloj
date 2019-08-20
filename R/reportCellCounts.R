#' Cell counts report.
#'
#' Generate plots and CSV files that report on cell counts over all cell
#' subsets. Report includes bar plots of cell counts by level.
#'
#' @param sample An Astrolabe samlpe.
#' @export
reportCellCounts <- function(sample) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  subset_cell_counts <- sample$aggregate_statistics$subset_cell_counts

  # Figure: Bar plot of subset cell counts, by level.
  lapply(orloj::nameVector(colnames(sample$cell_subsets)), function(level) {
    counts <- subset(subset_cell_counts, Parent == level)

    # Order cell subsets by name.
    counts$CellSubset <-
      factor(counts$CellSubset, levels = unique(sample$cell_subsets[[level]]))
    plt_list <-
      orloj::plotBarPlot(counts,
                         x = "CellSubset",
                         y = "N",
                         title = paste0(level, ": Cell Counts"),
                         theme = NULL)
    plt_list$plt <-
      plt_list$plt + ggplot2::scale_y_continuous(labels = scales::comma)
    plt_list
  })
}
