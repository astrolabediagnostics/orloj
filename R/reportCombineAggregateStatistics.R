#' Combined aggregate statistics reports.
#'
#' Generate plots and CSV files that report on aggregate statistics over all
#' experiment samples and all cell subsets. The report includes figures for the
#' "Assignment" level of the hierarchy only. It includes heatmaps of cell
#' subsets frequency and scaled frequency versus samples, and bar plots of
#' each subset frequency.
#'
#' @param experiment An Astrolabe experiment.
#' @export
reportCombineAggregateStatistics <- function(experiment) {
  samples <- experiment$samples
  # Load aggregate statistics.
  aggregate_statistics_filename <-
    file.path(experiment$analysis_path, "combine_aggregate_statistics.RDS")
  if (!file.exists(aggregate_statistics_filename)) {
    stop(paste0("unable to find ", aggregate_statistics_filename))
  }
  aggregate_statistics <- readRDS(aggregate_statistics_filename)

  analyses <- nameVector(c("Assignment", "Profiling"))
  lapply(analyses, function(analysis) {
    # Get cell counts for this analysis, match with sample names, and calculate
    # cell subset frequencies.
    cell_counts <- getCellCounts(aggregate_statistics, analysis)
    cell_counts <-
      dplyr::left_join(cell_counts, samples, by = "SampleId") %>%
      dplyr::rename(SampleName = Name)
    cell_counts <- cell_counts %>%
      dplyr::group_by(SampleName) %>%
      dplyr::mutate(Frequency = N / sum(N)) %>%
      dplyr::ungroup()

    # Add frequency of zero for missing subsets.
    cell_counts <-
      tidyr::complete(cell_counts, SampleName, CellSubset,
                      fill = list(Frequency = 0))
    # Calculate scaled frequency.
    cell_counts <- cell_counts %>%
      dplyr::group_by(CellSubset) %>%
      dplyr::mutate(ScaledFrequency = as.numeric(scale(Frequency))) %>%
      dplyr::ungroup()

    report <- list()

    # Figure: Heatmap of cell subset frequency by sample.
    report[["frequency_heatmap"]] <-
      plotHeatmap(cell_counts,
                  x = "CellSubset",
                  y = "SampleName",
                  value = "Frequency",
                  type = "frequency")
    # Figure: Heatmap of scaled cell subset frequency by sample.
    report[["scaled_frequency_heatmap"]] <-
      plotHeatmap(cell_counts,
                  x = "CellSubset",
                  y = "SampleName",
                  value = "ScaledFrequency",
                  type = "scaled_frequency")
    # Figure: Bar plot of each cell subset.
    report[["bar_plots"]] <- list()
    for (cell_subset in unique(cell_counts$CellSubset)) {
      cell_subset_counts <- dplyr::filter(cell_counts, CellSubset == cell_subset)
      plt_list <-
        plotBarPlot(cell_subset_counts,
                    x = "SampleName",
                    y = "Frequency",
                    title = cell_subset)
      plt_list$plt <- plt_list$plt +
        scale_y_continuous(labels = scales::percent)
      report[["bar_plots"]][[cell_subset]] <- plt_list
    }

    report
  })
}

#' Get cell counts for a given analysis type.
#'
#' Given an analysis type (such as "Assignment" or "Profiling") get the cell
#' counts for the respective cell subsets.
#'
#' @param aggregate_statistics List of data frames with aggregated statistics,
#' combined over all samples.
#' @param analysis Type of analysis for which to get count data.
#' @return A data fram with cell counts.
getCellCounts <- function(aggregate_statistics, analysis) {
  subset_cell_counts <- aggregate_statistics$subset_cell_counts

  if (analysis == "Assignment") {
    dplyr::filter(subset_cell_counts, Parent == "Assignment")
  } else if (analysis == "Profiling") {
    dplyr::filter(subset_cell_counts, Parent == "Profile")
  } else {
    stop("Unsupposed analysis type")
  }
}
