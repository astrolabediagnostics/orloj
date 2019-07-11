#' Combined aggregate statistics reports.
#'
#' Generate plots and CSV files that report on aggregate statistics over all
#' experiment samples and all cell subsets.
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

  analyses <- c()
  if ("Assignment" %in% aggregate_statistics$subset_cell_counts$Parent) {
    analyses <- c(analyses, "Assignment")
  }
  if ("Profile" %in% aggregate_statistics$subset_cell_counts$Parent) {
    analyses <- c(analyses, "Profiling")
  }
  analyses <- nameVector(analyses)

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
    # Add sample feature and reorganize columns.
    sample_features <- getExperimentSampleFeatures(experiment)
    cell_counts <-
      dplyr::left_join(cell_counts, sample_features, by = "SampleId")
    col_order <-
      c("SampleId", "Filename", "SampleName", experiment$features$FeatureName,
        "CellSubset", "N", "Frequency", "ScaledFrequency")
    cell_counts <- cell_counts[, col_order]

    # Figure: Heatmap of cell subset frequency by sample.
    frequency_heat_map <-
      plotHeatmap(cell_counts,
                  x = "CellSubset",
                  y = "SampleName",
                  value = "Frequency",
                  type = "frequency",
                  title = paste0(analysis, ": Frequency Heat Map"))
    # Figure: Heatmap of scaled cell subset frequency by sample.
    scaled_frequency_heat_map <-
      plotHeatmap(cell_counts,
                  x = "CellSubset",
                  y = "SampleName",
                  value = "ScaledFrequency",
                  type = "scaled_frequency",
                  title = paste0(analysis, ": Scaled Frequency Heat Map"))
    # Figures: Bar plot of each cell subset. 
    cell_subsets <- unique(cell_counts$CellSubset)
    bar_plots <- lapply(nameVector(cell_subsets), function(cell_subset) {
      cell_subset_counts <-
        dplyr::filter(cell_counts, CellSubset == cell_subset) 
      plt_list <- 
        plotBarPlot(cell_subset_counts, 
                    x = "SampleName", 
                    y = "Frequency",  
                    title = cell_subset)  
      plt_list$plt <- plt_list$plt +  
        scale_y_continuous(labels = scales::percent)  
      plt_list
    })

    # Generate report.
    channel_subset_statistics <-
      .getSubsetChannelStatistics(experiment, aggregate_statistics, analysis)
    list(
      # Spreadsheet: Count, frequency, and scaled frequency over all (sample,
      # subset) pairs.
      cell_subset_counts = list(data = cell_counts),
      frequency_heat_map = frequency_heat_map,
      scaled_frequency_heat_map = scaled_frequency_heat_map,
      # Spreadsheets: Median, mean, and CV of marker intensities over all
      # (sample, subset, marker) combinations.
      sample_subset_means = list(data = channel_subset_statistics$scs_mean),
      sample_subset_medians = list(data = channel_subset_statistics$scs_median),
      sample_subset_cvs = list(data = channel_subset_statistics$scs_cv),
      bar_plots = bar_plots
    )
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

.getSubsetChannelStatistics <- function(experiment, aggregate_statistics,
                                        analysis) {
  # Reorganize subset channel statistics into one CSVs: median and CV. Each CSV
  # includes statistics from all samples and cell subsets, including sample
  # features as columns.

  if (analysis == "Profiling") analysis <- "Profile"

  # Get same features.
  sample_features <- getExperimentSampleFeatures(experiment)
  if (nrow(sample_features) == 0) {
    sample_features <- experiment$samples
  } else {
    sample_features <-
      dplyr::left_join(experiment$samples, sample_features, by = "SampleId")
  }

  # Reshape channel aggregate statistics into wide (one for median, one for CV).
  scs <- aggregate_statistics$subset_channel_statistics
  scs <- scs[scs$Parent == analysis, ]
  scs$Cv <- scs$Sd / scs$Mean
  scs_mean <-
    reshape2::dcast(scs,
                    SampleId + ChannelName ~ CellSubset,
                    value.var = "Mean")
  scs_median <-
    reshape2::dcast(scs,
                    SampleId + ChannelName ~ CellSubset,
                    value.var = "Median")
  scs_cv <-
    reshape2::dcast(scs,
                    SampleId + ChannelName ~ CellSubset,
                    value.var = "Cv")

  # Combine each data frame with sample features and reorganize column order.
  scs_mean <- dplyr::left_join(scs_mean, sample_features, by = "SampleId")
  scs_median <- dplyr::left_join(scs_median, sample_features, by = "SampleId")
  scs_cv <- dplyr::left_join(scs_cv, sample_features, by = "SampleId")
  cols <- colnames(sample_features)
  cols <- c(cols, setdiff(colnames(scs_median), cols))
  scs_mean <- scs_mean[, cols]
  scs_median <- scs_median[, cols]
  scs_cv <- scs_cv[, cols]

  list(scs_mean = scs_mean, scs_median = scs_median, scs_cv = scs_cv)
}
