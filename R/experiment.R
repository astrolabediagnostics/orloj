# experiment.R
# Functions for interacting with Astrolabe experiments.

#' Load Astrolabe experiment.
#'
#' @param experiment_path Path to the Astrolabe experiment.
#' @return Astrolabe experiment.
#' @export
loadExperiment <- function(experiment_path) {
  experiment <- readRDS(file.path(experiment_path, "config.RDS"))

  experiment$analysis_path <- file.path(experiment_path, "analysis")
  experiment$reports_path <- file.path(experiment_path, "reports")
  experiment$samples_path <- file.path(experiment_path, "samples")

  experiment
}

#' Astrolabe experiment summary.
#'
#' @param experiment An Astrolabe experiment.
#' @export
experimentSummary <- function(experiment) {
  cat(paste0(
    "Astrolabe experiment with ",
    nrow(experiment$samples), " samples and ",
    nrow(experiment$channels), " channels\n\n"
  ))
  cat(paste0(
    "Classification channels:\n",
    paste(experiment$class_channels, collapse = ", "),
    "\n"
  ))
  if (nrow(experiment$features) > 0) {
    cat("\n")
    cat("Sample features:\n")
    for (row_idx in seq(nrow(experiment$features))) {
      row <- experiment$features[row_idx, ]
      feature_id <- paste0("feature_", row$FeatureId)
      feature_values <- unique(experiment$sample_features[[feature_id]])
      cat(paste0(
        row$FeatureName, ": ", paste(feature_values, collapse = ", "), "\n"))
    }
  }
}

.chooseLevel <- function(experiment) {
  # Choose default level for experiment.
  if (is.null(experiment$organism)) {
    return("Assignment")
  } else if (experiment$organism == "profiling_only") {
    return("Profiling")
  } else {
    return("Assignment")
  }
}

#' Experiment cell subset counts.
#'
#' Cell subset counts and frequencies for all of the samples in an experiment.
#'
#' @param experiment An Astrolabe experiment.
#' @param level Cell subset level. Currently supported levels are "Assignment"
#' and "Profiling". The default is "Profiling" for Profiling Only experiments
#' and "Assignment" otherwise.
#' @return Cell subset counts for that level.
#' @export
experimentCellSubsetCounts <- function(experiment,
                                       level = .chooseLevel(experiment)) {
  if (!(level %in% c("Assignment", "Profiling"))) {
    stop('level is not "Assignment" or "Profiling"')
  }
  if (level == "Profiling") level = "Profile"

  aggregate_statistics_filename <-
    file.path(experiment$analysis_path, "combine_aggregate_statistics.RDS")
  if (!file.exists(aggregate_statistics_filename)) {
    stop(paste0(aggregate_statistics_filename, "not found"))
  }
  aggregate_statistics <- readRDS(aggregate_statistics_filename)

  counts <- aggregate_statistics$subset_cell_counts
  if (!(level %in% counts$Parent)) {
    stop("level not found in cell subset counts")
  }

  counts <- counts %>%
    dplyr::filter(Parent == level) %>%
    dplyr::select(-Parent)
  counts <- dplyr::left_join(experiment$samples, counts, by = "SampleId")

  # Calculate frequencies.
  counts <- counts %>%
    dplyr::group_by(SampleId) %>%
    dplyr::mutate(Freq = N / sum(N)) %>%
    dplyr::ungroup()

  counts
}

#' Experiment cell subset channel statistics.
#'
#' Cell subset channel statistics for all of the samples in an experiment.
#'
#' @param experiment An Astrolabe experiment.
#' @param level Cell subset level. Currently supported levels are "Assignment"
#' and "Profiling". The default is "Profiling" for Profiling Only experiments
#' and "Assignment" otherwise.
#' @return Cell subset channel statistics for that level.
#' @export
experimentCellSubsetChannelStatistics <- function(experiment,
                                                  level = .chooseLevel(experiment)) {
  if (!(level %in% c("Assignment", "Profiling"))) {
    stop("level is not \"Assignment\" or \"Profiling\"")
  }
  if (level == "Profiling") level = "Profile"

  aggregate_statistics_filename <-
    file.path(experiment$analysis_path, "combine_aggregate_statistics.RDS")
  if (!file.exists(aggregate_statistics_filename)) {
    stop(paste0(aggregate_statistics_filename, "not found"))
  }
  aggregate_statistics <- readRDS(aggregate_statistics_filename)

  stats <- aggregate_statistics$subset_channel_statistics
  if (!(level %in% stats$Parent)) {
    stop("level not found in cell subset channel statistics")
  }

  stats <- stats %>%
    dplyr::filter(Parent == level) %>%
    dplyr::select(-Parent)
  stats <- dplyr::left_join(experiment$samples, stats, by = "SampleId")

  stats
}

#' Map experiment cell subsets to numerical values.
#'
#' Create a map from experiment cell subsets in a given level to unique
#' numerical values. This map is used when generating Assignment and Profiling
#' columns in exported FCS files.
#'
#' @param experiment An Astrolabe experiment.
#' @param level Cell subset level. Currently supported levels are "Assignment"
#' and "Profiling". The default is "Profiling" for Profiling Only experiments
#' and "Assignment" otherwise.
#' @return Map from cell subsets to numerical values.
#' @export
experimentCellSubsetMap <- function(experiment,
                                    level = .chooseLevel(experiment)) {
  if (!(level %in% c("Assignment", "Profiling"))) {
    stop("level is not \"Assignment\" or \"Profiling\"")
  }
  if (level == "Profiling") level = "Profile"

  aggregate_statistics_filename <-
    file.path(experiment$analysis_path, "combine_aggregate_statistics.RDS")
  if (!file.exists(aggregate_statistics_filename)) {
    stop(paste0(aggregate_statistics_filename, "not found"))
  }
  aggregate_statistics <- readRDS(aggregate_statistics_filename)

  stats <- aggregate_statistics$subset_channel_statistics
  if (!(level %in% stats$Parent)) {
    stop("level not found in cell subset channel statistics")
  }

  cell_subsets <-
    gtools::mixedsort(unique(stats$CellSubset[stats$Parent == level]))
  df <- data.frame(
    Value = seq(length(cell_subsets)),
    CellSubset = cell_subsets,
    stringsAsFactors = FALSE
  )

  # Add "bead", "dead", and "debris" at the end of the data frame.
  n <- nrow(df)
  df <-
    rbind(
      df,
      data.frame(
        Value = c(n + 1, n + 2, n + 3),
        CellSubset = c("Bead", "Debris", "Dead"),
        stringsAsFactors = FALSE
      ))
  df
}

#' Differential abundance analysis.
#'
#' Load the experiment differential abundance analysis, for a given cell subset
#' level.
#'
#' @param experiment An Astrolabe experiment.
#' @param level Cell subset level. Currently supported levels are "Assignment"
#' and "Profiling".
#' @param convert_ids Whether to convert Astrolabe IDs to feature names.
#' @return Differential abundance analysis list.
#' @export
differentialAbundanceAnalysis <- function(experiment,
                                          level = .chooseLevel(experiment),
                                          convert_ids = TRUE) {
  if (!(level %in% c("Assignment", "Profiling"))) {
    stop("level is not \"Assignment\" or \"Profiling\"")
  }

  # Value for NA feature value. Samples with this value should not be included
  # in the analysis.
  feature_na <- "__NA__"

  daa_filename <-
    file.path(experiment$analysis_path, "differential_abundance_analysis.RDS")
  if (!file.exists(daa_filename)) {
    stop(paste0(daa_filename, " not found"))
  }
  differential_abundance_analysis <- readRDS(daa_filename)

  daa <- differential_abundance_analysis$differential_abundance_analysis

  if (!(level %in% names(daa))) {
    return(NULL)
  }
  daa <- daa[[level]]

  # Find feature values.
  feature_values <- lapply(nameVector(names(daa)), function(feature_name) {
    values <- unique(experiment$sample_features[[feature_name]])
    values <- setdiff(values, feature_na)
    values
  })
  # Convert feature IDs to names.
  if (convert_ids) {
    m <-
      match(as.numeric(gsub("^feature_", "", names(daa))),
            experiment$features$FeatureId)
    names(daa) <- experiment$features$FeatureName[m]
    names(feature_values) <- experiment$features$FeatureName[m]
  }

  # Convert tables to tibbles and only keep p-value, FDR, and logFC.
  lapply(nameVector(names(daa)), function(feature_name) {
    if (is.null(daa[[feature_name]])) {
      NULL
    } else {
      tab <- as.data.frame(daa[[feature_name]]) %>%
        tibble::rownames_to_column("CellSubset")

      # Calculate max(logFC) for features with multiple values.
      log_fc_cols <- grep("logFC", colnames(tab))
      if (length(log_fc_cols) > 1) {
        log_fc <- tab[, log_fc_cols]
        max_log_fc <- apply(log_fc, 1, function(v) v[which.max(abs(v))])
        tab$logFC <- max_log_fc
      }

      cols <-
        c("CellSubset",
          "logFC",
          "PValue",
          "FDR",
          paste0("median_", feature_values[[feature_name]]))
      
      tab[, cols]
    }
  })
}

#' MDS map.
#'
#' Return the MDS map for the experiment, for a given level.
#'
#' @param experiment An Astrolabe experiment.
#' @param level Cell subset level. Currently supported levels are "Assignment"
#' and "Profiling".
#' @param convert_ids Whether to convert Astrolabe IDs to feature names.
#' @return MDS map.
#' @export
experimentMds <- function(experiment,
                          level = .chooseLevel(experiment),
                          convert_ids = TRUE) {
  if (!(level %in% c("Assignment", "Profiling"))) {
    stop("level is not \"Assignment\" or \"Profiling\"")
  }

  mds_filename <-
    file.path(experiment$analysis_path, "calculate_mds.RDS")
  if (!file.exists(mds_filename)) {
    stop(paste0(mds_filename, " not found"))
  }
  mds <- readRDS(mds_filename)

  mds <- mds[[level]]

  # Calculate and add mean frequency over all samples.
  cell_subset_counts <- experimentCellSubsetCounts(experiment, level = level)
  mds_freqs <- cell_subset_counts %>%
    dplyr::group_by(CellSubset) %>%
    dplyr::summarize(Freq = mean(Freq))
  mds <- dplyr::left_join(mds, mds_freqs, by = "CellSubset")

  # Calculate and add mean marker intensities over all samples.
  marker_stats <-
    experimentCellSubsetChannelStatistics(experiment, level = level)
  marker_stats <- marker_stats %>%
    dplyr::group_by(CellSubset, ChannelName) %>%
    dplyr::summarize(Median = mean(Median)) %>%
    reshape2::dcast(CellSubset ~ ChannelName, value.var = "Median")
  mds <- dplyr::left_join(mds, marker_stats, by = "CellSubset")

  # Add max(fold change) and -log10(FDR) for each feature.
  daa <- differentialAbundanceAnalysis(experiment, level = level, convert_ids)
  for (feature_name in names(daa)) {
    # Skip features for which we did not run DAA.
    if (is.null(daa[[feature_name]])) next

    tab <- daa[[feature_name]]
    tab <- tab[, c("CellSubset", "logFC", "FDR")]
    colnames(tab) <-
      c("CellSubset",
        paste0(feature_name, "_logFC"),
        paste0(feature_name, "_FDR"))
    mds <- dplyr::left_join(mds, tab, by = "CellSubset")
  }

  mds
}
