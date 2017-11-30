# experiment.R
# Functions for interacting with Astrolabe experiments.

#' Load Astrolabe experiment.
#'
#' @param experiment_path Path to the Astrolabe experiment.
#' @return Astrolabe experiment.
#' @export
loadExperiment <- function(experiment_path) {
  experiment <- readRDS(file.path(experiment_path, "config.RDS"))

  experiment$astrolabe_path <- NULL
  experiment$samples_path <- NULL
  experiment$analysis_path <- file.path(experiment_path, "analysis")
  experiment$reports_path <- file.path(experiment_path, "reports")

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

#' Experiment cell subset counts.
#'
#' Cell subset counts for all of the samples in an experiment.
#'
#' @param experiment An Astrolabe experiment.
#' @param level Cell subset level. Currently supported levels are "Assignment"
#' and "Profiling".
#' @return Cell subset counts for that level.
#' @export
experimentCellSubsetCounts <- function(experiment, level = "Assignment") {
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

  counts <- aggregate_statistics$subset_cell_counts
  counts <- counts %>%
    dplyr::filter(Parent == level) %>%
    dplyr::select(-Parent)
  counts <- dplyr::left_join(experiment$samples, counts, by = "SampleId")

  counts
}

#' Experiment cell subset channel statistics.
#'
#' Cell subset channel statistics for all of the samples in an experiment.
#'
#' @param experiment An Astrolabe experiment.
#' @param level Cell subset level. Currently supported levels are "Assignment"
#' and "Profiling".
#' @return Cell subset channel statistics for that level.
#' @export
experimentCellSubsetChannelStatistics <- function(experiment,
                                                  level = "Assignment") {
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
  stats <- stats %>%
    dplyr::filter(Parent == level) %>%
    dplyr::select(-Parent)
  stats <- dplyr::left_join(experiment$samples, stats, by = "SampleId")

  stats
}

#' Differential abundance analysis.
#'
#' Load the experiment differential abundance analysis, for a given cell subset
#' level.
#'
#' @param experiment An Astrolabe experiment.
#' @param level Cell subset level. Currently supported levels are "Assignment"
#' and "Profiling".
#' @return Differential abundance analysis list.
#' @export
differentialAbundanceAnalysis <- function(experiment, level = "Assignment") {
  if (!(level %in% c("Assignment", "Profiling"))) {
    stop("level is not \"Assignment\" or \"Profiling\"")
  }

  daa_filename <-
    file.path(experiment$analysis_path, "differential_abundance_analysis.RDS")
  if (!file.exists(daa_filename)) {
    stop(paste0(daa_filename, " not found"))
  }
  differential_abundance_analysis <- readRDS(daa_filename)

  daa <- differential_abundance_analysis$differential_abundance_analysis

  if (!(level %in% names(daa))) {
    stop("required level not found in differential abundance analysis")
  }
  daa <- daa[[level]]

  # Find feature values.
  feature_values <- lapply(nameVector(names(daa)), function(feature_name) {
    unique(experiment$sample_features[[feature_name]])
  })
  # Convert feature IDs to names.
  m <-
    match(as.numeric(gsub("^feature_", "", names(daa))),
          experiment$features$FeatureId)
  names(daa) <- experiment$features$FeatureName[m]
  names(feature_values) <- experiment$features$FeatureName[m]

  # Convert tables to tibbles and only keep p-value and FDR columns, and logFC for
  # features with two values.
  lapply(nameVector(names(daa)), function(feature_name) {
    if (is.null(daa[[feature_name]])) {
      NULL
    } else {
      tab <- as.data.frame(daa[[feature_name]]) %>%
        tibble::rownames_to_column("CellSubset")
      if ("logFC" %in% colnames(tab)) {
        cols <- c("CellSubset", "logFC", "PValue", "FDR")
      } else {
        cols <- c("CellSubset", "PValue", "FDR")
      }
      cols <- c(cols, paste0("median_", feature_values[[feature_name]]))
      
      tab[, cols]
    }
  })
}
