# sample.R
# Functions for interacting with Astrolabe samples.

# All possible Astrolabe analyses.
astrolabe_analyses <- c(
  "cell_assignments",
  "terminal_subset_profiling",
  "aggregate_statistics"
)

#' Load Astrolabe sample.
#'
#' Load FCS data and any additional analysis that has been exported by the
#' Astrolabe platform. Sample is identified through ID or name.
#'
#' @param experiment An Astrolabe experiment.
#' @param sample_id Astrolabe sample ID.
#' @param sample_name Astrolabe sample name.
#' @return Sample list with analysis data.
#' @export
loadSample <- function(experiment, sample_id = NULL, sample_name = NULL) {
  if (is.null(sample_id) & is.null(sample_name)) {
    stop("neither sample_id nor sample_name specified")
  }
  if (!is.null(sample_id) & !is.null(sample_name)) {
    stop("both sample_id and sample_name specified")
  }

  if (!is.null(sample_name)) {
    sample_id <-
      dplyr::filter(experiment$samples, Name == sample_name)$SampleId
    if (length(sample_id) == 0) {
      stop("sample name not found")
    }
  }

  analysis_path <- experiment$analysis_path

  fcs_data_filename <-
    file.path(analysis_path, paste0(sample_id, ".fcs_data.RDS"))
  if (!file.exists(fcs_data_filename)) {
    stop(paste0(fcs_data_filename, "not found"))
  }

  sample <- readRDS(fcs_data_filename)

  # Load any existing analysis files.
  for (analysis in astrolabe_analyses) {
    analysis_filename <-
      file.path(analysis_path, paste0(sample_id, ".", analysis, ".RDS"))
    if (file.exists(analysis_filename)) {
      sample[[analysis]] <- readRDS(analysis_filename)
    }
  }

  sample
}

#' Astrolabe sample summary.
#'
#' @param sample An Astrolabe sample.
#' @export
sampleSummary <- function(sample) {
  exprs <- fcsExprs(sample)
  cat(paste0(
    "An Astrolabe sample with ", nrow(exprs), " cells",
    " and ", ncol(sample$exprs), " channels\n"
  ))
  n_debris <-
    sum(sample$cell_assignments$cell_assignments$Assignment == "Debris")
  cat(paste0(
    nrow(sample$exprs) - length(sample$non_bead_indices), " bead events",
    " and ", n_debris, " debris events\n"
  ))
}

#' Get expression data from sample.
#'
#' Retrieves the exprs data frame from the sample, incorporating all relevant
#' pre-processing and Astrolabe analyses.
#'
#' @param sample An Astrolabe sample.
#' @param keep_debris Whether events in the Debris and Root/Unassigned labels
#' should be kept.
#' @return Expression (and analyses) data frame.
#' @export
fcsExprs <- function(sample, keep_debris = FALSE) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  exprs <- sample$exprs

  # Filter to non-bead events only.
  if (!is.null(sample$non_bead_indices)) {
    exprs <- exprs[sample$non_bead_indices, ]
  }

  # Incorporate cell assignments.
  if (!is.null(sample$cell_assignments)) {
    cell_assignments <- sample$cell_assignments
    exprs <- cbind(exprs, cell_assignments$cell_assignments)

    # Incorporate terminal subset profiling.
    if (!is.null(sample$terminal_subset_profiling)) {
      exprs$Profile <- exprs$Assignment
      exprs$Profile[exprs$Assignment != "Debris"] <-
        sample$terminal_subset_profiling$Profile
      exprs$Profile[is.na(exprs$Profile)] <-
        exprs$Assignment[is.na(exprs$Profile)]
    }

    if (!keep_debris) {
      exprs <- dplyr::filter(exprs, !(Assignment %in% astrolabe_debris_labels))
    }
  }

  exprs
}

#' Cell subset counts.
#'
#' Given an Astrolabe sample, return the cell subset counts for a labeling
#' level.
#'
#' @param sample An Astrolabe sample.
#' @param level Cell subset level. Currently supported levels are "Assignment"
#' and "Profiling".
#' @return Cell subset counts for that level.
sampleCellSubsetCounts <- function(sample, level = "Assignment") {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")
  if (!(level %in% c("Assignment", "Profile"))) {
    stop("level is not \"Assignment\" or \"Profiling\"")
  }

  sample$aggregate_statistics$subset_cell_counts %>%
    dplyr::filter(Parent == level) %>%
    dplyr::select(-Parent)
}

#' Cell subset channel statistics.
#'
#' Given an Astrolabe sample, return the cell subset channel intensity
#' statistics for a labeling level.
#'
#' @param sample An Astrolabe sample.
#' @param level Cell subset level. Currently supported levels are "Assignment"
#' and "Profiling".
#' @return Cell subset channel intensity statistics for that level.
sampleCellSubsetChannelStatistics <- function(sample, level = "Assignment") {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")
  if (!(level %in% c("Assignment", "Profile"))) {
    stop("level is not \"Assignment\" or \"Profiling\"")
  }

  sample$aggregate_statistics$subset_channel_statistics %>%
    dplyr::filter(Parent == level) %>%
    dplyr::select(-Parent)
}

#' Astrolabe cell subset levels.
#'
#' Return a data frame where each row corresponds to one level in the cell
#' subset classification. This includes levels that come from the hierarchy,
#' an "Assignment" level, and a terminal "Profile" level.
#'
#' @param sample An Astrolabe sample.
#' @return A data frame with cell subset levels.
#' @export
getCellSubsetLevels <- function(sample) {
  if (!isSample(sample)) stop("Expecting sample list object")

  exprs <- fcsExprs(sample)
  cols <- colnames(exprs)
  n_levels <- length(grep("Level_", colnames(exprs))) - 1

  levels <- tibble::tibble()

  if (n_levels > 0) {
    levels <- tibble::tibble(
      Parent = c("Level_0", paste0("Level_", seq(n_levels) - 1)),
      Level = c("Assignment", paste0("Level_", seq(n_levels)))
    )

    # Check for terminal subset profiling.
    if ("Profile" %in% cols) {
      levels <-
        rbind(
          levels,
          tibble::tibble(Parent = "Assignment", Level = "Profile")
        )
    }
  }

  levels
}
