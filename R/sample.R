# sample.R
# Functions for interacting with Astrolabe samples.

# All possible Astrolabe analyses.
astrolabe_analyses <- c(
  "cell_assignments",
  "subset_profiling_candidates",
  "subset_profiling_assignment",
  "aggregate_statistics",
  "calculate_qc_metrics"
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

#' Get sample name.
#'
#' @param experiment An Astrolabe experiment.
#' @param sample_id Astrolabe sample ID.
#' @return The sample's name, as supplied by user.
#' @export
getSampleName <- function(experiment, sample_id) {
  sample_row <- subset(experiment$samples, SampleId == sample_id)
  if (nrow(sample_row) != 1) stop("cannot find sample in experiment")
  sample_row$Name
}

#' Get sample filename.
#'
#' @param experiment An Astrolabe experiment.
#' @param sample_id Astrolabe sample ID.
#' @return The sample's filename.
#' @export
getSampleFilename <- function(experiment, sample_id) {
  sample_row <- subset(experiment$samples, SampleId == sample_id)
  if (nrow(sample_row) != 1) stop("cannot find sample in experiment")
  sample_row$Filename
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
  cat(paste0(
    sum(exprs$AstrolabeBead), " beads, ",
    sum(exprs$Debris), " debris and doublets, ",
    "and ",
    sum(exprs$Dead), " dead cells\n"))
}

#' Get expression data from sample.
#'
#' Retrieves the exprs data frame from the sample, incorporating all relevant
#' pre-processing and Astrolabe analyses.
#'
#' @param sample An Astrolabe sample.
#' @param keep_beads Whether bead events should be kept.
#' @param keep_debris Whether debris and doublet events should be kept.
#' @param keep_dead Whether dead events should be kept.
#' @return Expression data frame.
#' @export
fcsExprs <- function(sample,
                     keep_beads = FALSE,
                     keep_debris = FALSE,
                     keep_dead = FALSE) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  exprs <- sample$exprs

  exprs$AstrolabeBead <- FALSE
  exprs$Dead          <- FALSE
  exprs$Debris        <- FALSE

  # Mark beads.
  if (!is.null(sample$non_bead_indices)) {
    exprs$AstrolabeBead <- TRUE
    exprs$AstrolabeBead[sample$non_bead_indices] <- FALSE
  }

  # Mark debris and doublets.
  if (!is.null(sample$debris_indices)) {
    exprs$Debris[sample$non_bead_indices[sample$debris_indices]] <- TRUE
  }
  if (!is.null(sample$doublet_indices)) {
    exprs$Debris[sample$non_bead_indices[sample$doublet_indices]] <- TRUE
  }

  # Mark dead.
  if (!is.null(sample$live_indices)) {
    # Filter live_indices through beads and debris/doublets first.
    live_indices <- seq(nrow(exprs))
    live_indices <- live_indices[sample$non_bead_indices]
    live_indices <-
      live_indices[-c(sample$debris_indices, sample$doublet_indices)]
    live_indices <- live_indices[sample$live_indices]
    exprs$Dead <- TRUE
    exprs$Dead[exprs$AstrolabeBead | exprs$Debris] <- FALSE
    exprs$Dead[live_indices] <- FALSE
  }

  # Incorporate cell assignments.
  if (!is.null(sample$cell_assignments)) {
    cell_assignments <- sample$cell_assignments$cell_assignments
    ca_indices <- which(!exprs$AstrolabeBead & !exprs$Dead & !exprs$Debris)
    if (length(ca_indices) != nrow(cell_assignments)) {
      stop("number of rows in cell_assignments different than expected")
    }
    # Copy one column at a time.
    for (col_name in colnames(cell_assignments)) {
      exprs[[col_name]] <- NA
      exprs[ca_indices, col_name] <- cell_assignments[[col_name]]
    }
    # Update debris based on cell assignment.
    debris_labels <- astrolabeDebrisLabels()
    if (!is.null(sample$cm_is_debris)) {
      if (sample$cm_is_debris) {
        # Set CM to debris if required by experiment.
        debris_labels <- c(debris_labels, "CM-", "CM-_unassigned")
      }
    }
    if (!is.null(sample$granulocyte_is_debris)) {
      if (sample$granulocyte_is_debris) {
        # Set granulocyte to debris if required by experiment.
        debris_labels <- c(debris_labels, "Granulocyte")
      }
    }
    if (!is.null(sample$root_unassigned_is_debris)) {
      if (sample$root_unassigned_is_debris) {
        # Set granulocyte to debris if required by experiment.
        debris_labels <- c(debris_labels, "Root_unassigned")
      }
    }
    exprs$Debris[exprs$Assignment %in% debris_labels] <- TRUE

    # Update Assignment to AstrolabeBead/Debris/Dead values.
    exprs$Assignment[exprs$AstrolabeBead] <- "AstrolabeBead"
    exprs$Assignment[exprs$Debris] <- "Debris"
    exprs$Assignment[exprs$Dead] <- "Dead"

    # Incorporate profiling.
    if (!is.null(sample$subset_profiling_assignment)) {
      profile <- sample$subset_profiling_assignment$Profile
      profile_indices <-
        which(!exprs$AstrolabeBead & !exprs$Dead & !exprs$Debris)
      if (length(profile_indices) != length(profile)) {
        # Reverse compatibility: Length mismatch might be due to older version
        # of orloj treating Root_unassigned as debris.
        exprs$Debris[exprs$Assignment == "Root_unassigned"] <- TRUE
        profile_indices <-
          which(!exprs$AstrolabeBead & !exprs$Dead & !exprs$Debris)
        if (length(profile_indices) != length(profile)) {
          # Length still different, report error.
          stop("length of profile different than expected")
        }
      }
      exprs$Profile <- exprs$Assignment
      exprs$Profile[profile_indices] <- profile
    }
  }

  # Remove any unnecessary events.
  if (!keep_beads) exprs <- exprs[!exprs$AstrolabeBead, ]
  if (!keep_debris) exprs <- exprs[!exprs$Debris, ]
  if (!keep_dead) exprs <- exprs[!exprs$Dead, ]

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
#' @export
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
#' @export
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
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  exprs <- fcsExprs(sample)
  cols <- colnames(exprs)
  n_levels <- length(grep("Level_", colnames(exprs))) - 1

  levels <- tibble::tibble(Parent = character(0), Level = character(0))

  # Assignment levels.
  if (n_levels > 0) {
    levels <- tibble::tibble(
      Parent = c("Level_0", paste0("Level_", seq(n_levels) - 1)),
      Level = c("Assignment", paste0("Level_", seq(n_levels)))
    )
  }

  # Profiling level.
  if ("Profile" %in% cols) {
    levels <-
      rbind(
        levels,
        tibble::tibble(Parent = "Level_0", Level = "Profile")
      )
  }

  levels
}
