# sample.R
# Functions for interacting with Astrolabe samples.

# All possible Astrolabe analyses.
astrolabe_analyses <- c(
  "cell_assignments",
  "subset_profiling_candidates",
  "subset_profiling_assignment",
  "aggregate_statistics",
  "calculate_qc_metrics",
  "calculate_marker_thresholds"
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

  # Load the experiment-wide cell_subsets, if it exists.
  cell_subsets_filename <-
    file.path(experiment$analysis_path, "experiment_cell_subsets.RDS")
  if (file.exists(cell_subsets_filename)) {
    sample$cell_subsets <- readRDS(cell_subsets_filename)
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

  # Remove any columns from previous Astrolabe runs.
  if (!is.null(sample$cell_subsets)) {
    existing_cols <- intersect(colnames(exprs), colnames(sample$cell_subsets))
    exprs[, existing_cols] <- NULL
  }
  
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
    exprs$Dead <- NA
    exprs$Dead[sample$non_bead_indices] <- TRUE
    exprs$Dead[sample$non_bead_indices[sample$live_indices]] <- FALSE
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
      exprs[[col_name]] <- as.character(NA)
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
        debris_labels <-
          c(debris_labels, "Granulocyte", "Neutrophil", "Eosinophil")
      }
    }
    if (!is.null(sample$root_unassigned_is_debris)) {
      if (sample$root_unassigned_is_debris) {
        # Set granulocyte to debris if required by experiment.
        debris_labels <- c(debris_labels, "Root_unassigned")
      }
    }
    if (!is.null(sample$subsets_are_debris)) {
      debris_labels <- c(debris_labels, sample$subsets_are_debris)
    }
    exprs$Debris[exprs$Assignment %in% debris_labels] <- TRUE
    
    # Update Assignment to AstrolabeBead/Debris/Dead values.
    exprs$Assignment[exprs$AstrolabeBead] <- "AstrolabeBead"
    exprs$Assignment[exprs$Debris] <- "Debris"
    exprs$Assignment[exprs$Dead] <- "Dead"
    
    # Incorporate profiling.
    if (!is.null(sample$subset_profiling_assignment)) {
      profiling <- sample$subset_profiling_assignment$Profile
      profiling_indices <-
        which(!exprs$AstrolabeBead & !exprs$Dead & !exprs$Debris)
      if (length(profiling_indices) != length(profiling)) {
        # Reverse compatibility: Length mismatch might be due to older version
        # of orloj treating Root_unassigned as debris.
        exprs$Debris[exprs$Assignment == "Root_unassigned"] <- TRUE
        profiling_indices <-
          which(!exprs$AstrolabeBead & !exprs$Dead & !exprs$Debris)
        if (length(profiling_indices) != length(profiling)) {
          # Length still different, report error.
          stop("length of profiling different than expected")
        }
      }
      exprs$Profiling <- exprs$Assignment
      exprs$Profiling[profiling_indices] <- profiling
    }
  }
  
  # Incorporate Compartment.
  if (!is.null(sample$cell_subsets) &&
      "Assignment" %in% colnames(sample$cell_subsets) &&
      "Assignment" %in% colnames(exprs)) {
    ass_to_compartment_map <-
      unique(sample$cell_subsets[, c("Assignment", "Compartment")])
    exprs <- dplyr::left_join(exprs, ass_to_compartment_map, by = "Assignment")
    # Update beads, debris, and dead with correct values for Compartment.
    exprs$Compartment[exprs$AstrolabeBead] <- "AstrolabeBead"
    exprs$Compartment[exprs$Debris] <- "Debris"
    exprs$Compartment[exprs$Dead] <- "Dead"
    # Everything else is "other".
    exprs$Compartment[is.na(exprs$Compartment)] <- "Other"
  }
  
  if (nrow(exprs) != nrow(sample$exprs)) {
    # Make sure that we did not mess up the number of rows on exprs.
    stop("expression size changed from original")
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
#' @param level Cell subset level, such as "Assignment" or "Profiling".
#' @return Cell subset counts for that level.
#' @export
sampleCellSubsetCounts <- function(sample, level = "Assignment") {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")
  
  if (!level %in% sample$aggregate_statistics$subset_cell_counts) {
    stop("sample does not have this level")
  }
  
  subset(sample$aggregate_statistics$subset_cell_counts,
         Parent == level, select = c("CellSubset", "N"))
}

#' Cell subset channel statistics.
#'
#' Given an Astrolabe sample, return the cell subset channel intensity
#' statistics for a labeling level.
#'
#' @param sample An Astrolabe sample.
#' @param level Cell subset level, such as "Assignment" or "Profiling".
#' @return Cell subset channel intensity statistics for that level.
#' @export
sampleCellSubsetChannelStatistics <- function(sample, level = "Assignment") {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  scs <- sample$aggregate_statistics$subset_channel_statistics

  if (!level %in% sample$aggregate_statistics$subset_cell_counts) {
    stop("sample does not have this level")
  }

  scs %>% dplyr::filter(Parent == level) %>% dplyr::select(-Parent)
}
