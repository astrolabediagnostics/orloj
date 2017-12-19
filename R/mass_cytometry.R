# mass_cytometry.R
# Pre-processing mass cytometry data.

# The name of all of the bead channels.
orloj_bead_channel_names <-
  c(
    "Ce140Di",
    "Eu151Di",
    "Ho165Di",
    "Lu175Di"
  )

#' Preprocess mass cytometry sample.
#'
#' Quality control, cleaning, and transformation for mass cytometry sample.
#'
#' @inheritParams preprocess
#' @return Sample after the above steps are done.
massPreprocess <- function(sample, cofactor, bead_percentile) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  # Transform the mass channels.
  sample <- massTransformMassChannels(sample, cofactor)
  # Find non-bead events.
  sample <- massFindNonBeadEvents(sample, bead_percentile)

  sample
}

#' Transform expression data for mass cytometry sample.
#'
#' Choose which channels to transform and apply the asinh transformation. Mass
#' cytometry channels are identified using the naming standard from the Fluidigm
#' software.
#'
#' @param sample An Astrolabe sample.
#' @param cofactor Cofactor for asinh transformation.
#' @return Sample, with expression data transformed.
massTransformMassChannels <- function(sample, cofactor) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  name <- sample$parameter_name

  # Mass channel pattern for older CyTOF software. Include an exception for the
  # 110:114 CD3 trick that was used in some older experiments.
  mass_channel_pattern_old <-
    "\\([A-Z][a-z]\\d\\d\\d\\)|110:114"
  # Mass channel pattern for newer CyTOF software.
  mass_channel_pattern_new <- "[A-Z][a-z]?\\d\\d\\d?Di"

  mass_channel_indices_old <- grep(mass_channel_pattern_old, name)
  mass_channel_indices_new <- grep(mass_channel_pattern_new, name)

  if (length(mass_channel_indices_old) > length(name) / 2) {
    mass_channel_indices <- mass_channel_indices_old
  } else if (length(mass_channel_indices_new) > length(name) / 2) {
    mass_channel_indices <- mass_channel_indices_new
  } else {
    stop("Unable to identify mass cytometry channels by name pattern")
  }

  sample$exprs[, mass_channel_indices] <-
    asinh(sample$exprs[, mass_channel_indices] / cofactor)
  sample$mass_channel_indices <- mass_channel_indices

  sample
}

#' Find the indices of non-bead events.
#'
#' Non-bead event detection first clusters the data using the FlowSOM algorithm.
#' Then, for each cluster, the mean of all bead channels is calculated. Any
#' cluster that has at least one mean channel value below the bead_percentile
#' for that channel is a non-bead cluster. The indices of all cells in non-bead
#' clusters are non-bead events.
#'
#' @param sample An Astrolabe sample.
#' @param bead_percentile Clusters where at least one of the bead channels is
#' below this threshold are non-bead clusters.
#' @param min_n Minimum number of events in order to try and detect beads.
#' @return Sample, with additional non_bead_indices and non_bead_message.
#' fields.
massFindNonBeadEvents <- function(sample, bead_percentile, min_n = 200) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  n_bead_channels <- length(orloj_bead_channel_names)
  bead_channel_indices <- match(orloj_bead_channel_names, sample$parameter_name)

  # Default values (no bead removal).
  bead_channel_thresholds <- NA
  non_bead_indices <- seq(nrow(sample$exprs))

  if (any(is.na(bead_channel_indices))) {
    # We require all four channels in order to detect beads. If any are missing,
    # we assume no beads.
    non_bead_message <- "At least one bead channel is missing from data"
  } else if (nrow(fcsExprs(sample)) < min_n) {
    non_bead_message <- paste0("File has less than ", min_n, " events")
  } else {
    # Get bead expression data.
    bead_exprs <- fcsExprs(sample)[, bead_channel_indices]
    colnames(bead_exprs) <- orloj_bead_channel_names
    # Cluster with FlowSOM and compute channel means.
    som <- FlowSOM::SOM(as.matrix(bead_exprs))
    bead_exprs$Cluster <- som$mapping[, 1]
    bead_exprs_long <-
      reshape2::melt(bead_exprs,
                     id.vars = "Cluster",
                     variable.name = "BeadChannelName",
                     value.name = "Intensity")
    bead_means <- bead_exprs_long %>%
      dplyr::group_by(Cluster, BeadChannelName) %>%
      dplyr::summarize(MeanIntensity = mean(Intensity))
    # Tag each (cluster, channel) combination which is above the bead percentile
    # for that channel.
    bead_channel_thresholds <- bead_exprs_long %>%
      dplyr::group_by(BeadChannelName) %>%
      dplyr::summarize(Threshold =
                         as.numeric(quantile(Intensity, bead_percentile))) %>%
      tibble::deframe()
    bead_means <-
      dplyr::mutate(bead_means,
                    AboveThreshold = MeanIntensity >
                      bead_channel_thresholds[BeadChannelName])
    # A non-bead cluster will have at least one channel below threshold.
    non_bead_clusters <- bead_means %>%
      dplyr::group_by(Cluster) %>%
      dplyr::summarize(SumAboveThreshold = sum(AboveThreshold)) %>%
      dplyr::filter(SumAboveThreshold < n_bead_channels) %>%
      `$`(Cluster)
    # Set up variables for sample.
    non_bead_indices <- which(bead_exprs$Cluster %in% non_bead_clusters)
    n_beads_found <- nrow(bead_exprs) - length(non_bead_indices)
    if (n_beads_found == 0) {
      non_bead_message <- "No bead events found"
    } else {
      non_bead_message <- paste0(n_beads_found, " bead events found")
    }
  }

  sample$bead_channel_thresholds <- bead_channel_thresholds
  sample$non_bead_indices <- non_bead_indices
  sample$non_bead_message <- non_bead_message

  sample
}
