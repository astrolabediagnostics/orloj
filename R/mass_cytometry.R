# mass_cytometry.R
# Pre-processing mass cytometry data.

#' Bead channel names.
#' 
#' @return Vector of bead channel names.
#' @export
beadChannelNames <- function() {
  c(
    "Ce140Di",
    "Eu151Di",
    "Ho165Di",
    "Lu175Di"
  )
}

#' Preprocess mass cytometry sample.
#'
#' Quality control, cleaning, and transformation for mass cytometry sample.
#'
#' @inheritParams preprocess
#' @return Sample after the above steps are done.
massPreprocess <- function(sample, cofactor) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  # Transform the mass channels.
  sample <- massTransformMassChannels(sample, cofactor)

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
