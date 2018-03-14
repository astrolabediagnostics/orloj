# flow_cytometry.R
# Pre-processing flow cytometry data.

# Aurora ----------------------------------------------------------------------

.auroraTransformChannels <- function(sample, cofactor) {
  # Transform intensity values for Aurora samples.
  name <- sample$parameter_name

  # Channels to ignore.
  ignore_pattern <- "(^FSC-.$|^SSC-.$|^Time$)"
  transform_cols <- which(!grepl(ignore_pattern, name))
  for (col in transform_cols) {
    # Set minimum value to 0.
    sample$exprs[[col]] <- pmax(0, sample$exprs[[col]])
    # Transform with arcsinh.
    sample$exprs[[col]] <- asinh(sample$exprs[[col]] / cofactor)
  }

  sample$cofactor <- cofactor

  sample
}

.auroraScaleFscSsc <- function(sample) {
  # Scale SSC and FSC channels to the range in which we expect to find
  # transformed antibodies.
  fsc_ssc_pattern <- "(^FSC-.$|^SSC-.$)"
  fsc_ssc_cols <- which(grepl(fsc_ssc_pattern, sample$parameter_name))
  fsc_ssc_factor <-
    10 ^ floor(log10(quantile(unlist(sample$exprs[, fsc_ssc_cols]), 0.95)))
  sample$exprs[, fsc_ssc_cols] <- sample$exprs[, fsc_ssc_cols] / fsc_ssc_factor
  
  sample$fsc_ssc_factor <- fsc_ssc_factor

  sample
}

#' Preprocess Aurora sample.
#'
#' Quality control, cleaning, and transformation for Aurora sample.
#'
#' @param cofactor Cofactor for asinh transformation.
#' @inheritParams preprocess
#' @return Sample after the above steps are done.
auroraPreprocess <- function(sample, cofactor = 150) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  sample <- .auroraTransformChannels(sample, cofactor)
  sample <- .auroraScaleFscSsc(sample)

  sample
}
