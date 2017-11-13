# flow_cytometry.R
# Pre-processing flow cytometry data.

#' Preprocess flow cytometry sample.
#'
#' Quality control, cleaning, and transformation for flow cytometry sample.
#'
#' @inheritParams preprocess
#' @return Sample after the above steps are done.
flowPreprocess <- function(sample) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  stop("flowPreprocess not implemented yet")
}
