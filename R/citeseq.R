# citeseq.R
# Pre-processing citeseq.

#' Preprocess citeseq sample.
#'
#' Quality control, cleaning, and transformation for citeseq sample.
#'
#' @param cofactor Cofactor for asinh transformation.
#' @inheritParams preprocess
#' @return Sample after the above steps are done.
citeseqPreprocess <- function(sample, cofactor = 5) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  # Transform all channels.
  sample$exprs <- asinh(sample$exprs / cofactor)

  sample
}
