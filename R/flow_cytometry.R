# flow_cytometry.R
# Pre-processing flow cytometry data.

#' Preprocess flow cytometry sample.
#'
#' Quality control, cleaning, and transformation for flow cytometry sample.
#'
#' @param sample An Astrolabe sample.
#' @param cofactor Cofactor for asinh transformation.
#' @return Sample after the above steps are done.
#' @export
flowPreprocess <- function(sample, cofactor = 1000) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  sample <- .flowTransformChannels(sample, cofactor)

  sample
}


.flowTransformChannels <- function(sample, cofactor, new_max = 6) {
  # Transform expression data for flow cytometry sample.
  fsc_ssc_time_pattern <- "(^FSC|^SSC|^Time)"

  # Transform non-FSC/SSC/Time using asinh.
  name <- sample$parameter_name
  channel_indices <- grep(fsc_ssc_time_pattern, name, invert = TRUE)
  sample$exprs[, channel_indices] <-
    asinh(sample$exprs[, channel_indices] / cofactor)
  sample$transformed_channel_indices <- channel_indices
  sample$trans_function <- function(x) { asinh(x / cofactor) }
  sample$rev_trans_function <- function(x) { sinh(x) * cofactor }

  # Rescale FSC/SSC/Time to the range of transformed antibodies.
  fsc_ssc_cols <- grep(fsc_ssc_time_pattern, name)
  for (col in fsc_ssc_cols) {
    v <- sample$exprs[[col]]
    v[is.na(v)] <- 0 # set missing values to 0.
    v <- pmax(0, v) # set negative values to 0.
    sample$exprs[[col]] <- v / quantile(v, 0.99) * new_max
  }
  sample$fsc_ssc_cols <- fsc_ssc_cols

  sample
}
