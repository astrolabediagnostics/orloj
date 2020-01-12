# flow_cytometry.R
# Pre-processing flow cytometry data.

# Aurora ----------------------------------------------------------------------

#' Transform Aurora channels.
#'
#' Given a FlowCore flow frame, transform channel values using estimateLogicle.
#'
#' @param flow_frame FlowCore flow frame.
#' @return Transformed flow_frame.
#' @export
auroraTransformChannels <- function(flow_frame) {
  # Choose which channels to transform.
  ignore_pattern <- "(^FSC|^SSC|^Time)"
  channels <- 
    grep(ignore_pattern, flowCore::colnames(flow_frame), value = TRUE,
         invert = TRUE)

  # Calculate and apply expectedLogicle. The "^" is required  due to issues with
  # logicle's regular expressions.
  el <-
    flowCore::estimateLogicle(flow_frame,
                              channels = paste0("^", channels))
  # In extreme cases, estimateLogicle of a given channel might end up with
  # illegal parameters. Try to sidestep the situation by estimating on a subset
  # of the data.
  for (ch in channels) {
    el_apply <- try(el@transforms[[ch]]@f(1), silent = TRUE)
    if (class(el_apply) == "try-error") {
      message(paste0("auroraTransformChannels: ", ch, 
                     " leads to transform error"))
      q_thresh <- quantile(flow_frame@exprs[, ch], 0.2)
      el_ch <-
        flowCore::estimateLogicle(flow_frame[flow_frame[, ch] >= q_thresh, ],
                                  channels = paste0("^", ch))
      el@transforms[[ch]] <- el_ch@transforms[[ch]]
    }
  }
  
  flow_frame <- flowCore::transform(flow_frame, el)

  list(
    flow_frame = flow_frame,
    trans_function = el
  )
}

.auroraScaleFscSsc <- function(sample) {
  # Scale SSC and FSC channels to the range in which we expect to find
  # transformed antibodies.
  
  fsc_ssc_pattern <- "(^FSC|^SSC|^Time)"
  fsc_ssc_cols <- which(grepl(fsc_ssc_pattern, sample$parameter_name))
  # Set minimum value to 0.
  for (col in fsc_ssc_cols) {
    sample$exprs[[col]] <- pmax(0, sample$exprs[[col]])
  }
  # Calculate factor based on all FSC/SSC channels.
  fsc_ssc_factor <-
    10 ^ floor(log10(quantile(unlist(sample$exprs[, fsc_ssc_cols]), 0.95)))
  sample$exprs[, fsc_ssc_cols] <- sample$exprs[, fsc_ssc_cols] / fsc_ssc_factor
  
  sample$fsc_ssc_factor <- fsc_ssc_factor

  sample
}

#' Preprocess Aurora sample.
#'
#' Quality control and cleaning for Aurora sample.
#'
#' Transformation of Aurora data should be done before converting the flow frame
#' into an Astrolabe sample.
#'
#' @inheritParams preprocess
#' @seealso \code{\link{auroraTransformChannels}}
#' @return Sample after the above steps are done.
auroraPreprocess <- function(sample) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  sample <- .auroraScaleFscSsc(sample)

  sample
}
