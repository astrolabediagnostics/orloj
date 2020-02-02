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

.auroraScaleFscSsc <- function(sample, new_max = 7) {
  # Scale SSC and FSC channels to the range in which we expect to find
  # transformed antibodies.
  
  fsc_ssc_pattern <- "(^FSC|^SSC|^Time)"
  fsc_ssc_cols <- which(grepl(fsc_ssc_pattern, sample$parameter_name))
  for (col in fsc_ssc_cols) {
    sample$exprs[[col]] <- pmax(0, sample$exprs[[col]])
    v <- sample$exprs[[col]]
    v <- v[v < quantile(v, 0.99)]
    sample$exprs[[col]] <- sample$exprs[[col]] / max(v) * new_max
  }

  sample
}

#' Preprocess Aurora sample.
#'
#' Quality control and cleaning for Aurora sample.
#'
#' Transformation of Aurora data should be done before converting the flow frame
#' into an Astrolabe sample.
#'
#' @param sample An Astrolabe sample.
#' @seealso \code{\link{auroraTransformChannels}}
#' @return Sample after the above steps are done.
#' @export
auroraPreprocess <- function(sample) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  sample <- .auroraScaleFscSsc(sample)

  sample
}
