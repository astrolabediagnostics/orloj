#' Marker threshold report.
#'
#' Generate plot for the positive and negative cells of each marker.
#'
#' @param sample An Astrolabe sample.
#' @param max_n Subsample the data to a maximum of max_n cells.
#' @return An orloj report list with all of the required objects.
#' @export
reportMarkerThresholds <- function(sample, max_n = 100000) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")
  marker_thresholds <- sample$calculate_marker_thresholds
  if (is.null(marker_thresholds)) return(NULL)
  
  # Figure dimensions.
  fig_len <- 400
  # X/Y axis limit quantiles.
  limit_quantiles <- c(0.001, 0.999)
  
  exprs <- fcsExprs(sample)
  if (nrow(exprs) > max_n) exprs <- exprs[sample(seq(nrow(exprs)), max_n), ]
  
  # Decide on the X axis for the biaxial plots.
  x_axis_options <- c("CD45", "FSC_A", "Time")
  x_axis <- x_axis_options[x_axis_options %in% sample$parameter_desc][1]
  if (is.na(x_axis)) stop("cannot find channel for x_axis")
  x_axis_limit <- quantile(exprs[[x_axis]], limit_quantiles)
  
  # Decide on Y-axis limit, the limit that fits the widest channel.
  y_axis_limits <- lapply(names(marker_thresholds), function(channel_name) {
    quantile(exprs[[channel_name]], limit_quantiles)
  })
  y_axis_limit <- 
    c(min(unlist(lapply(y_axis_limits, function(l) l[1]))),
      max(unlist(lapply(y_axis_limits, function(l) l[2]))))
  
  # Figure: Biaxial plot of X-axis versus marker for each marker.
  lapply(orloj::nameVector(names(marker_thresholds)), function(channel_name) {
    pos_indices <- exprs[[channel_name]] >= marker_thresholds[[channel_name]]
    
    # Calculate bandwidth for geom_density2d. This is required in cases where
    # one of the columns could end with a bandwidth of 0.
    h <-
      c(MASS::bandwidth.nrd(exprs[[x_axis]]),
        MASS::bandwidth.nrd(exprs[[channel_name]]))
    h[h == 0] <- 0.1
    
    plt_title <- 
      paste0("Positive cells for ", channel_name,
             " (", round(mean(pos_indices) * 100, 1), "%)")
    
    # Generate figure.
    plt <- 
      ggplot(mapping = aes_string(x = x_axis, y = channel_name)) +
      geom_point(data = exprs, alpha = 0.1, color = "grey70", size = 0) +
      geom_point(
        data = exprs[pos_indices, ], alpha = 0.5, color = "#1C3C44", size = 1) +
      geom_density2d(data = exprs, h = h, color = "grey20") +
      labs(title = plt_title) +
      xlim(x_axis_limit) + ylim(y_axis_limit) +
      theme_linedraw() +
      theme(aspect.ratio = 1)
    list(plt = plt, width = fig_len, height = fig_len)
  })
}
