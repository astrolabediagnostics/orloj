# plotting.R
# Functions for plotting various figures of cytometry data and its analysis.

#' Export an orloj plot list.
#'
#' Exports an orloj plot list, applying any of the parameters that it includes.
#'
#' @param filename The name of the file to be exported. Note that the file
#' extension decides on the file type (such as png, svg, etc.).
#' @param plt_list An orloj plot list.
#' @param dpi Plot resolution. Only applies only to raster output types.
#' @import ggplot2
#' @export
exportPlot <- function(filename, plt_list, dpi = 100) {
  if (is.null(plt_list$plt)) {
    if (is.null(plt_list$data)) {
      stop("plt_list object must include either plt or data or both")
    }
  } else {
    # Plot included, check for width and height.
    if (is.null(plt_list$width)) {
      stop("plt_list object missing width")
    }
    if (is.null(plt_list$height)) {
      stop("plt_list object missing height")
    }
  }

  if (!is.null(plt_list$plt)) {
    # Export plot.
    ggsave(filename,
           plt_list$plt,
           width = plt_list$width / dpi,
           height = plt_list$height / dpi,
           limitsize = FALSE)
  }
  if (!is.null(plt_list$data)) {
    # Export data as CSV.
    write.csv(plt_list$data, paste0(filename, ".csv"))
  }
}

#' Export an orloj report.
#'
#' Exports an orloj report to a destination directory. Reports take orloj data
#' structures and generate a series of plots that summarize, illustrate, or
#' explain an analysis.
#'
#' @param dir The name of the dir where figures are to be exported.
#' @param report An orloj report.
#' @param file_format Exported plots file format (such as png, svg, etc.).
#' @inherit exportPlot
#' @export
exportReport <- function(dir, report, file_format = "png", dpi = 100) {
  for (plot_name in names(report)) {
    plot_contents <- names(report[[plot_name]])
    if (is.null(plot_contents)) {
      # Skip empty plots.
      next
    } else if ("plt" %in% plot_contents || "data" %in% plot_contents) {
      # A single plot.
      plot_filename <-
        file.path(dir, paste0(gsub("/", "_", plot_name), ".", file_format))
      exportPlot(plot_filename, report[[plot_name]], dpi)
    } else {
      # Another layer of plots, recursively export them.
      plot_path <- file.path(dir, plot_name)
      dir.create(plot_path)
      exportReport(plot_path, report[[plot_name]], file_format, dpi)
    }
  }
}

#' Bar plot object.
#'
#' Generates a bar plot (using geom_col) ggplot object.
#'
#' @param data Dataset to use for the plot. If sample, exprs will be used.
#' @param x,y Column names for the X-axis and Y-axis, respectively.
#' @param fill Column name for bar fill.
#' @param title Plot title.
#' @param scale_y_labels Scaling function for Y-axis tick labels.
#' @param theme Modifications to the default ggplot theme.
#' @return An orloj plot list with the plot object and any other parameters that
#' are required to export it.
#' @import ggplot2
#' @export
plotBarPlot <- function(data,
                        x,
                        y,
                        fill = NULL,
                        title = NULL,
                        scale_y_labels = NULL,
                        theme = NULL) {
  if (isSample(data)) {
    data <- fcsExprs(data)
  }

  # Add backticks to all variable names to avoid aes_string issues.
  x_bt      <- addBackticks(x)
  y_bt      <- addBackticks(y)
  fill_bt   <- addBackticks(fill)

  # Generate the plot.
  plt <- ggplot(data, aes_string(x = x_bt, y = y_bt))
  if (is.null(fill)) {
    plt <- plt + geom_col()
  } else {
    plt <- plt + geom_col(aes_string(fill = fill_bt))
  }
  plt <- plt +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.3))

  if (!is.null(title)) plt <- plt + labs(title = title)
  if (!is.null(scale_y_labels)) {
    plt <- plt + scale_y_continuous(labels = scale_y_labels)
  }
  if (!is.null(theme)) plt <- plt + theme

  # Set up width, height, data.
  width <-
    15 + # Y-axis title
    15 * max(nchar(as.character(data[[y]]))) + # Longest Y-axis label
    40 * length(unique(data[[x]])) # Number of X-axis values
  height <- 600
  if (is.null(fill)) {
    data <- data[, c(x, y)]
  } else {
    data <- data[, c(x, fill, y)]
  }

  list(
    plt = plt,
    width = width,
    height = height,
    data = data
  )
}

#' Scatter plot object.
#'
#' Generates a biaxial scatter plot ggplot object.
#'
#' @param data Dataset to use for the plot. If sample, exprs will be used.
#' @param x,y Column names for X-axis and Y-axis, respectively.
#' @param color Column name for the color aesthetic.
#' @param alpha geom_point alpha aesthetic.
#' @param title Plot title.
#' @param theme Modifications to the default ggplot theme.
#' @param force_data If false, plots with more than 100 points won't have plot
#' data included.
#' @return An orloj plot list with the plot object and any other parameters that
#' are required to export it.
#' @import ggplot2
#' @export
plotScatterPlot <- function(data,
                            x,
                            y,
                            color = NULL,
                            alpha = 0.5,
                            size = 1,
                            title = NULL,
                            theme = NULL,
                            force_data = FALSE) {
  if (isSample(data)) {
    data <- fcsExprs(data)
  }

  # Add backticks to all variable names to avoid aes_string issues.
  x_bt      <- addBackticks(x)
  y_bt      <- addBackticks(y)
  color_bt  <- addBackticks(color)

  # Generate the plot.
  plt <- ggplot(data, aes_string(x = x_bt, y = y_bt))
  if (is.null(color)) {
    plt <- plt + geom_point(alpha = alpha)
  } else {
    plt <-
      plt +
      geom_point(aes_string(color = color_bt),
                 alpha = alpha,
                 size = size) +
      guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
      theme(legend.position = "bottom")
  }

  if (!is.null(title)) plt <- plt + labs(title = title)
  if (!is.null(theme)) plt <- plt + theme

  # Default width and height for scatter plots.
  width <- 600
  height <- 600
  # Generate data.
  if (nrow(data) > 100 & !force_data) {
    data <- NULL
  } else {
    data <- data[, c(x, y)]
  }

  list(
    plt = plt,
    width = width,
    height = height,
    data = data
  )
}

#' Box plot object.
#'
#' @param data Dataset to use for the plot. If sample, exprs will be used.
#' @param x,y Column names for X-axis and Y-axis, respectively.
#' @param title Plot title.
#' @param scale_y_labels Scaling function for Y-axis tick labels.
#' @param theme Modifications to the default ggplot theme.
#' @return An orloj plot list with the plot object and any other parameters that
#' are required to export it.
#' @import ggplot2
#' @export
#'
plotBoxPlot <- function(data,
                        x,
                        y,
                        title = NULL,
                        scale_y_labels = NULL,
                        theme = NULL) {
  if (isSample(data)) {
    data <- fcsExprs(data)
  }

  # Add backticks to all variable names to avoid aes_string issues.
  x_bt <- addBackticks(x)
  y_bt <- addBackticks(y)

  plt <-
    ggplot(data, aes_string(x = x_bt, y = y_bt)) +
    geom_boxplot()

  if (!is.null(title)) plt <- plt + labs(title = title)
  if (!is.null(scale_y_labels)) {
    plt <- plt + scale_y_continuous(labels = scale_y_labels)
  }
  if (!is.null(theme)) plt <- plt + theme

  # Set up width, height, data.
  width <-
    15 + # Y-axis title
    15 * max(nchar(as.character(data[[y]]))) + # Longest Y-axis label
    40 * length(unique(data[[x]])) # Number of X-axis values
  height <- 600
  data <- data[, c(x, y)]

  list(
    plt = plt,
    width = width,
    height = height,
    data = data
  )
}

#' Heat map object.
#'
#' @param hm Data frame with heatmap data.
#' @param x,y Column names for X-axis and Y-axis, respectively.
#' @param value Column name for tile values.
#' @param type Heatmap type. Accepted values are NULL, "cluster_labels", and
#' "abundance". This will set some of the heatmap's formatting.
#' @param title Plot title.
#' @param x_axis_order,y_axis_order Order of X- and Y-axis tick labels. If none
#' specified, \code{\link[gtools]{mixedsort}} will be used.
#' @param theme Modifications to the default ggplot theme.
#' @return An orloj plot list with the plot object and any other parameters that
#' are required to export it.
#' @import ggplot2
#' @export
plotHeatmap <- function(hm,
                        x,
                        y,
                        value,
                        type = NULL,
                        title = NULL,
                        x_axis_order = NULL,
                        y_axis_order = NULL,
                        theme = NULL) {
  # Add backticks to all variable names to avoid aes_string issues.
  x_bt     <- addBackticks(x)
  y_bt     <- addBackticks(y)
  value_bt <- addBackticks(value)

  # Order according to mixedsort, unless otherwise instructed by parameter.
  if (is.null(x_axis_order)) {
    x_axis_order <- gtools::mixedsort(unique(as.character(hm[[x]])))
  }
  if (is.null(y_axis_order)) {
    y_axis_order <- gtools::mixedsort(unique(as.character(hm[[y]])))

  }
  hm[[x]] <- factor(hm[[x]], levels = x_axis_order)
  hm[[y]] <- factor(hm[[y]], levels = y_axis_order)

  # Generate the plot.
  plt <-
    ggplot(hm, aes_string(x = x_bt, y = y_bt)) +
    geom_tile(aes_string(fill = value_bt), color = "white") +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.4),
      legend.position = "bottom",
      panel.background = element_blank()
    )

  # Format according to type.
  if (!is.null(type)) {
    if (type == "cluster_labels") {
      plt <- plt +
        scale_fill_gradient(low = "white",
                            high = "#2F6148")
    } else if (type == "cluster_labels_cv") {
      plt <- plt +
        scale_fill_gradient(low = "white",
                            high = "red")
    } else if (type == "frequency") {
      plt <- plt +
        scale_fill_gradient(labels = scales::percent,
                            low = "white",
                            high = "purple")
    } else if (type == "scaled_frequency") {
      plt <- plt +
        scale_fill_gradient(low = "white",
                            high = "goldenrod4")
    } else {
      stop("Unknown heatmap type")
    }
  }

  # Add title and custom theme, if required.
  if (!is.null(title)) plt <- plt + labs(title = title)
  if (!is.null(theme)) plt <- plt + theme

  # Decide on width and height.
  width <-
    15 + # Y-axis title
    max(nchar(as.character(y_axis_order))) * 15 + # Longest Y-axis label
    15 * length(x_axis_order) # X-axis tiles
  height <-
    15 + # X-axis title
    30 + # Figure title
    max(nchar(as.character(x_axis_order))) * 15 + # Longest X-axis label
    15 * length(y_axis_order) + # Y-axis tiles
    60 # Legend

  # Generate the plt_list data object.
  data <-
    reshape2::dcast(hm,
                    as.formula(paste0(x, " ~ ", y)),
                    value.var = value)

  list(
    plt = plt,
    width = width,
    height = height,
    data = data
  )
}

#' Aggregate data and generate a heat map object.
#'
#' Given a data frame in long format and X- and Y-axes, calculate the mean value
#' of a value column for each (x, y) combination and plot as a heatmap.
#'
#' @param data Data frame to be aggregated and plotted.
#' @inheritParams plotHeatmap
#' @return An orloj plot list with the plot object and any other parameters that
#' are required to export it.
#' @import ggplot2
#' @export
plotHeatmapAggregate <- function(data,
                                 x,
                                 y,
                                 value,
                                 func = median,
                                 type = NULL,
                                 title = NULL,
                                 x_axis_order = NULL,
                                 y_axis_order = NULL,
                                 theme = NULL) {
  # Calculate mean values for each (x, y) combination.
  data$x <- data[[x]]
  data$y <- data[[y]]
  data$value <- data[[value]]
  hm <- data %>%
    dplyr::group_by(x, y) %>%
    dplyr::summarize(value = func(value)) %>%
    dplyr::ungroup()
  colnames(hm) <- c(x, y, value)

  # Generate heatmap.
  plotHeatmap(hm,
              x = x,
              y = y,
              value = value,
              type = type,
              title = title,
              x_axis_order = x_axis_order,
              y_axis_order = y_axis_order,
              theme = theme)
}
