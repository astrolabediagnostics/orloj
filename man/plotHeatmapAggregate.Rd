% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotHeatmapAggregate}
\alias{plotHeatmapAggregate}
\title{Aggregate data and generate a heat map object.}
\usage{
plotHeatmapAggregate(
  data,
  x,
  y,
  value,
  func = median,
  type = NULL,
  title = NULL,
  x_axis_order = NULL,
  y_axis_order = NULL,
  theme = NULL
)
}
\arguments{
\item{data}{Data frame to be aggregated and plotted.}

\item{x}{Column names for X-axis and Y-axis, respectively.}

\item{y}{Column names for X-axis and Y-axis, respectively.}

\item{value}{Column name for tile values.}

\item{type}{Heatmap type (NULL, cluster_labels, abundance, or change). This
will influence formatting.}

\item{title}{Plot title.}

\item{x_axis_order}{Order of X- and Y-axis tick labels. If none
specified, \code{\link[gtools]{mixedsort}} will be used.}

\item{y_axis_order}{Order of X- and Y-axis tick labels. If none
specified, \code{\link[gtools]{mixedsort}} will be used.}

\item{theme}{Modifications to the default ggplot theme.}
}
\value{
An orloj plot list with the plot object and any other parameters that
are required to export it.
}
\description{
Given a data frame in long format and X- and Y-axes, calculate the mean value
of a value column for each (x, y) combination and plot as a heatmap.
}
