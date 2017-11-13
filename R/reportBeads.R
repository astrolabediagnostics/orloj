#' Bead report.
#'
#' Generates plots and CSV files that report on the beads in a mass cytometry
#' file.
#'
#' @param sample An Astrolabe sample.
#' @return An orloj report list with all of the required objects.
#' @export
reportBeads <- function(sample) {
  # TODO statistics report
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  if (sample$source != "mass_cytometry") stop("Expecting mass cytometry data")

  # Get expression data and mark beads.
  bead_channel_indices <- match(orloj_bead_channel_names, sample$parameter_name)
  exprs <- sample$exprs[, bead_channel_indices]
  colnames(exprs) <- orloj_bead_channel_names
  exprs$Bead <- TRUE
  exprs$Bead[sample$non_bead_indices] <- FALSE

  # Plot: First bead channel versus all of the other bead channels
  x <- orloj_bead_channel_names[1]
  ys <- setdiff(orloj_bead_channel_names, x)
  report <- lapply(ys, function(y) {
    plotScatterPlot(exprs, x, y, "Bead")
  })
  names(report) <- paste0(x, "_vs_", ys)

  report
}
