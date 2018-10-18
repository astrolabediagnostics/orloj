#' Bead report.
#'
#' Generates plots that report on the beads in a mass cytometry file.
#'
#' @param sample An Astrolabe sample.
#' @return An orloj report list with all of the required objects.
#' @export
reportBeads <- function(sample) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  exprs <-
    fcsExprs(sample, keep_beads = TRUE, keep_dead = TRUE, keep_debris = TRUE)

  # No beads were found, nothing to report.
  if (sum(exprs$Bead) == 0) return(NULL)

  # Set up expression with standard bead channel names.
  bead_channel_indices <- match(beadChannelNames(), sample$parameter_name)
  existing_beads <- which(!is.na(bead_channel_indices))
  bead_channel_indices <- bead_channel_indices[existing_beads]
  bead_channel_names <- beadChannelNames()[existing_beads]

  bead_channel_idx <- which(colnames(exprs) == "Bead")
  exprs <- exprs[, c(bead_channel_indices, bead_channel_idx)]
  colnames(exprs) <- c(bead_channel_names, "Bead")

  # Figure: First bead channel versus all of the other bead channels
  x <- bead_channel_names[1]
  ys <- setdiff(bead_channel_names, x)
  report <- lapply(ys, function(y) {
    plotScatterPlot(exprs, x, y, "Bead")
  })
  names(report) <- paste0(x, "_vs_", ys)

  report
}
