#' Bead report.
#'
#' Generates plots that report on the beads in a mass cytometry file.
#'
#' @param sample An Astrolabe sample.
#' @return An orloj report list with all of the required objects.
#' @import patchwork
#' @export
reportBeads <- function(sample) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  exprs <-
    fcsExprs(sample, keep_beads = TRUE, keep_dead = TRUE, keep_debris = TRUE)
  
  # No beads were found, nothing to report.
  if (sum(exprs$AstrolabeBead) == 0) return(NULL)
  
  # Set up expression with standard bead channel names.
  bead_channel_indices <- match(beadChannelNames(), sample$parameter_name)
  existing_beads <- which(!is.na(bead_channel_indices))
  bead_channel_indices <- bead_channel_indices[existing_beads]
  bead_channel_names <- beadChannelNames()[existing_beads]
  
  bead_channel_idx <- which(colnames(exprs) == "AstrolabeBead")
  exprs <- exprs[, c(bead_channel_indices, bead_channel_idx)]
  colnames(exprs) <- c(bead_channel_names, "AstrolabeBead")
  
  # Rename True/False to Bead/Cell.
  exprs$AstrolabeBead <-
    factor(ifelse(exprs$AstrolabeBead, "AstrolabeBead", "Cell"),
           levels = c("Cell", "AstrolabeBead"))
  
  # Figure: First bead channel (on Y-axis) versus all of the other channels. Use
  # patchwork to combine all of them to single plot.
  y <- bead_channel_names[1]
  ylim <- c(0, ceiling(max(exprs[[y]]) / 0.25) * 0.25)
  xs <- setdiff(bead_channel_names, y)
  plt <- NULL
  width <- 0
  height <- 0
  for (x in xs) {
    x_plt <- plotScatterPlot(exprs, x, y, "AstrolabeBead", ylim = ylim)
    if (is.null(plt)) {
      plt <- x_plt$plt
    } else {
      plt <- plt + x_plt$plt
    }
    width <- width + x_plt$width
    height <- x_plt$height
  }
  plt <- plt + patchwork::plot_layout(nrow = 1)
  
  # Create a one-figure report.
  list(bead_report = list(plt = plt, width = width, height = height))
}
