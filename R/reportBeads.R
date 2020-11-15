#' Bead report.
#'
#' Generates plots that report on the bead detection.
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

  # Rename True/False to Bead/Cell.
  exprs$AstrolabeBead <-
    factor(ifelse(exprs$AstrolabeBead, "AstrolabeBead", "Cell"),
           levels = c("Cell", "AstrolabeBead"))

  if (sample$instrument == "mass_cytometry") {
    .reportMassBeads(sample, exprs)
  } else {
    .reportFlowBeads(sample, exprs)
  }
}

.reportMassBeads <- function(sample, exprs) {
  # Set up expression with standard bead channel names.
  bead_channel_indices <- match(beadChannelNames(), sample$parameter_name)
  existing_beads <- which(!is.na(bead_channel_indices))
  bead_channel_indices <- bead_channel_indices[existing_beads]
  bead_channel_names <- beadChannelNames()[existing_beads]
  
  bead_channel_idx <- which(colnames(exprs) == "AstrolabeBead")
  exprs <- exprs[, c(bead_channel_indices, bead_channel_idx)]
  colnames(exprs) <- c(bead_channel_names, "AstrolabeBead")
  
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

.reportFlowBeads <- function(sample, exprs) {
  # Reorganize exprs to include only relevant columns.
  astrolabe_bead_idx <- which(colnames(exprs) == "AstrolabeBead")
  bead_channel_idx <- sample$bead_channel_idx
  fsc_a_idx <- which(grepl("FSC_A", sample$parameter_desc, fixed = TRUE))
  if (length(fsc_a_idx) != 1) stop("unable to find FSC_A channel")
  exprs <- exprs[, c(fsc_a_idx, bead_channel_idx, astrolabe_bead_idx)]
  colnames(exprs)[2] <- "Bead"
  # Figure: Scatter plot of FSC versus bead channel, color-coded by Astrolabe's
  # bead identification.
  plotScatterPlot(exprs, x = "Bead", y = "FSC_A", color = "AstrolabeBead")
}
