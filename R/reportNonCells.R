#' Non-cells report.
#'
#' Generates plots that report on events with extreme FSC-A values, high FSC-W
#' values, or high LD values.
#'
#' @param sample An Astrolabe sample.
#' @return An orloj report list with all of the required objects.
#' @export
reportNonCells <- function(sample) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  # Only Aurora data has non-cell indices.
  if (sample$instrument != "aurora") return(NULL);

  if (nrow(sample$exprs) == length(sample$non_cell_indices)) {
    # No beads were found, nothing to report.
    return(NULL);
  }

  report <- list()

  # Get FSC-A, FSC-W, and LD, and set the various non-cell events.
  exprs <- orloj::fcsExprs(sample, keep_debris = TRUE)
  exprs <- exprs[, c("FSC_A", "FSC_W", "LD")]
  exprs$EventType <- c("Cell")
  exprs$EventType[sample$extreme_fsc_a_indices] <- "FSC-A High/Low"
  exprs$EventType[sample$high_fsc_w_indices] <- "FSC-W High (Doublet)"

  # Figure: FSC-A versus FSC-W.
  report$FSC_A_vs_FSC_W <-
    orloj::plotScatterPlot(exprs, "FSC_A", "FSC_W", "EventType",
                           size = 0, alpha = 0.1)
  # Table: Event type counts (after adding LD high).
  exprs$EventType[sample$high_ld_indices] <- "LD High (Dead)"
  counts <- exprs %>%
    dplyr::count(EventType) %>%
    dplyr::mutate(Per = n / sum(n))
  report$Counts <- list(data = counts)
  # Figure: FSC-A versus LD.
  exprs <- dplyr::filter(exprs, EventType %in% c("Cell", "LD High (Dead)"))
  report$FSC_A_vs_LD <- 
    orloj::plotScatterPlot(exprs, "FSC_A", "LD", "EventType",
                           size = 0, alpha = 0.1)

  report
}
