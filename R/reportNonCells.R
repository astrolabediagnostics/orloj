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

  exprs <- orloj::fcsExprs(sample, keep_debris = TRUE)
  # Collect the various non-cell events.
  df <- exprs[, c("FSC_A", "FSC_W", "LD")]
  df$EventType <- c("Cell")
  df$EventType[sample$extreme_fsc_a_indices] <- "FSC-A High/Low"
  df$EventType[sample$high_fsc_w_indices] <- "FSC-W High (Doublet)"
  df$EventType[sample$high_ld_indices] <- "LD High (Dead)"

  # Figure: FSC-A versus FSC-W.
  report$FSC_A_vs_FSC_W <-
    orloj::plotScatterPlot(df, "FSC_A", "FSC_W", "EventType")
  # Figure: FSC-A versus LD.
  report$FSC_A_vs_LD <- 
    orloj::plotScatterPlot(df, "FSC_A", "LD", "EventType")
  # Table: Event type counts.
  counts <- df %>%
    dplyr::count(EventType) %>%
    dplyr::mutate(Per = n / sum(n))
  report$Counts <- list(data = counts)

  report
}
