#' Non-cells report.
#'
#' For Aurora, generates plots that report on events with extreme FSC-A values,
# high FSC-W values, or high LD values.
#' For mass cytometry, generate plot for live/dead.
#'
#' @param sample An Astrolabe sample.
#' @return An orloj report list with all of the required objects.
#' @export
reportNonCells <- function(sample) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  if (sample$instrument == "aurora") {
    .reportNonCellsAurora(sample)
  } else if(sample$instrument == "mass_cytometry") {
    .reportNonCellsMassCytometry(sample)
  } else {
    NULL
  }
}

.reportNonCellsAurora <- function(sample) {
  if (nrow(sample$exprs) == length(sample$non_cell_indices)) {
    # No non-cell indices found, nothing to report.
    return(NULL)
  }

  exprs <- orloj::fcsExprs(sample, keep_debris = TRUE)

  if (!all(c("FSC_A", "FSC_W", "LD") %in% colnames(exprs))) {
    # One or more Aurora columns missing, nothing to report.
    return(NULL)
  }

  report <- list()

  # Get FSC-A, FSC-W, and LD, and set the various non-cell events.
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

.reportNonCellsMassCytometry <- function(sample) {
  exprs <- sample$exprs
  
  if (is.null(sample$live_indices) | is.null(sample$livedead_col)) {
    # Reverse compatibility: Older experiments did not have live_indices.
    sample$live_indices <- seq(nrow(exprs))
  }
  
  if (length(sample$live_indices) == nrow(exprs)) {
    # All cells are alive, nothing to report.
    return(NULL)
  }
  
  report <- list()
  
  exprs$EventType <- "Dead"
  exprs$EventType[sample$live_indices] <- "Alive"
  
  # Figure: DNA versus cisplatin.
  per_alive <- length(sample$live_indices) / nrow(exprs)
  exprs$DNA <- exprs[[sample$dna_col_idx]]
  exprs$LiveDeadStaining <- exprs[[sample$livedead_col]]
  plt <-
    ggplot(mapping = aes(x = DNA, y = LiveDeadStaining)) +
    geom_point(data = exprs, alpha = 0.1, color = "grey70", size = 0) +
    geom_point(data = exprs[sample$live_indices, ],
               alpha = 0.5, color = "#1C3C44", size = 1) +
    geom_density2d(data = exprs, color = "grey20") +
    labs(title =
           paste0(prettyNum(length(sample$live_indices), big.mark = ","), " (",
                  round(per_alive * 100, 1), "%) live events"),
         y = "Live/Dead Staining (Pt195)") +
    theme(aspect.ratio = 1)
  report$LiveDead <- list(plt = plt, width = 600, height = 600)
  
  report
}