#' Non-cells report.
#'
#' For Aurora, generates plots that report on events with extreme FSC-A values,
# high FSC-W values, or high LD values.
#' For mass cytometry, generate plot for debris and live/dead.
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

  stop("This code does not use the new fcsExprs structure")

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
  exprs <- fcsExprs(sample, keep_debris = TRUE, keep_dead = TRUE)

  # Set standard naming for Event Length and DNA columns.
  standard_channels <- findStandardMassCytometryChannels(sample)
  event_length_idx <- standard_channels$event_length_idx
  dna191_idx <- standard_channels$dna191_idx
  colnames(exprs)[event_length_idx] <- "Event Length"
  colnames(exprs)[dna191_idx] <- "DNA"

  # Decide on axis limits.
  dna_lim <- c(0, ceiling(max(exprs[[dna191_idx]])))
  event_length_lim <- c(0, ceiling(max(exprs[[event_length_idx]])))

  report <- list()

  if (any(exprs$Debris)) {
    # Generate debris and doublet figure.
    exprs$EventType <- "Cell"
    exprs$EventType[sample$debris_indices] <- "Debris"
    exprs$EventType[sample$doublet_indices] <- "Doublet"
    exprs$EventType <-
      factor(exprs$EventType, levels = c("Cell", "Debris", "Doublet"))

    # Figure: Event Length versus DNA.
    plt <- 
      ggplot(exprs, aes(x = `Event Length`, y = DNA)) +
      geom_point(alpha = 0.1, color = "grey70", size = 0) +
      geom_point(data = exprs[exprs$EventType != "Cell", ],
                 mapping = aes(color = EventType), alpha = 0.1, size = 1) +
      labs(title = "Preprocessing Debris and Doublets", x = "Event Length") +
      xlim(event_length_lim) + ylim(dna_lim) +
      theme(aspect.ratio = 1,
            legend.position = "bottom") +
      guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
      
    report$PreprocessingDebris <- list(plt = plt, width = 600, height = 600)
  }


  # Generate Ek'balam debris and Root_unassigned figure.
  exprs_clean <- exprs[!exprs$Debris & !exprs$Dead, ]
  exprs_clean$EventType <- "Cell"
  exprs_clean$EventType[exprs_clean$Assignment == "Debris"] <- "Debris"
  exprs_clean$EventType[exprs_clean$Assignment == "Root_unassigned"] <-
    "Unassigned"
  exprs_clean$EventType <-
    factor(exprs_clean$EventType, levels = c("Cell", "Debris", "Unassigned"))

  # Figure: Event Length versus DNA.
  plt <- 
    ggplot(exprs_clean, aes(x = `Event Length`, y = DNA)) +
    geom_point(alpha = 0.1, color = "grey70", size = 0) +
    geom_point(data = exprs_clean[exprs_clean$EventType != "Cell", ],
               mapping = aes(color = EventType), alpha = 0.1, size = 1) +
    labs(title = "Cell Labeling Debris and Unassigned", x = "Event Length") +
    xlim(event_length_lim) + ylim(dna_lim) +
    theme(aspect.ratio = 1,
          legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
  report$CellLabelingDebris <- list(plt = plt, width = 600, height = 600)

  
  if (any(exprs$Dead)) {
    # Generate live/dead report.
    exprs$EventType <- "Alive"
    exprs$EventType[exprs$Dead] <- "Dead"
    exprs$EventType <- factor(exprs$EventType, levels = c("Alive", "Dead"))
  
    # Figure: DNA versus cisplatin.
    n_alive <- sum(!exprs$Dead)
    per_alive <- mean(!exprs$Dead)
    exprs$DNA <- exprs[[sample$dna_col_idx]]
    exprs$LiveDeadStaining <- exprs[[sample$livedead_col]]
    plt <-
      ggplot(mapping = aes(x = DNA, y = LiveDeadStaining)) +
      geom_point(data = exprs, alpha = 0.1, color = "grey70", size = 0) +
      geom_point(data = exprs[sample$live_indices, ],
                 alpha = 0.5, color = "#1C3C44", size = 1) +
      geom_density2d(data = exprs, color = "grey20") +
      labs(title =
             paste0(prettyNum(n_alive, big.mark = ","), " (",
                    round(per_alive * 100, 1), "%) live events"),
           y = "Live/Dead Staining (Pt195)") +
      theme(aspect.ratio = 1)
    report$LiveDead <- list(plt = plt, width = 600, height = 600)
  }
    
  report
}
