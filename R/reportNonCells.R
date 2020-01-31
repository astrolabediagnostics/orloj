#' Non-cells report.
#'
#' Generate plots for debris, doublets, and live/dead events.
#'
#' @param sample An Astrolabe sample.
#' @param max_n Subsample the data to a maximum of max_n cells.
#' @return An orloj report list with all of the required objects.
#' @export
reportNonCells <- function(sample, max_n = Inf) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  # Figure dimensions.
  fig_len <- 400

  if (length(sample$debris_indices) == 0 &&
        length(sample$doublet_indices) == 0) return(NULL)
  
  # Instrument-specific report.
  if (sample$instrument == "mass_cytometry") {
    report <-
      .reportMassCytometryDebrisDoublets(sample, fig_len = fig_len,
                                         max_n = max_n)
  } else if (sample$instrument %in% c("aurora", "lsr_fortessa")) {
    report <-
      .reportFlowCytometryDebrisDoublets(sample, fig_len = fig_len,
                                         max_n = max_n)
  } else return(NULL)
  
  if (is.null(report)) return(report)

  # Live/Dead plot.
  report$LiveDead <- .plotLiveDead(sample, fig_len = fig_len, max_n = max_n)
    
  report
}

.reportMassCytometryDebrisDoublets <- function(sample, fig_len = 400,
                                               max_n = Inf) {
  # Generate debris and doublet report for mass cytometry data.
  
  # Set standard naming for Event Length and DNA columns.
  exprs <- fcsExprs(sample, keep_debris = TRUE, keep_dead = TRUE)
  standard_channels <- findStandardMassCytometryChannels(sample)
  event_length_idx <- standard_channels$event_length_idx
  dna191_idx <- standard_channels$dna191_idx
  # Don't return a report if DNA or event length are missing.
  if (is.null(dna191_idx) || is.null(event_length_idx)) return(NULL)
  colnames(exprs)[event_length_idx] <- "Event Length"
  colnames(exprs)[dna191_idx] <- "DNA"

  # Decide on axis limits.
  dna_lim <- c(0, ceiling(max(exprs[[dna191_idx]]) / 0.25) * 0.25)
  event_length_lim <-
    c(0, ceiling(max(exprs[[event_length_idx]]) / 0.25) * 0.25)

  report <- list()

  # Generate preprocessing debris and doublet figure.
  exprs$EventType <- "Cell"
  exprs$EventType[sample$debris_indices] <- "Debris"
  exprs$EventType[sample$doublet_indices] <- "Doublet"
  exprs$EventType <-
    factor(exprs$EventType, levels = c("Cell", "Debris", "Doublet"))

  if (nrow(exprs) > max_n) exprs <- exprs[sample(seq(nrow(exprs)), max_n), ]

  # Figure: Event Length versus DNA.
  plt <- 
    ggplot(exprs, aes(x = `Event Length`, y = DNA)) +
    geom_point(alpha = 0.1, color = "grey70", size = 1) +
    labs(title = "Preprocessing Debris and Doublets", x = "Event Length") +
    xlim(event_length_lim) + ylim(dna_lim) +
    facet_wrap(~ EventType) +
    theme_linedraw() +
    theme(aspect.ratio = 1)
  width <- length(unique(exprs$EventType)) * fig_len
    
  report$PreprocessingDebris <- list(plt = plt, width = width, height = fig_len)

  # Generate Ek'balam debris and Root_unassigned figure.
  exprs_clean <-
    exprs[setdiff(seq(nrow(exprs)),
                  c(sample$debris_indices, sample$doublet_indices)), ]
  exprs_clean$EventType <- "Cell"
  cm_neg_indices <-
    exprs_clean$Assignment == "CM-" |
    (grepl("CM-", exprs_clean$Assignment) &
      grepl("_unassigned", exprs_clean$Assignment))
  exprs_clean$EventType[cm_neg_indices] <- "Canonical Marker Negative"
  exprs_clean$EventType[exprs_clean$Assignment == "Root_unassigned"] <-
    "Unassigned"
  exprs_clean$EventType <-
    factor(exprs_clean$EventType,
           levels = c("Cell", "Canonical Marker Negative", "Unassigned"))

  # Figure: Event Length versus DNA.
  plt <- 
    ggplot(exprs_clean, aes(x = `Event Length`, y = DNA)) +
    geom_point(alpha = 0.1, color = "grey70", size = 1) +
    labs(title = "Cell Labeling Debris and Unassigned", x = "Event Length") +
    xlim(event_length_lim) + ylim(dna_lim) +
    facet_wrap(~ EventType) +
    theme_linedraw() +
    theme(aspect.ratio = 1)
  width <- length(unique(exprs_clean$EventType)) * fig_len
    
  report$CellLabelingDebris <- list(plt = plt, width = width, height = fig_len)
 
  report
}

.reportFlowCytometryDebrisDoublets <- function(sample, fig_len = 400,
                                               max_n = Inf) {
  # Generate debris and doublet report for flow cytometry data.
  exprs <- fcsExprs(sample, keep_debris = TRUE, keep_dead = TRUE)
  
  # Collect the FSC/SSC channels used for cleaning.
  fsc_a <- grep("FSC_A", sample$parameter_desc, fixed = TRUE, value = TRUE)
  fsc_h <- grep("FSC_H", sample$parameter_desc, fixed = TRUE, value = TRUE)
  ssc_a <- grep("SSC_A", sample$parameter_desc, fixed = TRUE, value = TRUE)

  if (length(fsc_a) != 1 || length(ssc_a) != 1) return(NULL)
  exprs$FSC_A <- exprs[[fsc_a]]
  if (length(fsc_h) == 1) exprs$FSC_H <- exprs[[fsc_h]]
  exprs$SSC_A <- exprs[[ssc_a]]
  
  # Populate with debris and doublets from cleaning step.
  exprs$Debris <- FALSE
  exprs$Debris[sample$debris_indices] <- TRUE
  exprs$Doublet <- FALSE
  exprs$Doublet[sample$doublet_indices] <- TRUE

  if (nrow(exprs) > max_n) exprs <- exprs[sample(seq(nrow(exprs)), max_n), ]

  # Figure: Debris, FSC_A versus SSC_A.
  fig_title <- 
    paste0(prettyNum(sum(exprs$Debris), big.mark = ","), " (",
           round(mean(exprs$Debris) * 100, 1), "%) debris events")
  plt_debris <- 
    ggplot(mapping = aes(x = FSC_A, y = SSC_A)) +
    geom_point(data = exprs, alpha = 0.1, color = "grey70", size = 0) +
    geom_point(data = subset(exprs, Debris),
               alpha = 0.5, color = "#1C3C44", size = 0) +
    geom_density2d(data = exprs, color = "grey20") +
    labs(title = fig_title) +
    theme_linedraw() +
    theme(aspect.ratio = 1)
  
  # Figure: Doublets, FSC_A versus FSC_H.
  if (length(fsc_h) != 1) {
    # No FSC_H, no doublet report.
    plt <- plt_debris
    width <- fig_len
    height <- fig_len
  } else {
    fig_title <- 
      paste0(prettyNum(sum(exprs$Doublet), big.mark = ","), " (",
             round(mean(exprs$Doublet) * 100, 1), "%) doublet events")
    plt_doublets <- 
      ggplot(mapping = aes(x = FSC_A, y = FSC_H)) +
      geom_point(data = subset(exprs, !Debris),
                 alpha = 0.1, color = "grey70", size = 0) +
      geom_point(data = subset(exprs, !Debris & Doublet),
                 alpha = 0.5, color = "#1C3C44", size = 0) +
      geom_density2d(data = exprs, color = "grey20") +
      labs(title = fig_title) +
      theme_linedraw() +
      theme(aspect.ratio = 1)
    plt <- plt_debris + plt_doublets
    width <- fig_len * 2
    height <- fig_len
  }
  
  list(PreprocessingDebris = list(plt = plt, width = width, height = height))
}

.plotLiveDead <- function(sample, fig_len = 400, max_n = Inf) {
  # Plot instrument-specific X axis marker versus Live/Dead.
  if (length(sample$live_dead_channel_name) == 0) return(NULL)
  
  # Reorganize data into a single data frame with generic column names.
  x_channel_idx <-
    which(sample$parameter_name == sample$live_dead_x_channel_name)
  live_dead_channel_idx <-
    which(sample$parameter_name == sample$live_dead_channel_name)
  exprs <- fcsExprs(sample, keep_dead = TRUE)
  if (nrow(exprs) > max_n) exprs <- exprs[sample(seq(nrow(exprs)), max_n), ]
  df <- data.frame(
    EventType =
      factor(ifelse(exprs$Dead, "Dead", "Alive"), levels = c("Alive", "Dead")),
    X = exprs[[x_channel_idx]],
    LiveDead = exprs[[live_dead_channel_idx]]
  )
  
  # Figure: Instrument-specific X versus Live/Dead.
  fig_title <- 
    paste0(prettyNum(sum(!exprs$Dead), big.mark = ","), " (",
           round(mean(!exprs$Dead) * 100, 1), "%) live events")
  plt <-
    ggplot(mapping = aes(x = X, y = LiveDead)) +
    geom_point(data = df, alpha = 0.1, color = "grey70", size = 0) +
    geom_point(data = subset(df, EventType == "Alive"),
               alpha = 0.5, color = "#1C3C44", size = 1) +
    geom_density2d(data = df, color = "grey20") +
    labs(title = fig_title,
         x = sample$live_dead_x_channel_name,
         y = sample$live_dead_channel_name) +
    theme_linedraw() +
    theme(aspect.ratio = 1)
  list(plt = plt, width = fig_len, height = fig_len)
}
