#' MDS Map report.
#'
#' Generate plots and CSV files that report on the MDS map.
#'
#' @param experiment An Astrolabe experiment.
#' @import ggplot2 viridis ggrepel patchwork
#' @export
reportMds <- function(experiment) {
  report <- list()

  analyses <- c("Assignment", "Profiling")
  if (experiment$organism == "profiling_only") {
    analyses <- c("Profiling")
  }

  # Iterate over all analyses and generate figures for each.
  for (level in analyses) {
    mds <- experimentMds(experiment, level = level)

    report[[level]] <- list()

    # Figure: Plain map, no color-coding.
    map_obj <-
      ggplot(mds, aes(x = V1, y = V2)) +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(aes(label = CellSubset)) +
      labs(title = paste0("MDS map for ", level)) +
      theme(aspect.ratio = 1,
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.background = element_blank())
    map_plt <- list(plt = map_obj, width = 900, height = 900, data = mds)
    report[[level]]$mds_map <- map_plt

    # Remove data for all further plots, no need to duplicate it.
    map_plt$data <- NULL

    # Figures: Color-coding by each channel.
    for (channel in experiment$analysis_channels) {
      channel_filename <- filenameify(channel)

      map_obj <-
        ggplot(mds, aes(x = V1, y = V2)) +
        geom_point(aes_string(color = channel), size = 3) +
        ggrepel::geom_text_repel(aes(label = CellSubset)) +
        viridis::scale_color_viridis(direction = -1) +
        labs(title = paste0("MDS map for ", level)) +
        theme(aspect.ratio = 1,
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              legend.position = "bottom",
              panel.background = element_blank())
      map_plt <- list(plt = map_obj, width = 900, height = 900)
      report[[level]][[paste0("channel_", channel_filename)]] <- map_plt
    }

    report[[level]]$shepard_plot <- .plotShepard(experiment, level, mds)
  }

  report
}

.plotShepard <- function(experiment, level, mds) {
  # Generate a Shepard plot for a given MDS map.

  # Load pairwise distances over all samples.
  sample_rds_files <-
    lapply(orloj::nameVector(experiment$samples$SampleId), function(sample_id) {
      readRDS(
        file.path(experiment$analysis_path,
                  paste0(sample_id,
                         ".calculate_cluster_pairwise_distances_mds.RDS")))
    })

  # Combine pairwise distances into a single data frame.
  level_df <- lapply(experiment$samples$SampleId, function(sample_id) {
    sample_rds_files[[sample_id]][[level]]
  }) %>% dplyr::bind_rows()
  # Calculate pairwise distances in MDS map.
  mtx <- as.matrix(mds[, c("V1", "V2")])
  rownames(mtx) <- mds$CellSubset
  dist_mtx <- as.matrix(dist(mtx))
  # Convert to data frame.
  dist_df <- dist_mtx %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column("From") %>%
    reshape2::melt(id.vars = "From",
                   variable.name = "To",
                   value.name = "MDSDist")
  # Combine with distances above.
  dist_df <- dplyr::left_join(level_df, dist_df, by = c("From", "To"))
  dist_df <- dplyr::filter(dist_df, !is.na(Dist), Dist > 0)

  pearson_rho <- round(with(dist_df, cor(Dist, MDSDist)), 2)

  plt_title <- paste0("MDS Map Shepard Plot (Pearson rho = ", pearson_rho, ")")
  plt <- orloj::plotScatterPlot(dist_df, x = "Dist", y = "MDSDist", alpha = 0.1)
  plt$plt <-
    plt$plt +
    labs(title = plt_title,
         x = "Original Sample Distance", y = "MDS Map Distance") +
    theme(aspect.ratio = 1)

  plt
}
