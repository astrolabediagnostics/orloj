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
  }

  report
}
