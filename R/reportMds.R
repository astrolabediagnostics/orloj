#' MDS Map report.
#'
#' Generate plots and CSV files that report on the MDS map.
#'
#' @param experiment An Astrolabe experiment.
#' @import ggplot2 viridis ggrepel patchwork
#' @export
reportMds <- function(experiment) {
  # Figure width and height. These values were chosen to match the aspect ratio
  # in the app.
  width <- 650
  height <- 500
  # Value for NA feature value. Samples with this value should not be included
  # in the MDS map for a given feature.
  feature_na <- "__NA__"
  # Consistent theme across all figures.
  figure_theme <-
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      panel.background = element_blank()
    )

  cell_subsets <-
    readRDS(file.path(experiment$analysis_path, "experiment_cell_subsets.RDS"))

  # Iterate over all cell subset levels and generate MDS maps for each.
  lapply(nameVector(colnames(cell_subsets)), function(level) {
    mds <- experimentMds(experiment, level = level)

    # Get cell subset counts for this level and join with sample features.
    cell_counts <- experimentCellSubsetCounts(experiment, level = level)
    cell_counts <- 
      dplyr::left_join(
        cell_counts,
        experiment$sample_features,
        by = "SampleId"
      )
    # Calculate maximum size of each dot (in figures where size is frequency).
    size_max <- ceiling(max(cell_counts$Freq) * 10) / 10

    report <- list()

    # Figure: Plain map, no color-coding.
    map_obj <-
      ggplot(mds, aes(x = V1, y = V2)) +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(aes(label = CellSubset)) +
      labs(title = paste0(level, ": MDS map"), x = "MDS 1", y = "MDS 2") +
      figure_theme
    map_plt <-
      list(plt = map_obj, width = width, height = height, data = mds)
    report$mds_map <- map_plt

    # Figures: Color-coding by each channel.
    for (channel in experiment$analysis_channels) {
      channel_filename <- filenameify(channel)
      channel_bt <- addBackticks(channel)
      channel_limits <-
        c(min(0, mds[[channel]]), ceiling(max(mds[[channel]]) * 10) / 10)

      map_obj <-
        ggplot(mds, aes(x = V1, y = V2)) +
        geom_point(aes_string(color = channel_bt), size = 3) +
        viridis::scale_color_viridis(direction = -1, limits = channel_limits) +
        labs(title = paste0(level, ": MDS map"), x = "MDS 1", y = "MDS 2") +
        figure_theme
      map_plt <- list(plt = map_obj, width = width, height = height)
      report[[paste0("channel_", channel_filename)]] <- map_plt
    }

    # Figures: MDS map for each sample, resize dot by subset frequency.
    report$by_sample <- list()
    for (sample_id in experiment$samples$SampleId) {
      sample_name <- subset(experiment$samples, SampleId == sample_id)$Name
      sample_mds <- mds[, c("CellSubset", "V1", "V2")]
      sample_mds <-
        dplyr::left_join(sample_mds,
                         subset(cell_counts, SampleId == sample_id),
                         by = "CellSubset")
      map_obj <-
        ggplot(sample_mds, aes(x = V1, y = V2)) +
        geom_point(aes(size = Freq), alpha = 0.75) +
        scale_size_continuous(
          name = "Frequency",
          labels = scales::percent,
          limits = c(0, size_max),
          range = c(1, 15)) +
        labs(title = paste0(level, ": MDS map"),
             subtitle = sample_name,
             x = "MDS 1", y = "MDS 2") +
        theme_linedraw() +
        figure_theme
      report$by_sample[[sample_name]] <- 
        list(plt = map_obj, width = width, height = height)
    }

    # Figures: MDS map for each feature, resize by median frequency.
    report$by_feature <- list()
    for (feature_id in experiment$features$FeatureId) {
      feature_id_long <- paste0("feature_", feature_id)
      feature_name <-
        subset(experiment$features, FeatureId == feature_id)$FeatureName
      
      # Get median frequency by subset for each value of this feature, excluding
      # samples whose value is NA, and merge with MDS map.
      feature_cell_counts <- cell_counts
      feature_cell_counts$Feature <- feature_cell_counts[[feature_id_long]]
      feature_cell_counts <- 
        feature_cell_counts %>%
        dplyr::group_by(CellSubset, Feature) %>%
        dplyr::summarize(Freq = mean(Freq))
      feature_cell_counts <- subset(feature_cell_counts, Feature != feature_na)
      feature_mds <- mds[, c("CellSubset", "V1", "V2")]
      feature_mds <- 
        dplyr::left_join(feature_mds, feature_cell_counts, by = "CellSubset")
      
      # Generate figure, make sure to facet on feature values.
      map_obj <-
        ggplot(feature_mds, aes(x = V1, y = V2)) +
        geom_point(aes(size = Freq), alpha = 0.75) +
        scale_size_continuous(
          name = "Frequency",
          labels = scales::percent,
          limits = c(0, size_max),
          range = c(1, 15)) +
        labs(title = paste0(level, ": MDS map"),
             subtitle = feature_name,
             x = "MDS 1", y = "MDS 2") +
        facet_wrap(~ Feature, ncol = 4) +
        theme_linedraw() +
        figure_theme
      
      # Calculate figure height and width based on the number of feature values.
      if (is.character(feature_mds$Feature)) {
        feature_mds$Feature <- factor(feature_mds$Feature)
      }
      n_features <- length(levels(feature_mds$Feature))
      feature_width <- min(n_features, 4) * width
      feature_height <- ceiling(n_features / 4) * height
      report$by_feature[[feature_name]] <- 
        list(plt = map_obj, width = feature_width, height = feature_height)
    }

    report$shepard_plot <- .plotShepard(experiment, level, mds)

    report
  })
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

  plt <- orloj::plotScatterPlot(dist_df, x = "Dist", y = "MDSDist", alpha = 0.1)
  plt$plt <-
    plt$plt +
    labs(title = paste0(level, ": MDS Map Shepard Plot"),
         subtitle = paste0("Pearson's rho = ", pearson_rho),
         x = "High-Dimensional Distance", y = "MDS Map Distance") +
    theme(aspect.ratio = 1)

  plt
}
