#' Differential abundance analysis report.
#'
#' Generate plots and CSV files that report on differential abundance analysis
#' over multiple samples, sample features, and cell subsets. For each sample
#' feature, the report includes box plots, bar plots, line plots, and volcano
#' plot of analysis results.
#'
#' @param experiment An Astrolabe experiment.
#' @import ggplot2
#' @export
reportDifferentialAbundanceAnalysis <- function(experiment) {
  samples <- experiment$samples
  features <- experiment$features
  sample_features <- experiment$sample_features

  # Load aggregate statistics and differential abundance analysis.
  aggregate_statistics_filename <-
    file.path(experiment$analysis_path, "combine_aggregate_statistics.RDS")
  if (!file.exists(aggregate_statistics_filename)) {
    stop(paste0("unable to find ", aggregate_statistics_filename))
  }
  aggregate_statistics <- readRDS(aggregate_statistics_filename)

  daa_filename <-
    file.path(experiment$analysis_path, "differential_abundance_analysis.RDS")
  if (!file.exists(daa_filename)) {
    stop(paste0("unable to find ", daa_filename))
  }
  differential_abundance_analysis <- readRDS(daa_filename)

  analyses <- nameVector(c("Assignment", "Profiling"))
  lapply(analyses, function(analysis) {
    cell_counts <- getCellCounts(aggregate_statistics, analysis)
    group_feature_label <- differential_abundance_analysis$group_feature_label
    daa <-
      differential_abundance_analysis$differential_abundance_analysis[[analysis]]

    # Join data with sample name and sample features.
    figure_data <- cell_counts %>%
      dplyr::left_join(sample_features, by = "SampleId") %>%
      dplyr::left_join(samples, by = "SampleId") %>%
      dplyr::rename(SampleName = Name)
    # Calculate frequencies.
    figure_data <- figure_data %>%
      dplyr::group_by(SampleId) %>%
      dplyr::mutate(Frequency = N / sum(N)) %>%
      dplyr::ungroup()

    # Generate figures for each feature separately.
    feature_names <- features$FeatureName
    lapply(nameVector(feature_names), function(feature_name) {
      feature_r_name <-
        paste0(
          "feature_",
          dplyr::filter(features, FeatureName == feature_name)$FeatureId
        )
      feature_report <- list()

      # Get analysis results for this feature.
      feature_top_tags <- daa[[feature_r_name]]$table
      feature_top_tags <-
        tibble::rownames_to_column(feature_top_tags, "CellSubset")
      cell_subset_labels <- feature_top_tags$CellSubset

      # Whether to include line plots (if sample features include patient as a
      # grouping feature).
      include_line_plots <- FALSE
      if (!is.null(group_feature_label)) {
        include_line_plots <- feature_r_name != group_feature_label
      }

      # Table: Result of differential analysis.
      feature_report[["analysis_results"]] <- list(data = feature_top_tags)

      # Figure: Volcano plot (if feature only has two values).
      if ("logFC" %in% colnames(feature_top_tags)) {
        feature_top_tags$negLog10Fdr = -log10(feature_top_tags$FDR)
        volcano_plt_list <-
          plotScatterPlot(feature_top_tags,
                          x = "logFC",
                          y = "negLog10Fdr")
        volcano_plt_list$plt <- volcano_plt_list$plt +
          ggrepel::geom_text_repel(aes(label = CellSubset))

        feature_report[["volcano"]] <- volcano_plt_list
      }

      # Figures: Line plots (for non-patient features in experiments that have
      # patient), box plots, and bar plots.
      feature_report[["box_plots"]] <- list()
      feature_report[["bar_plots"]] <- list()
      if (include_line_plots) feature_report[["line_plots"]] <- list()

      for (cell_subset_label in cell_subset_labels) {
        # Filter data down to this subset.
        cell_subset_data <-
          dplyr::filter(figure_data, CellSubset == cell_subset_label)

        # Format FDR nicely for figure title.
        fdr <-
          dplyr::filter(feature_top_tags, CellSubset == cell_subset_label)$FDR
        fdr_str <- paste0("(", formatPvalue(fdr, "FDR"), ")")
        fig_title <- paste0(cell_subset_label, " ", fdr_str)

        # Generate box plot.
        box_plot <-
          plotBoxPlot(cell_subset_data,
                      x = feature_r_name,
                      y = "Frequency",
                      title = fig_title,
                      scale_y_labels = scales::percent)
        box_plot$plt <- box_plot$plt + xlab(feature_name)
        feature_report[["box_plots"]][[cell_subset_label]] <- box_plot

        # Generate bar plot.
        sample_order <- dplyr::arrange(cell_subset_data, Frequency)$SampleName
        cell_subset_data$SampleName <-
          factor(cell_subset_data$SampleName, levels = sample_order)
        bar_plot <-
          plotBarPlot(cell_subset_data,
                      x = "SampleName",
                      y = "Frequency",
                      fill = feature_r_name,
                      title = fig_title,
                      scale_y_labels = scales::percent)
        feature_report[["bar_plots"]][[cell_subset_label]] <- bar_plot

        # Generate line plot, using the box_plot parameters as base.
        if (include_line_plots) {
          line_plot <- box_plot
          line_plot$plt <-
            ggplot(cell_subset_data,
                   aes_string(x = feature_r_name, y = "Frequency")) +
            geom_line(aes_string(group = group_feature_label)) +
            scale_y_continuous(labels = scales::percent) +
            labs(title = fig_title, x = feature_name)
          feature_report[["line_plots"]][[cell_subset_label]] <- line_plot
        }
      }

      feature_report
    })
  })
}
