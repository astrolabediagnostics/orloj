#' Differential abundance analysis report.
#'
#' Generate plots and CSV files that report on differential abundance analysis
#' over multiple samples, sample features, and cell subsets. For each sample
#' feature, the report includes spreadsheet and volcano plot of analysis
#' results.
#'
#' @param experiment An Astrolabe experiment.
#' @import ggplot2 ggrepel
#' @export
reportDifferentialAbundanceAnalysis <- function(experiment) {
  # Value for NA feature value. Samples with this value should not be included
  # in the analysis.
  feature_na <- "__NA__"

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
  # Get group feature label + name.
  group_feature_label <- differential_abundance_analysis$group_feature_label
  if (!is.null(group_feature_label)) {
    group_feature_name <-
      dplyr::filter(experiment$features,
                    FeatureId == gsub("feature_", "",
                                      group_feature_label))$FeatureName
  }

  analyses <-
    names(differential_abundance_analysis$differential_abundance_analysis)
  lapply(nameVector(analyses), function(analysis) {
    cell_counts <-
      subset(aggregate_statistics$subset_cell_counts, Parent == analysis)
    group_feature_label <- differential_abundance_analysis$group_feature_label
    daa <-
      differential_abundance_analysis$differential_abundance_analysis[[analysis]]

    # Get DAA in export-ready format.
    export_daa <- differentialAbundanceAnalysis(experiment, level = analysis)

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

      # Get analysis results for this feature.
      feature_top_tags <- daa[[feature_r_name]]$table
      if (is.null(feature_top_tags)) {
        return(feature_top_tags);
      }
      feature_top_tags <-
        tibble::rownames_to_column(feature_top_tags, "CellSubset")

      # Figure: Volcano plot
      feature_top_tags$negLog10Fdr = -log10(feature_top_tags$FDR)
      if ("logFC" %in% colnames(feature_top_tags)) {
        volcano_plot_x <- "logFC"
      } else if ("maxLogFc" %in% colnames(feature_top_tags)) {
        volcano_plot_x <- "maxLogFc"
      } else {
        stop("feature_top_tags does not include logFC or maxLogFc")
      }
      # Calculate symmetric X limit for volcano plot.
      volcano_x <- feature_top_tags[[volcano_plot_x]]
      volcano_xlim <- volcano_x[which.max(abs(volcano_x))]
      volcano_xlim <- ceiling(volcano_xlim / 0.25) * 0.25
      volcano_plt_list <-
        plotScatterPlot(feature_top_tags,
                        x = volcano_plot_x,
                        y = "negLog10Fdr",
                        xlim = c(-volcano_xlim, volcano_xlim),
                        title = paste0(feature_name, ": Volcano Plot"))
      volcano_plt_list$plt <- volcano_plt_list$plt +
        ggrepel::geom_text_repel(aes(label = CellSubset))

      list(
        # Spreadsheet: Statistical analysis results
        analysis_results = list(data = export_daa[[feature_name]]),
        volcano = volcano_plt_list
      )
    })
  })
}
