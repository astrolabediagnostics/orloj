#' Differential expression analysis report.
#'
#' Generate plots and CSV files that report on differential expression analysis
#' over multiple samples, sample features, and cell subsets.
#'
#' @param experiment An Astrolabe experiment.
#' @import ggplot2
#' @export
reportDifferentialExpressionAnalysis <- function(experiment) {
  # Value for NA feature value. Samples with this value should not be included
  # in the analysis.
  feature_na <- "__NA__"
  
  differential_expression_analysis <-
    experimentDifferentialExpressionAnalysis(experiment)
  
  # Separate directory for each cell subset level.
  lapply(nameVector(names(differential_expression_analysis)), function(level) {
    level_dea <- differential_expression_analysis[[level]]
    report <- list()
    
    # Report on each differential analysis kit.
    for (kit_name in names(level_dea)) {
      kit <- subset(experiment$de_analysis_kits, Name == kit_name)
      dea <- differential_expression_analysis[[level]][[kit_name]]$dea
      report[[kit_name]] <- list()
      
      # Get sample features for this kit and decide on sample order (sort by
      # primary and make sure that baseline is first).
      cell_subset_order <- rev(gtools::mixedsort(unique(dea$CellSubset)))
      sample_features <- experiment$sample_features
      sample_features$Primary <- sample_features[[kit$PrimaryFeatureId]]
      sample_features <- subset(sample_features, Primary != feature_na)
      sample_features$Primary <- factor(sample_features$Primary)
      sample_features$Primary <- 
        relevel(sample_features$Primary, ref = kit$PrimaryFeatureBaselineValue)
      sample_order <- sample_features$SampleId[order(sample_features$Primary)]
      sample_names <- 
        experiment$samples$Name[match(sample_order,
                                      experiment$samples$SampleId)]
      
      # Figure: Heat map of maximum fold changes.
      y_max <- ceiling(max(abs(dea$MaxFc), na.rm = TRUE) / 0.25) * 0.25
      fc_heat_map <- 
        plotHeatmap(dea,
                    x = "ChannelName",
                    y = "CellSubset",
                    value = "MaxFc",
                    y_axis_order = cell_subset_order)
      fc_heat_map$plt <- 
        fc_heat_map$plt + 
        scale_fill_gradientn(name = "max(log(fold change))",
                             colors = c("blue", "white", "red"),
                             limits = c(-y_max, y_max)) +
        labs(
          x = "Marker",
          y = level,
          title = paste0(kit_name, ": Fold Change")
        )
      report[[paste0(kit_name, "_fold_change")]] <- fc_heat_map
      
      # Figures: Heat map of marker intensity across samples for each marker.
      dea_long <-
        reshape2::melt(dea,
                       id.vars = c("CellSubset", "ChannelName"),
                       measure.vars = sample_order,
                       variable.name = "SampleId",
                       value.name = "Mean")
      for (channel_name in kit$ChannelNames[[1]]) {
        channel_dea_long <- subset(dea_long, ChannelName == channel_name)
        channel_heat_map <- 
          plotHeatmap(channel_dea_long,
                      x = "CellSubset",
                      y = "SampleId",
                      value = "Mean")
        channel_min <-
          floor(min(channel_dea_long$Mean, na.rm = TRUE) / 0.25) * 0.25
        channel_max <-
          ceiling(max(channel_dea_long$Mean, na.rm = TRUE) / 0.25) * 0.25
        channel_heat_map$plt <- 
          channel_heat_map$plt +
          scale_fill_gradient(low = "white", high = "#003300",
                              limits = c(channel_min, channel_max)) +
          scale_y_discrete(name = "Sample",
                           limits = sample_order,
                           labels = sample_names)
        report[[kit_name]][[paste0(channel_name, "_intensity")]] <- channel_heat_map
      }
    }
    
    report
  })
}
