#' Differential expression analysis report.
#'
#' Generate plots and CSV files that report on differential expression analysis
#' over multiple samples, sample features, and cell subsets.
#'
#' @param experiment An Astrolabe experiment.
#' @import ggplot2
#' @export
reportDifferentialExpressionAnalysis <- function(experiment, verbose = FALSE) {
  # Value for NA feature value. Samples with this value should not be included
  # in the analysis.
  FEATURE_NA <- "__NA__"
  # Minimum value for legends.
  LEGEND_MIN <- 1
  # Minimum max fold change to include ridge plot of marker.
  RIDGE_MIN_MAXFC <- 0.5
  
  differential_expression_analysis <-
    experimentDifferentialExpressionAnalysis(experiment)
  # Load channel densities across all samples.
    channel_dens <-
      lapply(nameVector(experiment$samples$SampleId), function(sample_id) {
        readRDS(file.path(
          experiment$analysis_path,
          paste0(sample_id,
                 ".aggregate_statistics.RDS")))$subset_channel_densities
      })
  
  # Separate directory for each cell subset level.
  lapply(nameVector(names(differential_expression_analysis)), function(level) {
    if (verbose) message(level)
    level_dea <- differential_expression_analysis[[level]]
    report <- list()
    
    # Report on each differential analysis kit.
    for (kit_name in names(level_dea)) {
      if (verbose) message(paste0("\t", kit_name))
      kit <- subset(experiment$de_analysis_kits, Name == kit_name)
      dea <- differential_expression_analysis[[level]][[kit_name]]$dea
      report[[kit_name]] <- list()
      
      # Get sample features for this kit and decide on sample order (sort by
      # primary and make sure that baseline is first).
      cell_subset_order <- rev(gtools::mixedsort(unique(dea$CellSubset)))
      sample_features <- experiment$sample_features
      sample_features$Primary <- sample_features[[kit$PrimaryFeatureId]]
      sample_features <- subset(sample_features, Primary != FEATURE_NA)
      sample_features$Primary <- factor(sample_features$Primary)
      sample_features$Primary <- 
        relevel(sample_features$Primary, ref = kit$PrimaryFeatureBaselineValue)
      sample_order <- sample_features$SampleId[order(sample_features$Primary)]
      sample_names <- 
        experiment$samples$Name[match(sample_order,
                                      experiment$samples$SampleId)]
      
      # Figure: Heat map of maximum fold changes.
      y_max <- ceiling(max(abs(dea$MaxFc), na.rm = TRUE))
      y_max <- max(y_max, LEGEND_MIN)
      fc_heat_map <- 
        plotHeatmap(dea,
                    x = "ChannelName",
                    y = "CellSubset",
                    value = "MaxFc",
                    type = "change",
                    y_axis_order = cell_subset_order,
                    fill_limits = c(-y_max, y_max))
      fc_heat_map$plt <- 
        fc_heat_map$plt + 
        labs(
          x = "Marker",
          y = level,
          title = paste0(kit_name, ": Fold Change")
        )
      # Transpose data to match heat map (rows are cell subsets).
      channel_names <- fc_heat_map$data$ChannelName
      fc_heat_map$data <- t(fc_heat_map$data[, -1])
      colnames(fc_heat_map$data) <- channel_names
      fc_heat_map$data <- fc_heat_map$data[rev(cell_subset_order), ]
      report[[paste0(kit_name, "_fold_change")]] <- fc_heat_map
      
      # Figures: Heat map of marker intensity across samples for each marker.
      if (verbose) message("\theat map of market intensity")
      dea_long <-
        reshape2::melt(dea,
                       id.vars = c("CellSubset", "ChannelName"),
                       measure.vars = sample_order,
                       variable.name = "SampleId",
                       value.name = "Mean")
      if (is.na(kit$ReferenceFeatureId)) {
        marker_legend_label <- "Marker Mean"
      } else {
        reference_feature_id <- 
          gsub("feature_", "", kit$ReferenceFeatureId, fixed = TRUE)
        reference_feature_name <-
          subset(experiment$features,
                 FeatureId == reference_feature_id)$FeatureName
        marker_legend_label <-
          paste0("Mean Intensity - (", reference_feature_name, ", ", 
                 kit$ReferenceFeatureBaselineValue, ")")
      }
      for (channel_name in kit$ChannelNames[[1]]) {
        channel_dea_long <- subset(dea_long, ChannelName == channel_name)
        channel_dea_long[[marker_legend_label]] <- channel_dea_long$Mean

        # Make the color coding friendlier for humans to look at by rounding and
        # setting reasonable min/max values.
        fill_min <-
          floor(min(channel_dea_long$Mean, na.rm = TRUE))
        # If there's no reference feature, minimum is 0 or close to 0, we don't
        # want to push it lower than that.
        if (!is.na(kit$ReferenceFeatureId)) {
          fill_min <- min(fill_min, -LEGEND_MIN)
        }
        fill_max <-
          ceiling(max(channel_dea_long$Mean, na.rm = TRUE))
        fill_max <- max(fill_max, LEGEND_MIN)

        channel_heat_map <- 
          plotHeatmap(channel_dea_long,
                      x = "CellSubset",
                      y = "SampleId",
                      value = marker_legend_label,
                      type = "change",
                      fill_limits = c(fill_min, fill_max))
        channel_heat_map$plt <- 
          channel_heat_map$plt +
          scale_y_discrete(name = "Sample",
                           limits = sample_order,
                           labels = sample_names) +
          guides(fill = guide_colorbar(title.position = "top"))
        # Transpose data to match heat map (rows are samples) and change sample
        # IDs to sample names.
        cell_subsets <- channel_heat_map$data$CellSubset
        channel_heat_map$data <- t(channel_heat_map$data[, -1])
        channel_heat_map$data <- channel_heat_map$data[rev(sample_order), ]
        colnames(channel_heat_map$data) <- cell_subsets
        rownames(channel_heat_map$data) <- rev(sample_names)

        report[[kit_name]][[paste0(channel_name, "_intensity")]] <- channel_heat_map
      }

      # Figures: Ridge plots of marker intensity distributions across samples
      # for each (Marker, Cell Subset) combination with MaxFC g.t.e
      # RIDGE_MIN_MAXFC.
      if (verbose) message("\tridge plots")
      reference_means <-
        differential_expression_analysis[[level]][[kit_name]]$reference_means
      ridge_dea <- subset(dea, abs(MaxFc) >= RIDGE_MIN_MAXFC)
      if (nrow(ridge_dea) > 0) {
        for (row_idx in seq(nrow(ridge_dea))) {
          channel_name <- ridge_dea$ChannelName[row_idx]
          if (is.null(report[[kit_name]][[channel_name]])) {
            report[[kit_name]][[channel_name]] <- list()
          }
          cell_subset <- ridge_dea$CellSubset[row_idx]
          if (verbose) message(paste0("\t\t", channel_name, ", ", cell_subset))

          # Find the X-axis values for each sample after subtracting
          # reference_means.
          x_vals <- 
            lapply(nameVector(sample_order), function(sample_id) {
              x <- channel_dens[[sample_id]][[level]][[cell_subset]][[channel_name]]$x
              if (is.null(reference_means)) {
                x
              } else {
                sample_rm <-
                  subset(reference_means,
                         SampleId == sample_id & CellSubset == cell_subset &
                           ChannelName == channel_name)$ReferenceMean
                x - sample_rm
              }
            })
          # Calculate global x-axis limits.
          x_min <- min(unlist(x_vals))
          x_max <- max(unlist(x_vals))
          
          # Generate plot for each sample.
          plts <- lapply(rev(sample_order), function(sample_id) {
            sample_name <- 
              subset(experiment$samples, SampleId == sample_id)$Name
            dens <-
              channel_dens[[sample_id]][[level]][[cell_subset]][[channel_name]]
            if (is.null(dens)) return(NULL)
            dens_df <- data.frame(X = x_vals[[sample_id]], Y = dens$y)
            ggplot(dens_df, aes(x = X, y = Y)) +
              geom_line() +
              xlim(c(x_min, x_max)) +
              labs(x = sample_name) +
              theme_linedraw() +
              theme(axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.x = element_text(size = 4),
                    axis.title.y = element_blank())
          })
          # Remove any missing samples from plts.
          plts <- plts[unlist(lapply(plts, function(l) !is.null(l)))]
          # Set up orloj plot and update report.
          plt <-
            Reduce(`+`, plts) +
            patchwork::plot_layout(ncol = 1) +
            patchwork::plot_annotation(
              subtitle = paste0(cell_subset, " (", channel_name, ")"))
          report[[kit_name]][[channel_name]][[cell_subset]] <- 
            list(plt = plt,
                 width = 300,
                 height = 50 + 50 * length(sample_order))
        }
      }
    }
    
    report
  })
}
