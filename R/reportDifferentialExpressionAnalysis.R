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
  # Minimum value for the legend range.
  LEGEND_MIN <- 1
  # Minimum max fold change to include ridge plot of marker.
  RIDGE_MIN_MAXFC <- 0.5
  # Minimum number of ridge plots to include.
  RIDGE_MIN_N <- 20
  # Maximum number of samples for ridge plots and dot plots.
  SAMPLE_COUNT_MAX_N <- 50
  
  differential_expression_analysis <-
    experimentDifferentialExpressionAnalysis(experiment)
  channel_dens <- .loadChannelDensities(experiment)
  
  # Generate the report for each cell subset level.
  lapply(nameVector(names(differential_expression_analysis)), function(level) {
    if (verbose) message(level)
    level_dea <- differential_expression_analysis[[level]]
    report <- list()
    
    # Report on each differential analysis kit.
    for (kit_name in names(level_dea)) {
      if (verbose) message(paste0("\t", kit_name))
      kit <- subset(experiment$de_analysis_kits, Name == kit_name)
      primary_feature_name <-
        subset(experiment$features,
               FeatureId == gsub("feature_", "",
                                 kit$PrimaryFeatureId))$FeatureName
      dea <- differential_expression_analysis[[level]][[kit_name]]$dea
      
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
      
      # Spreadsheet: DE analysis results, formatted nicely.
      report[[kit_name]] <-
        .generateDeResultsSpreadsheet(dea, sample_order, sample_names)
      
      # Figure: Heat map of maximum fold changes.
      report[[paste0(kit_name, "_fold_change")]] <-
        .generateFcHeatMap(kit_name, level, dea, cell_subset_order,
                           legend_min = LEGEND_MIN)
      
      marker_legend_label <- .getMarkerLegendLabel(experiment, kit)
      dea_long <-
        reshape2::melt(dea,
                       id.vars = c("CellSubset", "ChannelName"),
                       measure.vars = sample_order,
                       variable.name = "SampleId",
                       value.name = "Mean")
      
      # Figures: Heat map of marker intensity across samples for each marker.
      if (verbose) message("\theat map of marker intensity")
      report <- 
        .addMarkerIntensityHeatMap(
          report, dea_long, kit, kit_name, sample_order, sample_names,
          marker_legend_label, legend_min = LEGEND_MIN)
      
      # Figures: Dot and box plot of marker intensity across samples.
      report <-
        .addMarkerIntensityDotBoxPlots(
          report, kit, kit_name, dea_long, primary_feature_name,
          sample_features, cell_subset_order, sample_order, sample_names,
          marker_legend_label, sample_count_max_n = SAMPLE_COUNT_MAX_N)
      
      # Figures: Ridge plots of marker intensity distribution across samples.
      if (verbose) message("\tridge plots")
      reference_means <- 
        differential_expression_analysis[[level]][[kit_name]]$reference_means
      report <-
        .addMarkerIntensityRidgePlots(
          report, experiment, level, channel_dens, kit_name, dea,
          reference_means, sample_order,
          sample_count_max_n = SAMPLE_COUNT_MAX_N,
          ridge_min_maxfc = RIDGE_MIN_MAXFC, ridge_min_n = RIDGE_MIN_N)
    }
    
    report
  })
}

.loadChannelDensities <- function(experiment) {
  # Load channel densities across all samples.
  lapply(nameVector(experiment$samples$SampleId), function(sample_id) {
    readRDS(file.path(
      experiment$analysis_path,
      paste0(sample_id,
             ".aggregate_statistics.RDS")))$subset_channel_densities
  })
}

.generateDeResultsSpreadsheet <- function(dea, sample_order, sample_names) {
  # Generate the DE results spreadsheet, which is a nicer-looking version of the
  # dea data frame.
  first_cols <-
    c("CellSubset", "ChannelName", "MaxFc", "P.Value", "adj.P.Val")
  remove_cols <- c("AveExpr", "t", "B")
  contrast_cols <-
    setdiff(colnames(dea), c(first_cols, remove_cols, sample_order))
  dea_samples <- dea[, sample_order]
  colnames(dea_samples) <- sample_names
  dea_report <- cbind(dea[, c(first_cols, contrast_cols)], dea_samples)
  rownames(dea_report) <- NULL
  dea_report <- dea_report[order(dea_report$P.Value), ]
  list(data = dea_report)
}

.generateFcHeatMap <- function(kit_name, level, dea, cell_subset_order,
                               legend_min = 1) {
  # Generate a heat map and accommodating spreadsheet of maximum fold change
  # values for this differential expression analysis.
  y_max <- ceiling(max(abs(dea$MaxFc), na.rm = TRUE))
  y_max <- max(y_max, legend_min)
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
  fc_heat_map
}

.getMarkerLegendLabel <- function(experiment, kit) {
  # Build the marker legend label that will be used in later figures. The label
  # varies if the kit includes a reference feature.
  if (is.na(kit$ReferenceFeatureId)) return("Mean")
  
  reference_feature_id <- 
    gsub("feature_", "", kit$ReferenceFeatureId, fixed = TRUE)
  reference_feature_name <-
    subset(experiment$features,
           FeatureId == reference_feature_id)$FeatureName
  paste0("Mean - (", reference_feature_name, ", ", 
         kit$ReferenceFeatureBaselineValue, ")")
}

.shortenSampleNames <- function(v, max_nchar = 10) {
  # Shorten sample name vector, if it's longer than given character count.
  unlist(lapply(v, function(s) {
    if (nchar(s) <= max_nchar) {
      return(s)
    }
    paste0(substr(s, 1, max_nchar), "...")
  }))
}

.addMarkerIntensityHeatMap <- function(report, dea_long, kit, kit_name,
                                       sample_order, sample_names,
                                       marker_legend_label, legend_min = 1) {
  # Add the marker intensity heat maps to the report.
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
      fill_min <- min(fill_min, -legend_min)
    }
    fill_max <-
      ceiling(max(channel_dea_long$Mean, na.rm = TRUE))
    fill_max <- max(fill_max, legend_min)
    
    channel_heat_map <- 
      plotHeatmap(channel_dea_long,
                  x = "CellSubset",
                  y = "SampleId",
                  value = marker_legend_label,
                  type = "change",
                  title = channel_name,
                  fill_limits = c(fill_min, fill_max))
    channel_heat_map$plt <- 
      channel_heat_map$plt +
      scale_y_discrete(name = "Sample",
                       limits = sample_order,
                       labels = .shortenSampleNames(sample_names)) +
      guides(fill = guide_colorbar(title.position = "top"))
    # Transpose data to match heat map (rows are samples) and change sample
    # IDs to sample names.
    cell_subsets <- channel_heat_map$data$CellSubset
    channel_heat_map$data <- t(channel_heat_map$data[, -1])
    channel_heat_map$data <- channel_heat_map$data[rev(sample_order), ]
    colnames(channel_heat_map$data) <- cell_subsets
    rownames(channel_heat_map$data) <- rev(sample_names)
    
    report[[kit_name]][[paste0("heat_map.", channel_name)]] <- channel_heat_map
  }
  
  report
}

.addMarkerIntensityDotBoxPlots <- function(report,
                                           kit, kit_name,
                                           dea_long, primary_feature_name,
                                           sample_features, cell_subset_order,
                                           sample_order, sample_names,
                                           marker_legend_label,
                                           sample_count_max_n = 50) {
  # Add the marker intensity dot plots and box plots to the report.
  
  # Calculate figure size.
  width <- 900
  height <- ceiling(length(cell_subset_order) / 4) * 2.5 * 100
  
  for (channel_name in kit$ChannelNames[[1]]) {
    # Get DEA data for this channel and join with primary feature.
    channel_dea_long <- subset(dea_long, ChannelName == channel_name)
    channel_dea_long <-
      dplyr::left_join(channel_dea_long, sample_features, by = "SampleId")
    channel_dea_long[[kit$PrimaryFeatureId]] <- 
      factor(channel_dea_long[[kit$PrimaryFeatureId]])
    channel_dea_long[[kit$PrimaryFeatureId]] <- 
      relevel(channel_dea_long[[kit$PrimaryFeatureId]],
              ref = kit$PrimaryFeatureBaselineValue)
    
    # Dot plot. Only generate if sample count is l.t.e sample_count_max_n.
    if (nrow(sample_features) <= sample_count_max_n) {
      channel_dot_plot <-
        ggplot(channel_dea_long, aes(x = SampleId, y = Mean)) +
        geom_point(aes_string(color = kit$PrimaryFeatureId)) +
        scale_x_discrete(limits = sample_order,
                         labels = .shortenSampleNames(sample_names)) +
        scale_color_brewer(name = primary_feature_name, palette = "Dark2") +
        facet_wrap(~ CellSubset, scales = "free", ncol = 4,
                   labeller = label_wrap_gen()) +
        labs(y = paste0(channel_name, " ", marker_legend_label)) +
        theme_linedraw() +
        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0),
              axis.title.x = element_blank(),
              legend.position = "top")
      report[[kit_name]][[paste0("dot_plot.", channel_name)]] <-
        list(plt = channel_dot_plot, width = width, height = height)
    }
    
    # Box plot.
    set.seed(12345)
    channel_box_plot <- 
      ggplot(channel_dea_long,
             aes_string(x = kit$PrimaryFeatureId, y = "Mean")) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes_string(color = kit$PrimaryFeatureId),
                  width = 0.15, alpha = 0.4) +
      scale_color_brewer(name = primary_feature_name, palette = "Dark2") +
      facet_wrap(~ CellSubset, scales = "free", ncol = 4,
                 labeller = label_wrap_gen()) +
      labs(y = paste0(channel_name, " ", marker_legend_label)) +
      theme_linedraw() +
      theme(axis.title.x = element_blank(),
            legend.position = "top")
    report[[kit_name]][[paste0("box_plot.", channel_name)]] <-
      list(plt = channel_box_plot, width = width, height = height)
  }
  
  report
}

.addMarkerIntensityRidgePlots <- function(report,
                                          experiment, level, channel_dens,
                                          kit_name,
                                          dea, reference_means,
                                          sample_order,
                                          sample_count_max_n = 50,
                                          ridge_min_maxfc = 0.5,
                                          ridge_min_n = 20) {
  # Add the marker intensity ridge plots to the report. Ridge plots are
  # generated across samples for each (Marker, Cell Subset) combination with
  # MaxFC g.t.e ridge_min_maxfc. However, a minimum of ridge_min_n figures are
  # generated.
  
  # Don't generate ridge plots if we have a high sample count.
  if (length(sample_order) > sample_count_max_n) return(report)
  
  # Find combinations with MaxFC g.t.e threshold, and make sure we have the
  # required minimum.
  ridge_dea <- subset(dea, abs(MaxFc) >= ridge_min_maxfc)
  if (nrow(ridge_dea) < ridge_min_n) {
    # Add ridge plots to the count of ridge_min_n.
    ridge_dea <-
      rbind(
        ridge_dea,
        dea %>%
          dplyr::filter(!grepl("unassigned", CellSubset)) %>%
          dplyr::arrange(P.Value) %>%
          dplyr::slice(seq(ridge_min_n - nrow(ridge_dea)))
      )
  }
  
  if (nrow(ridge_dea) == 0) return(report)
  
  for (row_idx in seq(nrow(ridge_dea))) {
    channel_name <- ridge_dea$ChannelName[row_idx]
    if (is.null(report[[kit_name]][[channel_name]])) {
      report[[kit_name]][[channel_name]] <- list()
    }
    cell_subset <- ridge_dea$CellSubset[row_idx]
    
    # Find the X-axis values for each sample after subtracting reference_means.
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
  
  report
}
