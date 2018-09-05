#' Immport export.
#'
#' Export an Immport CyTOF_derived_data.tsv file for a single sample.
#'
#' @param file_path File path for the exported tsv file.
#' @param experiment An Astrolabe experiment.
#' @param sample_id An Astrolabe sample ID.
#' @export
exportImmport <- function(file_path, experiment, sample_id) {
  if (!(sample_id %in% experiment$samples$SampleId)) {
    stop("cannot find sample_id in experiment")
  }
  
  # Standard column names.
  column_name_col          <- "Column Name"
  sample_id_col            <- "Expsample ID"
  pop_name_reported_col    <- "Population Name Reported"
  gating_def_reported_col  <- "Gating Definition Reported"
  parent_pop_reported_col  <- "Parent Population Reported"
  value_col                <- "Population Statistic (count, percentile, etc)"
  unit_col                 <- "Population Stat Unit Reported" 
  workspace_file_col       <- "Workspace File"
  comments_col             <- "Comments"

  # Load Astrolabe hierarchy.
  hierarchy <- readRDS(file.path(experiment$analysis_path, "hierarchy.RDS"))
  hierarchy <- hierarchy$hierarchy

  # Get sample name, filename, and expression data.
  sample_name <-
    experiment$samples$Name[experiment$samples$SampleId == sample_id]
  sample_filename <-
    experiment$samples$Filename[experiment$samples$SampleId == sample_id]
  sample <- loadSample(experiment, sample_id = sample_id)
  exprs <- fcsExprs(sample)

  if (nrow(exprs) == 0) return()

  # Calculate count and percent over all profiling subsets.
  stats <- .calculateCountPercent(exprs, value_col, unit_col)
  # Calculate mean channel intensity over all (profile, channel) pairs.
  channel_stats <-
    .calculateMeanChannelIntensity(experiment$analysis_channels,
                                   exprs,
                                   value_col,
                                   unit_col)
  stats <-
    rbind(stats, channel_stats[, colnames(stats)])
  # Convert profile name to Immport convention and add gating definition.
  desc_df <- .convertHierarchyToDesc(hierarchy, experiment$class_channels)
  ass_profile_map <-
    .convertProfileToImmport(exprs,
                             desc_df,
                             pop_name_reported_col,
                             gating_def_reported_col)
  # Combine and fill in missing columns.
  stats <- dplyr::left_join(stats, ass_profile_map, by = "Profile")
  stats[[column_name_col]] <- ""
  stats[[sample_id_col]] <- sample_name
  stats[[parent_pop_reported_col]] <- "live cells"
  stats[[workspace_file_col]] <- sample_filename
  stats[[comments_col]] <- ""

  # Export to TSV file.
  .exportImmportStats(file_path, stats)
}

.calculateCountPercent <- function(exprs, value_col, unit_col) {
  stats <- exprs %>%
    dplyr::group_by(Profile) %>%
    dplyr::summarize(Count = n()) %>%
    dplyr::mutate(Percent = Count / sum(Count) * 100)
  stats <- stats %>%
    reshape2::melt(id.vars = "Profile",
                   value.name = value_col,
                   variable.name = unit_col)
  stats[[unit_col]] <- as.character(stats[[unit_col]])

  stats
}

.calculateMeanChannelIntensity <- function(channels,
                                           exprs,
                                           value_col,
                                           unit_col) {
  stats <- exprs %>%
    reshape2::melt(id.vars = "Profile",
                   measure.vars = channels,
                   variable.name = "Channel",
                   value.name = "Intensity")
  stats <- stats %>%
    dplyr::group_by(Profile, Channel) %>%
    dplyr::summarize(MeanIntensity = mean(Intensity)) %>%
    dplyr::mutate(Unit =
                    paste0("Mean asinh(", as.character(Channel), " / 5)")) %>%
    dplyr::ungroup()
  stats[[value_col]] <- stats$MeanIntensity
  stats[[unit_col]] <- stats$Unit

  stats
}

.convertHierarchyToDesc <- function(hierarchy, channel_names) {
  channel_names <- intersect(channel_names, colnames(hierarchy))
  
  lapply(hierarchy$CellSubset, function(cell_subset) {
    start_cell_subset <- cell_subset
    desc <- ""
    
    # Iterate backward across the hierarchy and get the markers in current
    # subset and each of its parents, stopping at "Root".
    while (cell_subset != "Root") {
      while (grepl("_unassigned", cell_subset)) {
        cell_subset <- gsub("_unassigned", "", cell_subset)
      }
      
      row <- dplyr::filter(hierarchy, CellSubset == cell_subset)
      rule <- as.matrix(row[, channel_names])
      rule <- rule[, !is.na(rule), drop = FALSE]
      
      # For each level, description includes positive markers (if any exist)
      # followed by negative markers (if any exist).
      desc <- 
        paste0(
          ifelse(
            any(rule),
            paste0(
              paste(
                paste0(colnames(rule)[rule], "+"),
                collapse = " "),
              ifelse(any(!rule), " ", "")
            ),
            ""
          ),
          ifelse(
            any(!rule),
            paste(paste0(colnames(rule)[!rule], "-"), collapse = " "),
            ""
          ),
          ifelse(nchar(desc) > 0, " ", ""),
          desc
        )
      
      cell_subset <- row$Parent
    }
    
    data.frame(
      Assignment = start_cell_subset,
      Desc = desc,
      stringsAsFactors = FALSE
    )
  }) %>% dplyr::bind_rows()
}

.convertProfileToImmport <- function(exprs,
                                     desc_df,
                                     pop_name_reported_col,
                                     gating_def_reported_col) {
  # Match Astrolabe assignment to Immport notation. "unassigned" cannot be
  # matched and keep Astrolabe assignments.
  ass_profile_map <- unique(exprs[, c("Assignment", "Profile")])
  ass_profile_map <-
    dplyr::left_join(ass_profile_map, astrolabe_immport_map, by = "Assignment")
  unassigned_indices <- grep("unassigned", ass_profile_map$Assignment)
  ass_profile_map[[pop_name_reported_col]][unassigned_indices] <- 
    ass_profile_map$Assignment[unassigned_indices]

  # For each profile, get the markers that were used for that profile.
  ass_profile_map$ProfileGate <- ""
  for (idx in seq(nrow(ass_profile_map))) {
    if (ass_profile_map$Assignment[idx] != ass_profile_map$Profile[idx]) {
      ass_profile_map$ProfileGate[idx] <-
        substr(ass_profile_map$Profile[idx],
               nchar(ass_profile_map$Assignment[idx]) + 1,
               nchar(ass_profile_map$Profile[idx]))
    }
  }
  # Immport profile names are of the format C & X, where C is taken from the
  # ontology and X from the profiling level.
  ass_profile_map[[pop_name_reported_col]] <- 
    paste0(ass_profile_map[[pop_name_reported_col]], " &",
           ass_profile_map$ProfileGate)

  # Create separate assignment column, replacing "unassigned" with blanks, and
  # match to assignment descriptions. Combine with profile gating description
  # to get gating definition for subset.
  ass_profile_map$AssignmentNU <- ass_profile_map$Assignment
  while (any(grepl("unassigned", ass_profile_map$AssignmentNU))) {
    ass_profile_map$AssignmentNU <-
      gsub("_unassigned", "", ass_profile_map$AssignmentNU)
  }
  ass_profile_map <-
    dplyr::left_join(ass_profile_map, desc_df,
                     by = c("AssignmentNU" = "Assignment"))
  ass_profile_map[[gating_def_reported_col]] <- 
    paste0(ass_profile_map$Desc, ass_profile_map$ProfileGate)

  ass_profile_map
}

.immportExportHeader <- function() {
  c(
    "cytof_derived_data\tSchema Version 3.18",
    "Please do not delete or edit this column",
    "Validation Level\tStandard"
  )
}

.immportExportColOrder <- function() {
  immport_export_col_order <- 
    c(
      "Column Name", 
      "Expsample ID", 
      "Population Name Reported", 
      "Gating Definition Reported", 
      "Parent Population Reported", 
      "Population Statistic (count, percentile, etc)", 
      "Population Stat Unit Reported", 
      "Workspace File", 
      "Comments"
    )
}

.exportImmportStats <- function(file_path, stats_long) {
  # Export the Immport stats data frame as a tsv, along with required
  # formatting.

  con <- file(file_path)
  writeLines(.immportExportHeader(), con)
  close(con)
  write.table(stats_long[, .immportExportColOrder()],
              file_path,
              append = TRUE,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE)
}
