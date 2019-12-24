# fcs.R
# Functions for the import, quality control, and processing of Flow Cytometry
# Standard (FCS) files.

#' Astrolabe debris labels.
#' 
#' @return Vector of Astrolabe debris labels.
#' @export
astrolabeDebrisLabels <- function() c("Debris")

# FCS File Interaction --------------------------------------------------------

#' Test if a given object is an Astrolabe sample.
#'
#' Check whether the given object is a list, and whether it includes all of the
#' fields that we would expect from the default importFcsFile.
#'
#' @param sample Object to be tested.
#' @return TRUE if the object is an Astrolabe sample, FALSE otherwise.
#' @export
isSample <- function(sample) {
  default_fields <- c("exprs", "parameter_name", "parameter_desc")

  if (!is.list(sample)) {
    FALSE
  } else {
    all(default_fields %in% names(sample))
  }
}

.massSuspectValues <- function() {
  # A list of suspect values. If too many of them are included in the
  # description field of an FCS file (and include additional values) we will
  # remove them.
  c(
    "89Y",   "Y89",   "102Pd", "Pd102", "103Rh", "Rh103", "104Pd", "Pd104",
    "105Pd", "Pd105", "106Pd", "Pd106", "108Pd", "Pd108", "110Pd", "Pd110",
    "113In", "In113", "115In", "In115", "118Sn", "Sn118", "120Sn", "Sn120",
    "127I",  "I127",  "131Xe", "Xe131", "133Cs", "Cs133", "138Ba", "Ba138",
    "139La", "La139", "140Ce", "Ce140", "141Pr", "Pr141", "142Nd", "Nd142",
    "142Ce", "Ce142", "143Nd", "Nd143", "144Nd", "Nd144", "145Nd", "Nd145",
    "146Nd", "Nd146", "147Sm", "Sm147", "148Nd", "Nd148", "148Sm", "Sm148",
    "149Sm", "Sm149", "150Nd", "Nd150", "151Eu", "Eu151", "152Gd", "Gd152",
    "152Sm", "Sm152", "153Eu", "Eu153", "154Gd", "Gd154", "154Sm", "Sm154",
    "155Gd", "Gd155", "156Gd", "Gd156", "157Gd", "Gd157", "158Gd", "Gd158",
    "159Tb", "Tb159", "160Dy", "Dy160", "160Gd", "Gd160", "161Dy", "Dy161",
    "162Dy", "Dy162", "162Er", "Er162", "163Dy", "Dy163", "164Er", "Er164",
    "164Dy", "Dy164", "165Ho", "Ho165", "166Er", "Er166", "167Er", "Er167",
    "168Er", "Er168", "169Tm", "Tm169", "170Er", "Er170", "171Yb", "Yb171",
    "172Yb", "Yb172", "173Yb", "Yb173", "174Yb", "Yb174", "175Lu", "Lu175",
    "176Lu", "Lu176", "176Yb", "Yb176", "195Pt", "Pt195", "208Pb", "Pb208",
    "209Bi", "Bi209"
  )
}

.removeMassFromDesc <- function(desc) {
  # Remove mass from channel descriptions.
  mass_str <- paste0("(", paste(.massSuspectValues(), collapse = "|"), ")")
  
  suspect_match <-
    unlist(lapply(desc, function(s) { grepl(mass_str, s) & grepl("_", s) }))
  if (sum(suspect_match) < 2) return(desc);
  
  # Try to remove them using regular expression.
  desc[suspect_match] <-
    unlist(lapply(desc[suspect_match], function(s) {
      s <- gsub(paste0("_", mass_str, "_"), "", s)
      s <- gsub(paste0(mass_str, "_"), "", s)
      s <- gsub(paste0("_", mass_str), "", s)
    }))
  
  return(desc);
}

.removeEqFromDesc <- function(desc) {
  # Remove "EQ" and other bead-related suffix from channel descriptions.
  desc <- unlist(lapply(desc, function(s) gsub(".EQ.BEADS", "", s)))
  desc <- unlist(lapply(desc, function(s) gsub(".BEADS", "", s)))
  desc <- unlist(lapply(desc, function(s) gsub(".EQ", "", s)))

  desc
}

.removeSpecialCasesFromDesc <- function(desc) {
  # Remove special cases which were encountered during our work with various
  # data sets.
  desc[grep("CD56.*CD56", desc)] <- "CD56"
  desc[grep("CD45RA.*CD45RA", desc)] <- "CD45RA"
  desc[grep("CD45RA.*Fluidigm", desc)] <- "CD45RA"
  desc[grep("TCR_gamma_delta", desc)] <- "gdTCR"
  desc[grep("NKp46__CD335", desc)] <- "NKp46"
  desc[desc == "FceR1a"] <- "FceRIa"
  desc[desc == "Singlec-8"] <- "Siglec8"

  # Remove the __v_ suffix from Cytobank viSNE analysis.
  desc <- gsub(" \\(v\\)$", "", desc)

  desc
}

#' Import FCS Channel Information.
#'
#' Import the TEXT section of an FCS file and extract channel names and
#' descriptions (the $PnN and $PnS fields, respectively).
#'
#' @param filename The name of the FCS file to import.
#' @return A dataframe with channel names and descriptions.
#' @export
importFcsChannels <- function(filename) {
  fcs_text <- flowCore::read.FCSheader(filename)
  fcs_text <- fcs_text[[1]]

  n_parameters <- as.integer(fcs_text["$PAR"])
  name_keywords <- paste("$P", seq(n_parameters), "N", sep = "")
  desc_keywords <- paste("$P", seq(n_parameters), "S", sep = "")
  channels <- tibble::tibble(
    Name = fcs_text[name_keywords],
    OrigDesc = fcs_text[desc_keywords]
  )
  channels$Desc <- channels$OrigDesc

  # Copy flow-specific channels from name to desc.
  flow_pattern <- "(^FSC-.$|^SSC-.$)"
  flow_names <- grepl(flow_pattern, channels$Name)
  channels$Desc[flow_names] <- channels$Name[flow_names]
  # Remove masses and EQ suffix from channel descriptions.
  channels$Desc <- .removeMassFromDesc(channels$Desc)
  channels$Desc <- .removeEqFromDesc(channels$Desc)
  channels$Desc <- .removeSpecialCasesFromDesc(channels$Desc)
  # NA/empty descriptions should copy from name.
  channels$Desc[is.na(channels$Desc) | channels$Desc == ""] <-
    channels$Name[is.na(channels$Desc) | channels$Desc == ""]

  # Remove all non-alphanumeric characters from description.
  channels$Desc <- gsub("[^[:alnum:]]", "_", channels$Desc)
  
  channels
}

#' Calculate parameter digest for sample.
#'
#' Apply the digest::digest function to the list of parameter descriptions and
#' parameter names for a given FCS data structure.
#'
#' @param sample An Astrolabe sample.
#' @param parameter_desc,parameter_name Parameter descriptions and names from
#' an FCS file.
#' @return Parameter digest for data.
#' @export
calculateFcsDigest <- function(sample, parameter_name = NULL) {
  if (is.null(parameter_name)) {
    # FCS data parameter.
    if (!isSample(sample)) stop("Expecting an Astrolabe sample")

    parameter_desc <- sample$parameter_desc
    parameter_name <- sample$parameter_name
  } else {
    # Parameters are explicit description and name.
    parameter_desc <- sample
  }

  digest::digest(list(desc = parameter_desc, name = parameter_name))
}

.identifyFcsInstrumentMissingCyt <- function(flow_frame) {
  # Identify instrument if the cyt field is missing.
  names <- flow_frame@parameters@data$name
  desc <- flow_frame@parameters@data$desc

  if (all(c("Ir191Di", "Ir193Di") %in% names)) {
    # Look for iridium channels to identify mass cytometry.
    return("mass_cytometry")
  }

  return("unknown")
}

#' Identify FCS file instrument.
#'
#' Read the CYT field from the FCS file header.
#'
#' @param flow_frame FlowCore flow frame.
#' @return Instrument for flow frame.
#' @export
identifyFcsInstrument <- function(flow_frame) {
  # Identify the instrument of a flow frame.

  cyt <- flow_frame@description$`$CYT`

  if (is.null(cyt)) {
    return(.identifyFcsInstrumentMissingCyt(flow_frame))
  }

  cyt <- tolower(cyt)
  if (grepl("cytof", cyt)) {
    return("mass_cytometry")
  } else if (grepl("aurora", cyt)) {
    return("aurora")
  } else if (grepl("lsrfortessa", cyt)) {
    return("lsr_fortessa")
  } else {
    return("unknown")
  }
}

#' Convert a FlowCore flow_frame to Astrolabe sample.
#'
#' Convert the flow_frame class into orloj's internal FCS list format.
#'
#' The orloj FCS list format will accumulate more fields as analyses are
#' applied to it. For example, pre-processing will add a mask to find the
#' non-bead indices. You can use \code{\link{fcsExprs}} the get the expression
#' data after applying all of the masks that are in the list.
#'
#' @param experiment An Astrolabe experiment.
#' @param filename The name of the FCS file to import.
#' @param flow_frame FlowCore flow frame.
#' @seealso \code{\link{isSample}}, \code{\link{fcsExprs}}
#' @return FCS data, in orloj internal FCS list format.
#' @export
convertFlowFrame <- function(experiment, filename, flow_frame) {
  # Import flow data and channels information.
  channels <- importFcsChannels(filename)
  # Update channel desc based on experiment desc.
  match_indices <- match(channels$Name, experiment$channels$Name)
  channels$Desc[!is.na(match_indices)] <-
    experiment$channels$Desc[match_indices[!is.na(match_indices)]]

  if (any(is.na(channels$Desc))) stop("desc cannot be NA")
  if (any(channels$Desc == "")) stop("desc cannot be empty")
  if (any(duplicated(channels$Desc))) stop("desc cannot have duplicates")

  desc <- channels$Desc
  name <- channels$Name
  exprs <- tibble::as_tibble(flow_frame@exprs)
  colnames(exprs) <- desc

  list(
    exprs = exprs,
    parameter_name = name,
    parameter_desc = desc
  )
}

#' Preprocess an Astrolabe sample.
#'
#' @param experiment An Astrolabe experiment.
#' @param sample An Astrolabe sample.
#' @return Sample after the above steps are done.
#' @export
preprocess <- function(experiment, sample) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  if (experiment$instrument == "mass_cytometry") {
    massPreprocess(sample)
  } else if (experiment$instrument %in% c("aurora", "lsr_fortessa")) {
    auroraPreprocess(sample)
  } else if (experiment$instrument == "citeseq") {
    citeseqPreprocess(sample)
  } else {
    stop("unknown sample instrument")
  }
}
