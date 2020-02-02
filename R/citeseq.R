# citeseq.R
# Pre-processing citeseq.

#' Verify File is CITEseq.
#'
#' Verify that a given file is a CITEseq file by examining the header and
#' looking for "cell_barcode", "adt_", and "gex_" values.
#'
#' @param filename The name of the file to test.
#' @param min_adt Minimum number of ADT columns.
#' @param min_gex Minimum number of GEX columns.
#' @return TRUE if file is CITEseq, FALSE otherwise.
#' @export
isFileCiteseq <- function(filename, min_adt = 5, min_gex = 500) {
  con <- file(filename, "r")
  line <- readLines(con, n = 1)
  close(con)
  line <- strsplit(line, "\t")[[1]]

  line[1] == "cell_barcode" &&
    length(grep("adt_", line)) > 5 &&
    length(grep("gex_", line)) > 500
}

#' Import CITEseq Channel Information.
#'
#' Import the antibody portion of CITEseq file.
#'
#' @param filename The name of the CITEseq file to import.
#' @return A dataframe with channel name and description.
#' @export
importCiteseqChannels <- function(filename) {
  if (!isFileCiteseq(filename)) stop(paste0(filename, " is not a CITEseq file"))

  # Import all "adt_" column headers from file.
  con <- file(filename, "r")
  line <- readLines(con, n = 1)
  close(con)
  line <- strsplit(line, "\t")[[1]]
  line <- grep("adt_", line, value = TRUE)

  # Set up data frame, remove adt_ prefix and all non-alphanumeric characters.
  channels <- data.frame(Name = line, OrigDesc = line, stringsAsFactors = FALSE)
  channels$Desc <- channels$OrigDesc
  channels$Desc <- gsub("adt_", "", channels$Desc)
  # Remove all non-alphanumeric characters from description.
  channels$Desc <- gsub("[^[:alnum:]]", "_", channels$Desc)
  # Remove double underscores and terminating underscores.
  channels$Desc <- gsub("__", "_", channels$Desc)
  channels$Desc <- gsub("_$", "", channels$Desc)

  channels
}

#' Preprocess citeseq sample.
#'
#' Quality control, cleaning, and transformation for citeseq sample.
#'
#' @param sample An Astrolabe sample.
#' @param cofactor Cofactor for asinh transformation.
#' @return Sample after the above steps are done.
#' @export
citeseqPreprocess <- function(sample, cofactor = 5) {
  if (!isSample(sample)) stop("Expecting an Astrolabe sample")

  # Transform all channels.
  sample$exprs <- asinh(sample$exprs / cofactor)

  sample
}
