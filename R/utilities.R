# utilities.R
# General R utilities.

`%>%` <- dplyr::`%>%`

#' Name elements of vector with their value
#'
#' Name all of the elements of a vector with their value. This is useful in
#' conjunction with lapply, as it allows you to operate on the values of the
#' vector and naming the list with these values.
#'
#' @param v A vector to be named.
#' @return Same vector, with names set to values.
#' @export
nameVector <- function(v) {
  names(v) <- v
  v
}

#' Format a p-value or FDR into a nice string.
#'
#' Unlike the R \code{\link{format.pval}}, this function refrains from using
#' scientific notation.
#'
#' @param v A p-value or FDR value to format.
#' @param label String describing the value being presented.
#' @param digits Number of decimal places to be used.
#' @return A string of a nicely formatted value.
#' @export
formatPvalue <- function(v, label = "p", digits = 4) {
  thresh <- 10 ^ (-(digits - 1))
  ifelse(v >= thresh,
         paste0(label, " = ", round(v, digits)),
         paste0(label, " < 10^", ceiling(log10(v))))
}

#' Add backticks to a string.
#'
#' Add backticks at the beginning and end of a string. This avoids situation
#' where lazy evaluation fails due to column names which are not legal as R
#' variable names.
#'
#' @param str A string.
#' @return String surrounded by backticks.
#' @export
addBackticks <- function(str) {
  if (is.null(str)) {
    NULL
  }
  else if (grepl("^`", str)) {
    str
  } else {
    paste0("`", str, "`")
  }
}

#' Convert a string into a filename-friendly format.
#'
#' Convert "+" to "pos" and "-" to "neg". Remove paranthesis. Convert any
#' non-alphanumeric character to underscore.
#'
#' @param str A string.
#' @param extension Optional: Filename extension.
#' @return A filename-friendly format of that string.
#' @export
filenameify <- function(str, extension = NULL) {
  str <- gsub("\\+$", "pos", str)
  str <- gsub("\\+ ", "pos ", str)
  str <- gsub("\\+)", "pos)", str)
  str <- gsub("\\+_", "pos_", str)

  str <- gsub("\\-$", "neg", str)
  str <- gsub("\\- ", "neg ", str)
  str <- gsub("\\-)", "neg)", str)
  str <- gsub("\\-_", "neg_", str)

  str <- gsub("\\)", "", str)
  str <- gsub("\\(", "", str)

  str <- gsub("[^[:alnum:]]", "_", str)

  if (!is.null(extension)) str <- paste0(str, ".", extension)

  str
}

#' Find indices of standard mass cytometry channels.
#'
#' Find the indices of a set of standard mass cytometry channels, including
#' Time, Event Length, Ir191/193 (DNA), and Live/Dead.
#'
#' @param sample An Astrolabe sample.
#' @return A list with indices for each of the standard channels.
#' @export
findStandardMassCytometryChannels <- function(sample) {
  # Allow parameter_name as input instead of sample.
  if (isSample(sample)) {
    parameter_name <- sample$parameter_name
  } else parameter_name <- sample

  time_idx <- grep("time", parameter_name, ignore.case = TRUE)

  # Event length or cell length.
  event_length_idx <-
    intersect(
      grep("event", parameter_name, ignore.case = TRUE),
      grep("length", parameter_name, ignore.case = TRUE)
    )
  if (length(event_length_idx) == 0) {
    event_length_idx <-
      intersect(
        grep("cell", parameter_name, ignore.case = TRUE),
        grep("length", parameter_name, ignore.case = TRUE)
      )
  }

  dna191_idx <- which(parameter_name == "Ir191Di")
  dna193_idx <- which(parameter_name == "Ir193Di")

  # Set missing values to NULL.
  if (length(time_idx) == 0)          time_idx <- NULL
  if (length(event_length_idx) == 0)  event_length_idx <- NULL
  if (length(dna191_idx) == 0)        dna191_idx <- NULL
  if (length(dna193_idx) == 0)        dna193_idx <- NULL

  list(
    time_idx = time_idx,
    event_length_idx = event_length_idx,
    dna191_idx = dna191_idx,
    dna193_idx = dna193_idx
  )
}
