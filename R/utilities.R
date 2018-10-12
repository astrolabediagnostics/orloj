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
#' @return A filename-friendly format of that string.
#' @export
filenameify <- function(str) {
  str <- gsub("\\+$", "pos", str)
  str <- gsub("\\+ ", "pos ", str)
  str <- gsub("\\+)", "pos)", str)

  str <- gsub("\\-$", "neg", str)
  str <- gsub("\\- ", "neg ", str)
  str <- gsub("\\-)", "neg)", str)

  str <- gsub("\\)", "", str)
  str <- gsub("\\(", "", str)

  str <- gsub("[^[:alnum:]]", "_", str)

  str
}

