% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{filenameify}
\alias{filenameify}
\title{Convert a string into a filename-friendly format.}
\usage{
filenameify(str, extension = NULL)
}
\arguments{
\item{str}{A string.}

\item{extension}{Optional: Filename extension.}
}
\value{
A filename-friendly format of that string.
}
\description{
Convert "+" to "pos" and "-" to "neg". Remove paranthesis. Convert any
non-alphanumeric character to underscore.
}
