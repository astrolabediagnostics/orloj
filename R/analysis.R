# analysis.R
# Functions for the analysis of cytometry data and of the output of clustering,
# dimensionality reduction, and other machine learning algorithms.

#' Jensen-Shannon Divergence
#'
#' Compute the Jensen-Shannon (JS) divergence between two one-dimensional
#' probability distributions. Original code adopted from:
#' http://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r
#'
#' @param p,q Density vectors of a one-dimensional probability distribution.
#' @return JS divergence between p and q.
#' @export
jsDivergence <- function(p, q) {
  # Verify p and q come from probability distributions.
  if (is.character(all.equal(sum(p), 1)) ||
      is.character(all.equal(sum(q), 1)) ||
      any(p < 0) ||
      any(q < 0)) {
    stop("p and q should be probability distributions")
  }

  m <- 0.5 * (p + q)
  0.5 * (sum(p * log(p / m)) + sum(q * log(q / m)))
}

#' Vectorized Two-dimensional Kernel Density Estimation
#'
#' A wrapper around \code{\link[MASS]{kde2d}} that converts the density
#' estimation into a one-dimensional probability distribution.
#'
#' @param mtx A numeric matrix which provides the data for the kernel density
#' estimation.
#' @param n Number of grid points in each direction. Can be a scalar or a
#' length-2 integer vector.
#' @param lims The limits of the rectangle covered by the grid as
#' c(xl, xu, yl, yu).
#' @return A vector for the one-dimensional probability distribution.
#' @export
kde2dVec <- function(mtx, n = 25, lims = c(range(mtx[, 1], mtx[, 2]))) {
  if (ncol(mtx) != 2) stop("mtx should be a two-dimensional matrix")

  density <- MASS::kde2d(mtx[, 1], mtx[, 2], n = n, lims = lims)
  density_vector <- as.vector(density$z)
  density_vector / sum(density_vector)
}

#' Jensen-Shannon Divergence Between Matrices
#'
#' Compute the Jensen-Shannon (JS) divergence between two matrices. The matrices
#' are each converted into a one-dimensional probability distribution which are
#' then feed into the JS equation.
#'
#' @param mtx1,mtx2 Numeric matrices which provide the data for the JS
#' divergence calculation.
#' @inheritParams kde2dVec
#' @return JS divergence between mtx1 and mtx2.
#' @export
jsDivergenceMtx <- function(mtx1, mtx2, n = 25) {
  if (ncol(mtx1) != 2 || ncol(mtx2) != 2) {
    stop("mtx1 and mtx2 should be two-dimensional matrices")
  }

  # Find lims for combined mtx1, mtx2.
  combined_mtx <- rbind(mtx1, mtx2)
  lims <- c(range(combined_mtx[, 1]), range(combined_mtx[, 2]))

  # Convert matrices to probability distributions.
  mtx1_prob <- kde2dVec(mtx1, n, lims)
  mtx2_prob <- kde2dVec(mtx2, n, lims)

  # Compute JS divergence
  jsDivergence(mtx1_prob, mtx2_prob)
}
