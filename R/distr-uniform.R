#' Coercion method for `unisph_sample` objects
#'
#' @param x A numeric matrix of shape \eqn{n \times 3}.
#'
#' @return An object of class `unisph_sample` and `sphere_sample`.
#'
#' @export
#' @examples
#' x <- runisph(100)
#' as_unisph_sample(x)
as_unisph_sample <- function(x) {
  if (!is.matrix(x) || ncol(x) != 3) {
    cli::cli_abort("{.arg x} must be a numeric matrix of shape n x 3")
  }
  class(x) <- c("unisph_sample", "sphere_sample", class(x))
  x
}

#' The Uniform distribution
#'
#' @param x A numeric matrix of shape \eqn{n \times 3} where \eqn{n} is the
#'   number of samples and the columns represent the x, y, z coordinates
#'   of the samples.
#' @param log A logical value indicating whether the log density should be
#'   returned. Defaults to `false`.
#' @param n An integer value indicating the number of samples to generate.
#'
#' @return
#' - `dunisph` returns a numeric vector of length \eqn{n} containing the
#' density values of the Uniform distribution.
#' - `punisph` returns a numeric vector of length \eqn{n} containing the
#' cumulative density values of the Uniform distribution.
#' - `qunisph` returns a numeric matrix of shape \eqn{n \times 3} containing
#' the quantile values of the Uniform distribution.
#' - `runisph` returns a numeric matrix of shape \eqn{n \times 3} containing
#' the random samples from the Uniform distribution.
#'
#' @examples
#' n <- 100
#' spl <- runisph(n)
#' dunisph(spl)
#' punisph(spl)
#'
#' @name Uniform
NULL

#' @rdname Uniform
#' @export
dunisph <- function(x, log = FALSE) {
  dunisph_impl(x, log)
}

#' @rdname Uniform
#' @export
punisph <- function(x) {
  punisph_impl(x)
}

#' @rdname Uniform
#' @export
runisph <- function(n) {
  as_unisph_sample(runisph_impl(n))
}
