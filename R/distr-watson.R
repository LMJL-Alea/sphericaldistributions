#' Coercion method for `watson_sample` objects
#'
#' @param x A numeric matrix of shape \eqn{n \times 3}.
#'
#' @return An object of class `watson_sample` and `sphere_sample`.
#'
#' @export
#' @examples
#' x <- rwatson(100, c(1, 0, 0), 10)
#' as_watson_sample(x)
as_watson_sample <- function(x) {
  if (!is.matrix(x) || ncol(x) != 3) {
    cli::cli_abort("{.arg x} must be a numeric matrix of shape n x 3")
  }
  class(x) <- c("watson_sample", "sphere_sample", class(x))
  x
}

#' The Watson distribution
#'
#' @param x A numeric matrix of shape \eqn{n \times 3} where \eqn{n} is the
#'   number of samples and the columns represent the x, y, z coordinates
#'   of the samples.
#' @param mu A numeric matrix of shape \eqn{1 \times 3} representing the mean
#'   axis of the Watson distribution.
#' @param kappa A positive numeric value representing the concentration
#'   parameter of the Watson distribution.
#' @param log A logical value indicating whether the log density should be
#'   returned. Defaults to `false`.
#' @param n An integer value indicating the number of samples to generate.
#'
#' @returns
#' - `dwatson` returns a numeric vector of length \eqn{n} containing the
#' density values of the Watson distribution.
#' - `pwatson` returns a numeric vector of length \eqn{n} containing the
#' cumulative density values of the Watson distribution.
#' - `qwatson` returns a numeric matrix of shape \eqn{n \times 3} containing
#' the quantile values of the Watson distribution.
#' - `rwatson` returns a numeric matrix of shape \eqn{n \times 3} containing
#' the random samples from the Watson distribution.
#'
#' @examples
#' mu <- c(1, 0, 0)
#' kappa <- 10
#' n <- 100
#' spl <- rwatson(n, mu, kappa)
#' dwatson(spl, mu, kappa)
#' pwatson(spl, mu, kappa)
#'
#' @name watson
NULL

#' @rdname watson
#' @export
dwatson <- function(x, mu, kappa, log = FALSE) {
  dwatson_impl(x, mu, kappa, log)
}

#' @rdname watson
#' @export
pwatson <- function(x, mu, kappa) {
  pwatson_impl(x, mu, kappa)
}

#' @rdname watson
#' @export
rwatson <- function(n, mu, kappa) {
  as_watson_sample(rwatson_impl(n, mu, kappa))
}

#' Compute the mean of a `watson_sample` object
#'
#' @param x A numeric matrix of shape \eqn{n \times 3}.
#' @param ... Additional arguments. Not used.
#'
#' @return A numeric matrix of shape \eqn{1 \times 3} storing the mean of the
#'   `watson_sample` object.
#'
#' @export
#' @examples
#' x <- rwatson(100, c(1, 0, 0), 10)
#' mean(x)
mean.watson_sample <- function(x, ...) {
  mean_watson_impl(x)
}
