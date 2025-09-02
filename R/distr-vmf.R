#' Coercion method for `vmf_sample` objects
#'
#' @param x A numeric matrix of shape \eqn{n \times 3}.
#'
#' @return An object of class `vmf_sample` and `sphere_sample`.
#'
#' @export
#' @examples
#' x <- rvmf(100, c(1, 0, 0), 10)
#' as_vmf_sample(x)
as_vmf_sample <- function(x) {
  if (!is.matrix(x) || ncol(x) != 3) {
    cli::cli_abort("{.arg x} must be a numeric matrix of shape n x 3")
  }
  if (inherits(x, "sphere_sample")) {
    class(x) <- c("vmf_sample", class(x)[-1])
  } else {
    class(x) <- c("vmf_sample", "sphere_sample", class(x))
  }
  x
}

#' The Von-Mises Fisher distribution
#'
#' @param x A numeric matrix of shape \eqn{n \times 3} where \eqn{n} is the
#'   number of samples and the columns represent the x, y, z coordinates
#'   of the samples.
#' @param mu A numeric matrix of shape \eqn{1 \times 3} representing the mean
#'   axis of the Von Mises-Fisher distribution.
#' @param kappa A positive numeric value representing the concentration
#'   parameter of the Von Mises-Fisher distribution.
#' @param log A logical value indicating whether the log density should be
#'   returned. Defaults to `false`.
#' @param n An integer value indicating the number of samples to generate.
#'
#' @returns
#' - `dvmf` returns a numeric vector of length \eqn{n} containing the
#' density values of the Von Mises-Fisher distribution.
#' - `pvmf` returns a numeric vector of length \eqn{n} containing the
#' cumulative density values of the Von Mises-Fisher distribution.
#' - `qvmf` returns a numeric matrix of shape \eqn{n \times 3} containing
#' the quantile values of the Von Mises-Fisher distribution.
#' - `rvmf` returns a numeric matrix of shape \eqn{n \times 3} containing
#' the random samples from the Von Mises-Fisher distribution.
#'
#' @examples
#' mu <- c(1, 0, 0)
#' kappa <- 10
#' n <- 100
#' spl <- rvmf(n, mu, kappa)
#' dvmf(spl, mu, kappa)
#' pvmf(spl, mu, kappa)
#'
#' @name VonMisesFisher
NULL

#' @rdname VonMisesFisher
#' @export
dvmf <- function(x, mu, kappa, log = FALSE) {
  dvmf_impl(x, mu, kappa, log)
}

#' @rdname VonMisesFisher
#' @export
pvmf <- function(x, mu, kappa) {
  pvmf_impl(x, mu, kappa)
}

#' @rdname VonMisesFisher
#' @export
rvmf <- function(n, mu, kappa) {
  as_vmf_sample(rvmf_impl(n, mu, kappa))
}

#' Fit a von Mises-Fisher distribution via MLE
#'
#' @param object An object of class `vmf_sample` specifying a sample of
#'   directions on the 3-sphere.
#' @param ... Additional arguments. Currently unused; must be empty.
#'
#' @return A list containing the fitted parameters `mu` (mean direction) and
#'   `kappa` (concentration parameter) of the von Mises-Fisher distribution.
#'
#' @export
#' @examples
#' x <- rvmf(100, c(1, 0, 0), 10)
#' fit(x)
fit.vmf_sample <- function(object, ...) {
  fit_vmf_impl(object)
}

#' Compute the mean of a `vmf_sample` object
#'
#' @param x A numeric matrix of shape \eqn{n \times 3}.
#' @param ... Additional arguments. Not used.
#'
#' @return A numeric matrix of shape \eqn{1 \times 3} storing the mean of the
#'   `vmf_sample` object.
#'
#' @export
#' @examples
#' x <- rvmf(100, c(1, 0, 0), 10)
#' mean(x)
mean.vmf_sample <- function(x, ...) {
  l <- fit_vmf_impl(x)
  l$mu
}
