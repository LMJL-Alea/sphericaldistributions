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
  class(x) <- c("vmf_sample", "sphere_sample", class(x))
  x
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
#' mean(as_vmf_sample(x))
mean.vmf_sample <- function(x, ...) {
  mean_vmf_impl(x)
}
