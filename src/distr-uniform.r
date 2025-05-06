#' Coercion method for `uni_sample` objects
#'
#' @param x A numeric matrix of shape \eqn{n \times 3}.
#'
#' @return An object of class `uni_sample` and `sphere_sample`.
#'
#' @export
#' @examples
#' x <- runi(100, c(1, 0, 0), 10)
#' as_uni_sample(x)
as_uni_sample <- function(x) {
  if (!is.matrix(x) || ncol(x) != 3) {
    cli::cli_abort("{.arg x} must be a numeric matrix of shape n x 3")
  }
  class(x) <- c("uni_sample", "sphere_sample", class(x))
  x
}

#' Compute the mean of a `uni_sample` object
#'
#' @param x A numeric matrix of shape \eqn{n \times 3}.
#' @param ... Additional arguments. Not used.
#'
#' @return A numeric matrix of shape \eqn{1 \times 3} storing the mean of the
#'   `uni_sample` object.
#'
#' @export
#' @examples
#' x <- runi(100, c(1, 0, 0), 10)
#' mean(as_uni_sample(x))
mean.uni_sample <- function(x, ...) {
  mean_uni_impl(x)
}