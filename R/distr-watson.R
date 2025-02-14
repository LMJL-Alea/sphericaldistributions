#' @export
as_watson_sample <- function(x) {
  if (!is.matrix(x) || ncol(x) != 3) {
    cli::cli_abort("{.arg x} must be a numeric matrix of shape n x 3")
  }
  class(x) <- c("watson_sample", "sphere_sample", class(x))
  x
}

#' @export
mean.watson_sample <- function(x, ...) {
  mean_watson_impl(x)
}
