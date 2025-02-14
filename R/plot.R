#' @importFrom graphics plot
#' @export
plot.sphere_sample <- function(x, radius = 0.95, res = 50, ...) {
  f <- \(x, y, z) { x^2 + y^2 + z^2 }
  x <- y <- z <- seq(-radius, radius, length.out = res)
  g <- expand.grid(x = x, y = y, z = z)
  voxel <- array(with(g, f(x, y, z)), dim = rep(res, 3))

  cont <- misc3d::computeContour3d(voxel, level = radius^2, x = x, y = y, z = z)
  idx <- matrix(0:(nrow(cont) - 1), ncol = 3, byrow = TRUE)

  plotly::plot_ly() |>
    plotly::add_trace(
      x = cont[, 1], y = cont[, 2], z = cont[, 3],
      i = idx[, 1], j = idx[, 2], k = idx[, 3],
      type = "mesh3d", opacity = 0.1, showlegend = FALSE, hoverinfo = "skip"
    ) |>
    plotly::add_markers(
      x = spl[, 1], y = spl[, 2], z = spl[, 3],
      type = "scatter3d", mode = "markers"
    )
}
