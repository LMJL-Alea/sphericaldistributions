create_sphere <- function(
  radius = 1,
  center = c(0, 0, 0),
  n = 50,
  opacity = 1,
  title = "3D Sphere Plot"
) {
  u <- seq(0, 2 * pi, length.out = n)
  v <- seq(0, pi, length.out = n)

  x <- outer(cos(u), sin(v)) * radius + center[1]
  y <- outer(sin(u), sin(v)) * radius + center[2]
  z <- outer(rep(1, n), cos(v)) * radius + center[3]

  p <- plotly::plot_ly(
    x = ~x,
    y = ~y,
    z = ~z,
    type = "surface",
    opacity = opacity
  )
  p <- plotly::layout(
    p,
    title = title,
    scene = list(
      xaxis = list(title = "X Axis"),
      yaxis = list(title = "Y Axis"),
      zaxis = list(title = "Z Axis")
    )
  )
  p <- plotly::style(p, hoverinfo = 'none')
  p <- plotly::hide_guides(p)
  p
}

#' Plot sphere sample
#'
#' @description This function creates a 3D scatter plot of points on a sphere
#'   using `plotly`. It visualizes a sample of points on a sphere, which can be
#'   generated using various functions like `rvmf()`, `rwatson()`, or similar
#'   functions that generate random samples on a sphere. The function allows
#'   customization of the sphere's resolution, opacity, and title of the plot.
#'
#' @param x An object of class `sphere_sample` specifying a sample of points on
#'   a sphere. It is typically produced by functions like `rvmf()`, `rwatson()`,
#'   or similar functions that generate random samples on a sphere. One can also
#'   add the `sphere_sample` class to a numeric matrix of shape \eqn{n \times 3}
#'   where \eqn{n} is the number of samples and the columns represent the x, y,
#'   z coordinates of the samples.
#' @param opacity A numeric value between 0 and 1 specifying the opacity of the
#'   sphere surface. Defaults to `0.1`.
#' @param title A string specifying the title of the plot. Defaults to `"3D
#'   Sphere Plot"`.
#' @param ... Additional arguments passed to `plotly::add_markers()` to
#'   customize the appearance of the points in the plot, such as `color`,
#'   `size`, etc.
#'
#' @returns
#' A plotly object representing a 3D scatter plot of points on a sphere.
#'
#' @importFrom graphics plot
#' @export
#'
#' @examples
#' x <- rwatson(100, c(1, 0, 0), 10)
#' as_watson_sample(x)
#' plot(x)
#' plot(x, color = I("darkblue"))
plot.sphere_sample <- function(
  x,
  opacity = 0.1,
  title = "3D Sphere Plot",
  ...
) {
  p <- create_sphere(radius = 0.95, n = 50L, opacity = opacity, title = title)

  p <- plotly::add_markers(
    p,
    x = x[, 1],
    y = x[, 2],
    z = x[, 3],
    type = "scatter3d",
    opacity = 1,
    ...
  )

  # use uniform color in trace 0
  p <- plotly::style(
    p,
    colorscale = list(c(0, 1), c("black", "black")),
    traces = 1
  )
  p
}
