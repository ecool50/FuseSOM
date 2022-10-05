# A function to generate the grid for the self organizing map
# function was obtained from https://rdrr.io/rforge/yasomi/ with some major modifications
#' @importFrom proxy dist
somGrid <- function(xDim, yDim, topo = c("rectangular", "hexagonal"), with.dist = TRUE) {
  topo <- match.arg(topo)
  if (xDim == 1) {
    tmp <- xDim
    xDim <- yDim
    yDim <- tmp
  }
  x <- seq(from = 1, by = 1, length.out = xDim)
  if (topo == "hexagonal" && yDim > 1) {
    y <- rev(seq(from = 1, by = sqrt(3) / 2, length.out = yDim))
  } else {
    y <- rev(seq(from = 1, by = 1, length.out = yDim))
  }
  pts <- as.matrix(expand.grid(x = x, y = y))
  if (topo == "hexagonal" && yDim > 1) {
    pts[, 1] <- pts[, 1] + rep(c(0, 0.5), each = xDim, length.out = nrow(pts))
  }
  if (topo == "rectangular") {
    diam <- sqrt(xDim^2 + yDim^2)
  } else {
    diam <- sqrt(sum((pts[1, ] - pts[nrow(pts), ])^2))
  }
  res <- list(
    pts = pts, xDim = xDim, yDim = yDim, topo = topo,
    size = xDim * yDim, diam = diam
  )
  if (with.dist) {
    res$dist <- as.matrix(dist(pts, method = "Euclidean"), diag = 0)
  }
  class(res) <- "somgrid"
  res
}
