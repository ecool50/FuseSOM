# A function to generate the radii values for the Gaussian kernel
# function was obtained from https://rdrr.io/rforge/yasomi/ with some major modifications

radiusExp <- function(min, max, steps) {
  max * (min / max)^(seq(0, 1, length.out = steps))
}

radiusLin <- function(min, max, steps) {
  seq(max, min, length.out = steps)
}

batchSomControl.default <- function(data, somGrid,
                                    mode = c("continuous", "stepwise"),
                                    minRadius, maxRadius, steps,
                                    decrease = c("power", "linear"), maxIter,
                                    kernel = c("gaussian", "linear"),
                                    normalised,
                                    assignment = c("single", "heskes"),
                                    cut = 1e-07, ...) {
  mode <- match.arg(mode, c("continuous", "stepwise"))
  decrease <- match.arg(decrease, c("power", "linear"))
  kernel <- match.arg(kernel, c("gaussian", "linear"))
  assignment <- match.arg(assignment, c("single", "heskes"))
  if (missing(maxRadius)) {
    maxRadius <- 2 / 3 * somGrid$diam + 1
  }
  if (missing(minRadius)) {
    minRadius <- switch(kernel,
      "gaussian" = 0.5,
      "linear" = 1
    )
  }
  if (maxRadius <= minRadius) {
    stop("maxRadius must be larger than minRadius")
  }
  if (minRadius <= 0) {
    stop("minRadius must be positive")
  }
  if (missing(steps)) {
    steps <- switch(mode,
      "stepwise" = max(2, ceiling(2 * maxRadius)),
      "continuous" = max(20, ceiling(5 * maxRadius))
    )
  }
  if (missing(maxIter)) {
    maxIter <- switch(mode,
      "stepwise" = 75,
      "continuous" = 1
    )
  }
  if (missing(normalised)) {
    normalised <- assignment == "heskes"
  }
  radii <- switch(decrease,
    "linear" = radiusLin(minRadius, maxRadius, steps),
    "power" = radiusExp(minRadius, maxRadius, steps)
  )
  list(
    mode = mode, radii = radii, maxIter = maxIter, kernel = kernel,
    normalised = normalised, assignment = assignment,
    cut = cut
  )
}
