

somInitPca <- function(data,somgrid,...) {
  UseMethod("somInitPca")
}

batchSom <- function(data,somGrid,init=c("pca"),prototypes,weights,
                     mode = c("continuous","stepwise"),
                     minRadius, maxRadius, steps,
                     decrease = c("power", "linear"), maxIter,
                     kernel = c("gaussian", "linear"), normalised,
                     assignment = c("single", "heskes"),
                     cut = 1e-07,
                     verbose=FALSE,keepdata=TRUE,...) {
  UseMethod("batchSom")
}

## in annealing.R

batchSomControl <- function(data,somGrid,
                             mode = c("continuous","stepwise"),
                             minRadius, maxRadius, steps,
                             decrease = c("power", "linear"), maxIter,
                             kernel = c("gaussian", "linear"),
                             normalised,
                             assignment = c("single", "heskes"),
                             cut = 1e-07,...) {
  UseMethod("batchSomControl")
}