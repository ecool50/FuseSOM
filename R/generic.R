# Definition for generic functions that computes the self organizing maps
# function was obtained from https://rdrr.io/rforge/yasomi/ with some major modifications

somInitPca <- function(data,somGrid,...) {
  UseMethod("somInitPca")
}

batchSom <- function(data,somGrid,init="pca",prototypes,weights,
                     mode="continuous",
                     minRadius, maxRadius, steps,
                     decrease="power", maxIter,
                     kernel="gaussian", normalised,
                     assignment="single",
                     cut=1e-07,
                     verbose=FALSE,keepdata=TRUE,...) {
  UseMethod("batchSom")
}

## in annealing.R

batchSomControl <- function(data,somGrid,
                             mode="continuous",
                             minRadius, maxRadius, steps,
                             decrease="power", maxIter,
                             kernel="gaussian",
                             assignment="single",
                             cut = 1e-07,...) {
  UseMethod("batchSomControl")
}
