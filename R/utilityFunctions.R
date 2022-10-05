
#' Function to do percentile normalizaton
#'
#' @param x Matrix to percentile normilse.
#' @importFrom stats quantile
.percentileNorm <- function(x){
  x <- as.matrix(x)
  percentiles <- quantile(x, probs = c(0.01, 0.99))
  q1 <- percentiles[[1]]
  q99 <- percentiles[[2]]
  
  x <- (x - q1)/(q99)
  return(x)
}

#' Function to do min max normalization
#'
#' @param x Matrix to min max nomalize.
#' @importFrom stats quantile
.minmaxNorm <- function(x){
  x <- as.matrix(x)
  percentiles <- quantile(x, probs = c(0.01, 0.99))
  minVal <- percentiles[[1]]
  maxVal <- percentiles[[2]]
  x <- (x - minVal)/(maxVal - minVal)
  return(x)
}

#' Function to do arsinh normalization
#' 
#' @param x A numeric or complex vector
#' @param cofactor Cofactor of the vector. Default is 5.
.arsinhNnorm <- function(x, cofactor=5){
  x <- asinh(x/cofactor)
  return(x)
}

#' Creates uniformly distributed data of same dimensionality as input data
#' this function was obtained from the Stab package
#'
#' @param data A data matrix.
.uniformData <- function(data) {
  
  # get dimensions of data
  n         <- nrow(data)
  unifDims  <- t(apply(data, 2, range))
  
  # sample data and combine to data matrix
  dataSynt  <- matrix(NA,nrow=nrow(data),ncol=ncol(data))
  for(i in seq_len(ncol(data))) {
    dataSynt[, i] <- stats::runif(n, unifDims[i, 1], unifDims[i, 2])
  }
  return(dataSynt)
}

#' A function to compute the elbow point given a set of points
#'
#' @param vals Values to compute the elbow point of.
.computeElbow <- function(vals){
  diffs <- diff(vals)
  diffs <- diffs[-1]
  optKb <- which.max(abs(diffs))+1
  return(optKb)
}