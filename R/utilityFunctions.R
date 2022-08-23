
# function to do percentile normalizaton
.percentileNorm <- function(x){
  x <- as.matrix(x)
  percentiles <- quantile(x, probs = c(0.01, 0.99))
  q1 <- percentiles[[1]]
  q99 <- percentiles[[2]]
  
  x <- (x - q1)/(q99)
  return(x)
}

# function to do min max normalization
.minmaxNorm <- function(x){
  x <- as.matrix(x)
  percentiles <- quantile(x, probs = c(0.01, 0.99))
  minVal <- percentiles[[1]]
  maxVal <- percentiles[[2]]
  x <- (x - minVal)/(maxVal - minVal)
  return(x)
}

# function to do arsinh normalization
.arsinhNnorm <- function(x, cofactor=5){
  x <- asinh(x/cofactor)
  return(x)
}

# Creates uniformly distributed data of same dimensionality as input data
# this function was obtained from the Stab package
.uniformData <- function(data) {
  
  # get dimensions of data
  n         <- nrow(data)
  unifDims  <- t(apply(data, 2, range))
  
  # sample data and combine to data matrix
  dataSynt  <- matrix(NA,nrow=nrow(data),ncol=ncol(data))
  for(i in 1:ncol(data)) dataSynt[, i] <- stats::runif(n, unifDims[i, 1], unifDims[i, 2])
  return(dataSynt)
}

# a function to compute the elbow point given a set of points
.computeElbow <- function(vals){
  diffs <- diff(vals)
  diffs <- diffs[-1]
  optKb <- which.max(abs(diffs))+1
  return(optKb)
}