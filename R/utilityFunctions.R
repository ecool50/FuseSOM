
# function to do percentile normalizaton
.percentileNorm <- function(x){
  x <- as.matrix(x)
  percentiles <- quantile(x, probs = c(0.01, 0.99))
  q_1 <- percentiles[[1]]
  q_99 <- percentiles[[2]]
  
  x <- (x - q_1)/(q_99)
  return(x)
}

# function to do min max normalization
.minmaxNorm <- function(x){
  x <- as.matrix(x)
  percentiles <- quantile(x, probs = c(0.01, 0.99))
  min_val <- percentiles[[1]]
  max_val <- percentiles[[2]]
  x <- (x - min_val)/(max_val - min_val)
  return(x)
}

# function to do arsinh normalization
.arsinhNnorm <- function(x, cofactor=5){
  x <- asinh(x/cofactor)
  return(x)
}

# this function was obtained from the Stab package
.instabLookup = function(x,y){
  lkup <- lookup()
  a = stabExp(x,lkup)
  b = stabExp(y,lkup)
  return(a * (1-b) + (1-a) * b)
}

# Creates uniformly distributed data of same dimensionality as input data
# this function was obtained from the Stab package
.uniformData <- function(data) {
  
  # get dimensions of data
  n         <- nrow(data)
  unifdims  <- t(apply(data, 2, range))
  
  # sample data and combine to data matrix
  data_synt  <- matrix(NA,nrow=nrow(data),ncol=ncol(data))
  for(i in 1:ncol(data)) data_synt[, i] <- stats::runif(n, unifdims[i, 1], unifdims[i, 2])
  return(data_synt)
}

# a function to compute the elbow point given a set of points
.computeElbow <- function(vals){
  diffs <- diff(vals)
  diffs <- diffs[-1]
  optkb <- which.max(abs(diffs))+1
  return(optkb)
}