
# function to estimate the number of clusters using discriminant analysis
# parts of this function is based on the sigclust2 package by Patrick Kimes
# see https://github.com/pkimes/sigclust2
.runDiscriminant <- function(distMat, minClusterSize, alpha=0.05){
  
  # do some house keeping
  n <- nrow(distMat)
  p <- ncol(distMat)
  ndType <- rep("", n-1)
  pEmp <- rep(0, n-1)
  
  numClusters <- 0
  
  # generate the initial tree
  initTree <- fastcluster::hclust(dist(distMat, method = 'maximum'), method = "average")
  
  # process the resulting dendrogram
  hcDat <- initTree
  idxHC <- .idxHc(initTree, n)
  cutoff <- .fwerCutoffmatrix(idxHC, alpha)
  pdMap <- .pdMap(initTree, n)
  ndType <- rep("", n-1)
  
  # run significance testing on each node
  for (k in seq_len(n-1)) {
    ## indices for subtree
    idxVals <- idxHC[k, ]
    idxSub <- unlist(idxHC[k, ])
    nSub <- length(idxSub)
    
    ## only calc p-values for branches w/ more than n_mini
    if (nSub < minClusterSize) {
      ndType[k] <- "n_small"
      next
    }
    
    if ((alpha < 1) && (k > 1) && (ndType[pdMap[k]] != "sig")) {
      ndType[k] <- "no_test"
      pEmp[k] <- 1
      next
    }
    
    # Generate initial assingments
    t <- c(idxVals[[1]], idxVals[[2]])
    xComb <- distMat[t,t]
    assignments <- kmeans(xComb, 2)$cluster
    
    # compute the discriminant projections
    xNew <- fpc::discrcoord(x=xComb, clvecd = assignments)$proj[, 1]
    resPval <- diptest::dip.test(xNew)$p.value
    
    # update results
    if(alpha < 1){
      if(resPval < alpha ){
        ndType[k] <- "sig"
        numClusters = numClusters + 1
      } else{
        ndType[k] <- "not_sig"
      }
      pEmp[k] <- resPval
    }
  }
  return(numClusters)
}

# identify parent node of each node in dendrogram
# parts of this function is based on the sigclust2 package by Patrick Kimes
# see https://github.com/pkimes/sigclust2
.pdMap <- function(hc, n) {
  ## determine parent branch node for all children nodes along dendrogram
  pdPairs <- rbind(cbind(hc$merge[, 1], seq_len(n-1)), 
                    cbind(hc$merge[, 2], seq_len(n-1)))
  pdMap <- data.frame(pdPairs[pdPairs[, 1] > 0, ])
  names(pdMap) <- c("dtr", "prt")
  pdMap <- pdMap$prt[order(pdMap$dtr)] #the parent of each daughter
  pdMap <- c(pdMap, n) #add final node without a parent
  
  ## flip index, hclust and shc use reversed ordering
  n - rev(pdMap)
}


## determine obs indices at each node of the dendrogram
# parts of this function is based on the sigclust2 package by Patrick Kimes
# see https://github.com/pkimes/sigclust2
.idxHc <- function(hc, n) {
  ## list array of cluster indices at each of the n-1 merges
  idxHC <- array(list(), c(2*n-1, 2))
  idxHC[1:n, 1] <- as.list(n:1)
  idxHC[(n+1):(2*n-1), ] <- hc$merge + n + (hc$merge<0)
  
  ## complete idxHC
  for (k in seq_len(n-1)) {
    idxHC[[n+k, 1]] <- unlist(idxHC[idxHC[[n+k, 1]], ])
    idxHC[[n+k, 2]] <- unlist(idxHC[idxHC[[n+k, 2]], ])
  }
  
  ## flip index, hclust and shc use revered ordering
  idxHC[(2*n-1):(n+1), ]
}

# parts of this function is based on the sigclust2 package by Patrick Kimes
# see https://github.com/pkimes/sigclust2
.fwerCutoffmatrix <- function(obj, alpha, ...) {
  alpha/(nrow(obj)+1) *
    apply(obj, 1, function(x) { length(unlist(x)) })
}
