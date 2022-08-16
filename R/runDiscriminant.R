
# function to estimate the number of clusters using discriminant analysis
# parts of this function is based on the sigclust2 package by Patrick Kimes
# see https://github.com/pkimes/sigclust2
.runDiscriminant <- function(dist_mat, minClusterSize, alpha=0.001){
  
  # do some house keeping
  n <- nrow(dist_mat)
  p <- ncol(dist_mat)
  nd_type <- rep("", n-1)
  p_emp <- rep(0, n-1)
  
  numClusters <- 0
  
  # generate the initial tree
  init_tree <- fastcluster::hclust(dist(dist_mat, method = 'maximum'), method = "average")
  
  # process the resulting dendrogram
  hc_dat <- init_tree
  idx_hc <- .idxHc(init_tree, n)
  cutoff <- .fwerCutoffmatrix(idx_hc, alpha)
  pd_map <- .pdMap(init_tree, n)
  nd_type <- rep("", n-1)
  
  # run significance testing on each node
  for (k in seq_len(n-1)) {
    ## indices for subtree
    idx_vals <- idx_hc[k, ]
    idx_sub <- unlist(idx_hc[k, ])
    n_sub <- length(idx_sub)
    
    ## only calc p-values for branches w/ more than n_mini
    if (n_sub < minClusterSize) {
      nd_type[k] <- "n_small"
      next
    }
    
    if ((alpha < 1) && (k > 1) && (nd_type[pd_map[k]] != "sig")) {
      nd_type[k] <- "no_test"
      p_emp[k] <- 1
      next
    }
    
    # Generate initial assingments
    t <- c(idx_vals[[1]], idx_vals[[2]])
    x_comb <- dist_mat[t,t]
    assignments <- kmeans(x_comb, 2)$cluster
    
    # compute the discriminant projections
    x_new <- fpc::discrcoord(x=x_comb, clvecd = assignments)$proj[, 1]
    res.pval <- diptest::dip.test(x_new)$p.value
    
    # update results
    if(alpha < 1){
      if(res.pval < alpha ){
        nd_type[k] <- "sig"
        numClusters = numClusters + 1
      } else{
        nd_type[k] <- "not_sig"
      }
      p_emp[k] <- res.pval
    }
  }
  return(numClusters)
}

# identify parent node of each node in dendrogram
# parts of this function is based on the sigclust2 package by Patrick Kimes
# see https://github.com/pkimes/sigclust2
.pdMap <- function(hc, n) {
  ## determine parent branch node for all children nodes along dendrogram
  pd_pairs <- rbind(cbind(hc$merge[, 1], seq_len(n-1)), 
                    cbind(hc$merge[, 2], seq_len(n-1)))
  pd_map <- data.frame(pd_pairs[pd_pairs[, 1] > 0, ])
  names(pd_map) <- c("dtr", "prt")
  pd_map <- pd_map$prt[order(pd_map$dtr)] #the parent of each daughter
  pd_map <- c(pd_map, n) #add final node without a parent
  
  ## flip index, hclust and shc use reversed ordering
  n - rev(pd_map)
}


## determine obs indices at each node of the dendrogram
# parts of this function is based on the sigclust2 package by Patrick Kimes
# see https://github.com/pkimes/sigclust2
.idxHc <- function(hc, n) {
  ## list array of cluster indices at each of the n-1 merges
  idx_hc <- array(list(), c(2*n-1, 2))
  idx_hc[1:n, 1] <- as.list(n:1)
  idx_hc[(n+1):(2*n-1), ] <- hc$merge + n + (hc$merge<0)
  
  ## complete idx_hc
  for (k in seq_len(n-1)) {
    idx_hc[[n+k, 1]] <- unlist(idx_hc[idx_hc[[n+k, 1]], ])
    idx_hc[[n+k, 2]] <- unlist(idx_hc[idx_hc[[n+k, 2]], ])
  }
  
  ## flip index, hclust and shc use revered ordering
  idx_hc[(2*n-1):(n+1), ]
}

# parts of this function is based on the sigclust2 package by Patrick Kimes
# see https://github.com/pkimes/sigclust2
.fwerCutoffmatrix <- function(obj, alpha, ...) {
  alpha/(nrow(obj)+1) *
    apply(obj, 1, function(x) { length(unlist(x)) })
}
