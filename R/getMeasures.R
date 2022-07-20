
# a function to compute the various measures when estimating the number of clusters
# see https://arxiv.org/abs/1608.07494
# function was obtained from https://github.com/cran/cstab with some minor modifications
.getMeasures <- function(data,
                        k,
                        linkage='average',
                        measures=c('wcd','sil','mse'))
  
{
  
  pear <- stats::cor(t(data), method = 'pearson')
  cosi <- coop::tcosine(data)
  spear <- stats::cor(t(data), method = 'spearman')
  
  dmat <- as.matrix(fuse(cor2dist(pear),cor2dist(cosi),cor2dist(spear)))
  
  hcobj <- fastcluster::hclust(stats::as.dist(dmat), method = linkage)
  cl <- stats::cutree(hcobj, k)
  
  # ------ WITHIN CLUSTER DISSIMILARITY
  
  # calc within cluster dissimilarity
  WCD = NULL
  if('wcd' %in% measures){
    norm_diss <- c()
    for(i in 1:k) norm_diss[i] <- sum(dmat[cl==i, cl==i]) / (2 * sum(cl == i))
    WCD <- sum(norm_diss)
  }
  
  
  # ------ SILHOUETTE
  
  Sil <- NULL
  if('sil' %in% measures) if(k > 1) Sil <- mean(cluster::silhouette(cl, stats::as.dist(dmat))[,3]) else Sil <- 0
  
  # ------ CLUSTER CENTERS
  
  # get centers
  centers = NULL
  if('mse' %in% measures | 'centers' %in% measures){
    if(k == 1)  centers <- colMeans(data)
    if(k > 1) {
        data_or <- as.data.frame(data)
        data_or$cl <- cl
        data_or <- data_or[order(data_or$cl),]
        centers <- matrix(NA, k, ncol(data))
        for(kk in 1:k) centers[kk,] <- colMeans(data_or[data_or$cl==kk, -(ncol(data)+1)])
      }
    }
  # ------ MSE FOR JUMP
  
  MSE = NULL
  if('mse' %in% measures) {
    if(k == 1) SE <- apply(data, 1, function(inst){diffs <- inst - centers; return(t(diffs) %*% diffs)})
    if(k > 1) {
      SE <- numeric()
      for(i in 1:nrow(data)) {
        diffs <- data[i,] - centers[cl[i],]
        se    <- t(diffs) %*% diffs
        SE[i] <- se
      }
    }
    MSE <- sum(SE) / nrow(data)
  }
  
  
  # ------ RETURN
  
  outlist <- list('WCD' = WCD, 'Sil'=Sil, 'MSE'=MSE, 'Centers'=centers)
  return(outlist)
}
