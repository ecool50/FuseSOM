percentile_norm <- function(x){
  x <- as.matrix(x)
  percentiles <- quantile(x, probs = c(0.01, 0.99))
  q_1 <- percentiles[[1]]
  q_99 <- percentiles[[2]]
  
  x <- (x - q_1)/(q_99)
  return(x)
}

minmax_norm <- function(x){
  x <- as.matrix(x)
  percentiles <- quantile(x, probs = c(0.01, 0.99))
  min_val <- percentiles[[1]]
  max_val <- percentiles[[2]]
  x <- (x - min_val)/(max_val - min_val)
  return(x)
}


arsinh_norm <- function(x, cofactor=5){
  x <- asinh(x/cofactor)
  return(x)
}


#' Estimate the optimal grid size for the Self Organizing Map
#' 
#' The function finds the eigenvalues of the sample covariance matrix. 
#' It will then return the number of significant eigenvalues according to 
#' the Tracy-Widom test.
#' The function is based on the estKW function fromt he SC3 package
#' @param dataset The optimal grid size.
#' @return the optimal grid size.
#' 
#' 
ComputeGridSize <- function(dataset) {
  
  p <- ncol(dataset)
  n <- nrow(dataset)
  
  # compute Tracy-Widom bound
  x <- scale(dataset)
  muTW <- (sqrt(n - 1) + sqrt(p))^2
  sigmaTW <- (sqrt(n - 1) + sqrt(p)) * (1/sqrt(n - 1) + 1/sqrt(p))^(1/3)
  sigmaHatNaive <- t(x)%*%x  # x left-multiplied by its transpose
  print(dim(sigmaHatNaive))
  bd <- 3.273 * sigmaTW + muTW  # 3.2730 is the p=0.001 percentile point for the Tracy-Widom distribution
  
  # compute eigenvalues and return the amount which falls above the bound
  evals <- eigen(sigmaHatNaive, symmetric = TRUE, only.values = TRUE)$values
  k <- 0
  for (i in 1:length(evals)) {
    if (evals[i] > bd) {
      k <- k + 1
    }
  }
  return(k+2)
}

#' Normalize Marker Intensities
#' The matrix of intensities is normalized based on one of four different method
#' These methods include Percentile, zscore, arsinh and minmax
#' 
#' @param data the raw intensity scores.
#' @param markers the markers of interest.
#' @param  method the normalizaton method
#' @return normalized matrix.
#' 
#' @export
NormalizeData <- function(data, markers, method = 'none'){
  if(!(method %in% c('none', 'percentile', 'zscore', 'arsinh', 'minmax'))){
    stop('Please provide a valid normalization method')
  }
  
  data <- as.matrix(data[, markers])
  
  if(method == 'none'){
    return(data)
  }else if(method == 'percentile'){
    return(percentile_norm(data))
  }else if(method == 'minmax'){
    return(minmax_norm(data))
  }else if(method == 'arsinh'){
    return(arsinh_norm(data))
  }
  return(scale(data))
}



#' Generate Self Organizing Maps
#' A self organizing map of the marker intensities is generated using the Yasomi package
#' The grid size is determined automatically
#' 
#' @param data the marker intensities
#' @return the self organizing map.
#' 
#' @import yasomi
#' @export
GeneratePrototypes <- function(data){
  
  print(paste('=====================Now generating SOM grid======================='))
  npcs <- ComputeGridSize(data)
  k = 10
  print(paste('Optimal Grid Size is: ', npcs))
  # 
  if((npcs*npcs) < k){
    npcs <- 10
  } else{
    npcs <- npcs + 2
  }
  
  # genearte the som grid based on the computed grid size
  sg <- somgrid(xdim=npcs,ydim=npcs,topo="hex")
  
  # generate the initial prototypes using the first two pcs
  init.res <- sominit.pca(as.matrix(data), somgrid = sg)
  print(paste('=====================Now tunning and running the SOM model====================='))
  
  # generate the self organizing map
  som_model <- batchsom(as.matrix(data), sg,
                                prototypes = init.res$prototypes,
                                verbose = T
  )
  return(som_model)
}


#' Generate Self Organizing Maps
#' A self organizing map of the marker intensities is generated using the Yasomi package
#' The grid size is determined automatically
#' 
#' @param som_model the self organizing map
#' @param numClusters the number of clusters to generate
#' @return the cluster labels
#' 
#' @import coop analogue psych FCPS
#' @export
ClusterPrototypes <- function(som_model, numClusters = NULL){
  if(is.null(numClusters)){
    stop('Please provide the number of clusters')
  }
  # get the prototypes
  prototypes <- som_model$prototypes
  print(paste('=====================Now clustering prototypes====================='))
  
  pear <- cor(t(prototypes), method = 'pearson')
  cosi <- tcosine(prototypes)
  spear <- cor(t(prototypes), method = 'spearman')
    
  S <- as.matrix(fuse(cor2dist(pear),cor2dist(cosi),cor2dist(spear)))
    
    
  res <- FCPS::HierarchicalClustering(S, ClusterNo = numClusters,
                                      Type = 'AverageL')$Cls
  print(paste('=====================Now mapping clusters to data ====================='))
  
  cluster_assignment <- res[som_model$classif]
  return(cluster_assignment)
}



#' Generate Heatmap of Clusters
#' A self organizing map of the marker intensities is generated using the Yasomi package
#' The grid size is determined automatically
#' 
#' @param data the marker intensities
#' @param markers the markers of interest
#' @param heatmap the heatmap of clusters by markers
#' @return the cluster label
#' 
#' @import pheatmap ggplotify
#' @export
generateHeatmap <- function(data, markers, labels){
  features <- data
  features_heatmap <- aggregate(.~as.character(labels),
                                features[,markers],
                                mean)
  rownames(features_heatmap) <- features_heatmap[,1]
  features_heatmap <- features_heatmap[,-1]
  
  features_heatmap <- sweep(features_heatmap,2, colMeans(features_heatmap), "-")
  features_heatmap <- sweep(features_heatmap,2, apply(features_heatmap,2,sd), "/")
  features_heatmap[features_heatmap>2] <- 2
  features_heatmap[features_heatmap< -2] <- -2
  
  annotation_row = data.frame(Clusters = rownames(features_heatmap))
  
  rn <- rownames(features_heatmap)
  features_heatmap <- as.matrix(features_heatmap)
  rownames(features_heatmap) <- rn
  rownames(annotation_row) <- rownames(features_heatmap)
  
  gaps_row <- which(!duplicated(substr(rownames(features_heatmap),1,2)))[-1]-1
  
  p <- as.ggplot(pheatmap(features_heatmap, gaps_row = gaps_row, 
                                     annotation_row = annotation_row, annotation_legend = FALSE, 
                                     cluster_rows = FALSE, cluster_cols = F))
  return(p)
}