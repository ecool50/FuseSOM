#' Estimate the optimal grid size for the Self Organizing Map
#' 
#' The function finds the eigenvalues of the sample covariance matrix. 
#' It will then return the number of significant eigenvalues according to 
#' the Tracy-Widom test.
#' The function is based on the estKW function from the SC3 package
#' @param dataset The optimal grid size.
#' @return the optimal grid size.
#' 
#' 
#' @author
#'   Elijah WIllie <ewil3501@uni.sydney.edu.au>
#'   
#' @export
#' 
computeGridsize <- function(dataset) {
  
  p <- ncol(dataset)
  n <- nrow(dataset)
  
  # compute Tracy-Widom bound
  x <- scale(dataset)
  muTW <- (sqrt(n - 1) + sqrt(p))^2
  sigmaTW <- (sqrt(n - 1) + sqrt(p)) * (1/sqrt(n - 1) + 1/sqrt(p))^(1/3)
  sigmaHatNaive <- t(x)%*%x  # x left-multiplied by its transpose
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
#' @author
#'   Elijah WIllie <ewil3501@uni.sydney.edu.au>
#' @export
normalizeData <- function(data, markers, method='none'){
  if(!(method %in% c('none', 'percentile', 'zscore', 'arsinh', 'minmax'))){
    stop('Please provide a valid normalization method')
  }
  
  data <- as.matrix(data[, markers])
  
  if(method == 'none'){
    return(data)
  }else if(method == 'percentile'){
    return(.percentileNorm(data))
  }else if(method == 'minmax'){
    return(.minaxNorm(data))
  }else if(method == 'arsinh'){
    return(.arsinhNnorm(data))
  }
  return(scale(data))
}



#' Generate a Self Organizing Map and return it's prototypes
#' A self organizing map of the marker intensities is generated using the Yasomi package
#' The grid size is determined automatically
#' 
#' @param data the marker intensities
#' @return the self organizing map object.
#' @importFrom yasomi somgrid 
#' @importFrom yasomi sominit.pca
#' @importFrom yasomi batchsom
#' @export
generatePrototypes <- function(data){
  
  message('Now Generating the Self Organizing Map Grid')
  npcs <- computeGridsize(data)
  message(paste('Optimal Grid Size is: ', npcs))
  
  # set a lower bound on the grid size
  if(npcs < 5){
    npcs <- 7
  } else{
    npcs <- npcs + 2
  }
  
  # genearte the som grid based on the computed grid size
  sg <- somgrid(xdim=npcs,ydim=npcs,topo="hex")
  
  # generate the initial prototypes using the first two pcs
  init.res <- sominit.pca(as.matrix(data), somgrid = sg)
  message('Now Running the Self Organizing Map Model')
  
  # generate the self organizing map
  som_model <- batchsom(as.matrix(data), sg,
                                prototypes = init.res$prototypes,
                                verbose = T
  )
  return(som_model)
}


#' Cluster the prototypes from the Self Organizing Map
#' 
#' Clustering is done using hierarchical clustering with 
#' the average linkage function
#' 
#' @param som_model the self organizing map
#' @param numClusters the number of clusters to generate
#' @return the cluster labels
#' 
#' @importFrom coop tcosine
#' @importFrom analogue fuse
#' @importFrom psych cor2dist
#' @importFrom FCPS HierarchicalClustering
#' @export
clusterPrototypes <- function(som_model, numClusters=NULL){
  if(is.null(numClusters)){
    stop('Please provide the number of clusters')
  }
  # get the prototypes
  prototypes <- som_model$prototypes
  message('Now Clustering the Prototypes')
  
  # compute the similarity matrices
  pear <- cor(t(prototypes), method = 'pearson')
  cosi <- tcosine(prototypes)
  spear <- cor(t(prototypes), method = 'spearman')
    
  # peform multiview integration
  fianl_dist <- as.matrix(fuse(cor2dist(pear),cor2dist(cosi),cor2dist(spear)))
    
  # cluster the final 
  clusters <- HierarchicalClustering(fianl_dist, ClusterNo = numClusters,
                                      Type = 'AverageL')$Cls
  message('Now Mapping Clusters to the Original Data')
  
  cluster_assignment <- clusters[som_model$classif]
  
  message('The Prototypes have been Clustered and Mapped Successfully')
  return(cluster_assignment)
}


#' A wrapper function to run the FuseSOM algorithm.
#' This function accepts a matrix, dataframe or a SingleCellExperiment object
#' 
#' 
#' 
#' @param data a matrix, dataframe or SingleCellExperiment object. For matrices
#' and dataframes, it is assumed that markers are the columns and samples rows
#' @param assay the assay of interest if SingleCellExperiment object is used
#' @param markers the markers of interest. If this is not provided, all columns will be used
#' @param  numCluster the number of clusters to be generated from the data
#' @return A list conataining the SOM model and the cluster labels if a dataframe 
#' or matrix is provided
#' @return A Slist containing the SOM model and a ingleCellExperiment object with 
#' labels in coldata if a SingleCellExperiment object is provided
#' 
#' @author
#'   Elijah WIllie <ewil3501@uni.sydney.edu.au>
#'
#' @importFrom SingleCellExperiment colData
#' @export
#' 
runFuseSOM <- function(data, assay=NULL, markers=NULL, numClusters=NULL){
  
  if(is.null(numClusters)){
    stop("Please provide the number of clusters")
  }
  
  flag = FALSE
  
  # if we have a dataframe or a matrix
  if(class(data) %in% c('data.frame', 'matrix')){
    message("You have provided a dataset of class data.frame or matrix")
    # if no markers are given, make sure all the columns are numeric
    if(is.null(markers)){
      num_numeric  <- sum(apply(data, 2, function(x) is.numeric(x)))
      if(num_numeric != ncol(data)){
        stop("If markers of interest are not provided, make sure the data contains all numeric columns")
      }
      data_new <- data
    }else{
      # extract the markers of interest
      data_new <- data[, markers]
    }
  }
  
  # if we have a single cell experiment object
  if(class(data) == "SingleCellExperiment"){
    flag = TRUE
    message("You have provided a dataset of class SingleCellExperiment")
    # make sure an assay is provided
    if(is.null(assay)){
      stop("If a SingleCellExperiment, make sure the appropriate assay is provided as well")
    }
    data_new <- t(assay(data, assay))
    
    # again if no markers are given, make sure all the columns are numeric
    if(is.null(markers)){
      num_numeric  <- sum(apply(data_new, 2, function(x) is.numeric(x)))
      if(num_numeric != ncol(data_new)){
        stop("If markers of interest are not provided, make sure the data contains all numeric columns")
      } 
    }else{
      # extract the markers of interest
      data_new <- data_new[, markers]
    }
    
    
  }
  
  # now we can run the FuseSOM algorithm
  message("Everything looks good. Now running the FuseSOM algorithm")
  data_new <- apply(data_new, 2, function(x) as.numeric(x))
  som_model <- generatePrototypes(data_new)
  clusters <- clusterPrototypes(som_model, numClusters = numClusters)
  
  message("The FuseSOM algorithm has completed successfully")
  
  if(flag){
    coldat <- cbind(colData(data), clusters)
    colData(data) <- coldat
    return(list(model=som_model, sce=data))
  }else{
    return(list(model=som_model, clusters=clusters))
  }
  
  
}



#' A function for estimating the number of clusters using various method
#' Methods available are: Discriminant, Distance (Gap, Silhouette, Slope, Jump, and Within Cluster Distance,) and Instability
#' 
#' 
#' 
#' @param som_model the Self Organizing Map object generated by generatePrototypes()
#' @param method one of Discriminant, Distance, Stability. By default, everything is run
#' @param kseq a sequence of the number of clusters to try. Default is 2:20 clusters
#' @return A list containing the optimal number of clusters for each class of method and the
#' values over the sequence of k provided
#' 
#' 
#' 
#' @author
#'   Elijah WIllie <ewil3501@uni.sydney.edu.au>
#' @import SingleCellExperiment
#' @importFrom coop tcosine
#' @importFrom analogue fuse
#' @importFrom psych cor2dist
#' @importFrom FCPS HierarchicalClustering
#' @export
#' 
estimateNumcluster <- function(som_model, 
                      method = c('Discriminant','Distance', 'Stability'),
                      kseq = 2:20 
)
{
  # extract the prototypes from the model
  prototypes <- som_model$prototypes
  
  k_discr <- NULL
  if('Discriminant' %in% method){
    
    message('Now Computing the Number of Clusters using Discriminant Analysis')
    # the minimum number of clusters
    nmin <- 20
    
    # compute the similarity matrices
    pear <- cor(t(prototypes), method='pearson')
    cosi <- tcosine(prototypes)
    spear <- cor(t(prototypes), method='spearman')
    
    # Get the multiview integration
    final_dist <- as.matrix(fuse(cor2dist(pear),cor2dist(cosi),cor2dist(spear)))
    
    # estimate the number of clusters uisng discriminant analysis
    k_discr = .runDiscriminant(final_dist,nmin)
  }
  
  # compute the number of clusters using distance methods
  k_dist <- NULL
  if('Distance' %in% method){
    message('Now Computing The Number Of Clusters Using Distance Analysis')
    
    # get the distances
    k_dist <- .cDistance(prototypes, kseq = kseq)
    
    # compute the optimal k values using the elbow method
    k_dist$k_Gap <- .computeElbow(k_dist$Gaps[1:max(kseq)])
    k_dist$k_Slope <- .computeElbow(k_dist$Slopes[1:max(kseq)])
    k_dist$k_Jump <- .computeElbow(k_dist$Jumps[1:max(kseq)])
    k_dist$k_WCD <- .computeElbow(k_dist$WCD[1:max(kseq)])
    k_dist$k_Sil <- .computeElbow(k_dist$Silhouettes[1:max(kseq)])
  }
  
  # compute the number of clusters using the instability metric
  k_stab <- NULL
  if('Stability' %in% method){
    message('Now Computing The Number Of Clusters Using Stability Analysis')
    
    # run the stability algorithm
    k_stab <- .cStability(prototypes,kseq = kseq, nB = 10)
  }
  
  outlist <- list('Discriminant'=k_discr,
                  'Distance'=k_dist,
                  'Stability'=k_stab)
  
  return(outlist)
  
}


