#' Estimate the optimal grid size for the Self Organizing Map
#' 
#' The function finds the eigenvalues of the sample covariance matrix. 
#' It will then return the number of significant eigenvalues according to 
#' the Tracy-Widom test.
#' The function is based on the estKW function from the SC3 package
#' @param dataset The optimal grid size.
#' @return the optimal grid size.
#' 
#' @examples
#' data("risom_dat")
#' risom_markers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' computeGridSize(risom_dat[, risom_markers])
#' 
#' 
#' @author
#'   Elijah WIllie <ewil3501@uni.sydney.edu.au>
#'   
#' @export
#' 
computeGridSize <- function(dataset) {
  
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

#' normalise Marker Intensities
#' The matrix of intensities is normalised based on one of four different method
#' These methods include Percentile, zscore, arsinh and minmax
#' 
#' @param data the raw intensity scores.
#' @param markers the markers of interest.
#' @param  method the normalizaton method
#' @param cofactor the cofactor for arsinh normalisation
#' @return normalised matrix.
#' 
#' @examples
#' data("risom_dat")
#' risom_markers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' normaliseData(risom_dat[, risom_markers])
#' 
#' @author
#'   Elijah WIllie <ewil3501@uni.sydney.edu.au>
#' @export
normaliseData <- function(data, markers, method='none', cofactor=5){
  if(!(method %in% c('none', 'percentile', 'zscore', 'arsinh', 'minmax'))){
    stop('Please provide a valid normalization method')
  }
  
  data <- as.matrix(data[, markers])
  
  if(method == 'none'){
    return(data)
  }else if(method == 'percentile'){
    return(.percentileNorm(data))
  }else if(method == 'minmax'){
    return(.minmaxNorm(data))
  }else if(method == 'arsinh'){
    return(.arsinhNnorm(data, cofactor))
  }
  return(scale(data))
}


#' normalize Marker Intensities
#' The matrix of intensities is normalised based on one of four different method
#' These methods include Percentile, zscore, arsinh and minmax
#' 
#' @param data the raw intensity scores.
#' @param markers the markers of interest.
#' @param method the normalizaton method
#' @param cofactor the cofactor for arsinh normalization
#' @return normalised matrix.
#' 
#' @examples
#' data("risom_dat")
#' risom_markers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' normalizeData(risom_dat[, risom_markers])
#' 
#' @author
#'   Elijah WIllie <ewil3501@uni.sydney.edu.au>
normalizeData <- function(data, markers, method='none', cofactor=5){
  if(!(method %in% c('none', 'percentile', 'zscore', 'arsinh', 'minmax'))){
    stop('Please provide a valid normalization method')
  }
  
  data <- as.matrix(data[, markers])
  
  if(method == 'none'){
    return(data)
  }else if(method == 'percentile'){
    return(.percentileNorm(data))
  }else if(method == 'minmax'){
    return(.minmaxNorm(data))
  }else if(method == 'arsinh'){
    return(.arsinhNnorm(data,cofactor))
  }
  return(scale(data))
}


#' Generate a Self Organizing Map and return it's prototypes
#' A self organizing map of the marker intensities is generated using the Yasomi package
#' The grid size is determined automatically
#' 
#' @param data the marker intensities
#' @param verbose should the progress be printed out
#' @return the self organizing map object
#' 
#' @examples 
#' data("risom_dat")
#' risom_markers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' generatePrototypes(risom_dat[, risom_markers])
#' 
#' @importFrom yasomi somgrid 
#' @importFrom yasomi sominit.pca
#' @importFrom yasomi batchsom
#' @export
generatePrototypes <- function(data, verbose=FALSE){
  
  message('Now Generating the Self Organizing Map Grid')
  size <- computeGridSize(data)
  message(paste('Optimal Grid Size is: ', size))
  
  # set a lower bound on the grid size
  if(size*size < 25){
    size <- 5
  } else{
    size <- size + 2
  }
  
  # genearte the som grid based on the computed grid size
  sg <- somgrid(xdim=size,ydim=size,topo="hex")
  
  # generate the initial prototypes using the first two pcs
  init.res <- sominit.pca(as.matrix(data), somgrid = sg)
  message('Now Running the Self Organizing Map Model')
  
  # generate the self organizing map
  som_model <- batchsom(as.matrix(data), sg,
                                prototypes = init.res$prototypes,
                                verbose = verbose
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
#' @examples  
#' data("risom_dat")
#' risom_markers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' prototypes <- generatePrototypes(risom_dat[, risom_markers])
#' clusters <- clusterPrototypes(prototypes, 23)
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
  cluster_assignment <- paste0('cluster_', cluster_assignment)
  
  message('The Prototypes have been Clustered and Mapped Successfully')
  return(cluster_assignment)
}


#' A wrapper function to run the FuseSOM algorithm.
#' This function accepts a matrix, dataframe or a SingleCellExperiment object
#' 
#' @param data a matrix, dataframe, SingleCellExperiment or SpatialExperiment object. 
#' For matrices and dataframes, it is assumed that markers are the columns and samples rows
#' @param assay the assay of interest if SingleCellExperiment object is used
#' @param markers the markers of interest. If this is not provided, all columns will be used
#' @param numCluster the number of clusters to be generated from the data
#' @param clusterCol the name of the column to store the clusters in
#' @param verbose should the generation of the Self Organising Map be printed
#' @return A list containing the SOM model and the cluster labels if a dataframe 
#' or matrix is provided
#' @return A SingleCellExperiment object with labels in coldata, and the SOM model 
#' in metadata if a SingleCellExperiment or SpatialExperiment object is provided
#' 
#' @examples 
#' data("risom_dat")
#' risom_markers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' res <- runFuseSOM(risom_dat, markers=risom_markers, numClusters=23)
#' 
#' @author
#'   Elijah WIllie <ewil3501@uni.sydney.edu.au>
#'
#' @importFrom SingleCellExperiment colData
#' @export
#' 
runFuseSOM <- function(data,markers=NULL, numClusters=NULL, assay=NULL,
                       clusterCol='clusters', verbose=FALSE){
  
  if(is.null(numClusters)){
    stop("Please provide the number of clusters")
  }
  
  flag = FALSE
  
  # if we have a dataframe or a matrix
  if(is(data, "data.frame") || is(data, "matrix")){
    message(paste("You have provided a dataset of class", class(data)[[1]]))
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
  } else if(is(data, "SingleCellExperiment") || is(data, "SpatialExperiment")){ # if we have a single cell experiment object
    flag = TRUE
    message(paste("You have provided a dataset of class", class(data)[[1]]))
    
    # make sure an assay is provided
    if(is.null(assay)){
      stop(paste("If a",class(data)[[1]],"make sure the appropriate assay is provided as well"))
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
    
    
    # check if the prototypes have already been generated
    if(!is.null(metadata(data)$SOM)){
      message("The prototypes have already been generated. Checking to see if the current markes were used to generate the prototypes")
      # check to see if the same markers were used to generate the previous SOM
      old_som <- metadata(data)$SOM
      old_markers <- colnames(old_som$prototypes)
      
      # if no markers were provided
      if(is.null(markers)){
        if(identical(old_markers,colnames(data_new))){
          message('The same markers were used to generate the prototypes. Will proceed to clustering the prototypes')
          clusters <- clusterPrototypes(old_som, numClusters = numClusters)
          colData(data)$clusters <- clusters
          # update the cluster name
          names(colData(data))[names(colData(data)) == "clusters"] <- clusterCol
          return(data)
        } else{
          message('Different markers were used to generate the prototypes. Will regenerate the self organizing map and then cluster prototypes')
        }
        
        # markers were provided
      } else{
        if(identical(old_markers,markers)){
          message('The same markers were used to generate the prototypes. Will proceed to clustering the prototypes')
          clusters <- clusterPrototypes(old_som, numClusters = numClusters)
          colData(data)$clusters <- clusters
          # update the cluster name
          names(colData(data))[names(colData(data)) == "clusters"] <- clusterCol
          return(data)
        } else{
          message('Different markers were used to generate the prototypes. Will regenerate the self organizing map and then cluster prototypes')
        }
        
      }
        
     
    }
  }else{
    stop("Please provide a dataset of type SingleCellExperiment, SpatialExperiment,
         data.frame or matrix")
  }
  
  # now we can run the FuseSOM algorithm
  message("Everything looks good. Now running the FuseSOM algorithm")
  data_new <- apply(data_new, 2, function(x) as.numeric(x))
  som_model <- generatePrototypes(data_new, verbose=verbose)
  clusters <- clusterPrototypes(som_model, numClusters = numClusters)
  
  message("The FuseSOM algorithm has completed successfully")
  
  if(flag){
    colData(data)$clusters <- clusters
    # update the cluster name
    names(colData(data))[names(colData(data)) == "clusters"] <- clusterCol
    metadata(data) <- append(metadata(data), list(SOM=som_model))
    return(data)
  }else{
    return(list(model=som_model, clusters=clusters))
  }
  
  
}



#' A function for estimating the number of clusters using various method
#' Methods available are: Discriminant, Distance 
#' (Gap, Silhouette, Slope, Jump, and Within Cluster Distance,) and Instability
#' 
#' @param data the Self Organizing Map object generated by generatePrototypes(), or
#' an object of class SingleCellExperiment or SpatialExperiment
#' @param method one of Discriminant, Distance, Stability. By default, everything is run
#' @param kseq a sequence of the number of clusters to try. Default is 2:20 clusters
#' @return A list containing the cluster estimations if a dataframe 
#' or matrix is provided
#' @return A SingleCellExperiment object with the cluster estimation in the metadata 
#' if a SingleCellExperiment or SpatialExperiment object is provided
#' 
#' @examples  
#' data("risom_dat")
#' risom_markers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' res <- runFuseSOM(risom_dat, markers=risom_markers, numClusters=23)
#' res.est.k <- estimateNumCluster(res$model, kseq=2:25)
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
estimateNumCluster <- function(data, 
                      method = c('Discriminant','Distance', 'Stability'),
                      kseq = 2:20 
)
{
  
  if(1 %in% kseq) stop('Please select a k sequence starting with 2: {2,3,...K}!')
  
  flag = FALSE
  
  # if we have a single cell experiment or spatial experiment object
  if(is(data, "SingleCellExperiment") || is(data, "SpatialExperiment")){
    flag = TRUE
    message(paste("You have provided a dataset of class:", class(data)[[1]]))
    som_model <- data@metadata$SOM
    
  } else{ # if we just have the som model
    som_model <- data
  }
  
  # extract the prototypes from the model
  prototypes <- som_model$prototypes
  
  k_discr <- NULL
  if('Discriminant' %in% method){
    
    message('Now Computing the Number of Clusters using Discriminant Analysis')
    # the minimum number of clusters
    nmin <- 10
    
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
    k_dist$k_Gap <- .computeElbow(k_dist$Gaps)
    k_dist$k_Slope <- .computeElbow(k_dist$Slopes)
    k_dist$k_Jump <- .computeElbow(k_dist$Jumps)
    k_dist$k_WCD <- .computeElbow(k_dist$WCD)
    k_dist$k_Sil <- .computeElbow(k_dist$Silhouettes)
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
  
  
  if(flag){
    metadata(data) <- append(metadata(data), list(clusterEstimation=outlist))
    return(data)
  }else{
    return(outlist)
  }
  
}



#' A function generating the elbow plot for the optimal number of clusters 
#' returned by the estimateNumcluster() function
#' Methods available are: 
#' Gap, Silhouette, Slope, Jump, and Within Cluster Distance(WCD)
#' 
#' @param k_est the estimation list returned by estimateNumcluster()
#' @param method one of 'jump', 'slope', 'wcd', 'gap', or 'silhouette'
#' @return an elbow plot object where the optimal number of clusters is marked
#' 
#' @examples 
#' data("risom_dat")
#' risom_markers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' res <- runFuseSOM(risom_dat, markers=risom_markers, numClusters=23)
#' res.est.k <- estimateNumCluster(res$model, kseq=2:25)
#' p <- optiPlot(res.est.k, method='jump')
#' 
#' @author
#'   Elijah WIllie <ewil3501@uni.sydney.edu.au>
#' @import ggplot2
#' @import ggpubr
#' @import stringr
#' @export
#' 

optiPlot <- function(data, method='jump'){
  
  # make sure a valid method is provided
  if(!(method %in% c('jump', 'slope', 'wcd', 'gap', 'silhouette'))){
    stop('Please provide a valid method')
  }
  
  # if we have a single cell experiment object
  if(is(data, "SingleCellExperiment") || is(data, "SpatialExperiment")){
    message(paste("You have provided a dataset of class:", class(data)[[1]]))
    k_est <- data@metadata$clusterEstimation
    
  } else { # if we just have the som model
    k_est <- data
  }
  
  # extract the relevant information for the method provided
  if(method == 'jump'){
    method <- stringr::str_to_title(method)
    jumps <- k_est$Distance$Jumps
    plot_data <- data.frame(Clusters=1:length(jumps), Jump=jumps)
    k_opti <- k_est$Distance$k_Jump
  } else if(method == 'slope'){
    method <- stringr::str_to_title(method)
    slopes <- k_est$Distance$Slopes
    plot_data <- data.frame(Clusters=1:length(slopes), Slope=slopes)
    k_opti <- k_est$Distance$k_Slope
  } else if(method == 'wcd'){
    method <- stringr::str_to_upper(method)
    wcds <- k_est$Distance$WCD
    plot_data <- data.frame(Clusters=1:length(wcds), WCD=wcds)
    k_opti <- k_est$Distance$k_WCD
  } else if(method == 'gap'){
    method <- stringr::str_to_title(method)
    gaps <- k_est$Distance$Gaps
    plot_data <- data.frame(Clusters=1:length(gaps), Gap=gaps)
    k_opti <- k_est$Distance$k_Gap
  } else{
    method <- stringr::str_to_title(method)
    silhouettes <- k_est$Distance$Silhouettes
    plot_data <- data.frame(Clusters=1:length(silhouettes), 
                           Silhouette=silhouettes)
    k_opti <- k_est$Distance$k_Sil
  } 
  
  # plot the data
  p <- ggpubr::ggline(plot_data, x = "Clusters", y = method, group = 1, color = "steelblue") +
    geom_vline(xintercept = k_opti, linetype = 2, color = 'steelblue') +
    labs(y = paste(method,"statistic (k)"), x = "Number of clusters k",
         title = "Optimal number of clusters using the elbow method") +
    theme(
      plot.title = element_text(color="Black", size=14, face="bold", hjust = 0.5),
      axis.title.x = element_text(color="Black", size=14, face="bold"),
      axis.title.y = element_text(color="Black", size=14, face="bold")
    )
  
  return(p)
  
  
}



#' A function for generating a heat map of marker expression across clusters 
#' 
#' @param data a matrix or dataframe where the rows are samples and columns are
#' markers
#' @param markers a list of markers of interest. If not provided, all columns
#' will be used
#' @param clusters a vector of cluster labels
#' @param threshold the value to threshold the marker expression at
#' @param clusterMarkers should the rows(markers) of the heatmap be clustered
#' @param fontSize the size of the text on the heatmap
#' @return a heatmap with the markers in the rows and clusters in the columns
#' 
#' @examples 
#' data("risom_dat")
#' risom_markers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' res <- runFuseSOM(risom_dat, markers=risom_markers, numClusters=23)
#' p.heat <- markerHeatmap(risom_dat, risom_markers, clusters=res$clusters)
#' 
#' @author
#'   Elijah WIllie <ewil3501@uni.sydney.edu.au>
#' @importFrom ggplotify as.ggplot
#' @importFrom pheatmap pheatmap
#' @export
#' 
markerHeatmap <- function(data, markers=NULL, clusters=NULL, threshold=2,
                          clusterMarkers=FALSE, fontSize=14){
  
  # do some house keeping
  if(is.null(clusters)){
    stop("Please provide a vector of cluster labels")
  }
  
  if(!((class(data) == 'matrix') || (class(data) == 'data.frame'))){
    stop("Make sure you are passing in a dataframe or matrix")
  }
  
  
  # again if no markers are given, make sure all the columns are numeric
  if(is.null(markers)){
    num_numeric  <- sum(apply(data, 2, function(x) is.numeric(x)))
    if(num_numeric != ncol(data)){
      stop("If markers of interest are not provided, make sure the data contains all numeric columns")
    } 
    message("No markers provided, will be using all columns as markers")
    markers <- colnames(data)
  }
  
  features <- data[, markers]
  # do some wrangling to get it in the proper format
  features_heatmap <- aggregate(.~as.character(clusters),
                                features[,markers],
                                mean)
  rownames(features_heatmap) <- features_heatmap[,1]
  features_heatmap <- features_heatmap[,-1]
  
  # compute the marker expression
  features_heatmap <- sweep(features_heatmap,2, colMeans(features_heatmap), "-")
  features_heatmap <- sweep(features_heatmap,2, apply(features_heatmap,2,sd), "/")
  features_heatmap[features_heatmap>threshold] <- threshold
  features_heatmap[features_heatmap< -threshold] <- -threshold
  
  # compute the heatmap annotations
  annotation_row = data.frame(Clusters = rownames(features_heatmap))
  
  rn <- rownames(features_heatmap)
  features_heatmap <- as.matrix(features_heatmap)
  rownames(features_heatmap) <- rn
  rownames(annotation_row) <- rownames(features_heatmap)
  
  gaps_row <- which(!duplicated(substr(rownames(features_heatmap),1,2)))[-1]-1
  
  p <- ggplotify::as.ggplot(pheatmap(features_heatmap, gaps_row = gaps_row, 
                                     annotation_row = annotation_row, annotation_legend = FALSE, 
                                     cluster_cols = clusterMarkers, cluster_rows = FALSE, 
                                     fontsize = fontSize))
  return(p)
}
