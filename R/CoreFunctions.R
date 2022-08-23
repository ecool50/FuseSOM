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
#' risomMarkers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' computeGridSize(risom_dat[, risomMarkers])
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
#' risomMarkers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' normaliseData(risom_dat[, risomMarkers])
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
#' risomMarkers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' normalizeData(risom_dat[, risomMarkers])
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
#' risomMarkers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' generatePrototypes(risom_dat[, risomMarkers])
#' 
#' @export
generatePrototypes <- function(data, verbose=FALSE){
  
  message('Now Generating the Self Organizing Map Grid')
  size <- computeGridSize(data)
  message(paste('Optimal Grid Size is: ', size))
  
  # set a lower bound on the grid size
  if(size*size < 25){
    size <- 5
  } else{
    size <- size
  }
  
  # genearte the som grid based on the computed grid size
  sg <- somGrid(xDim=size,yDim=size, topo = 'hexagonal')
  
  # generate the initial prototypes using the first two pcs
  initRes <- somInitPca(as.matrix(data), somGrid = sg)
  message('Now Running the Self Organizing Map Model')
  
  # generate the self organizing map
  somModel <- batchSom(as.matrix(data), sg,
                                prototypes = initRes$prototypes,
                                verbose = verbose
  )
  return(somModel)
}


#' Cluster the prototypes from the Self Organizing Map
#' 
#' Clustering is done using hierarchical clustering with 
#' the average linkage function
#' 
#' @param somModel the self organizing map
#' @param numClusters the number of clusters to generate
#' @return the cluster labels
#' 
#' @examples  
#' data("risom_dat")
#' risomMarkers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' prototypes <- generatePrototypes(risom_dat[, risomMarkers])
#' clusters <- clusterPrototypes(prototypes, 23)
#' 
#' @importFrom coop tcosine
#' @importFrom analogue fuse
#' @importFrom psych cor2dist
#' @importFrom FCPS HierarchicalClustering
#' @export
clusterPrototypes <- function(somModel, numClusters=NULL){
  if(is.null(numClusters)){
    stop('Please provide the number of clusters')
  }
  # get the prototypes
  prototypes <- somModel$prototypes
  message('Now Clustering the Prototypes')
  
  # compute the similarity matrices
  pear <- cor(t(prototypes), method = 'pearson')
  cosi <- tcosine(prototypes)
  spear <- cor(t(prototypes), method = 'spearman')
    
  # peform multiview integration
  finalDist <- as.matrix(fuse(cor2dist(pear),cor2dist(cosi),cor2dist(spear)))
    
  # cluster the final 
  clusters <- HierarchicalClustering(finalDist, ClusterNo = numClusters,
                                      Type = 'AverageL')$Cls
  message('Now Mapping Clusters to the Original Data')
  
  clusterAssignment <- clusters[somModel$classif]
  clusterAssignment <- paste0('cluster_', clusterAssignment)
  
  message('The Prototypes have been Clustered and Mapped Successfully')
  return(clusterAssignment)
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
#' 
#' @return A list containing the SOM model and the cluster labels if a dataframe 
#' or matrix is provided
#' @return A SingleCellExperiment object with labels in coldata, and the SOM model 
#' in metadata if a SingleCellExperiment or SpatialExperiment object is provided
#' 
#' @examples 
#' data("risom_dat")
#' risomMarkers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' res <- runFuseSOM(risom_dat, markers=risomMarkers, numClusters=23)
#' 
#' @author
#'   Elijah WIllie <ewil3501@uni.sydney.edu.au>
#'
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
      numNumeric  <- sum(apply(data, 2, function(x) is.numeric(x)))
      if(numNumeric != ncol(data)){
        stop("If markers of interest are not provided, make sure the data contains all numeric columns")
      }
      dataNew <- data
    }else{
      # extract the markers of interest
      dataNew <- data[, markers]
    }
  } else if(is(data, "SingleCellExperiment") || is(data, "SpatialExperiment")){ # if we have a single cell experiment object
    flag = TRUE
    message(paste("You have provided a dataset of class", class(data)[[1]]))
    
    # make sure an assay is provided
    if(is.null(assay)){
      stop(paste("If a",class(data)[[1]],"make sure the appropriate assay is provided as well"))
    }
    dataNew <- t(assay(data, assay))
    
    # again if no markers are given, make sure all the columns are numeric
    if(is.null(markers)){
      numNumeric  <- sum(apply(dataNew, 2, function(x) is.numeric(x)))
      if(numNumeric != ncol(dataNew)){
        stop("If markers of interest are not provided, make sure the data contains all numeric columns")
      } 
    }else{
      # extract the markers of interest
      dataNew <- dataNew[, markers]
    }
    
    
    # check if the prototypes have already been generated
    if(!is.null(metadata(data)$SOM)){
      message("The prototypes have already been generated. Checking to see if the current markes were used to generate the prototypes")
      # check to see if the same markers were used to generate the previous SOM
      oldSom <- metadata(data)$SOM
      oldMarkers <- colnames(oldSom$prototypes)
      
      # if no markers were provided
      if(is.null(markers)){
        if(identical(oldMarkers,colnames(dataNew))){
          message('The same markers were used to generate the prototypes. Will proceed to clustering the prototypes')
          clusters <- clusterPrototypes(oldSom, numClusters = numClusters)
          colData(data)$clusters <- clusters
          # update the cluster name
          names(colData(data))[names(colData(data)) == "clusters"] <- clusterCol
          return(data)
        } else{
          message('Different markers were used to generate the prototypes. Will regenerate the self organizing map and then cluster prototypes')
        }
        
        # markers were provided
      } else{
        if(identical(oldMarkers,markers)){
          message('The same markers were used to generate the prototypes. Will proceed to clustering the prototypes')
          clusters <- clusterPrototypes(oldSom, numClusters = numClusters)
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
  dataNew <- apply(dataNew, 2, function(x) as.numeric(x))
  somModel <- generatePrototypes(dataNew, verbose=verbose)
  clusters <- clusterPrototypes(somModel, numClusters = numClusters)
  
  message("The FuseSOM algorithm has completed successfully")
  
  if(flag){
    colData(data)$clusters <- clusters
    # update the cluster name
    names(colData(data))[names(colData(data)) == "clusters"] <- clusterCol
    metadata(data) <- append(metadata(data), list(SOM=somModel))
    return(data)
  }else{
    return(list(model=somModel, clusters=clusters))
  }
  
  
}



#' A function for estimating the number of clusters using various method
#' Methods available are: Discriminant, Distance 
#' (Gap, Silhouette, Slope, Jump, and Within Cluster Distance,) and Instability
#' 
#' @param data the Self Organizing Map object generated by generatePrototypes(), or
#' an object of class SingleCellExperiment or SpatialExperiment
#' @param method one of Discriminant, Distance, Stability. By default, everything is run
#' @param kSeq a sequence of the number of clusters to try. Default is 2:20 clusters
#' @return A list containing the cluster estimations if a dataframe 
#' or matrix is provided
#' @return A SingleCellExperiment object with the cluster estimation in the metadata 
#' if a SingleCellExperiment or SpatialExperiment object is provided
#' 
#' @examples  
#' data("risom_dat")
#' risomMarkers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' res <- runFuseSOM(risom_dat, markers=risomMarkers, numClusters=23)
#' res.est.k <- estimateNumCluster(res$model, kSeq=2:25)
#' 
#' @author
#'   Elijah WIllie <ewil3501@uni.sydney.edu.au>
#' @importFrom coop tcosine
#' @importFrom analogue fuse
#' @importFrom psych cor2dist
#' @importFrom FCPS HierarchicalClustering
#' @export
#' 
estimateNumCluster <- function(data, 
                      method = c('Discriminant','Distance'),
                      kSeq = 2:20 
)
{
  
  if(1 %in% kSeq) stop('Please select a k sequence starting with 2: {2,3,...K}!')
  
  flag = FALSE
  
  # if we have a single cell experiment or spatial experiment object
  if(is(data, "SingleCellExperiment") || is(data, "SpatialExperiment")){
    flag = TRUE
    message(paste("You have provided a dataset of class:", class(data)[[1]]))
    somModel <- data@metadata$SOM
    
  } else{ # if we just have the som model
    somModel <- data
  }
  
  # extract the prototypes from the model
  prototypes <- somModel$prototypes
  
  kDiscr <- NULL
  if('Discriminant' %in% method){
    
    message('Now Computing the Number of Clusters using Discriminant Analysis')
    # the minimum number of clusters
    nMin <- 10
    
    # compute the similarity matrices
    pear <- cor(t(prototypes), method='pearson')
    cosi <- tcosine(prototypes)
    spear <- cor(t(prototypes), method='spearman')
    
    # Get the multiview integration
    finalDist <- as.matrix(fuse(cor2dist(pear),cor2dist(cosi),cor2dist(spear)))
    
    # estimate the number of clusters uisng discriminant analysis
    kDiscr = .runDiscriminant(finalDist,nMin)
  }
  
  # compute the number of clusters using distance methods
  kDist <- NULL
  if('Distance' %in% method){
    message('Now Computing The Number Of Clusters Using Distance Analysis')
    
    # get the distances
    kDist <- .cDistance(prototypes, kSeq = kSeq)
    
    # compute the optimal k values using the elbow method
    kDist$kGap <- .computeElbow(kDist$Gaps)
    kDist$kSlope <- .computeElbow(kDist$Slopes)
    kDist$kJump <- .computeElbow(kDist$Jumps)
    kDist$kWCD <- .computeElbow(kDist$WCD)
    kDist$kSil <- .computeElbow(kDist$Silhouettes)
  }
  
  outList <- list('Discriminant'=kDiscr,
                  'Distance'=kDist)
  
  
  if(flag){
    metadata(data) <- append(metadata(data), list(clusterEstimation=outList))
    return(data)
  }else{
    return(outList)
  }
  
}



#' A function generating the elbow plot for the optimal number of clusters 
#' returned by the estimateNumcluster() function
#' Methods available are: 
#' Gap, Silhouette, Slope, Jump, and Within Cluster Distance(WCD)
#' 
#' @param kEst the estimation list returned by estimateNumcluster()
#' @param method one of 'jump', 'slope', 'wcd', 'gap', or 'silhouette'
#' @return an elbow plot object where the optimal number of clusters is marked
#' 
#' @examples 
#' data("risom_dat")
#' risomMarkers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' res <- runFuseSOM(risom_dat, markers=risomMarkers, numClusters=23)
#' resEstK <- estimateNumCluster(res$model, kSeq=2:25)
#' p <- optiPlot(resEstK, method='jump')
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
    kEst <- data@metadata$clusterEstimation
    
  } else { # if we just have the som model
    kEst <- data
  }
  
  # extract the relevant information for the method provided
  if(method == 'jump'){
    method <- stringr::str_to_title(method)
    jumps <- kEst$Distance$Jumps
    plotData <- data.frame(Clusters=1:length(jumps), Jump=jumps)
    kOpti <- kEst$Distance$kJump
  } else if(method == 'slope'){
    method <- stringr::str_to_title(method)
    slopes <- kEst$Distance$Slopes
    plotData <- data.frame(Clusters=1:length(slopes), Slope=slopes)
    kOpti <- kEst$Distance$kSlope
  } else if(method == 'wcd'){
    method <- stringr::str_to_upper(method)
    wcds <- kEst$Distance$WCD
    plotData <- data.frame(Clusters=1:length(wcds), WCD=wcds)
    kOpti <- kEst$Distance$kWCD
  } else if(method == 'gap'){
    method <- stringr::str_to_title(method)
    gaps <- kEst$Distance$Gaps
    plotData <- data.frame(Clusters=1:length(gaps), Gap=gaps)
    kOpti <- kEst$Distance$kGap
  } else{
    method <- stringr::str_to_title(method)
    silhouettes <- kEst$Distance$Silhouettes
    plotData <- data.frame(Clusters=1:length(silhouettes), 
                           Silhouette=silhouettes)
    kOpti <- kEst$Distance$kSil
  } 
  
  # plot the data
  pOpti <- ggpubr::ggline(plotData, x = "Clusters", y = method, group = 1, color = "steelblue") +
    geom_vline(xintercept = kOpti, linetype = 2, color = 'steelblue') +
    labs(y = paste(method,"statistic (k)"), x = "Number of clusters k",
         title = "Optimal number of clusters using the elbow method") +
    theme(
      plot.title = element_text(color="Black", size=14, face="bold", hjust = 0.5),
      axis.title.x = element_text(color="Black", size=14, face="bold"),
      axis.title.y = element_text(color="Black", size=14, face="bold")
    )
  
  return(pOpti)
  
  
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
#' risomMarkers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' res <- runFuseSOM(risom_dat, markers=risomMarkers, numClusters=23)
#' p.heat <- markerHeatmap(risom_dat, risomMarkers, clusters=res$clusters)
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
  
  if(!(is(data, "matrix") || is(data, "data.frame"))){
    stop("Make sure you are passing in a dataframe or matrix")
  }
  
  
  # again if no markers are given, make sure all the columns are numeric
  if(is.null(markers)){
    numNumeric  <- sum(apply(data, 2, function(x) is.numeric(x)))
    if(numNumeric != ncol(data)){
      stop("If markers of interest are not provided, make sure the data contains all numeric columns")
    } 
    message("No markers provided, will be using all columns as markers")
    markers <- colnames(data)
  }
  
  features <- data[, markers]
  # do some wrangling to get it in the proper format
  featuresHeatmap <- aggregate(.~as.character(clusters),
                                features[,markers],
                                mean)
  rownames(featuresHeatmap) <- featuresHeatmap[,1]
  featuresHeatmap <- featuresHeatmap[,-1]
  
  # compute the marker expression
  featuresHeatmap <- sweep(featuresHeatmap,2, colMeans(featuresHeatmap), "-")
  featuresHeatmap <- sweep(featuresHeatmap,2, apply(featuresHeatmap,2,sd), "/")
  featuresHeatmap[featuresHeatmap>threshold] <- threshold
  featuresHeatmap[featuresHeatmap< -threshold] <- -threshold
  
  # compute the heatmap annotations
  annotationRow = data.frame(Clusters = rownames(featuresHeatmap))
  
  rn <- rownames(featuresHeatmap)
  featuresHeatmap <- as.matrix(featuresHeatmap)
  rownames(featuresHeatmap) <- rn
  rownames(annotationRow) <- rownames(featuresHeatmap)
  
  gapRows <- which(!duplicated(substr(rownames(featuresHeatmap),1,2)))[-1]-1
  
  pHeat <- ggplotify::as.ggplot(pheatmap(featuresHeatmap, gaps_row = gapRows, 
                                     annotation_row = annotationRow, annotation_legend = FALSE, 
                                     cluster_cols = clusterMarkers, cluster_rows = FALSE, 
                                     fontsize = fontSize))
  return(pHeat)
}


#' A function for plotting the self organizing map
#' function was obtained from https://rdrr.io/rforge/yasomi/ with some minor modifications
#' 
#' @param som the self organizing map 
#' @param border a value for the plot border
#' @param withGrid option to add the som grid
#' @return None
#' @examples 
#' data("risom_dat")
#' risomMarkers <- c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD')
#' prototypes <- generatePrototypes(risom_dat[, risomMarkers])
#' plotSOM(prototypes)
#' 
#' @author
#'   Elijah WIllie <ewil3501@uni.sydney.edu.au>
#' @export
#' 
plotSOM <- function(som,border=NA,withGrid=TRUE,...) {
  args <- list(...)
  if(is.null(args$col)) {
    args$col <- "red"
  }
  if(withGrid) {
    add <- !is.null(args$add) && args$add
    plot(som$somGrid,add=add)
    args$add <- TRUE
  }
  args$border <- border
  sizes <- table(factor(som$classif,levels=1:nrow(som$prototypes)))
  sizes <- sizes/max(sizes)
  do.call("plot",c(list(x=som$somGrid,size=sizes),args))
}