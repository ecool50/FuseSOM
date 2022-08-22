# A function to initialize a somgrid using the first two pcs
# function was obtained from https://rdrr.io/rforge/yasomi/ with some major modifications

somInitPca.default <- function(data,somGrid,weights,withPrinComp=FALSE,...) {
  ### FIXME: data weights support
  ## we don't scale the data
  if(missing(weights) || is.null(weights)) {
    if(withPrinComp) {
      dataPca <- princomp(covmat=cov.wt(data))
    } else {
      dataPca <- prcomp(data)
    }
  } else {
    dataPca <- princomp(covmat=cov.wt(data,weights))
  }
  ## the more detailed axis is assigned to the first eigenvector
  if (somGrid$xDim>=somGrid$yDim) {
    xEv <- 1
    yEv <- 2
  } else {
    xEv <- 2
    yEv <- 1
  }

  xSpan <- somGrid$xDim - 1
  if(somGrid$yDim>1) {
    xSpan <- xSpan+0.5
  }
  x <- seq(from=-2*dataPca$sdev[xEv],by=4*dataPca$sdev[xEv]/xSpan,length.out=somGrid$xDim)
  y <- seq(from=2*dataPca$sdev[yEv],to=-2*dataPca$sdev[yEv],length.out=somGrid$yDim)
  base <- as.matrix(expand.grid(x = x, y = y))
  ## correction for hexagonal grids
  base[,1] <- base[,1]+rep(c(0,2*dataPca$sdev[xEv]/xSpan),each=somGrid$xDim,length.out=nrow(base))
  ## map back the grid to the original space
  if(inherits(dataPca,"prcomp")) {
    mapped <- tcrossprod(base,dataPca$rotation[,c(xEv,yEv)])
  } else {
    mapped <- tcrossprod(base,data.pca$loadings[,c(xEv,yEv)])
  }
  ## decentering
  prototypes <- sweep(mapped,2,dataPca$center,"+")
  list(prototypes=prototypes,dataPca=dataPca)
}

# A function to compute the best matching unit 
# function was obtained from https://rdrr.io/rforge/yasomi/ with some major modifications
bmu <- function(prototypes,data,weights) {
  if(ncol(prototypes)!=ncol(data)) {
    stop("'prototypes' and 'data' have different dimensions")
  }
  if(missing(weights)) {
    weights <- rep(1,nrow(data))
  } else if(length(weights)!=nrow(data)) {
    stop("'weights' and 'data' have different dimensions")
  }
  result <- .C("bmu",
               as.double(prototypes),
               as.integer(nrow(prototypes)),
               as.double(data),
               as.integer(nrow(data)),
               as.integer(ncol(prototypes)),
               as.double(weights),
               clusters=integer(nrow(data)),
               error=as.double(0),
               PACKAGE='FuseSOM')
  list(clusters=result$clusters+1,error=result$error)
}

# A function to compute second best matching unit
# function was obtained from https://rdrr.io/rforge/yasomi/ with some major modifications
secondBmu <- function(prototypes,data) {
  if(ncol(prototypes)!=ncol(data)) {
    stop("'prototypes' and 'data' have different dimensions")
  }
  matrix(.C("second_bmu",
            as.double(prototypes),
            as.integer(nrow(prototypes)),
            as.double(data),
            as.integer(nrow(data)),
            as.integer(ncol(prototypes)),
            clusters=integer(2*nrow(data)),
            PACKAGE='FuseSOM'
            )$clusters,ncol=2)
}

# A function to compute the self organinzing map
# function was obtained from https://rdrr.io/rforge/yasomi/ with some major modifications
batchSom.default <- function(data,somGrid,init="pca",prototypes,
                             weights,
                             mode="continuous",
                             minRadius, maxRadius, steps,
                             decrease="power", maxTter,
                             kernel="gaussian", normalised=0,
                             assignment="single",
                             cut=1e-07,
                             verbose=FALSE,keepdata=TRUE,...) {

  theCall <- match.call()
  if(verbose) {
    print(theCall)
  }
  theCall[[1]] <- batchSomControl
  control <- eval(theCall,envir = parent.frame())
  control$assignmentInt <- 0
  control$kernelInt <- 0
  control$normalised <- 0
  if(!missing(weights)) {
    if(length(weights)!=nrow(data)) {
      stop("'weights' and 'data' have different dimensions")
    }
  } else {
    ## keep weights to NULL for now to avoid princomp initialization
    weights <- NULL
  }

  if(ncol(prototypes)!=ncol(data)) {
    stop("'prototypes' and 'data' have different dimensions")
  }
  if(nrow(prototypes)!=somGrid$size) {
    stop("'prototypes' and 'somGrid' are not compatible")
  }
  
  if(is.null(weights)) {
    weights <- rep(1,nrow(data))
  }
  ## distances?
  if(is.null(somGrid$dist)) {
    somGrid$dist <- as.matrix(dist(somGrid$pts,method="euclidean"),diag=0)
  }
  pre <- batchSomLowLevelContinuous(somGrid,data,
                                     prototypes,weights,control,verbose)
  pre$control <- control
  if(keepdata) {
    pre$data  <- data
    pre$weights <- weights
  }
  return(pre)
}


# A wrapper function to call the C function that computes the self organizing map
# function was obtained from https://rdrr.io/rforge/yasomi/ with some major modifications
batchSomLowLevelContinuous <- function(somGrid,data,prototypes,weights,
                                        control,verbose) {
  result <- .C("batch_som_optim_continuous",
               proto=as.double(prototypes),
               as.integer(somGrid$size),
               as.double(data),
               as.integer(nrow(data)),
               as.integer(ncol(data)),
               as.double(weights),
               as.integer(control$assignmentInt),
               as.double(somGrid$dist),
               as.double(control$radii),
               as.integer(length(control$radii)),
               as.integer(control$kernelInt),
               as.integer(control$normalised),
               as.double(control$cut),
               as.integer(verbose),
               clusters=integer(nrow(data)),
               errors=as.double(rep(-1,1+length(control$radii))),
               PACKAGE='FuseSOM')
  prototypes <- matrix(result$proto,ncol=ncol(prototypes),
                       dimnames=list(NULL,dimnames(data)[[2]]))
  res <- list(somGrid=somGrid,
              prototypes=prototypes,
              classif=result$cluster+1,
              errors=result$errors[result$errors>=0])
  class(res) <- c("somnum","som")
  return(res)
}


