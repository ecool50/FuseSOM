# these functions were obtained from https://rdrr.io/rforge/yasomi/ with some major modifications

somInitPca.default <- function(data,somGrid,weights,with.princomp=FALSE,...) {
  ### FIXME: data weights support
  ## we don't scale the data
  if(missing(weights) || is.null(weights)) {
    if(with.princomp) {
      dataPca <- princomp(covmat=cov.wt(data))
    } else {
      dataPca <- prcomp(data)
    }
  } else {
    dataPca <- princomp(covmat=cov.wt(data,weights))
  }
  ## the more detailled axis is assigned to the first eigenvector
  if (somGrid$xDim>=somGrid$yDim) {
    xEv <- 1
    yEv <- 2
  } else {
    xEv <- 2
    yEv <- 1
  }
  if(somGrid$topo=="hexagonal") {
    xSpan <- somGrid$xDim - 1
    if(somGrid$yDim>1) {
      xSpan <- xSpan+0.5
    }
    x <- seq(from=-2*dataPca$sdev[xEv],by=4*dataPca$sdev[xEv]/xSpan,length.out=somGrid$xDim)
  } else {
    x <- seq(from=-2*dataPca$sdev[xEv],to=2*dataPca$sdev[xEv],length.out=somGrid$xDim)
  }
  y <- seq(from=2*dataPca$sdev[yEv],to=-2*dataPca$sdev[yEv],length.out=somGrid$yDim)
  base <- as.matrix(expand.grid(x = x, y = y))
  ## correction for hexagonal grids
  if(somGrid$topo=="hexagonal") {
    base[,1] <- base[,1]+rep(c(0,2*dataPca$sdev[xEv]/xSpan),each=somGrid$xDim,length.out=nrow(base))
  }
  ## map back the grid to the original space
  if(inherits(dataPca,"prcomp")) {
    mapped <- tcrossprod(base,dataPca$rotation[,c(xEv,yEv)])
  } else {
    mapped <- tcrossprod(base,dataPca$loadings[,c(xEv,yEv)])
  }
  ## decentering
  prototypes <- sweep(mapped,2,dataPca$center,"+")
  list(prototypes=prototypes,dataPca=dataPca)
}

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
               PACKAGE="FuseSOM")
  list(clusters=result$clusters+1,error=result$error)
}

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
            PACKAGE="FuseSOM")$clusters,ncol=2)
}

bmuHeskes <- function(prototypes,data,nv,weights) {
  if(ncol(prototypes)!=ncol(data)) {
    stop("'prototypes' and 'data' have different dimensions")
  }
  if(ncol(nv)!=nrow(nv)) {
    stop("'nv' is not a square matrix")
  }
  if(ncol(nv)!=nrow(prototypes)) {
    stop("'nv' and 'prototypes' have different dimensions")
  }
  if(missing(weights)) {
    weights <- rep(1,nrow(data))
  } else if(length(weights)!=nrow(data)) {
    stop("'weights' and 'data' have different dimensions")
  }
  ## nv must be in row major mode if normalised
  result <- .C("bmu_heskes",
               as.double(prototypes),
               as.double(t(nv)),
               as.integer(nrow(prototypes)),
               as.double(data),
               as.integer(nrow(data)),
               as.integer(ncol(prototypes)),
               as.double(weights),
               clusters=integer(nrow(data)),
               error=as.double(0),
               PACKAGE="FuseSOM")
  list(clusters=result$clusters+1,error=result$error)
}

batchSom.default <- function(data,somGrid,init=c("pca","random"),prototypes,
                             weights,
                             mode = c("continuous","stepwise"),
                             minRadius, maxRadius, steps,
                             decrease = c("power", "linear"), maxIter,
                             kernel = c("gaussian", "linear"), normalised,
                             assignment = c("single", "heskes"),
                             cut = 1e-07,
                             verbose=FALSE,keepdata=TRUE,...) {
  if(class(somGrid)!="somgrid") {
    stop("'somgrid' is not of somgrid class")
  }
  the.call <- match.call()
  if(verbose) {
    print(the.call)
  }
  the.call[[1]] <- batchSom.control
  control <- eval(the.call,envir = parent.frame())
  control$assignmentInt <- switch(control$assignment,"single"=0,"heskes"=1)
  control$kernelInt <- switch(control$kernel,"gaussian"=0,"linear"=1)
  if(!missing(weights)) {
    if(length(weights)!=nrow(data)) {
      stop("'weights' and 'data' have different dimensions")
    }
  } else {
    ## keep weights to NULL for now to avoid princomp initialization
    weights <- NULL
  }
  if(missing(prototypes)) {
    ## initialisation based on the value of init
    init <- match.arg(init)
    args <- list(...)
    params <- c(list(data=data,somGrid=somGrid,weights=weights),args)
    prototypes <- switch(init,
                         "pca"=do.call("somInitPca",params)$prototypes,
                         "random"=do.call("somInitRandom",params))
  } else {
    if(ncol(prototypes)!=ncol(data)) {
      stop("'prototypes' and 'data' have different dimensions")
    }
    if(nrow(prototypes)!=somGrid$size) {
      stop("'prototypes' and 'somgrid' are not compatible")
    }
  }
  if(is.null(weights)) {
    weights <- rep(1,nrow(data))
  }
  ## distances?
  if(is.null(somGrid$dist)) {
    somGrid$dist <- as.matrix(dist(somGrid$pts,method="Euclidean"),diag=0)
  }
  pre <- switch(control$mode,
                "stepwise"=batchsom.lowlevel(somGrid,data,prototypes,weights,
                                             control,verbose),
                "continuous"=batchsom.lowlevelcontinuous(somGrid,data,
                                                         prototypes,weights,control,verbose))
  pre$control <- control
  if(keepdata) {
    pre$data  <- data
    pre$weights <- weights
  }
  pre
}

batchSomLowLevel <- function(somGrid,data,prototypes,weights,control,verbose) {
  result <- .C("batch_som_optim",
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
               as.integer(control$maxIter),
               as.integer(control$kernelInt),
               as.integer(control$normalised),
               as.double(control$cut),
               as.integer(verbose),
               clusters=integer(nrow(data)),
               errors=as.double(rep(-1,1+length(control$radii)*control$maxIter)),
               PACKAGE="FuseSOM")
  prototypes <- matrix(result$proto,ncol=ncol(prototypes),
                       dimnames=list(NULL,dimnames(data)[[2]]))
  res <- list(somGrid=somGrid,
              prototypes=prototypes,
              classif=result$cluster+1,
              errors=result$errors[result$errors>=0])
  class(res) <- c("somnum","som")
  res
}

batchSomLowLeveLcontinuous <- function(somGrid,data,prototypes,weights,
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
               PACKAGE="yasomi")
  prototypes <- matrix(result$proto,ncol=ncol(prototypes),
                       dimnames=list(NULL,dimnames(data)[[2]]))
  res <- list(somGrid=somGrid,
              prototypes=prototypes,
              classif=result$cluster+1,
              errors=result$errors[result$errors>=0])
  class(res) <- c("somnum","som")
  res
}



