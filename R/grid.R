# A function to generate the grid for the self organizing map
# function was obtained from https://rdrr.io/rforge/yasomi/ with some major modifications

somGrid <- function(xDim,yDim,topo="hexagonal",withDist=TRUE) {
  topo <- match.arg(topo)
  if(xDim == 1) {
    tmp <- xDim
    xDim <- yDim
    yDim <- tmp
  }
  x <- seq(from=1,by=1,length.out=xDim)

  y <- rev(seq(from=1,by=sqrt(3)/2,length.out=yDim))
  
  pts <- as.matrix(expand.grid(x = x, y = y))
  pts[,1] <- pts[,1]+rep(c(0,0.5),each=xDim,length.out=nrow(pts))
  diam <- sqrt(sum((pts[1,]-pts[nrow(pts),])^2))
  
  res <- list(pts = pts, xDim = xDim, yDim = yDim, topo = topo,
              size = xDim*yDim, diam = diam)
  if(withDist) {
    res$dist <- as.matrix(dist(pts,method="euclidean"),diag=0)
  }
  class(res) <- "somgrid"
  return(res)
}

