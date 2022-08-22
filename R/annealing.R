# A function to generate the radii values for the Gaussian kernel
# function was obtained from https://rdrr.io/rforge/yasomi/ with some major modifications

radiusExp <- function(min,max,steps) {
  return(max*(min/max)^(seq(0,1,length.out=steps)))
}


# A function to set the default values for the self organizing map
# function was obtained from https://rdrr.io/rforge/yasomi/ with some major modifications
batchSomControl.default <- function(data,somGrid,
                                     mode="continuous",
                                     minRadius, maxRadius, steps,
                                     decrease="power", maxIter=1,
                                     kernel="gaussian",
                                     assignment="single",
                                     cut=1e-07,...)
{
  if(missing(maxRadius)) {
    maxRadius <- 2/3*somGrid$diam+1
  }
  if(missing(minRadius)) {
    minRadius <- 0.5
  }
  if(maxRadius<=minRadius) {
    stop("maxRadius must be larger than minRadius")
  }
  if(minRadius<=0) {
    stop("minRadius must be positive")
  }
  if(missing(steps)) {
    steps <- max(20,ceiling(5*maxRadius))
  }
  
  radii <- radiusExp(minRadius,maxRadius,steps)
  return(list(mode=mode,radii=radii,maxIter=maxIter,kernel=kernel,
       assignment=assignment,
       cut=cut))
}