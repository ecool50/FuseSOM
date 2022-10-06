# A function to estimate the number of clusters using the distance method
# see https://arxiv.org/abs/1608.07494
# function was obtained from https://github.com/cran/cstab with some
# minor modifications

.cDistance <- function(data, # n x p data matrix
                       kSeq, # sequence of ks to be checked
                       linkage = "average",
                       gapIter = 10) # n of simulated datasets in gap statistic
{

  # ----- INPUT TESTS

  # On Data
  if (sum(is.na(data)) > 0) stop("No missing values permitted!")

  # On k-sequence
  if (1 %in% kSeq) {
    stop("Please select a k sequence starting with 2: {2,3,...K}!")
  }


  # ----- HELPERS

  n <- nrow(data)
  dims <- ncol(data)
  if (!1 %in% kSeq) kSeq <- c(1, kSeq)


  # ----- EVALUATE REAL DATA

  WCD <- Sil <- MSE <- numeric()
  for (k in kSeq) {
    obj <- .getMeasures(data = data, k = k)
    WCD[k] <- obj$WCD
    Sil[k] <- obj$Sil
    MSE[k] <- obj$MSE
  }


  # ----- EVALUATE SYNTHETIC DATA (Gap-statistic)

  WCDRuns <- matrix(NA, nrow = gapIter, ncol = length(kSeq))
  for (i in seq_len(gapIter)) {
    dataSyn <- .uniformData(data)
    WCDs <- numeric()
    for (j in seq_along(kSeq)) {
      k <- kSeq[j]
      obj <- .getMeasures(data = dataSyn, k = k, measures = c("wcd"))
      WCDs[j] <- obj$WCD
    }
    WCDRuns[i, ] <- WCDs
  }
  WCDSyn <- colMeans(WCDRuns)

  # ----- COMPUTE MEASURES

  # Gap Statistic
  WCDDatLog <- log(WCD)
  WCDSynLog <- log(WCDSyn)
  WCDDatLog <- WCDDatLog - WCDDatLog[1]
  WCDSynLog <- WCDSynLog - WCDSynLog[1]
  gap <- WCDSynLog - WCDDatLog
  koptGap <- kSeq[gap == max(gap)]


  # Slope Statistic
  p <- 1
  slope <- -(Sil[-1] - Sil[-length(Sil)]) * Sil[-1]^p
  koptSlope <- kSeq[slope == max(slope)]

  ## Jump Statistic
  MSETr <- MSE^(-dims / 2)
  jump <- (MSETr - c(0, MSETr[-length(MSETr)])) # [-1]
  koptJump <- kSeq[jump == max(jump)]

  outlist <- list(
    "kGap" = koptGap,
    "kSlope" = koptSlope,
    "kJump" = koptJump,
    "WCD" = WCD,
    "WCDSyn" = WCDSyn,
    "Gaps" = gap,
    "Silhouettes" = Sil,
    "Slopes" = slope,
    "Jumps" = jump
  )

  return(outlist)
} # EoF
