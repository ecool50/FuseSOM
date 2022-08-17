# A fucntion to estimate the number of clusters using the instability method 
# see https://arxiv.org/abs/1608.07494
# function was obtained from https://github.com/cran/cstab with some minor modifications
.cStability <- function(data, # n x p data matrix
                        kSeq=2:20, # sequence of ks tested
                        nB=10, # number of bootstrap comparisons
                        norm=TRUE, # norm over pw equal assign,FALSE=as in Wang etal
                        predict=TRUE, # use prediction approach, if FALSE, use brute pair in equal cluster approach
                        linkage='average', # or average, or ...
                        pBar=TRUE
)
  
{
  
  
  # ---------- Input Checks ----------
  # On Data
  if(sum(is.na(data))>0) stop('No missing values permitted!')
  
  # On k-sequence
  if(1 %in% kSeq) stop('Please select a k sequence starting with 2: {2,3,...K}!')
  
  # On B
  if(round(nB)!=nB) stop('The number of bootstrap comparison has to be a positive integer value.')
  
  # ---------- Create Containers ----------
  
  nObj <- nrow(data)
  mInstab <- matrix(NA, nB, length(kSeq)) # storage: no normalization
  mInstabNorm <- matrix(NA, nB, length(kSeq)) # storage: normalization
  
  
  # ---------- Draw bootstrap samples ----------
  
  ind <- list()
  share <- list()
  for(b in seq_len(nB)) {
    tmp1 <- sample(seq_len(nObj), nObj, replace=TRUE)
    tmp2 <- sample(seq_len(nObj), nObj, replace=TRUE)
    tmp1 <- tmp1[order(tmp1)]
    tmp2 <- tmp2[order(tmp2)]
    ind[[b]] <- list(tmp1,tmp2)
    intersct <- intersect(tmp1,tmp2)
    share[[b]] <- list(tmp1 %in% intersct & !duplicated(tmp1), tmp2 %in% intersct & !duplicated(tmp2)) # leave duplicates in ?
  }
  
  
  # ---------- Calculate distance matrix ----------
  
  pear <- stats::cor(t(data), method = 'pearson')
  cosi <- coop::tcosine(data)
  spear <- stats::cor(t(data), method = 'spearman')
  
  distM <- as.matrix(fuse(cor2dist(pear),cor2dist(cosi),cor2dist(spear)))
  
  # ----- Loop over B comparisons -----
  
  if(pBar)  pB <- utils::txtProgressBar(min=0, max=nB, style = 2)
  
  for(b in seq_len(nB)) {
    for(k in kSeq) {
      #t = proc.time()[3]
      hc1 = fastcluster::hclust(stats::as.dist(distM[ind[[b]][[1]],ind[[b]][[1]]]), method = linkage)
      hc2 = fastcluster::hclust(stats::as.dist(distM[ind[[b]][[2]],ind[[b]][[2]]]), method = linkage)
      cl1 = stats::cutree(hc1, k)[share[[b]][[1]]]
      cl2 = stats::cutree(hc2, k)[share[[b]][[2]]]
      #print(proc.time()[3] - t)
      
      
      # check for equality of clusterings
      eq1 = equal(cl1)
      eq2 = equal(cl2)
      InStab <- mean(eq1 != eq2)
      
      # Normalize = FALSE
      mInstab[b, which(kSeq==k)] <- InStab
      
      # Normalize = TRUE
      normVal <- .instabLookup(table(cl1), table(cl2))
      mInstabNorm[b, which(kSeq==k)] <- InStab / normVal
      
    } # end for k
    
    if(pBar) utils::setTxtProgressBar(pB, b)
    
  } # end for B
  
  # taking the mean
  mInstabM      = colMeans(mInstab)
  mInstabNormM = colMeans(mInstabNorm)
  
  mInstabM[!is.finite(mInstabM)] <- Inf
  mInstabNormM[!is.finite(mInstabNormM)] <- Inf
  kOptInstab  = which(mInstabM == min(mInstabM))+(min(kSeq)-1)
  kOptInstabN = which(mInstabNormM == min(mInstabNormM))+(min(kSeq)-1)
  
  # replicate function call
  fCall <- list('kSeq'=kSeq,
                 'nB'=nB,
                 'norm'=norm,
                 'prediction'=predict,
                 'linkage'=linkage)
  
  outList <- list("kInstab"=kOptInstab,
                  "kInstabNorm"=kOptInstabN,
                  "instabPath"=mInstabM,
                  "instabPathNorm"=mInstabNormM,
                  "instabPathMatrix"=mInstab,
                  "instabPathNorm_matrix"=mInstabNorm,
                  'call'=fCall)
  
  return(outList)
  
} 