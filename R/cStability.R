# A fucntion to estimate the number of clusters using the instability method 
# see https://arxiv.org/abs/1608.07494
# function was obtained from https://github.com/cran/cstab with some minor modifications
.cStability <- function(data, # n x p data matrix
                       kseq = 2:20, # sequence of ks tested
                       nB   = 10, # number of bootstrap comparisons
                       norm = TRUE, # norm over pw equal assign,FALSE=as in Wang etal
                       predict = TRUE, # use prediction approach, if FALSE, use brute pair in equal cluster approach
                       linkage = 'average', # or average, or ...
                       pbar = T
                        )

{
  
  
  # ---------- Input Checks ----------
  set.seed(1994)
  # On Data
  if(sum(is.na(data))>0) stop('No missing values permitted!')
  
  # On k-sequence
  if(1 %in% kseq) stop('Please select a k sequence starting with 2: {2,3,...K}!')
  
  # On B
  if(round(nB)!=nB) stop('The number of bootstrap comparison has to be a positive integer value.')
  
    # ---------- Create Containers ----------
  
  n_obj <- nrow(data)
  m_instab <- matrix(NA, nB, length(kseq)) # storage: no normalization
  m_instab_norm <- matrix(NA, nB, length(kseq)) # storage: normalization
  
  
  # ---------- Draw bootstrap samples ----------
  
  ind <- list()
  share <- list()
  for(b in 1:nB) {
    tmp_1 <- sample(1:n_obj, n_obj, replace=T)
    tmp_2 <- sample(1:n_obj, n_obj, replace=T)
    tmp_1 <- tmp_1[order(tmp_1)]
    tmp_2 <- tmp_2[order(tmp_2)]
    ind[[b]] <- list(tmp_1,tmp_2)
    intersct <- intersect(tmp_1,tmp_2)
    share[[b]] <- list(tmp_1 %in% intersct & !duplicated(tmp_1), tmp_2 %in% intersct & !duplicated(tmp_2)) # leave duplicates in ?
  }
  
  
  # ---------- Calculate distance matrix ----------
  
  pear <- stats::cor(t(data), method = 'pearson')
  cosi <- coop::tcosine(data)
  spear <- stats::cor(t(data), method = 'spearman')
  
  distm <- as.matrix(fuse(cor2dist(pear),cor2dist(cosi),cor2dist(spear)))
  
  # ----- Loop over B comparisons -----
  
  if(pbar)  pb <- utils::txtProgressBar(min=0, max=nB, style = 2)
  
  for(b in 1:nB) {
    for(k in kseq) {
        #t = proc.time()[3]
      hc_1 = fastcluster::hclust(stats::as.dist(distm[ind[[b]][[1]],ind[[b]][[1]]]), method = linkage)
      hc_2 = fastcluster::hclust(stats::as.dist(distm[ind[[b]][[2]],ind[[b]][[2]]]), method = linkage)
      cl_1 = stats::cutree(hc_1, k)[share[[b]][[1]]]
      cl_2 = stats::cutree(hc_2, k)[share[[b]][[2]]]
        #print(proc.time()[3] - t)

      
      # check for equality of clusterings
      eq_1 = equal(cl_1)
      eq_2 = equal(cl_2)
      InStab <- mean(eq_1 != eq_2)
      
      # Normalize = FALSE
      m_instab[b, which(kseq==k)] <- InStab
      
      # Normalize = TRUE
      norm_val <- .instabLookup(table(cl_1), table(cl_2))
      m_instab_norm[b, which(kseq==k)] <- InStab / norm_val
      
    } # end for k
    
    if(pbar) utils::setTxtProgressBar(pb, b)
    
  } # end for B
  
  # taking the mean
  m_instab_M      = colMeans(m_instab)
  m_instab_norm_M = colMeans(m_instab_norm)
  
  m_instab_M[!is.finite(m_instab_M)] <- Inf
  m_instab_norm_M[!is.finite(m_instab_norm_M)] <- Inf
  kopt_instab  = which(m_instab_M == min(m_instab_M))+(min(kseq)-1)
  kopt_instabN = which(m_instab_norm_M == min(m_instab_norm_M))+(min(kseq)-1)
  
  # replicate function call
  f_call <- list('kseq'=kseq,
                 'nB'=nB,
                 'norm'=norm,
                 'prediction'=predict,
                 'linkage'=linkage)
  
  outlist <- list("k_instab"=kopt_instab,
                  "k_instab_norm"=kopt_instabN,
                  "instab_path"=m_instab_M,
                  "instab_path_norm"=m_instab_norm_M,
                  "instab_path_matrix"=m_instab,
                  "instab_path_nrom_matrix"=m_instab_norm,
                  'call'=f_call)
  
  return(outlist)
  
} 
