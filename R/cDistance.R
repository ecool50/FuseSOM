# A fucntion to estimate the number of clusters using the distance method 
# see https://arxiv.org/abs/1608.07494
# function was obtained from https://github.com/cran/cstab with some minor modifications

.cDistance <- function(data, # n x p data matrix
                      kseq, #sequence of ks to be checked
                      linkage = 'average',
                      gapIter = 10) # number of simulated datasets in gap statistic
{
  
  # ----- INPUT TESTS
  
  # On Data
  if(sum(is.na(data))>0) stop('No missing values permitted!')
  
  # On k-sequence
  if(1 %in% kseq) stop('Please select a k sequence starting with 2: {2,3,...K}!')
    
    
    # ----- HELPERS
    
  n    = nrow(data)
  dims = ncol(data)
  if(!1 %in% kseq) kseq = c(1,kseq)
  
  
  # ----- EVALUATE REAL DATA
  
  WCD <- Sil <- MSE <- numeric()
  for(k in kseq) {
    obj = .getMeasures(data = data, k = k)
    WCD[k] = obj$WCD
    Sil[k] = obj$Sil
    MSE[k] = obj$MSE
  }
  
  
  # ----- EVALUATE SYNTHETIC DATA (Gap-statistic)
  
  WCD_runs = matrix(NA,nrow=gapIter, ncol=length(kseq))
  for(i in 1:gapIter) {
    data_syn = .uniformData(data)
    WCDs = numeric()
    for(j in 1:length(kseq)) {
      k = kseq[j]
      obj = .getMeasures(data = data_syn, k = k, measures = c('wcd'))
      WCDs[j] = obj$WCD
    }
    WCD_runs[i,] = WCDs
  }
  WCD_syn = colMeans(WCD_runs)
  
  # ----- COMPUTE MEASURES
  
  # Gap Statistic
  WCD_dat_log = log(WCD)
  WCD_syn_log = log(WCD_syn)
  WCD_dat_log = WCD_dat_log - WCD_dat_log[1]
  WCD_syn_log = WCD_syn_log - WCD_syn_log[1]
  gap      = WCD_syn_log - WCD_dat_log
  kopt_gap = kseq[gap == max(gap)]
  
  
  # Slope Statistic
  p = 1
  slope = -(Sil[-1] - Sil[-length(Sil)]) * Sil[-1]^p
  kopt_slope = kseq[slope == max(slope)]
  
  ## Jump Statistic
  MSE_tr = MSE^(- dims/2)
  jump   = (MSE_tr - c(0, MSE_tr[-length(MSE_tr)])) #[-1]
  kopt_jump = kseq[jump == max(jump)]
  
  outlist <- list('k_Gap'=kopt_gap,
                  'k_Slope'=kopt_slope,
                  'k_Jump'=kopt_jump,
                  'WCD'=WCD,
                  'WCD_syn'=WCD_syn,
                  'Gaps'= gap,
                  'Silhouettes'=Sil,
                  'Slopes'=slope,
                  'Jumps'=jump)
  
  return(outlist)
  
} # EoF

