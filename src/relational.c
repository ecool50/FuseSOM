#include <R.h>

/* Partial sums computation.
 clusters(ndata): mapping between index and cluster
 ncluster: number of clusters
 diss(ndata,ndata): dissimilarity matrix
 sums(ndata,ncluster): result matrix; each line contains the sums over each
 cluster of the dissimilarities 
 bisums(ncluster,ncluster): second result matrix; contains the sums of all
 pairwise dissimilarities between elements of a given cluster
 */
void partial_sums(int *clusters,int *ndata,int *ncluster,double *diss,
                  double *sums,double *bisums) 
{
  int dataSize=*ndata;
  int nbCluster=*ncluster;
  int indiv,k,tmp,pos;
  
  /* should be useless */ 
  memset(sums,0,nbCluster*dataSize*sizeof(double));
  memset(bisums,0,nbCluster*nbCluster*sizeof(double));
  
  /* calculation of the partial sums */
  /* loop on individuals */
  for(indiv = 0 ; indiv < dataSize ; indiv++) {
    /* inner loop on individuals */
    for(k = 0 ; k < dataSize; k++) {
      sums[indiv * nbCluster + clusters[k]] += diss[indiv * dataSize + k];
    }
  }
  /* calculation of the partial sums sums */
  /* loop on individuals */
  for(indiv = 0 ; indiv < dataSize ; indiv++) {
    /* inner loop on clusters */
    tmp = indiv * nbCluster;
    pos = clusters[indiv] * nbCluster;
    for(k = 0 ; k < nbCluster; k++) {
      bisums[pos + k] += sums[tmp + k];
    }
  }
}

/* Weighted partial sums computation.
 clusters(ndata): mapping between index and cluster
 ncluster: number of clusters
 diss(ndata,ndata): dissimilarity matrix
 weights(ndata): weights
 sums(ndata,ncluster): result matrix
 bisums(ncluster,ncluster): second result matrix
 */
void weighted_partial_sums(int *clusters,int *ndata,int *ncluster,double *diss,
                           double *weights,double *sums,double *bisums) 
{
  int dataSize=*ndata;
  int nbCluster=*ncluster;
  int indiv,k,tmp,pos;
  double weight;
  
  /* should be useless */ 
  memset(sums,0,nbCluster*dataSize*sizeof(double));
  memset(bisums,0,nbCluster*nbCluster*sizeof(double));
  
  /* calculation of the partial sums */
  /* loop on individuals */
  for(indiv = 0 ; indiv < dataSize ; indiv++) {
    /* inner loop on individuals */
    for(k = 0 ; k < dataSize; k++) {
      sums[indiv * nbCluster + clusters[k]] += diss[indiv * dataSize + k] * weights[k];
    }
  }
  /* calculation of the partial sums sums */
  /* loop on individuals */
  for(indiv = 0 ; indiv < dataSize ; indiv++) {
    /* inner loop on clusters */
    tmp = indiv * nbCluster;
    weight = weights[indiv];
    pos = clusters[indiv] * nbCluster;
    for(k = 0 ; k < nbCluster; k++) {
      bisums[pos + k] += sums[tmp + k] * weight;
    }
  }
}



/* 
 Normalisation coefficients calculation for partial sums (t(h)BIPSh)
 bips(ncluster,ncluster): bi sums
 h(ncluster,ncluster): coefficients
 nf(ncluster): result vector
 */
void th_bips_h(double *bips,double *h,int *ncluster,double *nf) 
{
  int nbCluster=*ncluster;
  int j,c,d; /* cluster indexes */
double sum,interm;
int posC;

/* result loop */
for(j = 0; j<nbCluster; j++) {
  sum = 0;
  for(c = 0; c<nbCluster; c++) {
    interm = 0;
    posC = c*nbCluster;
    for(d = 0; d<nbCluster; d++) {
      interm += bips[posC+d]*h[j+d*nbCluster];
    }
    sum += interm*h[j+c*nbCluster];
  }
  nf[j] = sum;
}
}
