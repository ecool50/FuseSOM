#include <R.h>

int bmu(double *proto,int *nproto,double *data,int *ndata,int *dim,
        double *weights,int *winner,double *error) 
{
  double bestDist,dist,tmp;
  int i,j,k;
  int bestSoFar;
  int dataSize=*ndata,protoSize=*nproto,dimension=*dim;
  int changed = 0;
  double norm = 0;
  
  *error = 0;
  /* loop on the individuals */
  for(i = 0; i < dataSize; i++) {
    bestDist = R_PosInf;
    bestSoFar = -1;
    /* loop on prototypes */
    for(j = 0; j < protoSize; j++) {
      dist = 0;
      /* loop on dimensions */
      for(k = 0; k < dimension; k++) {
        tmp = data[i + k * dataSize] - proto[j + k * protoSize];
        dist += tmp * tmp;
      }
      if(dist < bestDist) {
        bestDist = dist;
        bestSoFar = j;
      }
    }
    *error += weights[i]*bestDist;
    norm += weights[i];
    if(bestSoFar != winner[i]) {
      winner[i] = bestSoFar;
      changed = 1;
    }
  }
  *error /= norm;
  return changed;
}

void second_bmu(double *proto,int *nproto,double *data,int *ndata,int *dim,
                int *winners) 
{
  double bestDist,secondBestDist,dist,tmp;
  int i,j,k;
  int bestSoFar,secondBestSoFar;
  int dataSize=*ndata,protoSize=*nproto,dimension=*dim;
  
  /* loop on the individuals */
  for(i = 0; i < dataSize; i++) {
    bestDist = R_PosInf;
    secondBestDist = R_PosInf;
    bestSoFar = -1;
    secondBestSoFar = -1;
    /* loop on prototypes */
    for(j = 0; j < protoSize; j++) {
      dist = 0;
      /* loop on dimensions */
      for(k = 0; k < dimension; k++) {
        tmp = data[i + k * dataSize] - proto[j + k * protoSize];
        dist += tmp * tmp;
      }
      if(dist < secondBestDist) {
        if(dist < bestDist) {
          secondBestDist = bestDist;
          secondBestSoFar = bestSoFar;
          bestDist = dist;
          bestSoFar = j;
        } else {
          secondBestDist = dist;
          secondBestSoFar = j;
        }
      }
    }
    winners[i] = bestSoFar + 1;
    winners[i+dataSize] = secondBestSoFar + 1;
  }
  
}

int bmu_heskes_ext_mem(double *proto,double *neigh,int *nproto,double *data,
                       int *ndata,int *dim,double *weights,int *winner,
                       double *error,double *distances) 
{
  double bestDist,dist,fullDist,tmp;
  int i,j,k,n;
  int base;
  int bestSoFar;
  int dataSize=*ndata,protoSize=*nproto,dimension=*dim;
  int changed = 0;
  double norm = 0;
  
  *error = 0;
  /* loop on the individuals */
  for(i = 0; i < dataSize; i++) {
    bestDist = R_PosInf;
    bestSoFar = -1;
    /* first compute all distances */
    for(j = 0; j < protoSize; j++) {
      dist = 0;
      /* loop on dimensions */
      for(k = 0; k < dimension; k++) {
        tmp = data[i + k * dataSize] - proto[j + k * protoSize];
        dist += tmp * tmp;
      }
      distances[j] = dist;
    }
    /* then compute Heskes distances */
    /* loop on prototypes */
    for(j = 0; j < protoSize; j++) {
      fullDist = 0;
      /* loop on "neighbors" */
      base = j * protoSize;
      for(n = 0; n < protoSize; n++) {
        fullDist += distances[n] * neigh[ base + n ];
      }
      if(fullDist < bestDist) {
        bestDist = fullDist;
        bestSoFar = j;
      }
    }
    *error += weights[i]*distances[bestSoFar];
    norm += weights[i];
    if(bestSoFar != winner[i]) {
      winner[i] = bestSoFar;
      changed = 1;
    }
  }
  *error /= norm;
  return changed;
}

int bmu_heskes(double *proto,double *neigh,int *nproto,double *data,
               int *ndata,int *dim,double *weights,int *winner,double *error) 
{
  /* storage for intermediate results */
  double *distances = (double *) R_alloc(*nproto, sizeof(double));
  
  return bmu_heskes_ext_mem(proto,neigh,nproto,data,ndata,dim,weights,winner,
                            error,distances);
}
