#include <R.h>
#include <Rmath.h>
#include "neighborhood.h"

/* Neighborhood function calculation. The kernel functions are defined such
   that neighborhood(x,y) is approximatly zero when the distance between x and
   y is higher or equal to the radius parameter plus one.
   When normalized, a row major approach is used to preserve memory locality. */

void neighborhood(double *distances,double *nv,int nbUnit,double radius,
		  int kernelType,int isNormalized) 
{
    int i,j,base;
    double tmp;
    int size = nbUnit*nbUnit;
    kernel_type eKernelType = (kernel_type) kernelType;
    
    switch(eKernelType) {
    case gaussian:
	/* Gaussian kernel */
	for(i = 0 ; i < size; i++) {
	    tmp = distances[i]/(radius+1);
	    nv[i] = exp (GAUSSIAN_COEFF*tmp*tmp);
	}
	break;
    case linear:
	/* Linear kernel */
	for(i = 0 ; i < size; i++) {
	    tmp = distances[i];
	    nv[i] = tmp < radius + 1 ? (1-tmp/(radius+1)) : 0;
	}
	break;
    case zeroone:
	/* Zero one kernel */
	for(i = 0 ; i < size; i++) {
	    tmp = distances[i];
	    nv[i] = tmp <= radius ? 1 : 0;
	}
	break;
    }
    if(isNormalized) {
	for(i = 0 ; i < nbUnit; i++) {
	    tmp = 0;
	    base = i*nbUnit;
	    for(j = 0 ; j < nbUnit; j++) {
		tmp += nv[base+j];
	    }
	    for(j = 0 ; j < nbUnit; j++) {
		nv[base+j] /= tmp;
	    }
	}
    }
}

void neighborhood_single(double *distances,double *nv,int *nbUnit,
			 double *radius,int *kernelType,int *isNormalized) 
{
    int i;
    double tmp;
    kernel_type eKernelType = (kernel_type) *kernelType;
    int nb=*nbUnit;
    double theRadius=*radius;
    
    switch(eKernelType) {
    case gaussian:
	/* Gaussian kernel */
	for(i = 0 ; i < nb; i++) {
	    tmp = distances[i]/(theRadius+1);
	    nv[i] = exp (GAUSSIAN_COEFF*tmp*tmp);
	}
	break;
    case linear:
	/* Linear kernel */
	for(i = 0 ; i < nb; i++) {
	    tmp = distances[i];
	    nv[i] = tmp < theRadius + 1 ? (1-tmp/(theRadius+1)) : 0;
	}
	break;
    case zeroone:
	/* Zero one kernel */
	for(i = 0 ; i < nb; i++) {
	    tmp = distances[i];
	    nv[i] = tmp <= theRadius ? 1 : 0;
	}
	break;
    }
    if(*isNormalized) {
	tmp = 0;
	for(i = 0 ; i < nb; i++) {
	    tmp += nv[i];
	}
	for(i = 0 ; i < nb; i++) {
	    nv[i] /= tmp;
	}
    }
}
