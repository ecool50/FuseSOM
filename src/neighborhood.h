#ifndef FUSESOM_NEIGHBORHOOD_H
#define FUSESOM_NEIGHBORHOOD__H

typedef enum {
    gaussian = 0,
    linear = 1,
    zeroone = 2
} kernel_type;

#define GAUSSIAN_COEFF -4.605170185988090914009

void neighborhood(double *distances,double *nv,int nbUnit,double radius,
		  int kernelType,int isNormalized);

void neighborhood_single(double *distances,double *nv,int *nbUnit,
			 double *radius,int *kernelType,int *isNormalized);

#endif /* !FUSESOM_NEIGHBORHOOD_H */
