#ifndef FUSESOM_BMU_H
#define FUSESOM_BMU_H

int bmu(double *proto,int *nproto,double *data,int *ndata,int *dim,
	double *weights,int *winner,double *error);

void second_bmu(double *proto,int *nproto,double *data,int *ndata,int *dim,
		int *winners);


int bmu_heskes(double *proto,double *neigh,int *nproto,double *data,
	       int *ndata,int *dim,double *weights,int *winner,double *error);

int bmu_heskes_ext_mem(double *proto,double *neigh,int *nproto,double *data,
		       int *ndata,int *dim,double *weights,int *winner,
		       double *error,double *distances);

void bmu_single(double *proto,int *nproto,double *one_data,int *dim,
		int *winner,double *error);

void bmu_single_full_data(double *proto,int *nproto,double *data,int *ndata,
			  int *dim,int *indiv,int *winner,double *error);


#endif /* !FUSESOM_BMU_H */
