#ifndef YASOMI_H
#define YASOMI_H

int bmu(double *proto,int *nproto,double *data,int *ndata,int *dim,
        double *weights,int *winner,double *error);

void second_bmu(double *proto,int *nproto,double *data,int *ndata,int *dim,
                int *winners);

int bmu_heskes(double *proto,double *neigh,int *nproto,double *data,
               int *ndata,int *dim,double *weights,int *winner,double *error);

int bmu_heskes_ext_mem(double *proto,double *neigh,int *nproto,double *data,
                       int *ndata,int *dim,double *weights,int *winner,
                       double *error,double *distances);

#endif /* !YASOMI_H */