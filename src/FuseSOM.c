#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void batch_som_optim(double *proto, int *nproto, double *data, int *ndata,
                            int *dim, double *weights, int *assign, double *grid,
                            double *radii, int *nradii, int *maxiter, int *kernel,
                            int *normalized, double *cut, int *verbose,
                            int *classif, double *errors);
extern void batch_som_optim_continuous(double *proto, int *nproto, double *data,
                                       int *ndata, int *dim, double *weights,
                                       int *assign, double *grid, double *radii,
                                       int *nradii, int *kernel,
                                       int *normalized, double *cut, int *verbose,
                                       int *classif, double *errors);
extern int bmu(double *proto, int *nproto, double *data, int *ndata, int *dim,
               double *weights, int *winner, double *error);
extern void second_bmu(double *proto, int *nproto, double *data, int *ndata, int *dim,
                       int *winners);

static const R_CMethodDef CEntries[] = {
    {"batch_som_optim", (DL_FUNC)&batch_som_optim, 17},
    {"batch_som_optim_continuous", (DL_FUNC)&batch_som_optim_continuous, 16},
    {"bmu", (DL_FUNC)&bmu, 8},
    {"second_bmu", (DL_FUNC)&second_bmu, 6},
    {NULL, NULL, 0}};

void R_init_FuseSOM(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}