#ifndef _ROOT_H_
#define _ROOT_H_

double find_root ( const double min, const double max, double (*f)(const double*,void*), void *info, double *fmin, double *fmax, const double tol, int * neval );

#endif
