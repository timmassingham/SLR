#ifndef _LINEMIN_H_
#define _LINEMIN_H_

double linemin_multid ( double (*fun)(const double *, void *), int dim, double * x, double *xnew, double * direct, void * info, const double min, const double max, const double tol, const int noisy, int * neval );
double linemin_1d ( double (*fun)(const double *, void *), double * x, void * info, const double min, const double max, const double tol, const int noisy, int * neval);
double linemin_backtrack ( double (*fun)(const double *, void *), int dim, double * x, double *xnew, double * direct, void * info, const double min, const double max, const double tol, const int noisy, int * neval );
#endif

