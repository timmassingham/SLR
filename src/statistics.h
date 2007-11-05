#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#define BONFERRONI	0
#define SIDAK		1

#ifndef _VEC_H_
#include "vec.h"
#endif

double * Pvalue_adjust_SingleStep ( const double * pval, const int n, const int method);
double * Pvalue_adjust_StepDown ( const double * pval, const int n, const int method);
double * Pvalue_adjust_StepUp ( const double * pval, const int n, const int method);
double mean ( const VEC v);
double median ( const VEC v);
double quantile (const VEC v, const double q);
VEC quantiles ( const VEC v, const VEC q);
double variance (const VEC v);
double sd(const VEC v);
double mad ( const VEC v);
#endif

