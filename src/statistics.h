#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#define BONFERRONI	0
#define SIDAK		1

double * Pvalue_adjust_SingleStep ( const double * pval, const int n, const int method);
double * Pvalue_adjust_StepDown ( const double * pval, const int n, const int method);
double * Pvalue_adjust_StepUp ( const double * pval, const int n, const int method);
#endif

