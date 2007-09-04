#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include "statistics.h"

#define MAXSTEP		1000
#define OOM(A) { if( (A) == NULL ) { fputs ("Out of memory",stderr); fflush(stderr); abort(); } }

static double pow1pm1 ( const double x, const double y);
static int IsPvals ( const double * pval, const int n);
static double CorrectPval ( const double p, const int n, const int method);
int CmpDoublePtr ( const void * dptr1, const void * dptr2);


static double pow1pm1 ( const double x, const double y){
  return expm1( y * log1p(x) );
}

static double CorrectPval ( const double p, const int n, const int method){
  double res;
  assert(p>=0. && p<=1.);
  assert(n>0);
  assert(method==BONFERRONI || method==SIDAK);

  switch(method){
    case BONFERRONI:	res = n * p; break;
    case SIDAK:		res = -pow1pm1(-p,n); break;
    default: puts("Unknown method in CorrectPval"); abort();
  }

  res = (res<=1.)?res:1.;

  assert(res>=0. && res<=1.);
  return res;
}

static int IsPvals ( const double * pval, const int n){
  int i;
  assert (NULL!=pval);
  assert (n>0);
  for ( i=0 ; i<n ; i++){
    if ( pval[i]<0. || pval[i]>1.){
      return 0;
    }
  }
  return 1;
}


double * Pvalue_adjust_SingleStep ( const double * pval, const int n, const int method){
  int i;
  double * pval_adj;
  assert(IsPvals(pval,n));
  assert(method==BONFERRONI || method==SIDAK);

  pval_adj = calloc(n,sizeof(double));
  OOM(pval_adj);

  for ( i=0 ; i<n ; i++){
    pval_adj[i] = CorrectPval(pval[i],n,method);
  }

  assert( IsPvals(pval_adj,n));
  return pval_adj;
}


double * Pvalue_adjust_StepDown ( const double * pval, const int n, const int method){
  int i;
  double * pval_adj,**work;
  assert(IsPvals(pval,n));
  assert(method==BONFERRONI || method==SIDAK);

  pval_adj = calloc(n,sizeof(double));
  OOM(pval_adj);
  work = calloc(n,sizeof(double *));
  OOM(work);
  for ( i=0 ; i<n ; i++){
    pval_adj[i] = pval[i];
    work[i] = &pval_adj[i];
  }

  qsort (work,n,sizeof(double *),CmpDoublePtr);

  *work[0] = CorrectPval(*work[0],n,method);
  for ( i=1 ; i<n ; i++){
    *work[i] = CorrectPval(*work[i],n-i,method);
    *work[i] = (*work[i]>*work[i-1])?(*work[i]):(*work[i-1]);
  }

  free(work);
  assert(IsPvals(pval_adj,n));

  return pval_adj;
}


double * Pvalue_adjust_StepUp ( const double * pval, const int n, const int method){
  int i;
  double * pval_adj,**work;
  assert(IsPvals(pval,n));

  if(method==SIDAK){
    fputs("Using Sidak adjustment with Hochberg\'s step-up proceedure.\nMay not be valid -- check theory!\n",stderr);
  }

  pval_adj = calloc(n,sizeof(double));
  OOM(pval_adj);
  work = calloc(n,sizeof(double *));
  OOM(work);
  for ( i=0 ; i<n ; i++){
    pval_adj[i] = pval[i];
    work[i] = &pval_adj[i];
  }

  qsort (work,n,sizeof(double *),CmpDoublePtr);

  for ( i=(n-2) ; i>=0 ; i--){
    *work[i] = CorrectPval( *work[i],n-i,method);
    *work[i] = (*work[i]<*work[i+1])?(*work[i]):(*work[i+1]);
  }

  free(work);
  assert(IsPvals(pval_adj,n));

  return pval_adj;
}

int CmpDoublePtr ( const void * dptr1, const void * dptr2){
  if ( **(const double **)dptr1 > **(const double **)dptr2)
    return 1;

  return -1;
}
