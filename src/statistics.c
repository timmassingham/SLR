/*
 *  Copyright 2003-2007 Tim Massingham (tim.massingham@ebi.ac.uk)
 *  Funded by EMBL - European Bioinformatics Institute
 */
/*
 *  This file is part of SLR ("Sitewise Likelihood Ratio")
 *
 *  SLR is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SLR is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SLR.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <err.h>
#include <errno.h>
#include "statistics.h"

#define MAXSTEP		1000
#define OOM(A) { if( (A) == NULL ) { fputs ("Out of memory",stderr); fflush(stderr); exit(EXIT_FAILURE); } }

static double pow1pm1 ( const double x, const double y);
static int IsPvals ( const double * pval, const int n);
static double CorrectPval ( const double p, const int n, const int method);
int CmpDoublePtr ( const void * dptr1, const void * dptr2);
int CmpDouble ( const void * dptr1, const void * dptr2);
double sample_variance_naive ( const VEC v);
double variance_naive ( const VEC v);

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

int CmpDouble ( const void * dptr1, const void * dptr2){
	if ( *(const double *)dptr1 > *(const double *)dptr2 ){
		return 1;
	}

	return -1;
}

/* Note: quick and dirty algorithm; not numerically stable.
 * Might need to be rewritten at some point
 */
double mean ( const VEC v){
	assert(NULL!=v);
	unsigned int len = v->n;
	double x = 0.;
	for ( unsigned int i=0 ; i<len ; i++){
		x += v->x[i];
	}
	return x/len;
}

double median ( const VEC v){
	return quantile(v,0.5);
}

double quantile_fromsorted ( const VEC v, double q){
	assert(NULL!=v);
	assert(q>=0. && q<=1.);

	const unsigned int idx = (unsigned int)(q*(vlen(v)-1));
	const double w = q*(vlen(v)-1) - idx;
	/* Deal with corner case when q exactly one. Due to
	 * (assumed) right continuity of indexing
	 */
	if ( idx == vlen(v) - 1){ return vget(v,idx);}
	return (1.-w)*vget(v,idx) + w*vget(v,idx+1);
}

double quantile (const VEC v, const double q){
	assert(NULL!=v);
	assert(q>=0. && q<=1.);
	VEC vcopy = copy_vec(v);
	qsort (vcopy->x,vlen(vcopy),sizeof(double),CmpDouble);
	const double med = quantile_fromsorted(vcopy,q);
	free_vec(vcopy);
	return med;
}

VEC quantiles ( const VEC v, const VEC q){
	assert(NULL!=q);
	assert(NULL!=v);

	VEC quant = create_vec(vlen(q));
	assert(NULL!=quant);
	VEC vcopy = copy_vec(v);
	qsort (vcopy->x,vlen(vcopy),sizeof(double),CmpDouble);

	for ( unsigned int i=0; i<vlen(q) ; i++){
		assert(vget(q,i)>=0. && vget(q,i)<=1.);
		vset(quant,i,quantile_fromsorted(vcopy,vget(q,i)));
	}
	free_vec(vcopy);
	return quant;
}

double sample_variance_naive ( const VEC v){
	assert(NULL!=v);
	const unsigned int len = vlen(v);
	const double vmean = mean(v);
	double sumsqr = 0.;
	for ( unsigned int i=0 ; i<len ; i++){
		sumsqr += (vget(v,i)-vmean)*(vget(v,i)-vmean);
	}
	return sumsqr;
}

double variance_naive ( const VEC v){
	return sample_variance_naive(v)/(vlen(v)-1);
}

/*  Variance of data by corrected two-pass algorithm  */
double variance (const VEC v){
	assert(NULL!=v);
	const unsigned int len = vlen(v);
	const double vmean = mean(v);
	const double correction = suma_vec(v,-vmean);

	return (sample_variance_naive(v)-correction*correction/len)/(len-1);
}

double sd(const VEC v){
	return sqrt(variance(v));
}

#define MADSCALE 1.4826 /* Constant to scale MAD so it is comparable to SD */
/*  Some unnecessary vector copies when this function is combined with median */
double mad ( const VEC v){
	assert(NULL!=v);
	const double vmedian = median(v);
	VEC vcopy = copy_vec (v);
	const unsigned int len = vlen(vcopy);
	for ( unsigned int i=0 ; i<len ; i++){
		vset(vcopy,i,vget(vcopy,i)-vmedian);
	}
	const double devmedian = median(vcopy);
	free_vec(vcopy);

	return devmedian * MADSCALE;
}
