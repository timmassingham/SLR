#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "linemin.h"
#include "spinner.h"

/* Prototype for Brent*/
double fminbr ( double a, double b, double (*f)(double), double tol, int * neval);

static double (*linemin_function)(const double *, void *) = NULL;
static int linemin_dim = 0;
static double * linemin_x;
static double * linemin_xnew;
static double * linemin_direct;
static double * linemin_igrad;
static void * linemin_info;
static int linemin_noisy;
static SPINNER * linemin_spinner;

static void setf ( double (*f)(const double *, void *), int dim, double * x, double *xnew, double *direct, double * igrad, void * info, int noisy );
static void unsetf ( void);
static double fun_wrapper ( double x);
static double qfun_wrapper (double x);
static int CheckLinemin ( void);
static void setf1d ( double (*f)(const double *, void *), void * info, int noisy);
static void unsetf1d (void);
static double fun_wrapper1d ( double x);
static int CheckLinemin1D (void);
static void ScaleDirection( double * direct, const double min, const double max, const int n);




double linemin_multid ( double (*fun)(const double *, void *), int dim, double * x, double *xnew, double * direct, void * info, const double min, const double max, const double tol, const int noisy, int * neval ){
  int i;
	int localneval=0;
  double res,fx;
  assert(NULL!=fun);
  assert(dim>1);
  assert(NULL!=x);
  assert(NULL!=xnew);
  assert(NULL!=direct);
  assert(NULL!=info);
  assert(min<max);
  assert(tol>=0.);
  assert(noisy==0 || noisy==1);

  for ( i=0 ; i<dim ; i++){
    x[i] += min * direct[i];
  }
  ScaleDirection(direct,min,max,dim);
  setf (fun, dim, x, xnew, direct, NULL, info, noisy);
  res = fminbr(0.,1.,fun_wrapper,tol,&localneval);
  fx = fun_wrapper (res); localneval++;
  unsetf();
  for ( i=0 ; i<dim ; i++){
    x[i] += res * direct[i];
  }
	*neval += localneval;
//printf("Linesearch res = %e (%e,%e)\t%e\t%d evals\n",min+res*(max-min),min,max,res,localneval);
  return fx;
}

double linemin_backtrack ( double (*fun)(const double *, void *), int dim, double * x, double *xnew, double * direct, void * info, const double min, const double max, const double tol, const int noisy, int * neval ){
  int i;
  double res,fx;
  assert(NULL!=fun);
  assert(dim>1);
  assert(NULL!=x);
  assert(NULL!=xnew);
  assert(NULL!=direct);
  assert(NULL!=info);
  assert(min<max);
  assert(tol>=0.);
  assert(noisy==0 || noisy==1);

	double range = max-min;	
	for ( int i=0 ; i<dim ; i++){
		xnew[i] = x[i] + min*direct[i];
	}
  double fmin = fun(xnew,info);
  *neval = *neval + 1;

  double fact = 2.;  
  double f;
  do {
		fact /= 2.;
		for ( int i=0 ; i<dim ; i++){
			xnew[i] = x[i] + (min+fact*range) * direct[i];
		}
		f = fun(xnew,info);
		*neval = *neval + 1;
  } while ( f>fmin && fact>tol);
  
  if ( fact>tol){
  	/*  Found acceptable point. Keep bisecting until no
  	 * further improvement
  	 */
  	do{
  		fmin = f;
  		fact /= 2.;
			for ( int i=0 ; i<dim ; i++){
				xnew[i] = x[i] + (min+fact*range) * direct[i];
			}
			f = fun(xnew,info);
       *neval = *neval + 1;
  	} while ( f<fmin && fact>tol);
  	/*  Exited because point found is worse */
  	if ( fact>tol){
			for ( int i=0 ; i<dim ; i++){
				xnew[i] = x[i] + (min+2.*fact*range) * direct[i];
			}
			f = fmin;
  	}	
  }

	//printf("Backtrack factor = %e\t%e\t%e\n",fact,min,max);

  for ( int i=0 ; i<dim ; i++){
		x[i] = xnew[i];
  }

  return f;
}


void ScaleDirection( double * direct, const double min, const double max, const int n){
  int i;
  assert(NULL!=direct);
  assert(max>min);
  assert(n>0);

  for ( i=0 ; i<n ; i++){
    direct[i] = direct[i] * (max-min);
  };
}

double linemin_1d ( double (*fun)(const double *, void *), double * x, void * info, const double min, const double max, const double tol, const int noisy, int * neval ){
  double res,fx;
  assert(NULL!=fun);
  assert(NULL!=x);
  assert(NULL!=info);
  assert(min<max);
  assert(tol>0.);
  assert(noisy==0 || noisy==1);

  setf1d (fun,info,noisy);
  res = fminbr ( min,max,fun_wrapper1d,tol,neval);
  fx = fun_wrapper1d(res); *neval = *neval + 1;
  x[0] = res;
  unsetf1d();

  return fx;
}

static void setf1d ( double (*fun)(const double *, void *), void * info, int noisy){
  assert (NULL!=fun);
  assert (NULL!=info);

  linemin_info = info;
  linemin_function = fun;
  linemin_noisy = noisy;
  if ( noisy){
    linemin_spinner = CreateSpinner(2);
  }

  assert(CheckLinemin1D());
}

static void unsetf1d ( void ){
  assert(CheckLinemin1D());

  linemin_info = NULL;
  linemin_function = NULL;
  if ( linemin_noisy){
    DeleteSpinner (linemin_spinner);
    linemin_noisy = 0;
  }
  linemin_spinner = NULL;
}

static double fun_wrapper1d ( double x){
  double res;
  assert (CheckLinemin1D());

  if ( linemin_noisy){
    UpdateSpinner(linemin_spinner);
  }
  res = linemin_function ( &x, linemin_info);

  return res;
}


static void setf ( double (*f)(const double *, void *), int dim, double * x, double *xnew, double *direct, double * igrad, void * info, int noisy ){
  assert(NULL!=f);
  assert(dim>0);
  assert(NULL!=x);
  assert(NULL!=xnew);
  assert(NULL!=direct);
  assert(NULL!=info);
  assert(NULL==linemin_function);

  linemin_function = f;
  linemin_dim = dim;
  linemin_x = x;
  linemin_xnew = xnew;
  linemin_direct = direct;
  linemin_igrad = igrad;
  linemin_info = info;
  linemin_noisy = noisy;
  if (noisy){
    linemin_spinner = CreateSpinner(2);
  }

  assert(CheckLinemin());
}

static void unsetf ( void){
  assert(CheckLinemin());
  linemin_function = NULL;
  linemin_dim = 0;
  linemin_x = NULL;
  linemin_xnew = NULL;
  linemin_direct = NULL;
  linemin_igrad = NULL;
  linemin_info = NULL;
  if (linemin_noisy){
    DeleteSpinner (linemin_spinner);
  }
  linemin_noisy = 0;
  linemin_spinner = NULL;
}

static double fun_wrapper ( double x){
  int i;
  double res;
  assert(CheckLinemin());

  if ( linemin_noisy){
    UpdateSpinner(linemin_spinner);
  }
  for ( i=0 ; i<linemin_dim ; i++){
    linemin_xnew[i] = linemin_x[i] + x * linemin_direct[i];
  }
  res = linemin_function (linemin_xnew, linemin_info);
  //printf ("f(%e) = %e\n",x,res);
    
  return res;
}

static double qfun_wrapper ( double x){
	if ( linemin_noisy){ UpdateSpinner(linemin_spinner); }
	for ( int i=0 ; i<linemin_dim ; i++){
		linemin_xnew[i] = linemin_x[i] + x * linemin_igrad[i] + x*x * (linemin_direct[i]-linemin_igrad[i]);
  	}
  	return linemin_function (linemin_xnew, linemin_info);
}

static int CheckLinemin ( void){
  assert (linemin_function!=NULL);
  assert (linemin_dim>0);
  assert (linemin_xnew!=NULL);
  assert (linemin_direct!=NULL);
  assert (linemin_info!=NULL);
  assert (linemin_x!=NULL);
  assert (linemin_noisy==0 || linemin_noisy==1);
  if ( linemin_noisy){
    assert (linemin_spinner!=NULL);
  }

  return 1;
}

static int CheckLinemin1D ( void){
  assert (NULL!=linemin_function);
  assert (NULL!=linemin_info);
  assert (linemin_noisy==0 || linemin_noisy==1);
  if ( linemin_noisy){
    assert (NULL!=linemin_spinner);
  }

  return 1;
}

