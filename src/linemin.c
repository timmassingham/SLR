#include <stdlib.h>
#include <assert.h>
#include "linemin.h"
#include "spinner.h"

/* Prototype for Brent*/
double brentmin ( double lb, const double * flbp, double ub, const double * fubp, double x, double * fxp, double (*fun)(const double, void *), const double tol, void * info, int * neval);

static double (*linemin_function)(const double *, void *) = NULL;
static void * linemin_info;
static int linemin_noisy;
static SPINNER * linemin_spinner;

static void setf1d ( double (*f)(const double *, void *), void * info, int noisy);
static void unsetf1d (void);
static double fun_wrapper1d ( double x, void * info);
static int CheckLinemin1D (void);



double linemin_backtrack ( double (*fun)(const double *, void *), int dim, double * x, double *xnew, double * direct, void * info, const double min, const double max, const double tol, const int noisy, int * neval ){
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


double linemin_1d ( double (*fun)(const double *, void *), double * x, void * info, const double min, const double max, const double tol, const int noisy, int * neval ){
  double res,fx;
  assert(NULL!=fun);
  assert(NULL!=x);
  assert(NULL!=info);
  assert(min<max);
  assert(tol>0.);
  assert(noisy==0 || noisy==1);

  setf1d (fun,info,noisy);
  res = brentmin (min,NULL,max,NULL,x[0],NULL,fun_wrapper1d,1e-5,info,neval);
  fx = fun_wrapper1d(res,info); *neval = *neval + 1;
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

static double fun_wrapper1d ( double x, void * info){
	return linemin_function(&x,info);
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

