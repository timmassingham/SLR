#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include "gamma.h"
#ifdef WINDOWS
  #include "math_win.h"
#endif

static double gamma_series ( const double x, const double a);
static double gamma_contfrac ( const double x, const double a);

/*void main ( int argc, char ** argv){
  double x,a;
  double p,q;
  if ( argc!=3){
    exit(EXIT_FAILURE);
  }
  sscanf(argv[1],"%le",&x);
  sscanf(argv[2],"%le",&a);
  

  p = pchisq (x,a,0);
  q = pchisq (x,a,1);
  
  printf ("pchisq(%e,%e) = %e\n",x,a,p);
  printf ("1-pchisq(%e,%e) = %e\n",x,a,q);
  printf ("sum = %e\n",p+q);
}*/

  

static double gamma_series ( const double x, const double a){
  double res,term;
  int i;
  assert(x>=0.);
  assert(a>0.);

  if ( x==0.){
    return 0.;
  }
  if ( x==HUGE_VAL){
    return 1.;
  }

  term = 1. / tgamma (a);
  term /= a;
  term *= exp(-x) * pow(x,a);

  res = term;
  for ( i=1 ; ; i++){
    term *= x;
    term /= a+i;
    res += term;
    if ( fabs(term)<=DBL_EPSILON*fabs(res))
      break;
  }

  return res;
}

static double gamma_contfrac ( const double x, const double a){
  int i;
  double f,c,d,dn,cn,term;
  assert(x>=0.);
  assert(a>0.);

  if ( x==0.){
    return 1.;
  }
  if ( x==HUGE_VAL){
    return 0.;
  }

  /* Manually unrolled first iteration of Lentz's algorithm
   */
  f = 1e-300 + 1./(x+1.-a);
  c = (x+1-a) + 1e300;
  d = 1. / (x+1-a);

  for ( i=1 ; ; i++){
    dn = (x + 2*i+1 - a) - (i*(i-a)) * d;
    if ( dn==0.)
      dn = 1e-300;
    cn = (x + 2*i+1 - a) - (i*(i-a)) / c;
    if ( cn==0.)
      cn = 1e-300;
    dn = 1. / dn;
    term = cn * dn;
    f *= term;
    if ( fabs(1.-term)<DBL_EPSILON )
      break;
    c = cn;
    d = dn;
  }

  f *= exp(-x) * pow(x,a) / tgamma(a);

  return f;
}


double pgamma ( const double x, const double a, const double b, const int tail){
  double res,p,q;
  assert(x>=0.);
  assert(a>0.);
  assert(b>0);
  assert(tail==0 || tail==1);

  /* Use continued fraction approx
   */
  if ( 1.+a<x*b){
    q = gamma_contfrac(x*b,a);
    if ( tail==1)
      res = q;
    else
      res = 1.-q;
  } else {
    p = gamma_series(x*b,a);
    if ( tail==1)
      res = 1.-p;
    else
      res = p;
  }

  return res;
}

double pchisq ( double x, double deg, int tail){
  return pgamma (x/2.,deg/2.,1.,tail);
}
