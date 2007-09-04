#include <assert.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "math_win.h"

unsigned long __ul_nan[2] = {0xffffffff,0x7fffffff};



/* Calculates log(1+x) accurately when x close to 0.
 */
double log1p ( const double x){
  double res;

  assert(x>-1.);
  
  if ( fabs(x)>sqrt(DBL_EPSILON)){
    res = log(1.+x);
  } else {
    res = x * ( 1. - x/2.);
  }

  return res;   
} 


/*  Calculates exp(x)-1 accurately when x is close to zero
 */
double expm1 ( const double x){
  double res;

  if ( fabs(x)>sqrt(DBL_EPSILON)){
    res = exp(x)-1.;
  } else {
    res = x * ( 1. + x/2.);
  }

  assert (res>=0.);
  return res;
}


int finite ( const double x){
  if ( -HUGE_VAL==x ||HUGE_VAL == x || NAN==x){
    return 0;
  }
  
  return 1;
}

double tgamma ( double x) {
	static double c[6] = {76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
	double cum,tmp,y; 
	int j;
	
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	y=x+1.;
	cum = 1.000000000190015;
	for ( j=0 ; j<6 ; j++){
		cum += c[j]/y;
		y++;
	}
	
	return exp(-tmp+log(2.5066282746310005*cum/x));
}
