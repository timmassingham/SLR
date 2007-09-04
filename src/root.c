#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>



int sign (double x){
   if ( x<0 ) return -1;
   if ( x>0 ) return 1;
   return 0;
}
/*  Root finding via Ridders' method. 
 */
double find_root ( double min, double max, double (*f)(const double*,void*), void *info, double *fmin, double *fmax, const double tol, int * neval ){
   double fl,fu,fc,fx;
   double c,x,disc;
   double delta;

   assert (max>=min);
   assert (NULL!=f);
   assert (NULL!=info);
   assert (tol>0.);

   if ( NULL!=fmin){ fl = *fmin;}
   else { fl = f(&min,info); *neval = *neval + 1;}
   if ( NULL!=fmax){ fu = *fmax;}
   else { fu = f(&max,info); *neval = *neval + 1;}

   if ( fabs(fl)<tol || fabs(fu)<tol) {
      return ( (fabs(fl)<fabs(fu))?min:max );
   }

   while ( (fabs(fu)>tol || fabs(fl)>tol) && fabs(max-min)>0.5*3e-8*fabs(min+max)){
      assert (sign(fl)!=sign(fu)); /* Ensure that [min,max] brackets a root*/
      c = (min+max)/2.;
      fc = f(&c,info);	*neval = *neval + 1;
      if ( 0. == fabs(fc)) return c;
      disc = sqrt (fc*fc-fl*fu);
      if ( disc == 0.){ return c;} /* Should only occur if fc=fl=fu=0 since fu*fl strictly positive*/
      x = c + (c-min) * fc*sign(fl-fu)/disc;
      fx = f(&x,info);  *neval = *neval + 1;
      if ( 0. == fabs(fx)) return x;

      /*  Sort out new brackets -- not necessarily respecting order of min and max */
      /*  This deviates from the algorithm as described in numerical recipes, since both
       * newly calculated values are taken into account when choosing new bracket. This
       * guarantees that the method will converge at worst half the speed of bisection
       */
      if ( sign(fx) != sign(fc) ){		/* Root is between new points, both of which */
	 min = x;	fl = fx;		/* lie between min and max                   */
	 max = c;	fu = fc;
      } else if ( sign(fx) == sign(fl) ){			/* New lower bound. Choose x */
	 if ( fabs(max-x)<fabs(max-c) ) { min = x; fl = fx;} 	/* or c to minimise bracket. */
	 else 				{ min = c; fl = fc;}
      } else if ( sign(fx) == sign(fu) ){			/*  New upper bound. Choose x */
	 if ( fabs(min-x)<fabs(min-c) ) { max = x; fu = fx;}	/* or c to minimise bracket   */
	 else				{ max = c; fu = fc;}
      } else {
	 assert(0.==fabs(fx));  /*  Ensured that fl!=0 and fu!=0, so if here fx=0*/
	 return x;
      }
   }

   /*  Final linear interpolation*/
   delta = (fabs(fu-fl)>fabs(fu)*DBL_EPSILON)?fl/(fu-fl):0.5;
   x = min - (max-min) * delta;
   return x;
}

