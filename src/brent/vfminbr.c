/*
 *=======================================================================
 *			Verify FMINBR routine
 */

#include "math.h"

double fminbr(double a, double b, double (*f)(), double tol);

static int counter;			/* Iteration counter	*/

test(a,b,f,msg)			/* Run a test			*/
double a,b;				/* Range the root is seeked for */
double (*f)(double x);			/* Functiom under examination	*/
char * msg;				/* Explanation message		*/
{
  double minloc;
  counter = 0;
  printf("\nFor function %s\nin [%g,%g] min found is at\t%.9e\n",msg,a,b,
	 (minloc=fminbr(a,b,f,EPSILON)) );
  printf("Min function value found\t%.4e\nNo. of iterations\t\t%d\n",
	 (*f)(minloc), counter);
}

double f1(x)				/* Test of the Forsythe book	*/
double x;
{
  counter++;
  return (powi(x,2)-2.0)*x - 5.0;
}

double f2(x)				/* Modified test 1            	*/
double x;
{
  counter++;
  return powi( (powi(x,2)-2.0)*x - 5.0, 2 );
}

double f3(x)				/* My test                  	*/
double x;
{
  counter++;
  return powi( cos(x) - x,2 ) - 2;
}

double f4(x)				/* My test                  	*/
double x;
{
  counter++;
  return powi( sin(x) - x,2 ) + 1;
}


main()
{
  test(0.0,1.0,f1,"x^3 - 2*x - 5");
  printf("Exact min is at\t\t0.81650\n");

  test(2.0,3.0,f2,"(x^3 - 2*x - 5)^2");
  printf("Exact root is \t\t2.0945514815\n");

  test(2.0,3.0,f3,"(cos(x)-x)^2 - 2");
  test(-1.0,3.0,f3,"(cos(x)-x)^2 - 2");
  test(-1.0,3.0,f4,"(sin(x)-x)^2 + 1");
}
