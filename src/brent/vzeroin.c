/*
 *=======================================================================
 *			Verify ZEROIN routine
 */

#include "math.h"

double zeroin(double ax, double bx, double (*f)(), double tol);

static int counter;			/* Iteration counter	*/

test(a,b,f,msg)			/* Run a test			*/
double a,b;				/* Range the root is seeked for */
double (*f)(double x);			/* Functiom under examination	*/
char * msg;				/* Explanation message		*/
{
  double root;
  counter = 0;
  printf("\nFor function %s\nin [%g,%g] root is\t%.9e\n",msg,a,b,
	 (root=zeroin(a,b,f,0.0)) );
  printf("Function value at the root found\t%.4e\nNo. of iterations\t\t%d\n",
	 (*f)(root), counter);
}

double f1(x)				/* Test of the Forsythe book	*/
double x;
{
  counter++;
  return (pow(x,2)-2.0)*x - 5.0;
}

double f2(x)				/* My test                  	*/
double x;
{
  counter++;
  return cos(x) - x;
}

double f3(x)				/* My test                  	*/
double x;
{
  counter++;
  return sin(x) - x;
}


main()
{
  test(2.0,3.0,f1,"x^3 - 2*x - 5");
  printf("Exact root is \t\t2.0945514815\n");

  test(2.0,3.0,f2,"cos(x)-x");
  test(-1.0,3.0,f2,"cos(x)-x");
  test(-1.0,3.0,f3,"sin(x)-x");
}
