#ifndef _OPTIMIZE_H_
#define _OPTIMIZE_H_
/* Constants needed for optimizer:
 * 	opt->trust	Initial trust region, default 1
 * 	tol		Tolerance, default 1e-8
 * 	RESTART		Restart optimizer at optima.
 * 	max_step	Maximum number of steps to take.
 */

/* Suggested improvements:
 * 	Need testing in more extreme cases.
 * 	Termination criteria need to be worked on. Current look at last step
 * and norm of current gradient.
 * 	Possibility of incorrect termination in a saddle-point? May only fail
 * if we land exactly in the saddle point as Hessian is guaranteed to be
 * positive definite
 */

void Optimize (  double * x, int n, void (*df)(const double *,double *, void *), double (*f)(const double *, void*), double * fx, void * data, double *bd, int noisy);

#endif
