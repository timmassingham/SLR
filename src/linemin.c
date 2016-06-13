/*
 *  Copyright 2003-2008 Tim Massingham (tim.massingham@ebi.ac.uk)
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
#include <assert.h>

/* Prototype for Brent*/
double brentmin(double lb, const double *flbp, double ub, const double *fubp,
                double x, double *fxp, double (*fun) (const double, void *),
                const double tol, void *info, int *neval);

static double (*linemin_function) (const double *, void *) = NULL;
static void *linemin_info;

static void setf1d(double (*f) (const double *, void *), void *info);
static void unsetf1d(void);
static double fun_wrapper1d(double x, void *info);
static int CheckLinemin1D(void);

/**  Back-tracking approximate line search

Finds a new point where the sufficient decrease criterion is met

@param fun     Pointer to function to be evalutated.
@param finit   Initial function value.
@param ndim    Dimension of input arrays.
@param x       Pointer to initial position.
@param xnew    Pointer to new position (working memory).
@param grad    Pointer to gradient.
@param direct  Pointer to direction to be searched along.
@param info    Pointer passed through to function.
@param step    Initial step size (may be negative).
@param neval   Pointer to nuber of evaluations used.
*/
double linemin_backtrack(double (*fun) (const double *, void *), double finit,
                         int ndim, double *x, double *xnew, const double *grad,
                         const double *direct, void *info, double step,
                         int *neval)
{
    const int niteration = 8;
    const double factor = 0.5;
    const double tol = 1e-4;
    assert(NULL != fun);
    assert(ndim >= 1);
    assert(NULL != x);
    assert(NULL != xnew);
    assert(NULL != direct);
    assert(NULL != info);

    // Criterion for sufficient decrease
    double sufficient = 0.0;
    for (int i = 0; i < ndim; i++) {
        sufficient += direct[i] * grad[i];
    }
    sufficient *= tol;

    double f;
    int it = 0;
    for (it = 0; it < niteration; it++) {
        // Evaluate new point
        for (int i = 0; i < ndim; i++) {
            xnew[i] = x[i] + step * direct[i];
        }
        f = fun(xnew, info);
        *neval = *neval + 1;

        // Check for sufficient decrease
        double deltaf = f - finit;
	if (deltaf < step * sufficient) {
            break;
        }
	step *= factor;
    }

    // Store new value of x
    if (it < niteration){
        for (int i = 0; i < ndim; i++) {
            x[i] = xnew[i];
	}
    }

    return f;
}

double linemin_1d(double (*fun) (const double *, void *), double *x, void *info,
                  const double min, const double max, const double tol,
                  int *neval)
{
    double res, fx;
    assert(NULL != fun);
    assert(NULL != x);
    assert(NULL != info);
    assert(min < max);
    assert(tol > 0.);

    setf1d(fun, info);
    res =
        brentmin(min, NULL, max, NULL, x[0], NULL, fun_wrapper1d, 1e-5, info,
                 neval);
    fx = fun_wrapper1d(res, info);
    *neval = *neval + 1;
    x[0] = res;
    unsetf1d();

    return fx;
}

static void setf1d(double (*fun) (const double *, void *), void *info)
{
    assert(NULL != fun);
    assert(NULL != info);

    linemin_info = info;
    linemin_function = fun;

    assert(CheckLinemin1D());
}

static void unsetf1d(void)
{
    assert(CheckLinemin1D());

    linemin_info = NULL;
    linemin_function = NULL;
}

static double fun_wrapper1d(double x, void *info)
{
    return linemin_function(&x, info);
}

static int CheckLinemin1D(void)
{
    assert(NULL != linemin_function);
    assert(NULL != linemin_info);
    return 1;
}
