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

#include <assert.h>
#include <cblas.h>
#include <err.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "linemin.h"
#include "matrix.h"
#include "utility.h"

#define RESTART		1
#define RESET		100
#define MAX_TRUST	10.0
#define MIN_TRUST	1e-4
#define FACTOR_TRUST    2.0
#define MINSCALE	0.001
#define MAXSCALE	1000.0

#define HESSIAN_NONPD		  1
#define PARAM_BOUND		  2
#define NEWTON			  4
#define TRUNC_BOUND		  8
#define INVALID_STEP		 16
#define BAD_STEP		 32
#define REARRANGED		 64
#define DOG_LEG			128

#define BOUND_TOL	1e-5

#define TEMPFILE 	"parameters.tmp"

typedef struct {
    double *x;
    double *xn;
    double *dx;
    double *dxn;
    double *H;
    double *space;
    double *lb;
    double *ub;
    void *state;
    double fc, fn;
    double (*f) (const double *, void *);
    void (*df) (const double *, double *, void *);

    int *onbound;
    int n;
    int neval;
    double trust;
} OPTOBJ;

struct scaleinfo {
    double *sx;
    int dim;
    void (*df) (const double *, double *, void *);
    double (*f) (const double *, void *);
    double *scale;
    void *state;
};

void
UpdateH_BFGS(double *H, const double *x, double *xn, const double *dx,
             double *dxn, double *scale, const int n, double *space,
             const int *onbound);
double TakeStep(OPTOBJ * opt, const double tol, double *factor, int *newbound);

void
InitializeOpt(OPTOBJ * opt, double *x, int n,
              void (*df) (const double *, double *, void *),
              double (*f) (const double *, void *), double fx, void *data,
              double *bd);
OPTOBJ *NewOpt(int n);
void FreeOpt(OPTOBJ * opt);
void MakeErrString(char **string, int errn);
void InitializeH(OPTOBJ * opt);
void TestIdentity(double *A, double *B, int n);
double TrimAtBoundaries(const double *x, const double *direct,
                        const double *scale, const int n, const double *lb,
                        const double *ub, const int *onbound, int *idx);
int UpdateActiveSet(const double *x, double *grad, const double *scale,
                    double *InvHess, const double *lb, const double *ub,
                    int *onbound, const int n, int *newbound);
double GetNewtonStep(double *direct, const double *InvHess,
                     const double *grad, const int n, const int *onbound);
void ScaledStep(const double factor, const double *x, double *xn,
                const double *direct, const int n);
void MakeMatrixDiagonal(double *A, int n);
int CheckScaleInfo(struct scaleinfo *sinfo);
void dfWrap(const double *x, double *grad, void *info);
double fWrap(const double *x, void *info);
void Rescale(double *x, double *dx, double *H, int n, double *scale);
void AnalyseOptima(double *x, double *dx, int n, int *onbound, double *lb,
                   double *ub);
void check_grad(const char *str, OPTOBJ * opt);
int step, reset;
int errn = 0;

/** Read optimisation parameters from a file

@param filename A filename to read from.
@param opt An OPTOBJ to write the values to.
*/
void read_opt_parameters(const char *filename, OPTOBJ * opt)
{
    assert(NULL != filename);
    assert(NULL != opt);

    FILE *fp = fopen(filename, "r");
    if (NULL == fp) {
        errx(EXIT_FAILURE, "Failed from read from temporary file");
    }
    if ((fread(opt->x, sizeof(double), opt->n, fp) != opt->n) ||
        (fread(opt->dx, sizeof(double), opt->n, fp) != opt->n) ||
        (fread
         (((struct scaleinfo *)opt->state)->scale, sizeof(double), opt->n,
          fp) != opt->n) || (fread(&opt->fc, sizeof(double), 1, fp) != 1)) {
        errx(EXIT_FAILURE, "Failed from read from temporary file");
    }
    fclose(fp);
}

/** Read optimisation parameters from a file

@param filename A filename to write to.
@param opt An OPTOBJ to write the values from.
@returns success or failure
*/
bool write_opt_parameters(const char *filename, const OPTOBJ * opt)
{
    assert(NULL != filename);
    assert(NULL != opt);

    FILE *fp = fopen(filename, "w");
    if (NULL == fp) {
        warnx("Failed to open temporary file");
        return false;
    }

    fwrite(opt->x, sizeof(double), opt->n, fp);
    fwrite(opt->dx, sizeof(double), opt->n, fp);
    fwrite(((struct scaleinfo *)opt->state)->scale, sizeof(double), opt->n, fp);
    fwrite(&opt->fc, sizeof(double), 1, fp);
    fclose(fp);

    return true;
}


double calcerr(const double x, const double y)
{
    return x - y;               /* For absolute errors */
}

void
Optimize(double *x, int n, void (*df) (const double *, double *, void *),
         double (*f) (const double *, void *), double *fx, void *data,
         double *bd, const bool writeTemp, const bool readTemp)
{
    OPTOBJ *opt;
    double fo, fn, tol, md;
    int max_restart, restarts;
    char *errstring = NULL;
    double fact;
    double *scale;
    int newbound = 1;
    bool tempOk = true;

    opt = NewOpt(n);
    if (NULL == opt) {
        return;
    }
    InitializeOpt(opt, x, n, df, f, *fx, data, bd);
    if (readTemp) {
        read_opt_parameters(TEMPFILE, opt);
    }

    tol = 3e-8;
    max_restart = 20;
    restarts = -1;

    /* Do optimization, allowing restarts so don't get bogged down. */

    printf("Initial\tf: %8.6f\nStep     f(x)      delta\n", opt->fc);
    do {
        fo = opt->fc;
        InitializeH(opt);
        fact = 1.;
        do {
            fn = opt->fc;
            errn = 0;
            md = TakeStep(opt, tol, &fact, &newbound);
            MakeErrString(&errstring, errn);
            step++;
            printf("%3d: %9f %10.5e %4d %s\t%9.3f\n", step,
                   opt->fc, fabs(opt->fc - fn), opt->neval, errstring, md);

            // Write temporary values to file
            if (writeTemp && tempOk) {
                tempOk = write_opt_parameters(TEMPFILE, opt);
            }

        }
        while ((calcerr(fn, opt->fc) > tol) || newbound);

        printf("***\n");
        restarts++;
    }
    while (restarts < max_restart && (calcerr(opt->fc, fo) > tol)
           && RESTART);

    if (restarts == max_restart) {
        printf
            ("Didn't converge after %d restarts. Returning best value.\n",
             restarts);
    }
    scale = ((struct scaleinfo *)opt->state)->scale;
    for (int i = 0; i < n; i++) {
        x[i] = opt->x[i] * scale[i];
        opt->dx[i] *= scale[i];
    }
    *fx = opt->fc;

    FreeOpt(opt);
}

void MakeErrString(char **str, int errn)
{
    int nerr = 0;
    char *string;

    string = *str;

    if (NULL != string)
        free(string);

    if ((errn & HESSIAN_NONPD) == HESSIAN_NONPD)
        nerr++;
    if ((errn & PARAM_BOUND) == PARAM_BOUND)
        nerr++;
    if ((errn & NEWTON) == NEWTON)
        nerr++;
    if ((errn & TRUNC_BOUND) == TRUNC_BOUND)
        nerr++;
    if ((errn & INVALID_STEP) == INVALID_STEP)
        nerr++;
    if ((errn & BAD_STEP) == BAD_STEP)
        nerr++;
    if ((errn & REARRANGED) == REARRANGED)
        nerr++;
    if ((errn & DOG_LEG) == DOG_LEG)
        nerr++;

    string = malloc((nerr + 1) * sizeof(char));

    nerr = 0;
    if ((errn & HESSIAN_NONPD) == HESSIAN_NONPD)
        string[nerr++] = '-';
    if ((errn & PARAM_BOUND) == PARAM_BOUND)
        string[nerr++] = 'B';
    if ((errn & NEWTON) == NEWTON)
        string[nerr++] = 'N';
    if ((errn & TRUNC_BOUND) == TRUNC_BOUND)
        string[nerr++] = 'T';
    if ((errn & INVALID_STEP) == INVALID_STEP)
        string[nerr++] = 'V';
    if ((errn & BAD_STEP) == BAD_STEP)
        string[nerr++] = 'W';
    if ((errn & REARRANGED) == REARRANGED)
        string[nerr++] = 'R';
    if ((errn & DOG_LEG) == DOG_LEG)
        string[nerr++] = 'D';
    string[nerr] = '\0';

    *str = string;
}

OPTOBJ *NewOpt(int n)
{
    OPTOBJ *opt;

    if (n < 1)
        return NULL;

    opt = malloc(sizeof(OPTOBJ));
    if (NULL == opt)
        return NULL;

    opt->n = n;
    opt->x = malloc(n * sizeof(double));
    opt->xn = malloc(n * sizeof(double));
    opt->dx = malloc(n * sizeof(double));
    opt->dxn = malloc(n * sizeof(double));
    opt->lb = malloc(n * sizeof(double));
    opt->ub = malloc(n * sizeof(double));
    opt->H = malloc(n * n * sizeof(double));
    opt->space = malloc(4 * n * sizeof(double));
    opt->onbound = malloc(n * sizeof(int));
    opt->neval = 0;
    opt->trust = 0.1;

    if (NULL == opt->x || NULL == opt->xn || NULL == opt->dx
        || NULL == opt->dxn || NULL == opt->lb || NULL == opt->ub
        || NULL == opt->H || NULL == opt->space || NULL == opt->onbound) {
        FreeOpt(opt);
        opt = NULL;
    }
    return opt;
}

void FreeOpt(OPTOBJ * opt)
{
    struct scaleinfo *sinfo;
    free(opt->x);
    free(opt->xn);
    free(opt->dx);
    free(opt->dxn);
    free(opt->lb);
    free(opt->ub);
    free(opt->H);
    free(opt->space);
    free(opt->onbound);
    sinfo = (struct scaleinfo *)opt->state;
    free(sinfo->sx);
    free(sinfo->scale);
    free(sinfo);
    free(opt);
}

void
InitializeOpt(OPTOBJ * opt, double *x, int n,
              void (*df) (const double *, double *, void *),
              double (*f) (const double *, void *), double fx, void *data,
              double *bd)
{
    struct scaleinfo *sinfo;

    step = 0;

    for (int i = 0; i < n; i++) {
        opt->x[i] = x[i];
        opt->lb[i] = bd[i];
        opt->ub[i] = bd[n + i];
        opt->onbound[i] = 0;
    }

    sinfo = malloc(sizeof(struct scaleinfo));
    sinfo->dim = n;
    sinfo->f = f;
    sinfo->df = df;
    sinfo->state = data;
    sinfo->sx = malloc(n * sizeof(double));
    sinfo->scale = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++)
        sinfo->scale[i] = 1.;

    opt->df = dfWrap;
    opt->f = fWrap;
    opt->state = sinfo;
    opt->f(opt->x, opt->state);
    opt->neval++;
    opt->df(opt->x, opt->dx, opt->state);
    opt->neval++;

    for (int i = 0; i < n; i++)
        if ((opt->x[i] <= opt->lb[i] && opt->dx[i] >= 0.)
            || (opt->x[i] >= opt->ub[i] && opt->dx[i] <= 0.))
            opt->onbound[i] = 1;
    opt->fc = fx;

}


double TakeStep(OPTOBJ * opt, const double tol, double *factor, int *newbound)
{
    int idx, neval=0;

    double * direct = opt->space;
    double * space = opt->space + opt->n;

    *newbound = 0;
    double norm = GetNewtonStep(direct, opt->H, opt->dx, opt->n, opt->onbound);
    if (norm > opt->trust) {
        scale_vector(direct, opt->n, opt->trust / norm);
    }

    *factor = 1.;
    double maxfactor =
        TrimAtBoundaries(opt->x, direct,
                         ((struct scaleinfo *)opt->state)->scale, opt->n,
                         opt->lb, opt->ub, opt->onbound, &idx);
    maxfactor = (maxfactor > 1.0) ? 1.0 : maxfactor;
    
    for (int i = 0; i < opt->n; i++) {
        opt->xn[i] = opt->x[i];
    }
    opt->fn = linemin_backtrack(opt->f, opt->fc, opt->n, opt->xn, space, opt->dx, direct, opt->state, maxfactor, &neval);
    opt->neval += neval;

    if( neval == 1){
        opt->trust *= FACTOR_TRUST;
        opt->trust = (opt->trust >= MAX_TRUST) ? MAX_TRUST : opt->trust;
    } else if( neval > 2){
        opt->trust /= FACTOR_TRUST;
        opt->trust = (opt->trust <= MIN_TRUST) ? MIN_TRUST : opt->trust;
    }

    opt->df(opt->xn, opt->dxn, opt->state);
    opt->neval++;
    for (int i = 0; i < opt->n; i++) {
        direct[i] = -opt->dxn[i];
    }
    UpdateH_BFGS(opt->H, opt->x, opt->xn, opt->dx, opt->dxn,
                 ((struct scaleinfo *)opt->state)->scale, opt->n, space,
                 opt->onbound);
    UpdateActiveSet(opt->xn, direct,
                    ((struct scaleinfo *)opt->state)->scale, opt->H,
                    opt->lb, opt->ub, opt->onbound, opt->n, newbound);

    opt->fc = opt->fn;
    for (int i = 0; i < opt->n; i++) {
        opt->x[i] = opt->xn[i];
        opt->dx[i] = opt->dxn[i];
    }

    norm = 0.;
    for (int i = 0; i < opt->n; i++)
        if (!opt->onbound[i])
            norm += opt->dx[i] * opt->dx[i];
    return sqrt(norm);
}

double
TrimAtBoundaries(const double *x, const double *direct,
                 const double *scale, const int n, const double *lb,
                 const double *ub, const int *onbound, int *idx)
{
    double bound, maxerr, epserr;
    double maxfact, fact;

    assert(NULL != x);
    assert(NULL != direct);
    assert(NULL != scale);
    assert(n > 0);
    assert(NULL != lb);
    assert(NULL != ub);
    assert(NULL != onbound);

    maxfact = DBL_MAX;
    maxerr = 0.;
    for (int i = 0; i < n; i++) {
        if (!onbound[i] && fabs(direct[i]) > DBL_EPSILON) {
            bound = ((direct[i] > 0.) ? ub[i] : lb[i]) / scale[i];
            fact = (bound - x[i]) / direct[i];
            epserr = (fabs(bound) + fabs(x[i])) / fabs(direct[i]);
            assert(fact > 0.);
            if (fact < maxfact) {
                *idx = i;
                maxfact = fact;
                maxerr = epserr;
            }
        }
    }

    /* Correction for rounding to ensure bounds are never exceeded */
    maxfact -= (maxerr + maxfact) * DBL_EPSILON;

    for (int i = 0; i < n; i++) {
        assert(onbound[i]
               || (((x[i] + maxfact * direct[i]) * scale[i] - lb[i] >= 0.)
                   && ((x[i] + maxfact * direct[i]) * scale[i] - ub[i] <= 0.)));
    }

    return maxfact;
}

int
UpdateActiveSet(const double *x, double *direct, const double *scale,
                double *InvHess, const double *lb, const double *ub,
                int *onbound, const int n, int *newbound)
{
    int nremoved = 0;
    bool inverted = false;

    assert(NULL != x);
    assert(NULL != direct);
    assert(NULL != scale);
    assert(NULL != lb);
    assert(NULL != ub);
    assert(NULL != onbound);
    assert(NULL != InvHess);
    assert(n > 0);
    for (int i = 0; i < n; i++) {
        assert(finite(x[i]));
        assert(finite(direct[i]));
        for (int j = 0; j < n; j++)
            assert(finite(InvHess[i * n + j]));
    }

    /*Check boundaries */
    for (int i = 0; i < n; i++) {
        if ((x[i] * scale[i] - lb[i] < BOUND_TOL)
            || (ub[i] - x[i] * scale[i] < BOUND_TOL)) {
            if (onbound[i] == 0) {
                /* Newly on boundary. Modify Hessian */
                nremoved++;
                errn = errn | PARAM_BOUND;
		if( ! inverted){
		    InvertMatrix(InvHess, n);
		    inverted = true;
		}
                double diag = InvHess[i * n + i];
                for (int j = 0; j < n; j++) {
                    InvHess[i * n + j] = 0.;
                    InvHess[j * n + i] = 0.;
                }
                InvHess[i * n + i] = 1.0; //fabs(diag);
	    }
	}
	onbound[i] = 0;
        if ((x[i] * scale[i] - lb[i] < BOUND_TOL && direct[i] <= 0.)
            || (ub[i] - x[i] * scale[i] < BOUND_TOL && direct[i] >= 0.)) {
            onbound[i] = 1;
            direct[i] = 0.;
	}
    }

    if (inverted) {
        InvertMatrix(InvHess, n);
    }
    *newbound += nremoved;
    return nremoved;
}

double
GetNewtonStep(double *direct, const double *InvHess,
              const double *grad, const int n, const int *onbound)
{
    assert(NULL != direct);
    assert(NULL != InvHess);
    assert(NULL != grad);
    assert(NULL != onbound);
    assert(n > 0);

    double * g = malloc(n * sizeof(double));
    memcpy(g, grad, n * sizeof(double));
    for ( int i = 0 ; i < n ; i++){
        if( onbound[i] ){
	    g[i] = 0.0;
	}
    }
    cblas_dsymv(CblasColMajor, CblasLower, n, -1.0, InvHess, n, g, 1, 0.0, direct, 1);

    double norm = 0.0;
    for (int i = 0; i < n; i++) {
	if(onbound[i]){
	    direct[i] = 0.0;
	}
        norm += direct[i] * direct[i];
    }

    free(g);

    return sqrt(norm);
}

void
ScaledStep(const double factor, const double *x, double *xn,
           const double *direct, const int n)
{
    assert(NULL != x);
    assert(NULL != xn);
    assert(NULL != direct);
    assert(n > 0);

    for (int i = 0; i < n; i++)
        xn[i] = x[i] + direct[i] * factor;
}

/* Space n+n+n */
void 
UpdateH_BFGS(double *H, const double *x, double *xn, const double *dx,
             double *dxn, double *scale, const int n, double *space,
             const int *onbound)
{
    const double update_tol = 1e-4;
    double gd = 0., *Hg, gHg = 0.;
    double *g, *d;

    g = space;
    space += n;
    d = space;
    space += n;
    Hg = space;
    space += n;

    gd = 0.;
    memset(d, 0, n * sizeof(double));
    memset(g, 0, n * sizeof(double));
    for (int i = 0; i < n; i++) {
        if (!onbound[i]) {
            d[i] = xn[i] - x[i];
            g[i] = dxn[i] - dx[i];
            gd += g[i] * d[i];
        }
    }

    if (gd <= update_tol) {
        errn = errn | HESSIAN_NONPD;
        MakeMatrixIdentity(H, n);
        return;
    }


    cblas_dsymv(CblasColMajor, CblasLower, n, 1.0, H, n, g, 1, 0.0, Hg, 1);
    gHg = 0.;
    for (int i = 0; i < n; i++) {
        gHg += g[i] * Hg[i];
    }

    /*
     * Update inverse hessian using BFGS. Using fact that InvHessian is
     * symmetric, so only need to do half the matrix.
     */
    const double f = 1. + gHg / gd;
    cblas_dsyr(CblasColMajor, CblasLower, n, f / gd, d, 1, H, n);
    cblas_dsyr2(CblasColMajor, CblasLower, n, -1.0 / gd, d, 1, Hg, 1, H, n); 

    Rescale(xn, dxn, H, n, scale);
}

void InitializeH(OPTOBJ * opt)
{
    MakeMatrixIdentity(opt->H, opt->n);
    reset = 1;
}

void TestIdentity(double *A, double *B, int n)
{

    double tmp;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            tmp = 0.;
            for (int k = 0; k < n; k++)
                tmp += A[i * n + k] * B[k * n + j];
            if (i == j && fabs(1. - tmp) > 1e-8)
                printf("%d %d = %e\n", i, j, tmp);
            else if (i != j && fabs(tmp) > 1e-8)
                printf("%d %d = %e\n", i, j, tmp);
        }
}

double fWrap(const double *x, void *info)
{
    struct scaleinfo *sinfo;
    double fx;
    assert(NULL != x);
    assert(NULL != info);

    sinfo = (struct scaleinfo *)info;
    assert(CheckScaleInfo(sinfo));

    for (int i = 0; i < sinfo->dim; i++) {
        sinfo->sx[i] = x[i] * sinfo->scale[i];
    }

    fx = sinfo->f(sinfo->sx, sinfo->state);

    return fx;
}

void dfWrap(const double *x, double *grad, void *info)
{
    struct scaleinfo *sinfo;
    assert(NULL != x);
    assert(NULL != info);

    sinfo = (struct scaleinfo *)info;
    assert(CheckScaleInfo(sinfo));

    for (int i = 0; i < sinfo->dim; i++) {
        sinfo->sx[i] = x[i] * sinfo->scale[i];
    }

    sinfo->df(sinfo->sx, grad, sinfo->state);

    for (int i = 0; i < sinfo->dim; i++) {
        grad[i] *= sinfo->scale[i];
    }

    return;
}

int CheckScaleInfo(struct scaleinfo *sinfo)
{
    assert(NULL != sinfo);
    assert(sinfo->dim > 0);
    assert(NULL != sinfo->sx);
    assert(NULL != sinfo->df);
    assert(NULL != sinfo->f);
    assert(NULL != sinfo->scale);
    assert(NULL != sinfo->state);

    return 1;
}

void Rescale(double *x, double *dx, double *H, int n, double *scale)
{
    double scalefact;

    assert(NULL != x);
    assert(NULL != dx);
    assert(NULL != H);
    assert(n > 0);
    assert(NULL != scale);

    for (int i = 0; i < n; i++) {
        scalefact = sqrt(H[i * n + i]);
        scalefact = (scalefact > MINSCALE) ? scalefact : MINSCALE;
        scalefact = (scalefact < MAXSCALE) ? scalefact : MAXSCALE;
        for (int j = 0; j < n; j++) {
            H[i * n + j] /= scalefact;
            H[j * n + i] /= scalefact;
        }
        x[i] /= scalefact;
        dx[i] *= scalefact;
        scale[i] *= scalefact;
    }
}

void
AnalyseOptima(double *x, double *dx, int n, int *onbound, double *lb,
              double *ub)
{
    assert(NULL != x);
    assert(NULL != dx);
    assert(n > 0);
    assert(NULL != onbound);
    assert(NULL != lb);
    assert(NULL != ub);

    for (int i = 0; i < n; i++) {
        printf("%4d: ", i);
        if (onbound[i]) {
            printf("Boundary. lb: %e, ub: %e. Grad %e\n",
                   x[i] - lb[i], ub[i] - x[i], dx[i]);
        } else {
            printf("Not on boundary. lb: %e, ub: %e, grad = %e\n",
                   x[i] - lb[i], ub[i] - x[i], dx[i]);
        }
    }
}

void check_grad(const char *str, OPTOBJ * opt)
{
    assert(NULL != opt);
    puts(str);
    for (int i = 0; i < opt->n; i++) {
        double oldx = opt->x[i];
        opt->x[i] += 1e-5;
        double fp = opt->f(opt->x, opt->state);
        opt->x[i] -= 2e-5;
        double fm = opt->f(opt->x, opt->state);
        opt->x[i] = oldx;
        printf("%d:\t%e\t%e\t%e\n", i, (fp - fm) / 2e-5, opt->dx[i],
               fabs((fp - fm) / 2e-5 - opt->dx[i]));
    }
}
