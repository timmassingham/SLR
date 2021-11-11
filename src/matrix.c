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
#include <float.h>
#include <lapacke.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "utility.h"

void Matrix_Matrix_Mult ( const double * A, const int nr1, const int nc1, const double * B, const int nr2, const int nc2, double * C){
  assert(NULL!=A);
  assert(NULL!=B);
  assert(NULL!=C);
  assert(nc1==nr2);

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nc2, nr1, nc1, 1.0, B, nc2, A, nc1, 0.0, C, nc2);
}

void Matrix_MatrixT_Mult ( const double * A, const int nr1, const int nc1, const double * B, const int nr2, const int nc2, double *C){
  assert(NULL!=A);
  assert(NULL!=B);
  assert(NULL!=C);
  assert(nc1==nr2);

  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, nc2, nr1, nc1 , 1.0, B, nr2, A, nc1, 0.0, C, nc2);
}

void MatrixT_Matrix_Mult ( const double * A, const int nr1, const int nc1, const double * B, const int nr2, const int nc2, double *C){
  assert(NULL!=A);
  assert(NULL!=B);
  assert(NULL!=C);
  assert(nc1==nr2);

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, nc2, nr1, nc1, 1.0, B, nc2, A, nr1 , 0.0, C, nc2);
}



void MakeMatrixIdentity (double * restrict mat, const size_t n){
	assert(NULL!=mat);
	assert(n>0);

	memset(mat,0,n*n*sizeof(double));
	for ( size_t i=0 ; i<n ; i++){
		mat[i*n+i] = 1.;
	}
}

void
MakeMatrixDiagonal(double * restrict A, const size_t n)
{
        assert(NULL != A);
        assert(n > 0);

        for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                        if (i != j)
                                A[i * n + j] = 0.;
                }
                A[i * n + i] = fabs(A[i * n + i]);
        }
}


void CopyMatrix ( const double * restrict A, double * restrict B, size_t n){
	for (size_t i=0 ; i<n ; i++)
		for (size_t j=0 ; j<n ; j++)
			B[i*n+j] = A[i*n+j];
}


void NormalizeRows ( double * restrict mat, size_t n){
	for (size_t i=0 ; i<n ; i++){
		double a = 0.;
		for (size_t j=0 ; j<n ; j++)
			a += mat[i*n+j] * mat[i*n+j];
		a = sqrt(a);
		for (size_t j=0 ; j<n ; j++)
			mat[i*n+j] /= a;
	}
}


void TransposeMatrix (const double * restrict a, double * restrict b, size_t n){
	for (size_t i=0 ; i<n ; i++)
		for (size_t j=0 ; j<n ; j++)
			b[j*n+i] = a[i*n+j];
}

void NormalizeColumns ( double * restrict mat, size_t n){
	for (size_t i=0 ; i<n ; i++){
		double a = 0.;
		for (size_t j=0 ; j<n ; j++)
			a += mat[j*n+i] * mat[j*n+i];
		a = sqrt(a);
		for (size_t j=0 ; j<n ; j++)
			mat[j*n+i] /= a;
	}
}

double MatrixFMax ( double * restrict A, size_t n){
	double a=0.;

	for (size_t i=0 ;i<n ; i++)
		for (size_t j=0 ;j<n ; j++)
			if (i!=j)
			a = (a>fabs(A[i*n+j]))?a:fabs(A[i*n+j]);
	return a;
}

double VectorNorm ( double * restrict A, size_t n){
	double a=0.;

	for (size_t i=0 ;i<n; i++)
		a += A[i]*A[i];
	return sqrt(a);
}

double VectorDotProduct ( const double * restrict A, const double * restrict B, const size_t n){
	double a=0.;

	for (size_t i=0 ; i<n ; i++)
		a += A[i] * B[i];
	return a;
}

void GramSchmidtTranspose ( double * restrict A, size_t n){
	double max=0.;

	
	NormalizeRows (A, n);
	for (size_t i=0 ; i<n ; i++){
		double a = VectorNorm(&A[i*n],n);
		for (size_t j=0 ; j<n ; j++)
			A[i*n+j] /= a;
		for (size_t j=i+1 ; j<n ; j++){
			a = VectorDotProduct(&A[j*n],&A[i*n],n);
			max = (max>fabs(a))?max:fabs(a);
			for (size_t k=0 ; k<n ; k++)
				A[j*n+k] -= a * A[i*n+k];
		}
	}
	if(max>0.1)
		printf ("Large change in GS\n");
}



int Factorize ( double * A, double * val, size_t n){
	int INFO = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'L', n, A, n, val);
	
	return INFO;
}


int InvertMatrix ( double * A, size_t n){
	int INFO;
	lapack_int * ipiv = malloc(n * sizeof(lapack_int));
	INFO = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, A, n, ipiv);
        if (INFO==0){
		LAPACKE_dgetri(LAPACK_COL_MAJOR, n, A, n, ipiv);
        }
	free(ipiv);

        if (INFO!=0)
                return -1;
        return 0;
}


/**	Hadamard multiplication A.B of two square matrices

   Result is stored in B

@param A  Input array
@param B  Input / output array
@param n  Size of matrices
**/
void HadamardMult (const double * restrict A, double * restrict B, size_t n){
    assert(NULL!=A);
    assert(NULL!=B);
    assert(n>0);

    for (size_t i=0 ; i<n*n ; i++){
        B[i] *= A[i];
    }
}

double MatrixMaxElt ( double * restrict A, size_t n){
	double max = -DBL_MAX;
	
	assert (NULL!=A);

	for (size_t i=0 ; i<n ; i++)
		for (size_t j=0 ; j<n ; j++){
			if ( A[i*n+j] >max)
				max = A[i*n+j];
		}

	return max;
}



double MatrixMinElt (double * restrict A, size_t n){
	double min=DBL_MAX;
	
	assert (NULL!=A);

	for (size_t i=0 ; i<n ; i++)
		for (size_t j=0 ; j<n ; j++){
			if ( A[i*n+j] >min)
				min = A[i*n+j];
		}

	return min;
}



bool IsFiniteVector ( const double * restrict a, const size_t n){
  assert(NULL!=a);
  assert(n>0);

  for (size_t i=0 ; i<n ; i++){
    if(!isfinite(a[i]))
      return false;
  }
  return true;
}


bool IsZeroVector ( const double * restrict a, const size_t n){
  bool has_zero = true;
  
  assert(NULL!=a);
  assert(n>0);

  for (size_t i=0 ; i<n ; i++){
    has_zero &= (a[i]==0.);
  }

  return has_zero;
}
