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

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdbool.h>

void Matrix_Matrix_Mult ( const double * A, const int nr1, const int nc1, const double * B, const int nr2, const int nc2, double *C);
void Matrix_MatrixT_Mult ( const double * A, const int nr1, const int nc1, const double * B, const int nr2, const int nc2, double *C);
void MatrixT_Matrix_Mult ( const double * A, const int nr1, const int nc1, const double * B, const int nr2, const int nc2, double *C);

void NormalizeColumns ( double * restrict mat, size_t n);
void NormalizeRows ( double * restrict mat, size_t n);
void TransposeMatrix (const double * restrict a, double * restrict b, size_t n);
double MatrixFMax ( double * restrict A, size_t n);
double VectorNorm ( double * restrict A, size_t n);
double VectorDotProduct ( const double * restrict A, const double * restrict B, const size_t n);
void GramSchmidtTranspose ( double * restrict A, size_t n);
void CopyMatrix ( const double * restrict A, double * restrict B, size_t n);
int Factorize ( double * A, double * val, size_t n);
void HadamardMult (const double * restrict A, double * restrict B, size_t n);
void MakeMatrixIdentity (double * restrict mat, const size_t n);
void MakeMatrixDiagonal(double * restrict A, const size_t n);
double MatrixMaxElt ( double * restrict A, size_t n);
double MatrixMinElt (double * restrict A, size_t n);
int InvertMatrix (double * A, size_t n);

bool IsFiniteVector (const double * restrict a, const size_t n);
bool IsZeroVector ( const double * restrict a , const size_t n);
#endif
