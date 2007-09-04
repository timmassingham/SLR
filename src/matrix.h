#ifndef _MATRIX_H_
#define _MATRIX_H_

void Matrix_Matrix_Mult ( const double * A, const int nr1, const int nc1, const double * B, const int nr2, const int nc2, double *C);
void Matrix_MatrixT_Mult ( const double * A, const int nr1, const int nc1, const double * B, const int nr2, const int nc2, double *C);
void MatrixT_Matrix_Mult ( const double * A, const int nr1, const int nc1, const double * B, const int nr2, const int nc2, double *C);

void NormalizeColumns ( double * mat, int n);
void NormalizeRows ( double * mat, int n);
void TransposeMatrix (double *a, double *b, int n);
double MatrixFMax ( double * A, int n);
double VectorNorm ( double * A, int n);
double VectorDotProduct ( const double * A, const double * B, const int n);
void GramSchmidtTranspose ( double * A,int n);
void CopyMatrix ( double *A, double *B, int n);
int Factorize ( double * A, double * val, int n);
void HadamardMult ( double * A, double * B, int n);
void MakeMatrixIdentity (double * mat, const int n);
void MakeMatrixDiagonal(double *A, const int n);
double MatrixMaxElt ( double * A, int n);
double MatrixMinElt (double * A, int n);
int InvertMatrix (double * A, int n);

int IsFiniteVector (const double * a, const int n);
int IsZeroVector ( const double *a , const int n);
#endif
