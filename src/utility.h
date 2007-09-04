#ifndef _UTILITY_H_
#define _UTILITY_H_

#ifndef _STDIO_H_
#include <stdio.h>
#endif

void slrwarn ( int i, char * s);
void PrintMatrix ( const double * m, const int n);
void PrintVector ( const double * x, const int n);
void fprint_vector ( FILE * fp, const char sep, const double * x, const int n);
void fprint_ivector ( FILE * fp, const char sep, const int * x, const int n);
void Free ( void * mem);
int NumberPairs (int i);
void PrintMatrixAsBinary ( double *m, int n);
void PrintMatrixAsSign (double *m, int n);
int UpperTriangularCoordinate ( int i, int j, int n);
int LowerTriangularCoordinate ( int i, int j, int n);
int TriangularCoordinate ( int i, int j, int n);
int ReadVectorFromOpenFile (double * x, int n, FILE * fp);

char * ReadString ( FILE * fp);
char gchar ( FILE * fp);

char * itoa ( int i);

double * norm_vector ( double * vec, const int n, const double sum);
double * scale_vector ( double * vec, const int n, const double fact);

int teeint ( FILE * fp, const char * s, const int i);

#endif
