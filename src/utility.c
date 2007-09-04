#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <assert.h>

#ifndef OOM
#define OOM(A) if ( A == NULL){ \
                  printf ("Out of Memory! %s:%d\n",__FILE__,__LINE__); \
                  exit (EXIT_FAILURE); }
#endif


void 
Free(void *mem)
{
	if (NULL != mem)
		free(mem);
}

void 
slrwarn(int i, char *s)
{
#ifdef WARNINGS
	assert(NULL != s);
	if (0 == i) {
		fputs(s, stdout);
	}
#endif
}





void 
PrintMatrix(double *m, int n)
{
	int             i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			printf("%3.2e ", m[i * n + j]);
		printf("\n");
	}
	printf("\n");
}

void 
PrintVector(double *x, int n)
{
	int             i;

	for (i = 0; i < n; i++)
		printf("%6.5e ", x[i]);
	printf("\n");
}


int 
NumberPairs(int n)
{
	if (n < 1)
		return -1;
	return (((n - 1) * n) / 2);
}


void 
PrintMatrixAsBinary(double *m, int n)
{
	int             i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			printf("%d ", (fabs(m[i * n + j]) <= DBL_EPSILON) ? 0 : 1);
		printf("\n");
	}
}

void 
PrintMatrixAsSign(double *m, int n)
{
	int             i, j;
	char            c;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (m[i * n + j] < -DBL_EPSILON)
				c = '-';
			else if (m[i * n + j] > DBL_EPSILON)
				c = '+';
			else
				c = '0';
			printf("%c ", c);
		}
		printf("\n");
	}
}



int 
UpperTriangularCoordinate(int i, int j, int n)
{
	int             a;

	if (i == j)
		return -1;
	/* Set i to be the smallest of the two */
	a = i;
	i = (i > j) ? j : i;
	j = (i == a) ? j : a;
	a = ((n - 1) * (n - 2) - (n - 1 - i) * (n - 2 - i)) / 2 + (j - 1);

	return a;
}

int 
LowerTriangularCoordinate(int i, int j, int n)
{
	int             a;

	if (i == j)
		return -1;
	/* Set i to be the smallest of the two */
	a = i;
	i = (i > j) ? j : i;
	j = (i == a) ? j : a;
	a = (j * (j - 1)) / 2 + i;

	return a;
}

int 
TriangularCoordinate(int i, int j, int n)
{
	int             a;

	if (i == j)
		return -1;
	/* Set i to be the smallest of the two */
	a = i;
	i = (i > j) ? j : i;
	j = (i == a) ? j : a;
	a = ((n) * (n - 1) - (n - i) * (n - 1 - i)) / 2 + j;

	return a;
}


int 
ReadVectorFromOpenFile(double *x, int n, FILE * fp)
{
	int             i, e;

	if (NULL == x || NULL == fp)
		return EOF;

	for (i = 0; i < n; i++) {
		e = fscanf(fp, "%le", &x[i]);
		if (EOF == e)
			return e;
	}

	return i;
}

char 
gchar(FILE * fp)
{
	unsigned char   c;

	while (isalpha(c = getc(fp)) == 0 && c != '-');

	return c;
}



#define STRLEN 20
char           *
ReadString(FILE * fp)
{
	char           *string, *tmp;
	char            c;
	int             a = 0, b;

	do
		c = getc(fp);
	while (isspace(c) != 0);

	string = calloc((size_t) (STRLEN + 1), sizeof(char));
	OOM(string);
	do {
		string[a++] = c;
		if (a % 20 == 0) {
			tmp = calloc((size_t) (a + STRLEN), sizeof(char));
			OOM(tmp);
			for (b = 0; b < a; b++)
				tmp[b] = string[b];
			free(string);
			string = tmp;
		}
	} while ((c = getc(fp)) != '\n' && c != EOF);
	string[a] = '\0';

	tmp = calloc((size_t) (a + 1), sizeof(char));
	OOM(tmp);
	for (b = 0; b < a + 1; b++)
		tmp[b] = string[b];
	free(string);

	return tmp;
}

char           *
itoa(int i)
{
	int             len;
	char           *a;

	len = (int) ceil(log10((double) abs(i + 1)));
	if (i < 0) {
		len++;
	}
	len++;			/* space for null  */
	a = malloc(len * sizeof(char));
	a[len - 1] = '\0';
	if (i < 0) {
		a[0] = '-';
	}
	i = abs(i);
	for (len -= 2; i > 0; i /= 10) {
		a[len--] = (i % 10) + 48;
	}

	return a;
}

#ifdef TEST
int 
main(int argc, char *argv[])
{
	int             i, a;
	char           *b;

	for (a = 1; a < argc; a++) {
		sscanf(argv[a], "%d", &i);
		b = itoa(i);
		printf("itoa(%d) = %s\n", i, b);
		free(b);
	}
}
#endif

/* Scale elements of vector so their sum is sum */
double *
norm_vector(double *vec, const int n, const double sum)
{
	double          total = 0, fact;

	assert(NULL != vec);
	assert(n >= 0);

	for (int i = 0; i < n; i++) {
		total += vec[i];
	}
	fact = sum / total;
	for (int i = 0; i < n; i++) {
		vec[i] *= fact;
	}

	return vec;
}

double * 
scale_vector ( double * vec, const int n, const double fact){
	assert(NULL!=vec);
	assert(n>0);

	for ( int i=0 ; i<n ; i++){
		vec[i] *= fact;
	}

	return vec;
}

