#ifndef _VEC_H_
#define _VEC_H_

struct __vec {
	double * x;
	unsigned int n;
};

typedef struct __vec * VEC;

VEC create_vec ( const unsigned int n);
VEC create_zerovec ( const unsigned int n);
void initialize_vec ( VEC v, const double val);
void free_vec ( const VEC v); 
double dotproduct_vec ( const VEC a, const VEC B);
double norm_vec ( const VEC v);
double minelt_vec ( const VEC v);
double maxelt_vec ( const VEC v);
VEC copy_vec ( const VEC v);
void fprint_vec (FILE * fp, const char * prefix, const char * sep, const char * suffix, const VEC v);
void fprint_rvec(FILE * fp, const char * name, const VEC v);

#define vget(VEC,I) VEC->x[(I)]
#define vset(VEC,I,VAL) VEC->x[(I)] = (VAL)
#define vlen(VEC) VEC->n

#endif
