#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include "vec.h"

VEC create_vec ( const unsigned int n){
	assert(n>0);

	VEC v = malloc(sizeof(struct __vec));
	if ( NULL==v){return NULL;}

	v->n = n;
	v->x = malloc(n*sizeof(double));
	if ( NULL==v){ free(v); return NULL;}

	return v;
}

VEC create_zerovec ( const unsigned int n){
	assert(n>0);

        VEC v = malloc(sizeof(struct __vec));
        if ( NULL==v){return NULL;}

        v->n = n;
        v->x = calloc(n,sizeof(double));
        if ( NULL==v){ free(v); return NULL;}

        return v;
}

void initialize_vec ( VEC v, const double val){
	assert(NULL!=v);
	const int len = v->n;
	for ( unsigned int i=0 ; i<len ; i++){
		v->x[i] = val;
	}
}

void free_vec ( const VEC v){
	assert(NULL!=v);
	assert(NULL!=v->x);

	free(v->x);
	free(v);
}
	
double dotproduct_vec ( const VEC a, const VEC b){
	assert(NULL!=a);
	assert(NULL!=b);
	assert(a->n == b->n);

	const int len = a->n;
	double sum = 0.;
	for ( unsigned int i=0 ; i<len ; i++){
		sum += a->x[i] * b->x[i];
	}
	return sum;
}

double norm_vec ( const VEC v){
	double norm = dotproduct_vec(v,v);
	assert(norm>=0.);
	return norm;
}

double minelt_vec ( const VEC v){
	assert(NULL!=v);
	double minelt = DBL_MAX;
	const unsigned int len = v->n;
	for ( unsigned int i=0 ; i<len ; i++){
		if ( minelt>v->x[i]){ minelt = v->x[i];}
	}
	return minelt;
}

double maxelt_vec ( const VEC v){
	assert(NULL!=v);
        double maxelt = -DBL_MAX;
        const unsigned int len = v->n;
        for ( unsigned int i=0 ; i<len ; i++){
                if ( maxelt<v->x[i]){ maxelt = v->x[i];}
        }
	return maxelt;
}

VEC copy_vec ( const VEC v){
	assert(NULL!=v);
	const int len = v->n;
	VEC vcopy = create_vec(len);
	assert(NULL!=vcopy);
	for ( unsigned int i=0 ; i<len ; i++){
		vcopy->x[i] = v->x[i];
	}
	return vcopy;
}

void fprint_vec (FILE * fp, const char * prefix, const char * sep, const char * suffix, const VEC v){
	assert(NULL!=fp);
	assert(NULL!=sep);
	assert(NULL!=v);

	fprintf (fp,"%s%e",prefix,vget(v,0));
	const unsigned int len = vlen(v);
	for ( unsigned int i=1 ; i<len ; i++){
		fprintf (fp, "%s%e",sep,vget(v,i));
	}
	fputs(suffix,fp);
}

void fprint_rvec(FILE * fp, const char * name, const VEC v){
	assert(NULL!=fp);
	assert(NULL!=name);
	assert(NULL!=v);
	char * prefix = malloc((strlen(name)+5)*sizeof(char));
	memcpy(prefix,name,strlen(name)*sizeof(char));
	memcpy(prefix,"<-c(",5*sizeof(char));
	fprint_vec(fp,prefix,",",");\n",v);
	free(prefix);
}

double suma_vec ( const VEC v, const double a){
	assert(NULL!=v);
	const unsigned int len = vlen(v);
	double sum = 0.;
	for ( unsigned int i=0 ; i<len ; i++){
		sum += vget(v,i) + a;
	}

	return sum;
}


double sum_vec ( const VEC v){
        assert(NULL!=v);
        const unsigned int len = vlen(v);
        double sum = 0.;
        for ( unsigned int i=0 ; i<len ; i++){
                sum += vget(v,i);
        }

        return sum;
}

