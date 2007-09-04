double log1p ( const double x);
double expm1 ( const double x);
int finite ( const double x);
double tgamma (double x);

extern unsigned long __ul_nan[];

#define NAN (*(double *) __ul_nan)