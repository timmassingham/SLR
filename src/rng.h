#ifndef _RNG_H_
#define _RNG_H_

#define RL_LAGGED	0
#define RL_LINEAR	1
void RL_Init(const unsigned int seed);
void RL_Close(void);

double RandomStandardUniform ( void);
double RandomExp (double mean);
double RandomGamma ( double shape, double scale);
int RandomDirichlet ( double *a, double *p, int n);
int RandomDirichlet_rejection ( double *a, double *p, int n);
void SetRandomGenerator (int gen);
#endif

