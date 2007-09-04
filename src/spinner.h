#ifndef _SPINNER_H_
#define _SPINNER_H_

typedef struct {
	int step;
	int nstates;
	int type;
} SPINNER;

SPINNER * CreateSpinner (int type);
void UpdateSpinner (SPINNER * spin);
void DeleteSpinner (SPINNER * spin);
void CreateIndependentSpinner (int type);
void * IndependentSpinner (void * ptr);
#endif

