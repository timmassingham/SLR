#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "spinner.h"


int maxtype = 3;
int SpinNType[3] = {4,4,5};
char spin1states[4] = {'\\','|','/','-'};
char spin2states[4] = {'.','o','O','o'};
char spin3states[5] = {'.',':','+','*','0'};
char * SpinStates[3] = {spin1states,spin2states,spin3states};

SPINNER * CreateSpinner ( int type){
	SPINNER * spin;

	if(type<0 || type>=maxtype)
		return NULL;

	spin = malloc (sizeof(SPINNER));
	if (NULL!=spin){
		spin->type = type;
		spin->nstates = SpinNType[type];
		spin->step = 0;
	}

	return spin;
}

void UpdateSpinner ( SPINNER * spin){
	if(NULL==spin)
		return;

	if(0!=spin->step)
		printf ("\b");
	putchar (SpinStates[spin->type][(spin->step)%(spin->nstates)]);
	fflush(stdout);
	spin->step++;
}

void DeleteSpinner (SPINNER * spin){
	if(NULL==spin)
		return;
	if(0!=spin->step){
		putchar('\b');
		fflush(stdout);
	}
	free(spin);
}

