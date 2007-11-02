#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>

void warn ( const char * fmt, ...){
	va_list args;

	assert(NULL!=fmt);
	va_start(args,fmt);

	fputs("Warning: ",stderr);
	vfprintf(stderr,fmt,args);
	
	va_end(args);
}


void err ( const int errnum, const char * fmt, ...){
        va_list args;

	assert(NULL!=fmt);
        va_start(args,fmt);

        fprintf(stderr, "Error: ");
        vfprintf(stderr,fmt,args);
        
        va_end(args);
	abort();
}


