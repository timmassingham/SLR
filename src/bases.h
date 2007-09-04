#ifndef _BASES_H_
#define _BASES_H_

int NumberPossibleBases ( const int seq_type, const int gencode);
int CodonAsQcoord (int base, int seqtype, int gencode);
int GapChar(int seqtype);
char NucleoAsChar ( int a);
char AminoAsChar ( int a);
int ToAmino ( char c);
int ToNucleo ( char c);
int IsSeqtype ( const int seqtype);
int IsValidBase ( const int base, const int seqtype, const int gencode);


double *  ConvertCodonFreqsToQcoord ( const double * freqs, const int gencode);


#define	SEQTYPE_NUCLEO		0
#define SEQTYPE_AMINO		1
#define SEQTYPE_CODON		2
#define SEQTYPE_CODONQ		3

#endif
