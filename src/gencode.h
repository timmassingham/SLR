#ifndef _GENCODE_H_
#define _GENCODE_H_
int IsValidGencode ( const int gencode);
int CodonToAmino (const int codon, const int gencode);
int IsStop (const int codon, const int gencode);
int CodonToQcoord ( const int codon, const int gencode);
int NumberSenseCodonsInGenCode ( const int gencode);
int IsNonSynonymous ( const int codon1, const int codon2, const int gencode);
int HasTransition ( const int codon1, const int codon2);
int NumberNucChanges ( int codon1, int codon2);
int QcoordToCodon ( const int qcoord, const int gencode);
int QcoordToAmino ( const int qcoord, const int gencode);
int GetGeneticCode ( const char * gencode_str);
int FourfoldDegenerate ( const int codon, const int gencode);
int Degeneracy ( const int codon, const int gencode);

double * GetAminoFrequencies ( const double * codons, const int gencode);


#define GENCODE_UNIVERSAL               0
#define GENCODE_MAMMALIAN_MITOCHONDRIAL 1
#endif
