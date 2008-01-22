/*
 *  Copyright 2003-2007 Tim Massingham (tim.massingham@ebi.ac.uk)
 *  Funded by EMBL - European Bioinformatics Institute
 */
/*
 *  This file is part of SLR ("Sitewise Likelihood Ratio")
 *
 *  SLR is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SLR is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SLR.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _BASES_H_
#define _BASES_H_

#ifndef _STDBOOL_H_
#include <stdbool.h>
#endif

enum SEQTYPE { SEQTYPE_NULL, SEQTYPE_NUCLEO, SEQTYPE_NUCLEOAMBIG, SEQTYPE_AMINO, SEQTYPE_AMINOAMBIG, SEQTYPE_CODON, SEQTYPE_CODONAMBIG, SEQTYPE_CODONQ};
struct triple_nuc { int fst,snd,trd;};

int NumberPossibleBases ( const enum SEQTYPE seqtype, const int gencode);
int GapChar(const enum SEQTYPE seqtype);
char char_from_base (const int base, const enum SEQTYPE seqtype);
int base_from_char (const char c, const enum SEQTYPE seqtype);
int aminoambig_from_char ( char c);
int amino_from_char( char c);
int nucleo_from_char( char c);
int nucleoambig_from_char ( char c);
bool IsSeqtype ( const enum SEQTYPE seqtype);
bool IsValidBase ( const int base, const enum SEQTYPE seqtype, const int gencode);
bool is_ambiguous ( const int base, const enum SEQTYPE seqtype );
bool is_ambiguous_seqtype (const enum SEQTYPE seqtype);
unsigned int ambiguity_char (const enum SEQTYPE seqtype);
unsigned int nonambig_codon_from_codon (const unsigned int codon);
enum SEQTYPE ambig_seqtype_from_seqtype(const enum SEQTYPE seqtype);
enum SEQTYPE nonambig_seqtype_from_seqtype(const enum SEQTYPE seqtype);
int codon_from_nucs ( const int base1, const int base2, const int base3);
int codon_from_triplenuc (struct triple_nuc nucs);
struct triple_nuc triplenuc_from_codon ( const int codon );
struct triple_nuc triplenuc_from_codonambig ( const int codon);
unsigned int codon_from_codonambig ( const unsigned int codon);
int nonambig_base_from_base (const int base, const enum SEQTYPE seqtype);


double *  ConvertCodonFreqsToQcoord ( const double * freqs, const int gencode);

void (*pp_base(const enum SEQTYPE seqtype))(FILE *, const unsigned int, const enum SEQTYPE);
void fprint_base ( FILE * fp, const unsigned int base, const enum SEQTYPE seqtype);
void fprint_codon ( FILE * fp, const unsigned int codon, const enum SEQTYPE seqtype);
void fprint_codonq (FILE * fp, const unsigned int codonq, const enum SEQTYPE seqtype);

#endif
