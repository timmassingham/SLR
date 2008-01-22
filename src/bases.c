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

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <err.h>
#include "gencode.h"
#include "bases.h"

char AminoAmbig[23] = "ARNDCQEGHILKMFPSTWYVX-";
char NucleoAmbig[7] = "ACGTN-";
char Amino[22] = "ARNDCQEGHILKMFPSTWYV-";
char Nucleo[6] = "ACGT-";


void fprint_base ( FILE * fp, const unsigned int base, const enum SEQTYPE seqtype){
	fputc(char_from_base(base,seqtype),stdout);
}

void fprint_codon ( FILE * fp, const unsigned int codon, const enum SEQTYPE seqtype){
	struct triple_nuc nucs;
	enum SEQTYPE nuc_seqtype;
	switch(seqtype){
	case SEQTYPE_CODON:
		nucs = triplenuc_from_codon (codon);
		nuc_seqtype = SEQTYPE_NUCLEO;
		break;
	case SEQTYPE_CODONAMBIG:
		nucs = triplenuc_from_codonambig (codon);
		nuc_seqtype = SEQTYPE_NUCLEOAMBIG;
		break;
	default:
		err(EXIT_FAILURE,"Called print_codon with non-codon seqtype %d\n");
	}
	fprint_base(stdout,nucs.fst, nuc_seqtype);
	fprint_base(stdout,nucs.snd, nuc_seqtype);
	fprint_base(stdout,nucs.trd, nuc_seqtype);
}

void fprint_codonq (FILE * fp, const unsigned int codonq, const enum SEQTYPE seqtype){
	assert(SEQTYPE_CODONQ==seqtype);
	fprintf(stdout,"%3d",codonq);
}

void (*pp_base(const enum SEQTYPE seqtype))(FILE *, const unsigned int, const enum SEQTYPE){
	switch (seqtype){
	case SEQTYPE_NUCLEO:
	case SEQTYPE_NUCLEOAMBIG:
	case SEQTYPE_AMINO:
	case SEQTYPE_AMINOAMBIG:
		return fprint_base; break;
	case SEQTYPE_CODON:
		return fprint_codon; break;
	case SEQTYPE_CODONQ:
		return fprint_codonq; break;
	default:
		err(EXIT_FAILURE,"Unrecognised seqtype %d at %s:%d\n",seqtype,__FILE__,__LINE__);
	}
	return NULL;
}

bool is_readable_seqtype (enum SEQTYPE seqtype){
	switch(seqtype){
	case SEQTYPE_NUCLEO:
	case SEQTYPE_NUCLEOAMBIG:
	case SEQTYPE_AMINO:
	case SEQTYPE_AMINOAMBIG:        return true;
	default:
		return false;
	}
	return false;
}


bool is_printable_seqtype (enum SEQTYPE seqtype){
	switch(seqtype){
	case SEQTYPE_NUCLEO:
	case SEQTYPE_NUCLEOAMBIG:
	case SEQTYPE_AMINO:
	case SEQTYPE_AMINOAMBIG:	return true;
	default:
		return false;
	}
	return false;
}
	

char char_from_base ( const int base, const enum SEQTYPE seqtype){
	/*  Genetic code argument in IsValidBase is irrelevant since no
         * codon seqtype is (directly) printable
	 */
	assert(is_printable_seqtype(seqtype) && IsValidBase(base,seqtype,-1));
	switch(seqtype){
	case SEQTYPE_AMINO:		return Amino[base];
	case SEQTYPE_NUCLEO:		return Nucleo[base];
	case SEQTYPE_AMINOAMBIG:	return AminoAmbig[base];
	case SEQTYPE_NUCLEOAMBIG:	return NucleoAmbig[base];
	default:
		err(EXIT_FAILURE,"Called char_from_base with nonprintable SEQTYPE %d.\n",seqtype);
	}
	return '\0';
}

int base_from_char ( const char c, const enum SEQTYPE seqtype){
	assert(is_readable_seqtype(seqtype));
	switch(seqtype){
	case SEQTYPE_NUCLEO:		return nucleo_from_char(c);
	case SEQTYPE_NUCLEOAMBIG:	return nucleoambig_from_char(c);
	case SEQTYPE_AMINO:		return amino_from_char(c);
	case SEQTYPE_AMINOAMBIG:	return aminoambig_from_char(c);
	default:
		err(EXIT_FAILURE,"Trying to read non-readable seqtype %d (base was %d). %s:%d\n",seqtype,c,__FILE__,__LINE__);
	}
	return -1;
}


int codonambig_from_triplenuc (struct triple_nuc nucs){
	const int nuc_gap = GapChar(SEQTYPE_NUCLEOAMBIG);
	/*  Case: have gap */
	if ( nuc_gap==nucs.fst && nuc_gap==nucs.snd && nuc_gap==nucs.trd)
		return GapChar(SEQTYPE_CODONAMBIG);
	/*  Case: some but not all of codon positions are gaps. Not legal */
	if ( nuc_gap==nucs.fst || nuc_gap==nucs.snd || nuc_gap==nucs.trd)
		err(EXIT_FAILURE,"Illegal gapping pattern in codon in triple_to_codon: %c%c%c\n", char_from_base(nucs.fst,SEQTYPE_NUCLEOAMBIG), char_from_base(nucs.snd,SEQTYPE_NUCLEOAMBIG), char_from_base(nucs.trd,SEQTYPE_NUCLEOAMBIG));

	/* Hard coded numbers should be removed */
	return 25*nucs.fst + 5*nucs.snd + nucs.trd;
}

int codon_from_triplenuc (struct triple_nuc nucs){
	const int nuc_gap = GapChar(SEQTYPE_NUCLEO);
	/*  Case: have gap */
	if ( nuc_gap==nucs.fst && nuc_gap==nucs.snd && nuc_gap==nucs.trd)
		return GapChar(SEQTYPE_CODON);
	/*  Case: some but not all of codon positions are gaps. Not legal */
	if ( nuc_gap==nucs.fst || nuc_gap==nucs.snd || nuc_gap==nucs.trd)
		err(EXIT_FAILURE,"Illegal gapping pattern in codon in triple_to_codon: %c%c%c\n", char_from_base(nucs.fst,SEQTYPE_NUCLEO), char_from_base(nucs.snd,SEQTYPE_NUCLEO), char_from_base(nucs.trd,SEQTYPE_NUCLEO));

        /* Hard coded numbers should be removed */
        return 16*nucs.fst + 4*nucs.snd + nucs.trd;
}


struct triple_nuc triplenuc_from_codonambig ( const int codon ){
	struct triple_nuc nucs;
	const int nuc_gap = GapChar(SEQTYPE_NUCLEOAMBIG);
	const int codon_gap = GapChar(SEQTYPE_CODONAMBIG);
	/* Deal with all gaps */
	if ( codon == codon_gap){
		nucs.fst = nucs.snd = nucs.trd = nuc_gap;
		return nucs;
	}
	nucs.fst = codon / 25;
	int rem = codon % 25;
	nucs.snd = rem / 5;
	nucs.trd = codon % 5;
	return nucs;
}

struct triple_nuc triplenuc_from_codon ( const int codon ){
	struct triple_nuc nucs;
	const int nuc_gap = GapChar(SEQTYPE_NUCLEO);
	const int codon_gap = GapChar(SEQTYPE_CODON);
	/* Deal with all gaps */
	if ( codon == codon_gap){
		nucs.fst = nucs.snd = nucs.trd = nuc_gap;
		return nucs;
	}
	nucs.fst = codon / 16;
	int rem = codon % 16;
	nucs.snd = rem / 4;
	nucs.trd = codon % 4;
	return nucs;
}


int codon_from_nucs ( const int base1, const int base2, const int base3){
	struct triple_nuc nucs;
	nucs.fst=base1; nucs.snd=base2; nucs.trd=base3;
	return codon_from_triplenuc (nucs);
}

int codonambig_from_nucs ( const int base1, const int base2, const int base3){
	struct triple_nuc nucs;
	nucs.fst=base1; nucs.snd=base2; nucs.trd=base3;
	return codonambig_from_triplenuc (nucs);
}


unsigned int ambiguity_char (const enum SEQTYPE seqtype){
	switch(seqtype){
	case SEQTYPE_NUCLEOAMBIG: return 4;
	case SEQTYPE_AMINOAMBIG: return 20;
	case SEQTYPE_CODONAMBIG: err(EXIT_FAILURE,"Called \"ambiguity_character\" with argument SEQTYPE_CODONAMBIG\n");
	case SEQTYPE_NULL:
	case SEQTYPE_NUCLEO:
	case SEQTYPE_AMINO:
	case SEQTYPE_CODON:
	case SEQTYPE_CODONQ: err(EXIT_FAILURE,"No \"ambiguity_character\" for seqtype %d\n",seqtype);
	}
	return -1; /* Should never reach here */
}

unsigned int codon_from_codonambig ( const unsigned int codon){
	if (is_ambiguous(codon,SEQTYPE_CODONAMBIG))
		err(EXIT_FAILURE,"Called codon_from_codonambig with ambiguous codon\n");
	if(codon==GapChar(SEQTYPE_CODONAMBIG)){ return GapChar(SEQTYPE_CODON);}
	struct triple_nuc nucs = triplenuc_from_codonambig(codon);
	return 16*nucs.fst + 4*nucs.snd + nucs.trd;
}

bool is_ambiguous_seqtype (const enum SEQTYPE seqtype){
	switch(seqtype){
	case SEQTYPE_NUCLEOAMBIG:
	case SEQTYPE_AMINOAMBIG:
	case SEQTYPE_CODONAMBIG:	return true;
	default:			return false;
	}
	return false;
}

enum SEQTYPE ambig_seqtype_from_seqtype(const enum SEQTYPE seqtype){
	switch(seqtype){
	case SEQTYPE_NUCLEO:	return SEQTYPE_NUCLEOAMBIG;
	case SEQTYPE_AMINO:	return SEQTYPE_AMINOAMBIG;
	case SEQTYPE_CODON:	return SEQTYPE_CODONAMBIG;
	case SEQTYPE_NUCLEOAMBIG:
	case SEQTYPE_AMINOAMBIG:
	case SEQTYPE_CODONAMBIG:
		return seqtype;
	case SEQTYPE_CODONQ:
	case SEQTYPE_NULL:
		err(EXIT_FAILURE,"No ambiguous SEQTYPE for %d. %s:%d\n",seqtype,__FILE__,__LINE__);
	}
	return SEQTYPE_NULL;
}

enum SEQTYPE nonambig_seqtype_from_seqtype(const enum SEQTYPE seqtype){
	switch(seqtype){
	case SEQTYPE_NUCLEOAMBIG:	return SEQTYPE_NUCLEO;
	case SEQTYPE_AMINOAMBIG:	return SEQTYPE_AMINO;
	case SEQTYPE_CODONAMBIG:	return SEQTYPE_CODON;
	default:
		return seqtype;
	}
	return SEQTYPE_NULL;
}


bool is_ambiguous ( const int base, const enum SEQTYPE seqtype ){
	int ambig_nucleo;
	struct triple_nuc codon;
	switch(seqtype){
	case SEQTYPE_NULL:
	case SEQTYPE_NUCLEO:
	case SEQTYPE_AMINO:
	case SEQTYPE_CODON:
	case SEQTYPE_CODONQ: return false;
	case SEQTYPE_AMINOAMBIG:
	case SEQTYPE_NUCLEOAMBIG:
		if ( ambiguity_char(seqtype)==base) return true;
		else return false;
	case SEQTYPE_CODONAMBIG:
		ambig_nucleo = ambiguity_char(SEQTYPE_NUCLEOAMBIG);
		codon = triplenuc_from_codonambig (base);
		if (ambig_nucleo==codon.fst) return true;
		if (ambig_nucleo==codon.snd) return true;
		if (ambig_nucleo==codon.trd) return true;
		return false;
	}
	return false;
}

/*  Convert a single character alphabetical representation of an amino acid to
 * integer*/
int aminoambig_from_char (char c)
{

  c = toupper (c);
  switch (c) {
  case 'A': return 0;
  case 'R': return 1;
  case 'N': return 2;
  case 'D': return 3;
  case 'C': return 4;
  case 'Q': return 5;
  case 'E': return 6;
  case 'G': return 7;
  case 'H': return 8;
  case 'I': return 9;
  case 'L': return 10;
  case 'K': return 11;
  case 'M': return 12;
  case 'F': return 13;
  case 'P': return 14;
  case 'S': return 15;
  case 'T': return 16;
  case 'W': return 17;
  case 'Y': return 18;
  case 'V': return 19;
  case 'X': return 20;
  case '-': return 21;
  }

  return -1;
}

int amino_from_char (char c)
{

  c = toupper (c);
  switch (c) {
  case 'A': return 0;
  case 'R': return 1;
  case 'N': return 2;
  case 'D': return 3;
  case 'C': return 4;
  case 'Q': return 5;
  case 'E': return 6;
  case 'G': return 7;
  case 'H': return 8;
  case 'I': return 9;
  case 'L': return 10;
  case 'K': return 11;
  case 'M': return 12;
  case 'F': return 13;
  case 'P': return 14;
  case 'S': return 15;
  case 'T': return 16;
  case 'W': return 17;
  case 'Y': return 18;
  case 'V': return 19;
  case '-': return 20;
  }

  return -1;
}



/*  Convert single character alphabetic representation of nucleotide to integer i
 */
int nucleoambig_from_char (char c)
{
  c = toupper (c);
  switch (c) {
  case 'A': return 0;
  case 'C': return 1;
  case 'G': return 2;
  case 'T': return 3;
  case 'N': return 4;
  case '-': return 5;
  }

  return -1;
}

int nucleo_from_char (char c)
{
  c = toupper (c);
  switch (c) {
  case 'A': return 0;
  case 'C': return 1;
  case 'G': return 2;
  case 'T': return 3;
  case '-': return 4;
  }

  return -1;
}



/*  Returns number of possible bases for each sequence type. Bases are numbered
 * zero to NumberPossibleBases(). The number representing gap is equal to
 * NumberPossibleBases(), except for SEQTYPE_CODONQ.
 * Includes ambiguity characters.
 */
int NumberPossibleBases (const enum SEQTYPE seqtype, const int gencode)
{
  assert (IsSeqtype (seqtype));
  assert ( seqtype!=SEQTYPE_CODONQ || IsValidGencode(gencode));


  switch (seqtype) {
  case SEQTYPE_NULL:		return 0;
  case SEQTYPE_NUCLEOAMBIG: 	return 5;
  case SEQTYPE_NUCLEO:		return 4;
  case SEQTYPE_AMINOAMBIG:	return 21;
  case SEQTYPE_AMINO:		return 20;
  case SEQTYPE_CODONAMBIG:	return 125;
  case SEQTYPE_CODON:		return 64;
  case SEQTYPE_CODONQ:
    return NumberSenseCodonsInGenCode (gencode);
  }

  return -1;
}


/*  Return integer corresponding to gapchar for sequence type
 */
int GapChar (const enum SEQTYPE seqtype)
{
  assert(IsSeqtype(seqtype));

  switch (seqtype) {
  case SEQTYPE_NULL:	err(EXIT_FAILURE,"No gap character for SEQTYPE_NULL %s:%d\n",__FILE__,__LINE__);
  case SEQTYPE_NUCLEO:		return 4;
  case SEQTYPE_NUCLEOAMBIG:	return 5;
  case SEQTYPE_AMINO:		return 20;
  case SEQTYPE_AMINOAMBIG:	return 21;
  case SEQTYPE_CODONAMBIG:	return 125;
  case SEQTYPE_CODON:
  case SEQTYPE_CODONQ:		return 64;
  }

  return -1;
}


/*  Is sequence type valid?
 */
bool IsSeqtype (const enum SEQTYPE seqtype)
{
  switch (seqtype) {
  case SEQTYPE_NULL:
  case SEQTYPE_NUCLEO:
  case SEQTYPE_NUCLEOAMBIG:
  case SEQTYPE_AMINO:
  case SEQTYPE_AMINOAMBIG:
  case SEQTYPE_CODON:
  case SEQTYPE_CODONAMBIG:
  case SEQTYPE_CODONQ:
    return true;
  }

  return false;
}

/*  Is base valid for sequence type? For non-codon data, or codons in Q coord's,
 * the gencode argument in not needed and can be set to any value
 */ 
bool IsValidBase ( const int base, const enum SEQTYPE seqtype, const int gencode){
  int gapchar;

  gapchar = GapChar(seqtype);
  
  if ( 0<=base  &&  (base<NumberPossibleBases(seqtype,gencode)  || base==gapchar) )
    return true;

  return false;
}


double *  ConvertCodonFreqsToQcoord ( const double * freqs, const int gencode){
  double * qfreqs;
  int max_codon,max_qcodon,codon,qcodon;

  assert(NULL!=freqs);
  assert(IsValidGencode(gencode));
  max_codon = NumberPossibleBases(SEQTYPE_CODON,gencode);
  max_qcodon = NumberPossibleBases(SEQTYPE_CODONQ,gencode);

  #ifndef NDEBUG
  {
    double sum = 0.;
    for ( codon=0 ; codon<max_codon ; codon++){
      assert(freqs[codon]>=0. && freqs[codon]<=1.);
      sum += freqs[codon];
    }
    assert(fabs(1.-sum)<DBL_EPSILON*max_codon);
  }
  #endif

  qfreqs = calloc ( max_qcodon,sizeof(double));
  if(NULL==qfreqs)
    return NULL;

  for ( codon=0 ; codon<max_codon; codon++){
    qcodon = CodonToQcoord (codon,SEQTYPE_CODON,gencode);
    if ( -1 != qcodon){
      qfreqs[qcodon] = freqs[codon];
    }
  }

  #ifndef NDEBUG
  {
    double sum = 0.;
    for ( qcodon=0 ; qcodon<max_qcodon ; qcodon++){
      assert(qfreqs[qcodon]>=0. && qfreqs[qcodon]<=1.);
      sum += qfreqs[qcodon];
    }
    assert(fabs(1.-sum)<DBL_EPSILON*max_qcodon);
  }
  #endif

  return qfreqs;
}

int nonambig_base_from_base (const int base, const enum SEQTYPE seqtype){
	const int gapchar = GapChar(seqtype);
	if (is_ambiguous(base,seqtype)){
		err(EXIT_FAILURE,"Called nonambigous_base_from_base with ambiguous base. %s:%d\n",__FILE__,__LINE__);
	}
	switch(seqtype){
	case SEQTYPE_NUCLEOAMBIG:
	case SEQTYPE_AMINOAMBIG:
		if (gapchar==base){ return GapChar(nonambig_seqtype_from_seqtype(seqtype));}
		return base;
	case SEQTYPE_CODONAMBIG:
		codon_from_codonambig (base);
	default:
		return base;
	}
}
