#include <stdio.h>
#include "model.h"
#include "options.h"
#include "rng.h"
#include "data.h"

#define VERSIONSTRING   "ensembl-development"

/*   Strings describing options and defaults */
int n_options = 23;
char *options[] =
  { "seqfile", "treefile", "outfile", "kappa", "omega", "codonf",
"nucleof", "aminof", "reoptimise", "nucfile", "aminofile", "positive_only",
"gencode","timemem","ldiff", "paramin", "paramout", "skipsitewise", "seed",
"saveseed", "freqtype", "cleandata", "branopt" };
char *optiondefault[] =
  { "incodon", "intree", "slr.res", "2.0", "0.1", "0", "0", "0", "1",
"nuc.dat", "amino.dat", "0", "universal","0", "3.841459", "", "", "0", "0", "1", "0", "0","1" };
char optiontype[] =
  { 's', 's', 's', 'f', 'f', 'd', 'd', 'd', 'd', 's', 's', 'd', 's', 'd', 'f', 
's', 's', 'd', 'd', 'd', 'd','d','d'};
int optionlength[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1 };
char *default_optionfile = "slr.ctl";


int main ( const int argc, const char * argv[]){
	/*  Sort out options */
	/* TODO: replace with implementation based on getopt_long */
	ReadOptions (argc, argv);

	double kappa = *(double *)   GetOption ("kappa");
	double omega = *(double *)   GetOption ("omega");
	char * seqfile = (char *)    GetOption ("seqfile");
	char * treefile = (char *)   GetOption ("treefile");
	char * outfile = (char *)    GetOption ("outfile");
	int codonf = *(int *)     GetOption ("codonf");
	int nucleof = *(int *)    GetOption ("nucleof");
	int aminof = *(int *)     GetOption ("aminof");
	char * nucfile = (char *)    GetOption ("nucfile");
	char * aminofile = (char *)  GetOption ("aminofile");
	int reoptimise = *(int *) GetOption ("reoptimise");
	int positive = *(int *)   GetOption ("positive_only");
	char * gencode_str = (char *)GetOption ("gencode");
	int timemem = *(int *)    GetOption ("timemem");
	double ldiff = *(double *)   GetOption ("ldiff");
	char * paramin = (char *)    GetOption ("paramin");
	char * paramout = (char *)   GetOption ("paramout");
	int skipsitewise = *(int *) GetOption ("skipsitewise");
	unsigned int seed = *(unsigned int *) GetOption("seed");
	unsigned int saveseed = *(unsigned int *) GetOption("saveseed");
	unsigned int freqtype = *(unsigned int *) GetOption("freqtype");
	unsigned int cleandata = *(unsigned int *) GetOption("cleandata");
	enum model_branches branopt = *(enum model_branches *) GetOption("branopt");

	PrintOptions();

	RL_Init();

	fputs("# SLR \"Sitewise Likelihood Ratio\" selection detection program. Version ",stdout);
	fputs(VERSIONSTRING,stdout);
	fputc('\n',stdout);
	
	SetAminoAndCodonFuncs (nucleof, aminof, nucfile, aminofile);
	gencode = GetGeneticCode (gencode_str);

	/*  Read in sequence data and convert to (unique) codons */
	DATA_SET * nuc_data = read_data (seqfile,SEQTYPE_NUCLEO);
	DATA_SET * codon_data = ConvertNucToCodon (nuc_data,gencode);
	DATA_SET * zcodon_data = compress_data(sort_data(codon_data));
	DATA_SET * data_zn = RemoveTrivialObs(zcodon_data);
	FreeDataSet(nuc_data);
	FreeDataSet(codon_data);
	FreeDataSet(zcodon_data);
	if(NULL==zncodon_data){
		fputs("Error reading in sequences\n",stderr);
		exit(EXIT_FAILURE);
	}
	printf ("# Read seqfile file %s. %d species, %d sites.\n", seqfile, data_zn->n_sp, data_zn->n_pts);

	double * freqs = GetBaseFreqs (data, 0); // Needs to be calculated before conversion to Qcoord

	DATA_SET * data_znq = ConvertCodonToQcoord(data_zn);
	FreeDataSet(data_zn);
	if (data_znq->n_pts != data_znq->n_unique_pts){
		printf ("# Redundancy. Reduced sites from %d to %d\n", data->n_pts, data->n_unique_pts);
	}

	/*  Read in tree */
	TREE ** trees = read_tree_strings (treefile);
	create_tree(trees[0]);
	if ( data_znq->n_sp != trees[0]->nsp ){ 
		fprintf (stderr,"Tree and sequence file are incompatible. %d sequences in tree, %d in sequence file\n", trees[0]->nsp, data_znq->n_sp);
		exit(EXIT_FAILURE);
	}
	printf ("# Read tree from %s.\n", treefile);
	print_tree (stdout, trees[0]->tree, NULL, trees[0]);


	

	return EXIT_SUCCESS;
}
