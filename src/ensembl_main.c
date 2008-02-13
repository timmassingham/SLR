/*
 *  Copyright 2008 Tim Massingham (tim.massingham@ebi.ac.uk)
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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include "model.h"
#include "options.h"
#include "rng.h"
#include "data.h"
#include "tree.h"
#include "bases.h"
#include "codonmodel.h"
#include "gencode.h"
#include "statistics.h"
#include "tree_data.h"
#include "optimize.h"
#include "like.h"
#include "root.h"
#include "gamma.h"
#include "linemin.h"
#include "matrix.h"
#include "utility.h"

#define VERSIONSTRING   "ensembl-development"
#define GRIDSIZE        50


/*   Strings describing options and defaults */
int n_options = 21;
char *options[] =
  { "seqfile", "treefile", "outfile", "kappa", "omega", "codonf",
"nucleof", "aminof", "nucfile", "aminofile", "positive_only",
"gencode","ldiff", "seed",
"saveseed", "freqtype", "cleandata", "pvalthresh", "branchpval", "saturated", "subtreeout" };
char *optiondefault[] =
  { "incodon", "intree", "slr.res", "2.0", "0.1", "0", "0", "0",
"nuc.dat", "amino.dat", "0", "universal", "3.841459", "0", "1", "0", "0", "0.01","-0.", "1.",
"subtrees.out" };
char optiontype[] =
  { 's', 's', 's', 'f', 'f', 'd', 'd', 'd', 's', 's', 'd', 's', 'f', 
 'd', 'd', 'd','d', 'f', 'f', 'f' ,'s'};
int optionlength[] = {  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
char *default_optionfile = "slr.ctl";

struct selectioninfo {
	double * llike_neu;
	double * llike_max;
	double * omega_max;
	double * lbound, *ubound;
	int * type;
};

struct ml_res {
	double loglike, * score, * info, * var;
	double * param;
	unsigned int nparam;
	double scale_neu,scale_synon;
};

struct midpoint {
	unsigned int branch;
	double left,right,length;
};


const char *OutString[5] = { "All gaps", "Single char", "Synonymous", "", "Constant" };

VEC create_grid ( const unsigned int len, const int positive);
int FindBestX (const double *grid, const int site, const int n);
DATA_SET *ReadData (const char *name, const int gencode);
struct ml_res * OptimizeTree ( const DATA_SET * data, TREE * tree, double * freqs, double * x, const unsigned int freqtype, const int codonf, const enum model_branches branopt);
struct selectioninfo *  CalculateSelection ( TREE * tree, DATA_SET * data, double kappa, double omega, double * freqs, const double ldiff, const unsigned int freqtype, const int codonf);
void PrintResults ( char * outfile, struct selectioninfo * selinfo, const double * pval, const double * pval_adj, const int nsites);
double * CalculatePvals ( const double * lmax, const double * lneu, const int n, const int positive_only);
double * AdjustPvals ( const double * pval, DATA_SET * data);
double * CalculateEntropy ( const DATA_SET * data, const double * freqs);
void PrintSummary ( FILE * file, const struct selectioninfo * selinfo, const double * pval, const double * pval_adj,const int n_pts);
double CalcLike_Wrapper ( const double * x, void * info);
void Set_CalcLike_Wrapper ( double (*f)(const double *,void *), double diff);
void GradLike_Single (const double *param, double *grad, void *data);
double CalcLike_Single (const double *param, void *data);
double * GradLike_Full (const double *param, double *grad, void *data);
double * InfoLike_Full (const double * param, double * info, void *data);
double ztest_twotail ( const double stat, const double mean, const double var);
double * change_variance_basis ( const double * m, const double * var, const unsigned int n);
VEC matrix_vec_mult ( const double * m, const unsigned int nrow, const unsigned int ncol, VEC v);
unsigned int get_scaling_basis_sub (double * basis, const NODE * node, const NODE * parent, const unsigned int nparam, const bool from_root);
double * get_scaling_basis ( const unsigned int root, const TREE * tree, const unsigned int nparam);
bool * analyse_branches ( const TREE * tree, const VEC blength, const double ds_scale, const double thresh);
TREE * remove_branches_from_tree ( TREE * tree, const bool * keep );
TREE * remove_orphan_branches ( TREE * tree );
unsigned int remove_orphan_branches_sub ( const NODE * node );
TREE * delete_orphan_nodes ( TREE * tree);
unsigned int count_edges ( const NODE * node);
void assess_orphan_node ( const NODE * node);
unsigned int print_subtrees ( FILE * fp, const TREE * tree);
void prune_into_subtrees ( TREE * tree, const bool * keep, const char * subtreefile);
double * max_tipdistance ( const TREE * tree);
void max_tipdistance_pass1 (const NODE * node, const NODE * parent, double * LRmax);
void max_tipdistance_pass2 (const NODE * node, const NODE * parent, const NODE * gparent, double *left, double * right);
struct midpoint * midpoint_tree ( const TREE * tree );
double * subsqrmatrix ( const double * m, const unsigned int n, const bool * insub);
unsigned int branch_number ( const NODE * node, const NODE * parent);





int main ( int argc, char * argv[]){
	bool skip_lineage = false;
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
	int positive = *(int *)   GetOption ("positive_only");
	char * gencode_str = (char *)GetOption ("gencode");
	double ldiff = *(double *)   GetOption ("ldiff");
	unsigned int seed = *(unsigned int *) GetOption("seed");
	unsigned int saveseed = *(unsigned int *) GetOption("saveseed");
	unsigned int freqtype = *(unsigned int *) GetOption("freqtype");
	unsigned int cleandata = *(unsigned int *) GetOption("cleandata");
	const double pval_thresh = *(double *) GetOption("pvalthresh");
	double bran_pval_thresh = *(double *) GetOption("branchpval");
	double saturated = *(double *) GetOption("saturated");
	char * subtreeout = (char *) GetOption("subtreeout");
	if ( -0. == bran_pval_thresh){
		bran_pval_thresh = pval_thresh;
	}

	PrintOptions();

	RL_Init(seed);

	fputs("# SLR \"Sitewise Likelihood Ratio\" selection detection program. Version ",stdout);
	fputs(VERSIONSTRING,stdout);
	fputc('\n',stdout);
	
	SetAminoAndCodonFuncs (nucleof, aminof, nucfile, aminofile);
	const unsigned int gencode = GetGeneticCode (gencode_str);

	/*  Read in sequence data and convert to (unique) codons */
	DATA_SET * nuc_data = read_data (seqfile,SEQTYPE_NUCLEO);
	DATA_SET * codon_data = ConvertNucToCodon (nuc_data,gencode);
	DATA_SET * zcodon_data = compress_data(sort_data(codon_data));
	DATA_SET * data_zn = RemoveTrivialObs(zcodon_data);
	FreeDataSet(nuc_data);
	FreeDataSet(codon_data);
	FreeDataSet(zcodon_data);
	if(NULL==data_zn){
		fputs("Error reading in sequences\n",stderr);
		exit(EXIT_FAILURE);
	}
	printf ("# Read seqfile file %s. %d species, %d sites.\n", seqfile, data_zn->n_sp, data_zn->n_pts);

	double * freqs = GetBaseFreqs (data_zn, 0); // Needs to be calculated before conversion to Qcoord

	DATA_SET * data_znq = ConvertCodonToQcoord(data_zn);
	if (data_znq->n_pts != data_znq->n_unique_pts){
		printf ("# Redundancy. Reduced sites from %d to %d\n", data_znq->n_pts, data_znq->n_unique_pts);
	}

	/*  Read in tree */
	TREE ** trees = read_tree_strings (treefile);
	create_tree(trees[0]);
	if ( data_znq->n_sp != trees[0]->n_sp ){ 
		fprintf (stderr,"Tree and sequence file are incompatible. %d sequences in tree, %d in sequence file\n", trees[0]->n_sp, data_znq->n_sp);
		exit(EXIT_FAILURE);
	}
	printf ("# Read tree from %s.\n", treefile);
	print_tree (stdout, trees[0]->tree, NULL, trees[0], NULL);

	fputs("# Start of checks. Holding branches proportional and optimising other parameters.\n",stdout);
	/*  Check that tree does have branch lengths */
	for ( int bran=0 ; bran<trees[0]->n_br ; bran++){
		NODE * node = trees[0]->branches[bran];
		assert(NULL!=node);
		if ( node->blength[0] < 0.){
			if(!skip_lineage){fputs("# At least one branch length missing from tree, Skipping lineage-specific test.\n",stderr);}
			skip_lineage = true;
			node->blength[0] = RandomExp(0.1);
			const unsigned int a = find_connection(node->branch[0],node);
			assert(a>=0 && a<node->branch[0]->nbran);
			node->branch[0]->blength[ a ] = node->blength[0];
		}
	}
skip_lineage=false;
	struct ml_res * loglike_prop = NULL;
	if ( ! skip_lineage ){	
		int nparam = 3; /* kappa, omega and scale */
		double * x = calloc(nparam,sizeof(double));
		x[0] = 1.; x[1] = kappa; x[2] = omega;
		assert(kappa>=0.); assert(omega>=0.);
		loglike_prop = OptimizeTree (data_znq,trees[0],freqs,x,freqtype,codonf,Branches_Proportional);
		double scale = x[0]; kappa = x[1]; omega = x[2];
		fprintf(stdout,"#\tlnL = %e\tkappa = %e\tomega = %e\tscale = %e\n",loglike_prop->loglike,kappa,omega,scale);
		ScaleTree(trees[0],scale);
		free(x);
	
		print_tree(stdout,trees[0]->tree,NULL,trees[0], NULL);
		VEC old_branches = branchlengths_from_tree (trees[0]);
		free_summary(fprint_summary(stdout,"branch lengths",summarise_vec(copy_vec(old_branches))));
		bool * keep = analyse_branches (trees[0],old_branches,loglike_prop->scale_synon,saturated);
		/* factor following into function */
		if ( trees[0]->n_br != sum_bool(keep,trees[0]->n_br) ){
			prune_into_subtrees(trees[0],keep,subtreeout);
		}
		free(keep);
		free_vec(old_branches);
	}
	
	/*  Full optimisation using previous as initial values  */
	unsigned int nparam = 2 + trees[0]->n_br;
	double * x = calloc(nparam,sizeof(double));
	VEC old_branches = branchlengths_from_tree (trees[0]);
	for( unsigned int i=0 ; i<trees[0]->n_br ; i++){
		x[i] = vget(old_branches,i);
	}
	x[0+trees[0]->n_br] = kappa; x[1+trees[0]->n_br] = omega;
	fprintf(stdout, "#\tkappa = %e\tomega = %e\n",kappa,omega);
	add_lengths_to_tree(trees[0],x);
	print_tree(stdout,trees[0]->tree,NULL,trees[0], NULL);
	struct ml_res * loglike_full = OptimizeTree (data_znq,trees[0],freqs,x,freqtype,codonf,Branches_Variable);
	fprintf(stdout, "#\tlnL = %e\tkappa = %e\tomega = %e\n",loglike_full->loglike,kappa,omega);
	print_tree(stdout,trees[0]->tree,NULL,trees[0], NULL);
	VEC new_branches = branchlengths_from_tree (trees[0]);
	kappa = x[0+trees[0]->n_br]; omega = x[1+trees[0]->n_br];
	free_summary(fprint_summary(stdout,"branch lengths",summarise_vec(copy_vec(new_branches))));
	free(x);
	bool * keep = analyse_branches (trees[0],new_branches,loglike_full->scale_synon,saturated);
	/* factor following into function */
	if ( trees[0]->n_br != sum_bool(keep,trees[0]->n_br) ){
		prune_into_subtrees(trees[0],keep,subtreeout);
	}
	free(keep);

	/*  Check whether there has been a significant change in any branches  */
	if ( ! skip_lineage ){
		double lrt = -2.*(loglike_full->loglike - loglike_prop->loglike);
		double lrt_pval = pchisq(lrt,trees[0]->n_br - 1,1);
		fputs("# Test statistic and pvalue for lineage-specific changes in branch length.\n",stdout);
		fputs("# Note: Test may be inaccurate (but conservative) when some branch lengths\n# are zero under the null model (see Self and Liang, 1984).\n",stdout);
		fprintf(stdout,"#\tLRT = %e\tpval = %e\tdeg = %d\n",lrt,lrt_pval, trees[0]->n_br - 1);
fprintf(stdout,"#\tNeutral scale = %e\tSynonymous scale = %e\n",loglike_full->scale_neu,loglike_full->scale_synon);
		if ( lrt_pval < pval_thresh){
			fputs("# Warning: significant lineage-specific changes in branch length detected\n",stdout);
			fputs("# Checking branches for significance (branch-wise z-test)\n",stdout);
			const unsigned int nbr = trees[0]->n_br;
			const unsigned int nparam = loglike_full->nparam;
			unsigned int nbran_signif = 0;
			double * pvals = calloc(nbr,sizeof(double));
			for ( int br=0 ; br<nbr ; br++){
				pvals[br] = ztest_twotail(vget(old_branches,br),vget(new_branches,br),loglike_full->var[br*nparam+br]);
			}
			double * corrected_pvals = Pvalue_adjust_StepUp (pvals,nbr,BONFERRONI);
			for ( int br=0 ; br<nbr ; br++){
				if (pvals[br]<bran_pval_thresh){
					nbran_signif++;
					if ( NULL!=trees[0]->branches[br]->name){
						fprintf(stdout,"# Branch %d (%s) significant. pval = %e (corrected %e)\n",br,trees[0]->branches[br]->name,pvals[br],corrected_pvals[br]); 
					} else {
						fprintf(stdout,"# Branch %d significant. pval = %e (corrected %e)\n",br,pvals[br],corrected_pvals[br]);
					}
				}
			}
			free(corrected_pvals);
			if ( 0==nbran_signif){
				fputs("# No branches individually significant\n",stdout);
			}

			/* Start of code for checking bipartitions. Need functions for midpt rooting */
			fputs("# Checking treewise bipartitions of data for significance (z-test)\n",stdout);
			unsigned int npart_signif = 0;
			struct midpoint * midpt = midpoint_tree(trees[0]);
			if (NULL!=trees[0]->branches[midpt->branch]->name){
				fprintf(stdout,"# Midpoint rooting. Root on branch %d (%s)\n",midpt->branch,trees[0]->branches[midpt->branch]->name);
			} else {
				fprintf(stdout,"# Midpoint rooting. Root on branch %d\n",midpt->branch);
			}
			double * scaling_basis = get_scaling_basis(midpt->branch,trees[0],nparam);

			bool * isbranch = calloc(nparam,sizeof(bool));
			for ( unsigned int i=0 ; i<nbr ; i++){ isbranch[i] = true;}
			double * branch_transform = subsqrmatrix(scaling_basis,nparam,isbranch);

			free(isbranch);
			VEC new_scale = matrix_vec_mult (branch_transform,nbr,nbr,new_branches);
			VEC old_scale = matrix_vec_mult (branch_transform,nbr,nbr,old_branches);
			double * scaling_var = change_variance_basis (scaling_basis,loglike_full->var,nparam);
			for ( int i=0 ; i<nbr ; i++){
				pvals[i] = ztest_twotail(vget(old_scale,i),vget(new_scale,i),scaling_var[i*nparam+i]);
			}
			corrected_pvals = Pvalue_adjust_StepUp (pvals,nbr,BONFERRONI);
			for ( int i=0 ; i<nbr ; i++){
				if(pvals[i]<bran_pval_thresh){
					npart_signif++;
					if ( NULL!=trees[0]->branches[i]->name){
						fprintf(stdout,"# Partition %d (%s) significant. pval = %e (corrected %e)\n",i,trees[0]->branches[i]->name,pvals[i],corrected_pvals[i]);
					} else {
						fprintf(stdout,"# Partition %d significant. pval = %e (corrected %e)\n",i,pvals[i],corrected_pvals[i]);
					}
				}
			}
			free(pvals);
			free(corrected_pvals);
			if ( 0==npart_signif ){
				fputs("# No partitions individually significant\n",stdout);
			}
			if ( 0==npart_signif && 0==nbran_signif){
				fputs("# No partitions or branches found significant.\n# Not a problem if original result was marginally significant\n",stdout);
			}
			/* Worth looking at eigen factors? -- major deformations */
		}	
	}

	/*  Do normal analysis here  */
	struct selectioninfo * selinfo = CalculateSelection ( trees[0], data_znq, kappa,omega, freqs, ldiff,freqtype,codonf);
	double * pval = CalculatePvals ( selinfo->llike_max, selinfo->llike_neu, data_znq->n_pts, positive);
	double * pval_adj = AdjustPvals ( pval,data_znq);

	PrintResults ( outfile, selinfo, pval, pval_adj, data_znq->n_pts);
	PrintSummary ( stdout, selinfo, pval, pval_adj, data_znq->n_pts);

	if(saveseed){RL_Close();}
	return EXIT_SUCCESS;
}


int FindBestX (const double *grid, const int site, const int n)
{
  assert(NULL!=grid);
  assert(site>=0);
  assert(n>=0);

  double min = DBL_MAX;
  unsigned int bestidx = 0;
  for (unsigned int a = 0; a < n; a++) {
    if (grid[site * n + a] < min) {
      min = grid[site * n + a];
      bestidx = a;
    }
  }

  return bestidx;
}



DATA_SET *ReadData (const char *name, const int gencode)
{
  DATA_SET *tmp, *data;

  assert (NULL != name);

  /* Read nucleotides and convert into codons
   */
  tmp = read_data (name, SEQTYPE_NUCLEO);
  if (NULL == tmp)
    return NULL;
  data = ConvertNucToCodon (tmp,gencode);
  if (NULL == data) {
    puts
      ("Error converting nucleotides to codons. Returning uncompressed sequence.");
    return tmp;
  }
  FreeDataSet (tmp);

	/*  Check if sequence contains stop codons  */
	int nstop = count_alignment_stops(data);
	if (0!=nstop){
		fputs("Alignment contains stop codons. Cannot continue.\n",stderr);
		exit(EXIT_FAILURE);
	}
  /*  Sort and compress sequence to remove redundency
   */
  sort_data (data);
  tmp = compress_data (data);
  if (NULL == tmp) {
    puts ("Error compressing sequence! Returning uncompressed set");
    return data;
  }
  FreeDataSet (data);

  /*  Find and mask trivial observations in data. Likelihood for these
   * observations can be calculated trivially without using the pruning
   * algorithm
   */
  data = RemoveTrivialObs (tmp);
  if (NULL == data) {
    puts
      ("Error removing trivial observations (single chars and all gaps).\nReturning compressed sequence.\n");
    return tmp;
  }
  FreeDataSet (tmp);
  if ( data->n_pts != data->n_unique_pts){
    printf ("# Redundency. Reduced sites from %d to %d\n", data->n_pts, data->n_unique_pts);	}

  return data;
}

struct ml_res *  OptimizeTree ( const DATA_SET * data, TREE * tree, double * freqs, double * x, const unsigned int freqtype, const int codonf, const enum model_branches branopt){
  struct single_fun *info;
  double *bd,fx;
  int i;
  MODEL * model;

  CheckIsDataSet (data);
  CheckIsTree (tree);
  assert(NULL!=x);

  printf("# Reoptimising parameters, branches %s\n",model_branches_string[branopt]);

  const unsigned int nbr = tree->n_br;
  model = NewCodonModel_full ( data->gencode, x[nbr+0], x[nbr+1], freqs, codonf ,freqtype, branopt);
  OOM(model);
  model->exact_obs = 1;

  const unsigned int nparam = model->nparam + ((Branches_Variable==branopt)?nbr:0);
  bd = calloc ( 2*nparam,sizeof(double)); OOM(bd);

  /* Set boundaries
   */
  for ( i=0 ; i<nparam; i++){
    bd[i] = 1e-8;
    bd[i+nparam] = 99.;
  }

  /*  Check that initial estimates are within boundaries
   */
  for ( i=0 ; i<nparam ; i++){
    if (x[i]<=bd[i]) 
      x[i] = bd[i]+1e-5;
    if (x[i]>=bd[nparam+i])
      x[i] = bd[nparam+i]-1e-5;
  }

  info = calloc (1,sizeof (struct single_fun));
  OOM(info);
  info->tree = tree;
  info->p = calloc (data->n_pts * 2,sizeof (double));
  OOM(info->p);
  info->model = model;

  add_data_to_tree (data, tree, model);
  //x[nbr-1] = 1.;
  //CheckModelDerivatives(model,0.5,x+nbr,1e-5);
  fx = CalcLike_Single ( x, info);

  Optimize (x, nparam, GradLike_Full, CalcLike_Single, &fx, (void *) info, bd, 2);

  struct ml_res * result = calloc(1,sizeof(struct ml_res));
  result->loglike = fx;
  result->score = GradLike_Full (x, calloc(nparam,sizeof(double)), (void *)info);
  result->info = InfoLike_Full (x, calloc(nparam*nparam,sizeof(double)), (void *)info);
  result->var = InvertMatrix(CopyMatrix(result->info,calloc(nparam*nparam,sizeof(double)),nparam),nparam);
  result->param = CopyVector(x,calloc(nparam,sizeof(double)),nparam);
  result->nparam = nparam;
  {
	const double scale = GetScale_single(model,x[nparam-1]);
	result->scale_neu = scale/GetScale_single(model,1.);
	result->scale_synon = scale/GetScale_single(model,0.);
  }
  FreeModel (model);
  free(bd);
  free(info->p);
  free(info);

  return result;
}




struct selectioninfo * CalculateSelection ( TREE * tree, DATA_SET * data, double kappa, double omega, double * freqs, const double ldiff, const unsigned int freqtype, const int codonf){
  double x[1];
  struct selectioninfo * selinfo;
  int positive;
  double factor;
  MODEL * model;
  DATA_SET * data_single;
  struct single_fun * info;
  double bd[2];
  double * likelihood_grid;
  int col;
  int * done_usite;
  int species,bufflen;

  CheckIsTree(tree);
  CheckIsDataSet(data);
  assert(kappa>=0.);
  assert(omega>=0.);
  assert(freqs!=NULL);

  const int dosupport = (0.0==ldiff)?0:1;

  selinfo = malloc(sizeof(struct selectioninfo));		OOM(selinfo);
  selinfo->llike_neu = calloc(data->n_pts,sizeof(double));	OOM(selinfo->llike_neu);
  selinfo->llike_max = calloc(data->n_pts,sizeof(double));      OOM(selinfo->llike_max);
  selinfo->omega_max = calloc(data->n_pts,sizeof(double));      OOM(selinfo->omega_max);
  if ( dosupport ){
  	selinfo->lbound = calloc(data->n_pts,sizeof(double));         OOM(selinfo->lbound);
  	selinfo->ubound = calloc(data->n_pts,sizeof(double));         OOM(selinfo->ubound);
  } else {
	selinfo->lbound = NULL;
	selinfo->ubound = NULL;
  }
  selinfo->type = calloc ( data->n_pts,sizeof(int));		OOM(selinfo->type);

  positive = *(int *) GetOption("positive_only");

  model = NewCodonModel_single ( data->gencode, kappa,omega,freqs,codonf,freqtype);
  OOM(model);
  model->exact_obs = 1;

  /* Calculate scale factor relative to neutral evolution and scale
   * tree appropriately.
   */
  {
    factor = GetScale_single (model,omega);
    factor /= GetScale_single (model, 1.);
    ScaleTree (tree,factor);
    printf ("# Scaling tree to neutral evolution. Factor = %3.2f\n", factor);
  }


  //  One site data set to be used in all optimizations
  data_single = CreateDataSet (1, data->n_sp);
  OOM (data_single);
  for ( species=0 ; species<data->n_sp ; species++){
     bufflen = 1 + strlen(data->sp_name[species]);
     data_single->sp_name[species] = malloc(bufflen*sizeof(char));
     strncpy(data_single->sp_name[species],data->sp_name[species],bufflen);
   }

  info = calloc(1,sizeof(struct single_fun));
  info->tree = tree;
  info->p = calloc(2*data->n_unique_pts,sizeof(double));

  //  Set boundaries. Lower bound is 1. if only positive selection is
  // is of interest.
  if (positive == 1)
    bd[0] = 1.;
  else
    bd[0] = 0.;
  bd[1] = 99.;
  info->model = model;

  

  /*  Calculate grid of sitewise likelihoods for many omega, use to
   * provide good starting values for each sitewise observation.
   *  Do this since it is relative quick to calculate the sitewise for all
   * sites for a single omega (due to memory effects rather than algorithms).
   */
  puts ("# Calculating initial estimates of sitewise conservation");
  add_data_to_tree (data, tree, model);
  const VEC omega_grid = create_grid(GRIDSIZE,positive);
  /*  Fill out sitewise likelihoods for grid  */
  likelihood_grid = calloc (data->n_unique_pts * GRIDSIZE, sizeof (double));
  OOM(likelihood_grid);
  for ( unsigned int row=0; row< GRIDSIZE ; row++){
    x[0] = vget(omega_grid,row);
    CalcLike_Single (x,info);
    for ( unsigned int pt=0 ; pt<data->n_unique_pts ; pt++){
      likelihood_grid[pt*GRIDSIZE+row] = -(tree->tree)->scalefactor - log (info->p[pt]);
    }
  }
  /*  Fill out vector of likelihoods for neutral evolution */
  double * likelihood_neutral = calloc(data->n_unique_pts,sizeof(double));
  x[0] = 1.;
  CalcLike_Single (x,info);
  for ( unsigned int pt=0 ; pt<data->n_unique_pts ; pt++){
    likelihood_neutral[pt] = -(tree->tree)->scalefactor - log (info->p[pt]);
  }

  
  puts ("# Calculating conservation at each site. This may take a while.");
  col = 0;
  done_usite = calloc ( data->n_unique_pts, sizeof(int));
  for ( unsigned int site=0 ; site<data->n_unique_pts ; site++){
    done_usite[site] = -1;
  }

  for ( unsigned int site=0 ; site<data->n_pts ; site++){
    double fm,fn;
    double lb,ub;
    double omegam;
    int type;

    if ( col%50 == 0){
      printf ("\n%4d:  ",col+1);
    }
    col++;

    // Is site all gaps?`
    if ( data->index[site] == -INT_MAX){
      omegam = 1.;
      fm = 0.;
      fn = 0.;
      type = 0;
      if ( dosupport ){
      	lb = 0.;
      	ub = HUGE_VAL;
      }
    }
    // Does site only exist in one sequence
    else if ( data->index[site] < 0){
      omegam = 1.;
      fm = -log(model->pi[-data->index[site]-1]);
      fn = fm;
      type = 1;
      if ( dosupport ){
      	lb = 0.;
      	ub = HUGE_VAL;
      }
      //printf ("%5d recent insert\n",site);
    }
    else if (done_usite[data->index[site]]!=-1){
      int usite = done_usite[data->index[site]];
      fn = selinfo->llike_neu[usite];
      fm = selinfo->llike_max[usite];
      omegam = selinfo->omega_max[usite];
      if ( dosupport ){
      	lb = selinfo->lbound[usite];
      	ub = selinfo->ubound[usite];
      }
      type = selinfo->type[usite];
      //printf ("%5d same as site %d\n",site,usite);
    } else {
      int start;
      // General case
      CopySiteToDataSet (data, data_single, site);
      add_data_to_tree (data_single, tree, model);
      start = FindBestX (likelihood_grid, data->index[site], GRIDSIZE);
      fn = likelihood_neutral[data->index[site]];

      bd[0] = (start>0)?vget(omega_grid,start-1):0.;
      if(positive){ bd[0]=1.;}
      bd[1] = (start<GRIDSIZE-1)?vget(omega_grid,start+1):99.;
      x[0] = vget(omega_grid,start);
      
      int neval=0;
      fm = linemin_1d ( CalcLike_Single, x, (void *)info, bd[0], bd[1], 1e-5,0,&neval);
      omegam = model->param[1];
      if ( IsConserved(data,site) ){
	type = 4;
      } else if ( IsSiteSynonymous(data, site, data->gencode)){
	type = 2;
      } else {
	type = 3;
      }

      /*  Find confidence interval for omega (actually "support") */
      if ( dosupport ){
      	neval=0;
      	if ( likelihood_grid[data->index[site]*GRIDSIZE]-fm<=ldiff/2. ){
        	lb = (double)positive;
      	} else {
                double initial_lb = (double)positive;
        	Set_CalcLike_Wrapper (CalcLike_Single,fm+ldiff/2.);
        	lb = find_root (initial_lb,omegam,CalcLike_Wrapper,(void*)info,NULL,NULL,1e-3,&neval);
      	}

        if ( likelihood_grid[data->index[site]*GRIDSIZE + GRIDSIZE - 1]-fm<=ldiff/2. ){
          ub = 99.;
	} else {
          Set_CalcLike_Wrapper (CalcLike_Single,fm+ldiff/2.);
          neval = 0;
          ub = find_root (omegam,99.,CalcLike_Wrapper,(void*)info,NULL,NULL,1e-3,&neval);
      	}
      }


      assert(data->index[site]>=0);
      done_usite[data->index[site]] = site;
    }

    selinfo->llike_neu[site] = fn;
    selinfo->llike_max[site] = fm;
    selinfo->omega_max[site] = omegam;
    if ( dosupport ){
    	selinfo->lbound[site] = lb;
    	selinfo->ubound[site] = ub;
    }
    selinfo->type[site] = type;
    putchar('.');
    fflush(stdout);
  }
  free(done_usite);
  free(likelihood_grid);
  free(likelihood_neutral);
  free_vec(omega_grid);
  putchar('\n');

  return selinfo;
}

void PrintResults ( char * outfile, struct selectioninfo * selinfo, const double * pval, const double * pval_adj, const int nsites){
  int site,i;
  FILE * out_fp;
  double stat;
  char result[7],sign;

  assert(NULL!=outfile);
  assert(NULL!=selinfo);
  assert(NULL!=pval);
  assert(NULL!=pval_adj);
  assert(nsites>0);

  const int dosupport = (NULL==selinfo->lbound)?0:1;
  assert ( !dosupport || (NULL!=selinfo->ubound && NULL!=selinfo->lbound) );
  assert (  dosupport || (NULL==selinfo->ubound && NULL==selinfo->lbound) );
  

  out_fp = fopen (outfile,"w");
  if ( NULL==out_fp){
    fprintf (stderr,"Could not open file \"%s\" for output\n",outfile); 
    return;
  }

  if ( dosupport ){
  	fputs ("# Site  Neutral  Optimal   Omega    lower    upper LRT_Stat    Pval     Adj.Pval    Q-value Result Note\n", out_fp);
  } else {
	fputs ("# Site  Neutral  Optimal   upper LRT_Stat    Pval     Adj.Pval    Q-value Result Note\n", out_fp);
  }

  for ( site=0 ; site<nsites ; site++){
    stat = 2. * (selinfo->llike_neu[site] - selinfo->llike_max[site]);
    for ( i=0 ; i<6 ; i++){
      result[i] = ' ';
    }
    result[6] = '\0';
    if (selinfo->omega_max[site]>1.){
      sign = '+';
    } else {
      sign = '-';
    }
    if ( pval[site]<=0.05)	result[0] = sign;
    if ( pval[site]<=0.01)	result[1] = sign;
    if ( pval_adj[site]<=0.05)	result[2] = sign;
    if ( pval_adj[site]<=0.01)	result[3] = sign;
   
    if ( dosupport){ 
    	fprintf (out_fp, " %4d %8.2f %8.2f %8.4f %8.4f %8.4f %8.4f %6.4e %6.4e %6.4e %s %s\n",
		site + 1, selinfo->llike_neu[site], selinfo->llike_max[site], 
		selinfo->omega_max[site], selinfo->lbound[site], selinfo->ubound[site], 
		stat, pval[site], pval_adj[site],pval_adj[site+nsites],result,OutString[selinfo->type[site]]);
    } else {
        fprintf (out_fp, " %4d %8.2f %8.2f %8.4f %8.4f %6.4e %6.4e %6.4e %s %s\n",
                site + 1, selinfo->llike_neu[site], selinfo->llike_max[site],
                selinfo->omega_max[site], stat, pval[site], pval_adj[site], pval_adj[site+nsites], result,
		OutString[selinfo->type[site]]);
    }
  }

  fclose(out_fp);
}



double * CalculatePvals ( const double * lmax, const double * lneu, const int n, const int positive_only){
  int site;
  double * pval;
  double x;
  assert(NULL!=lmax);
  assert(NULL!=lneu);
  assert(n>1);

  pval = malloc(n*sizeof(double));
  OOM(pval);
  for ( site=0 ; site<n ; site++){
    x = -2.*(lmax[site]-lneu[site]);
    if ( x<0.)
      x = 0.;
    pval[site] = pchisq(x,1.,1);
  }

  /*  If positive only, then pvals are from chisq-bar distribution
   * (exactly half those from chisq_1).
   */
  if ( positive_only ){
    for ( site=0 ; site<n ; site++){
      pval[site] /= 2.;
    }
  }

  return pval;
}

double * AdjustPvals ( const double * pval, DATA_SET * data){
  int site,idx;
  int n_idx;
  double * adj_tmp, *adj, *qval_tmp;
  assert(NULL!=pval);
  CheckIsDataSet (data);

  adj = malloc(2*data->n_pts*sizeof(double));
  OOM(adj);
 
  n_idx = 0; 
  for ( site=0 ; site<data->n_pts ; site++){
    idx = data->index[site];
    if ( idx>=0){
      adj[n_idx] = pval[site];
      adj[n_idx+data->n_pts] = pval[site];
      n_idx++;
    }
  }
  
  adj_tmp = Pvalue_adjust_StepUp ( adj, n_idx, BONFERRONI);
  qval_tmp = qvals_storey02(adj+data->n_pts,n_idx,DELTA);
  
  n_idx = 0;
  for ( site=0 ; site<data->n_pts ; site++){
    idx = data->index[site];
    if ( idx>=0){
      adj[site] = adj_tmp[n_idx];
      adj[site+data->n_pts] = qval_tmp[n_idx];
      n_idx++;
    } else {
      adj[site] = 1.;
      adj[site+data->n_pts] = 1.;
    }
  }

  free(adj_tmp);
  free(qval_tmp);
  return adj;
}


void PrintSummary ( FILE * file, const struct selectioninfo * selinfo, const double * pval, const double * pval_adj,const int n_pts){
  const double * omegam, * lliken, * llikem;
  int site, npos[4]={0,0,0,0}, ncons[4]={0,0,0,0};
  assert(NULL!=selinfo);
  assert(NULL!=selinfo->llike_neu);
  assert(NULL!=selinfo->llike_max);
  assert(NULL!=selinfo->omega_max);
  assert(NULL!=pval);
  assert(NULL!=pval_adj);
  assert(n_pts>0);

  omegam = selinfo->omega_max;
  lliken = selinfo->llike_neu;
  llikem = selinfo->llike_max;

  for ( site=0 ; site<n_pts ; site++){
    if ( omegam[site]>1.){
      if (pval_adj[site]<0.01){
	npos[0]++;
      }
      if (pval_adj[site]<0.05){
	npos[1]++;
      }
      if (pval[site]<0.01){
	npos[2]++;
      }
      if (pval[site]<0.05){
	npos[3]++;
      }
    } else if ( omegam[site]<1.){
      if (pval_adj[site]<0.01){
        ncons[0]++;
      }
      if (pval_adj[site]<0.05){
        ncons[1]++;
      }
      if (pval[site]<0.01){
        ncons[2]++;
      }
      if (pval[site]<0.05){
        ncons[3]++;
      }
    }
  }

  fprintf (file, "# Positively selected sites (cumulative)\n");
  fprintf (file, "# Significance  Number sites\n");
  fprintf (file, "# 99%% corrected  %5d\n",npos[0]);
  fprintf (file, "# 95%% corrected  %5d\n",npos[1]);
  fprintf (file, "# 99%%            %5d\n",npos[2]);
  fprintf (file, "# 95%%            %5d\n",npos[3]);
  fputc ('\n',file);
  fprintf (file, "# Conserved sites (cumulative)\n");
  fprintf (file, "# Significance  Number sites\n");
  fprintf (file, "# 99%% corrected  %5d\n",ncons[0]);
  fprintf (file, "# 95%% corrected  %5d\n",ncons[1]);
  fprintf (file, "# 99%%            %5d\n",ncons[2]);
  fprintf (file, "# 95%%            %5d\n",ncons[3]);

}



double (*CalcLike_Wrapper_fun)(const double *,void *);
double CalcLikeWrapper_const;

double CalcLike_Wrapper ( const double * x, void * info){
  return (CalcLike_Wrapper_fun (x,info) - CalcLikeWrapper_const);
}

void Set_CalcLike_Wrapper ( double (*f)(const double *,void *), double diff){
  CalcLike_Wrapper_fun = f;
  CalcLikeWrapper_const = diff;
}


#define OMEGAMAX	50.0
#define OMEGAEXPCONST	0.5

VEC create_grid ( const unsigned int len, const int positive){
	assert(len>1);
	assert(0==positive || 1==positive);
	
	VEC grid = create_vec (len);
	const double expconst = (positive==0) ? (OMEGAMAX / expm1(OMEGAEXPCONST*(double)(len-1)))
                                              : ((OMEGAMAX-1.) / expm1(OMEGAEXPCONST*(double)(len-1)));
	for ( unsigned int i=0 ; i<len ; i++){
		vset(grid,i,expconst*expm1(OMEGAEXPCONST*(double)i)+((positive==0)?0:1));
	}
	return grid;
}

/* Two-tailed z-test. Returns p-value based on normal distribution */
double ztest_twotail ( const double stat, const double mean, const double var){
	double zstat = (stat-mean)/sqrt(var);
	return erfc(fabs(zstat)/sqrt(2.));
}


double * change_variance_basis ( const double * m, const double * var, const unsigned int n){
	assert (NULL!=m);
	assert (NULL!=var);
	double * tmp_mat = calloc(n*n,sizeof(double));
	double * result = calloc(n*n,sizeof(double));
	Matrix_MatrixT_Mult (var,n,n,m,n,n,tmp_mat);
	Matrix_Matrix_Mult (m,n,n,tmp_mat,n,n,result);
	free(tmp_mat);
	return result;
}

VEC matrix_vec_mult ( const double * m, const unsigned int nrow, const unsigned int ncol, VEC v){
	assert(NULL!=m);
	assert(NULL!=v);
	assert(vlen(v)==ncol);
	VEC res = create_vec(nrow);
	for ( unsigned int row=0 ; row<nrow ; row++){
		double dot = 0.;
		for ( unsigned int col=0 ; col<ncol ; col++){
			dot += m[row*ncol+col] * vget(v,col);
		}
		vset(res,row,dot);
	}
	return res;
}

double * subsqrmatrix ( const double * m, const unsigned int n, const bool * insub){
	assert(NULL!=m);
	assert(NULL!=insub);

	const unsigned int dimsub = sum_bool(insub,n);
	double * subm = malloc(dimsub*dimsub*sizeof(double));
	for ( unsigned int i=0,si=0 ; i<n ; i++){
		if ( true==insub[i]){
			for ( unsigned int j=0,sj=0 ; j<n ; j++){
				if ( true==insub[j] ){
					subm[si*dimsub+sj] = m[i*n+j];
					sj++;
				}
			}
			si++;
		}
	}

	return subm;
}

double * get_scaling_basis ( const unsigned int root, const TREE * tree, const unsigned int nparam){
	assert(NULL!=tree);
	assert(root<tree->n_br);
	assert(tree->n_br<=nparam);

	const unsigned int nbr = tree->n_br;
	const NODE * rnode = tree->branches[root];
	double * basis = calloc(nparam*nparam,sizeof(double));

	/*  Left and right subtrees from root */
	get_scaling_basis_sub(basis,rnode->branch[0],rnode,nparam,true);
	get_scaling_basis_sub(basis,rnode,rnode->branch[0],nparam,true);
	/* Last element of basis is scaling of root node */
	basis[root*nparam+root] = 1.;
	/* Parameters other than branch lengths are left the same */
	for ( unsigned int param=nbr ; param<nparam ; param++){
		basis[param*nparam+param] = 1.;
	}
	return basis;
}

unsigned int get_scaling_basis_sub (double * basis, const NODE * node, const NODE * parent, const unsigned int nparam,const bool from_root){
	assert(NULL!=basis);
	assert(NULL!=node);
	/* Depth first recursion */
	const unsigned int mybran = branch_number(node,parent);
	for ( unsigned int a=0 ; a<node->nbran ; a++){
		if( NULL==CHILD(node,a) || parent==CHILD(node,a) ){continue;}
		const unsigned int childbran = get_scaling_basis_sub(basis,CHILD(node,a),node,nparam,false);
		if (!from_root){
			for ( unsigned int i=0 ; i<nparam ; i++){
				/*  For checking. Should be just a copy*/
				basis[mybran*nparam+i] += basis[childbran*nparam+i];
			}
		}
	}

	if ( !from_root){ /* Root is special case */
		basis[mybran*nparam+mybran] = 1.;
	}

	return mybran;
}

bool * analyse_branches ( const TREE * tree, const VEC blength, const double ds_scale, const double thresh){
	assert(NULL!=tree);
	assert(NULL!=blength);
	assert(ds_scale>0.);
	assert(thresh>0.);
	const unsigned int nbr = tree->n_br;

	fputs("# Analysing branches for synonymous saturation\n",stdout);
	bool * unsaturated = calloc(nbr,sizeof(bool));
	unsigned int br_count = 0;
	for ( unsigned int br=0 ; br<nbr ; br++){
		if ( vget(blength,br) > thresh / ds_scale ){
			fprintf (stdout,"Branch %d is saturated. Synonymous length = %e (%e)\n",br,vget(blength,br)*ds_scale,vget(blength,br));
			br_count++;
		} else {
			unsaturated[br] = true;
		}
	}
	if ( 0==br_count){
		fputs("# No branches found to be saturated\n",stdout);
	}
	return unsaturated;
}

/*  Blanks branches that should not be kept. Note that branch arrays
 * are no longer correctly NULL terminated (rely on node->nbran
 * instead).
 */
TREE * remove_branches_from_tree ( TREE * tree, const bool * keep ){
	assert(NULL!=tree);
	assert(NULL!=keep);
	const unsigned int nbr = tree->n_br;

	for ( unsigned int br=0 ; br<nbr ; br++){
		if ( false==keep[br] ){
			const NODE * node = tree->branches[br];
			const unsigned int a = find_connection(node->branch[0],node);
			assert(a>=0 && a<node->branch[0]->nbran);
			node->branch[0]->branch[a] = NULL;
			node->branch[0] = NULL;
		}
	}
	return tree;
}


unsigned int print_subtrees ( FILE * fp, const TREE * tree){
	assert(NULL!=tree);

	bool * seen = calloc(tree->n_br,sizeof(bool));
	unsigned int nsubtree = 0;
	for ( RBITER iter=iter_rbtree(tree->leaves) ; next_rbtree(iter) ; ){
		NODE * leaf = (NODE *) itervalue_rbtree(iter);
		/*  Has leaf already been printed as part of tree, or is it
		 * an orphan (branch leading to it deleted).
		 */
		if ( (! seen[leaf->bnumber]) && NULL!=leaf->branch[0]){
			print_tree (fp,leaf,NULL,tree, seen);
			nsubtree++;
		}
	}
	free(seen);
	return nsubtree;
}

/*  Internal nodes with only edge leading from them lead nowhere.
 * Delete.
 */
TREE * remove_orphan_branches ( TREE * tree ){
	assert(NULL!=tree);

	unsigned int ndel=0;
	do {
		ndel = 0;
		for ( unsigned int br=0 ; br<tree->n_br ; br++){
			ndel += remove_orphan_branches_sub(tree->branches[br]);
		}
		ndel += remove_orphan_branches_sub(tree->tree);
	} while (ndel!=0);

	return tree;
}

unsigned int remove_orphan_branches_sub ( const NODE * node ){
	assert(NULL!=node);
	const unsigned int nedge = node->nbran;
	unsigned int ndel = 0;
	/* Case: is leaf */
	if(1==nedge){ return ndel;}
	/* Count number of edges left after deletion */
	/* Note: is it enough to consider branches that only have parent left? */
	unsigned int edges_left=0;
	for ( unsigned int edge=0 ; edge<nedge ; edge++){
		if ( NULL!=node->branch[edge] ){ edges_left++;}
	}
	/*  One edge left but not leaf means the branch leads
	 * nowhere and so is an orphan */
	if ( 1==edges_left ){
		unsigned int deledge = 0;
		for ( unsigned int edge=0 ; edge<nedge ; edge++){
			if ( NULL!=node->branch[edge] ){
				deledge = edge;
				break;
			}
		}
		unsigned int conn = find_connection(node->branch[deledge],node);
		assert(conn>=0 && conn<node->branch[deledge]->nbran);
		node->branch[deledge]->branch[conn] = NULL;
		node->branch[deledge] = NULL;
		ndel++;
	}
	return ndel;
}


/*  Remove internal nodes that only have two edges */
TREE * delete_orphan_nodes ( TREE * tree){
	assert(NULL!=tree);

	const unsigned int nbr = tree->n_br;
	for ( unsigned int br=0 ; br<nbr ; br++){
		assess_orphan_node ( tree->branches[br] );
	}
	assess_orphan_node(tree->tree);

	return tree;
}

unsigned int count_edges ( const NODE * node){
	assert(NULL!=node);
	unsigned int nedges = 0;
	for ( unsigned int edge=0 ; edge<node->nbran ; edge++){
		if (NULL!=node->branch[edge]){nedges++;}
	}
	return nedges;
}

void assess_orphan_node ( const NODE * node){
	assert(NULL!=node);
	/*  Case: node is leaf */
	if ( node->nbran == 1){return;}
	unsigned int nedges = count_edges(node);
	if ( nedges==2){
		unsigned int edge=0;
		for ( ; edge<node->nbran ; edge++){ if(node->branch[edge]!=NULL) break;}
		const unsigned int fst_edge = edge++;
		for ( ; edge<node->nbran ; edge++){ if(node->branch[edge]!=NULL) break;}
		const unsigned int snd_edge = edge;
		assert(fst_edge!=snd_edge);
		assert(edge!=node->nbran);
		const double length = node->blength[fst_edge] + node->blength[snd_edge];
		const unsigned int fst_connect = find_connection(node->branch[fst_edge],node);
		assert(fst_connect>=0 && fst_connect<node->branch[fst_edge]->nbran);
		const unsigned int snd_connect = find_connection(node->branch[snd_edge],node);
		assert(snd_connect>=0 && snd_connect<node->branch[snd_edge]->nbran);
		node->branch[fst_edge]->branch[fst_connect] = node->branch[snd_edge];
		node->branch[fst_edge]->blength[fst_connect] = length;
		node->branch[snd_edge]->branch[snd_connect] = node->branch[fst_edge];
		node->branch[snd_edge]->blength[snd_connect] = length;
		node->branch[fst_edge] = NULL;
		node->branch[snd_edge] = NULL;
	}
}


void prune_into_subtrees ( TREE * tree, const bool * keep, const char * subtreefile){
	assert(NULL!=tree);
	assert(NULL!=keep);
	assert(NULL!=subtreefile);

	FILE * fp = fopen(subtreefile,"w");
	if ( NULL==fp){
		fprintf(stderr,"Failed to open subtree file \"%s\" for output!\n", subtreefile);
		fp = stdout;
	}
	remove_branches_from_tree (tree,keep);
	remove_orphan_branches (tree);
	delete_orphan_nodes (tree);
	fputs("# Remaining subtrees\n",stdout);
	const unsigned int nsubtree = print_subtrees(fp,tree);
	fprintf(stdout,"#\tSubtrees found = %d\n",nsubtree);
	if ( fp!=stdout){
		fprintf(stdout,"# Subtrees saved in file \"%s\"\n",subtreefile);
		fclose(fp);
	}
	exit(EXIT_SUCCESS);
}

double * max_tipdistance ( const TREE * tree){
	assert(NULL!=tree);
	const unsigned int nbr = tree->n_br;
	#ifndef NDEBUG
		/*  Code to check that tree has branch lengths  */
		for ( unsigned int br=0 ; br<nbr ; br++){
			assert(tree->branches[br]->blength[0]>=0.);
		}
	#endif

	double * LRmax = calloc(2*nbr,sizeof(double));
	for ( unsigned int br=0 ; br<2*nbr ; br++){
		LRmax[br] = -1.;
	}

	/*  Note: want to start on leaf? */
	max_tipdistance_pass1 ( tree->tree, NULL, LRmax );
	max_tipdistance_pass2 ( tree->tree, NULL, NULL, LRmax, LRmax+nbr );
	return LRmax;
}

void max_tipdistance_pass1 (const NODE * node, const NODE * parent, double * LRmax){
	assert(NULL!=node);
	assert(NULL!=LRmax);

	/*  Depth first recursion through tree */
	double maxlen = 0.;
	for ( unsigned int a=0 ; a<node->nbran ; a++){
		if ( NULL==CHILD(node,a) || parent==CHILD(node,a) ){ continue; } 
		max_tipdistance_pass1 ( CHILD(node,a), node, LRmax );
		const double len = LRmax[CHILD(node,a)->bnumber] + node->blength[a];
		maxlen = (maxlen>len)?maxlen:len;
	}
	if ( NULL!=parent) LRmax[node->bnumber] = maxlen;
}

void max_tipdistance_pass2(const NODE * node, const NODE * parent, const NODE * gparent, double *left, double * right){
	assert(NULL!=node);
	assert(NULL!=left && NULL!=right);

	if ( NULL!=parent ){
		double maxlen = 0.;
		for ( unsigned int a=0 ; a<parent->nbran ; a++){
			const NODE * cnode = CHILD(parent,a);
			if ( NULL==cnode || node==cnode){ continue;}
			const double len = parent->blength[a] + ((cnode!=gparent)? left[cnode->bnumber] : right[parent->bnumber]);
			maxlen = (maxlen>len)?maxlen:len;
		}
		right[node->bnumber] = maxlen;
	}
	for ( unsigned int a=0 ; a<node->nbran ; a++){
		if ( NULL==CHILD(node,a) || parent==CHILD(node,a) ){ continue; }
		max_tipdistance_pass2 ( CHILD(node,a), node, parent, left, right );
	}
}

struct midpoint * midpoint_tree ( const TREE * tree){
	assert(NULL!=tree);

	const unsigned int nbr = tree->n_br;
	double * left = max_tipdistance(tree);
	const double * right = left + nbr;

	struct midpoint * midpt = malloc(sizeof(struct midpoint));	
	for ( unsigned int br=0 ; br<nbr ; br++){
		const double len = tree->branches[br]->blength[0];
		/*  Midpt satifies max(left,right)<=0.5*(left+right+len)
		 * which reduces to below
		 */
		if ( fabs(left[br]-right[br]) < len ){
			midpt->branch = br;
			midpt->left = left[br];
			midpt->right = right[br];
			midpt->length = len;
			break;
		}
	}

	free(left);
	return midpt;
}

unsigned int branch_number ( const NODE * node, const NODE * parent){
	assert(NULL!=node);

	/*  Traversal consistent with original rooting */
	if ( parent==node->branch[0]){ 
		if ( parent->branch[0]==node){
			// Either node or parent is root
			return (node->bnumber<parent->bnumber)?node->bnumber:parent->bnumber;
		}
		return node->bnumber;
	}

	return parent->bnumber;
}
