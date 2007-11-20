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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tree.h"
#include "data.h"
#include "model.h"
#include "bases.h"
#include "tree_data.h"
#include "utility.h"


static int memadd_plik_tree ( TREE * tree, const int size, const int exact_obs, const int nbase);
static int memfree_plik_tree ( TREE * tree);
static int memadd_seq_tree ( TREE * tree, const int size);
static int memfree_seq_tree ( TREE * tree);


void add_single_site_to_tree ( TREE * tree, const DATA_SET * data, const MODEL * model, const int a){
        int i,j,leaf;
        double b;

        CheckIsTree(tree);

        for ( i=0 ; i<data->n_sp ; i++){
	    leaf = find_leaf(i,tree,data);
	    (tree->leaves[leaf])->seq[a] = data->seq[i][a];
	}
        if ( model->exact_obs != 1)
                for ( i=0 ; i<data->n_sp ; i++){
			leaf = find_leaf(i,tree,data);
                        if(data->seq[i][a] != GapChar(model->seqtype))
                                b = 0.;
                        else
                                b = 1.;
                        for ( j=0 ; j<model->nbase ; j++)
                                (tree->leaves)[leaf]->plik[j] = b;
                        if (data->seq[i][a] != GapChar(model->seqtype))
                                (tree->leaves)[leaf]->plik[a] = 1.;
                }

        CheckIsTree(tree);
}

int add_data_to_tree (const DATA_SET * data_old, TREE * tree, MODEL * model)
{
  int a, b, c,leaf;
  double *tmp, *param;
  const DATA_SET *data;

  CheckIsTree (tree);
  CheckIsDataSet (data_old);

  data = data_old;

  model->n_unique_pts = data->n_unique_pts;
  model->n_pts = data->n_pts;
  if (model->pt_freq != NULL)
    free (model->pt_freq);
  model->pt_freq = calloc ((size_t) model->n_unique_pts, sizeof (double));
  OOM (model->pt_freq);
  for (b = 0; b < model->n_unique_pts; b++)
    model->pt_freq[b] = data->freq[b];
  if (model->tmp_plik != NULL)
    free (model->tmp_plik);
  model->tmp_plik =
    malloc ((size_t) (model->n_unique_pts * model->nbase) * sizeof (double));
  if (model->index != NULL)
    free (model->index);
  model->index = malloc (model->n_pts * sizeof (int));
  OOM (model->index);
  for (b = 0; b < model->n_pts; b++)
    model->index[b] = data->index[b];



  if (data->n_sp != tree->n_sp) {
    printf ("Error, data and tree have different number of species!\n");
    exit(EXIT_FAILURE);
  }

  (void) memadd_plik_tree (tree, data->n_unique_pts * model->nbase,
			   model->exact_obs, model->nbase);
  (void) memadd_seq_tree (tree, data->n_unique_pts);

  param = model->param;
  for (a = 0; a < data->n_sp; a++) {
    leaf = find_leaf(a,tree,data);
    for (b = 0; b < data->n_unique_pts; b++)
      (tree->leaves[leaf])->seq[b] = data->seq[a][b];

    if (model->exact_obs != 1) {
      for (b = 0, tmp = (tree->leaves[leaf])->plik;
	   b < data->n_unique_pts * model->nbase; b++, tmp++)
	*tmp = 0.;
      tmp = (tree->leaves[leaf])->plik;
      for (b = 0; b < data->n_unique_pts; b++) {
	/*  Deal with non-gap */
	if ((tree->leaves[leaf])->seq[b] != GapChar (model->seqtype))
	  tmp[b * model->nbase + (tree->leaves[leaf])->seq[b]] = 1.0;
	else
	  for (c = 0; c < model->nbase; c++)
	    tmp[b * model->nbase + c] = 1.0;
      }
    }
  }


  CheckIsTree (tree);
  return 0;
}



static int memadd_plik_tree ( TREE * tree, const int size, const int exact_obs, const int nbase){
        int a;

        if ( (tree->tree)->plik != NULL)
                (void)memfree_plik_tree ( tree);

        (tree->tree)->plik = calloc ( (size_t)size, sizeof(double));

        for ( a=0 ; a<tree->n_br ; a++){
                (tree->branches[a])->plik = calloc ( size,sizeof(double));
                (tree->branches[a])->back = calloc ( size,sizeof(double));
                (tree->branches[a])->mid = calloc (size,sizeof(double));
                (tree->branches[a])->dback = calloc (size,sizeof(double));
                (tree->branches[a])->mat = calloc (nbase*nbase,sizeof(double));
                (tree->branches[a])->bmat = calloc (nbase*nbase,sizeof(double));
        }

        return 0;
}

#define Free(A) if(*A != NULL){ free(*A); *A=NULL;}

static int memfree_plik_tree ( TREE * tree){
        int a;
        NODE * node;

        CheckIsTree(tree);

        node = tree->tree;
        Free(&node->plik);
        Free(&node->back);
        Free(&node->mid);
        Free(&node->mat);
        Free(&node->bmat);
        Free(&node->dback);

        for ( a=0 ; a<tree->n_br ; a++){
                node = tree->branches[a];
                Free(&node->plik);
                Free(&node->back);
                Free(&node->mid);
                Free(&node->mat);
                Free(&node->bmat);
                Free(&node->dback);
        }

        return 0;

        free ( (tree->tree)->plik);
        (tree->tree)->plik = NULL;

        for ( a=0 ; a<tree->n_br ; a++){
                if ( (tree->branches[a])->back != NULL)
                        free ( (tree->branches[a])->back );
                if ( (tree->branches[a])->mid !=  NULL)
                        free ( (tree->branches[a])->mid );
                if ( (tree->branches[a])->mat != NULL)
                        free ( (tree->branches[a])->mat );
                if ( (tree->branches[a])->bmat != NULL)
                        free ( (tree->branches[a])->bmat );
                if ( (tree->branches[a])->plik != NULL)
                        free ( (tree->branches[a])->plik );
                if ( (tree->branches[a])->dback != NULL);
                        free ( (tree->branches[a])->dback );
        }

        CheckIsTree(tree);
        return 0;
}


static int memadd_seq_tree ( TREE * tree, const int size){
        int a,b;

        CheckIsTree(tree);

        if ( (tree->tree)->seq != NULL)
                (void)memfree_seq_tree (tree);

        (tree->tree)->seq = calloc ( (size_t)size, sizeof(int));
        OOM ( (tree->tree)->seq );
        for ( a=0 ; a<tree->n_br ; a++){
                (tree->branches[a])->seq = calloc ( (size_t)size, sizeof(int));
                OOM ( (tree->branches[a])->seq );
                for ( b=0 ; b<size ; b++)
                        (tree->branches[a])->seq[b] = 0;
        }

        CheckIsTree(tree);
        return 0;
}


static int memfree_seq_tree ( TREE * tree){
        int a;
        CheckIsTree(tree);

        free ( (tree->tree)->seq);
        (tree->tree)->seq = NULL;
        for ( a=0 ; a<tree->n_br ; a++){
                free ( (tree->branches[a])->seq);
                (tree->branches[a])->seq = NULL;
        }

        CheckIsTree(tree);
        return 0;
}



int find_leaf ( const int i, const TREE * tree, const DATA_SET * data){
   int a;
   char * num;
   CheckIsTree(tree);
   CheckIsDataSet(data);
   assert(tree->n_sp == data->n_sp);
   assert(i>=0 && i<tree->n_sp);
   
   a = find_leaf_by_name(data->sp_name[i],tree);
   if ( -1 == a){ /*  Name does not exist -- leaves may be numbered instead */
      num = itoa((i+1));
      a = find_leaf_by_name(num,tree);
      if( -1 == a){
	 fprintf(stderr,"Cannot find species %s in tree (also tried %s for species num %d).\n",data->sp_name[i],num,i+1);
	 exit(EXIT_FAILURE);
      }
      free(num);
   }

   assert(a>=0 && a<tree->n_sp);
   return a;
}

