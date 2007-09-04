#ifndef _TREE_H_
#define _TREE_H_

#define MAX_SP          1000
#define MAX_BR          2010

#define ISLEAF(A)       ( A->branch[1] == NULL )
#define CHILD(A,B)      A->branch[B]
#define CHILDP(A,B)     ((A->branch[B])->plik);
#define LENGTH(A)       A->length[B]

#define DBL_EQUALS(A,B) ( fabs((A)-(B))<=2*DBL_EPSILON*(fabs(A)+fabs(B)) )
#ifndef OOM
#define OOM(A) if ( A == NULL){ \
                  printf ("Out of Memory!\n"); \
                  exit (EXIT_FAILURE); }
#endif

#include "../libdict-0.2.1/dict.h"



struct node {
        struct node **   branch;
        double          * blength;
        int             bnumber;
	int		nbran;
	int 		maxbran;
        int             *seq;
	char 		*name;
        /*  Part_lik is array, length N_BASES*N_PTS, indexed by a*N_BASES+b*/
        double          *plik;
        double          *mid,*back,*dback;
        double          *mat, *bmat;
        double          scalefactor,bscalefactor;
        int scale,bscale;
};

typedef struct node NODE;

typedef struct {
        int n_sp;
        int n_br;
        char * tstring;
        NODE * tree;
        NODE * branches[MAX_BR];
        NODE * leaves[MAX_SP];
	dict * leaf_names;
} TREE;


void Recurse_forward ( const TREE * tree, void (*fun)(void *, int,int), void * info);
void Recurse_backward ( const TREE * tree, void (*fun)(void *, int,int), void *info);
void CheckIsTree ( const TREE * tree);
void create_tree (TREE * tree);
void print_tree ( FILE * out, const NODE * node, const NODE * parent, const TREE * tree);
int find_leaf_number ( const NODE * leaf, const TREE * tree);
int find_branch_number ( const NODE * branch, const TREE * tree);
int find_connection ( const NODE * from, const NODE * to);
int add_lengths_to_tree ( TREE * tree, double *lengths);
void PrintBranchLengths (FILE * fp, const TREE * tree);
void ScaleTree ( TREE * tree, const double f);
TREE * CopyTree ( const TREE * tree);
TREE * CloneTree ( TREE * tree);
void FreeTree ( TREE * tree);
void FreeNode( NODE * node, NODE * parent);

TREE ** read_tree_strings ( char * filename);
TREE * copy_tree_strings ( const TREE * tree );
int save_tree_strings ( char * filename, TREE ** trees);
int add_lengths_to_tree ( TREE * tree, double *lengths);

int find_leaf_by_name ( const char * name, const TREE * tree);

#endif

