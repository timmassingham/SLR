#ifndef _TREE_DATA_H_
#define _TREE_DATA_H_

#ifndef _TREE_H_
#include "tree.h"
#endif
#ifndef _DATA_H_
#include "data.h"
#endif
#ifndef _MODEL_H_
#include "model.h"
#endif

struct single_fun {
        TREE * tree;
        MODEL * model;
        double * p;
};


int add_data_to_tree ( const DATA_SET * data, TREE * tree, MODEL * model);
void add_single_site_to_tree ( TREE * tree, const DATA_SET * data, const MODEL * model, const int a);
int find_leaf ( const int i, const TREE * tree, const DATA_SET * data);
#endif

