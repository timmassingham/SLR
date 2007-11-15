#ifndef _RBTREE_H_
#define _RBTREE_H_

#ifndef _STDBOOL_H_
#include <stdbool.h>
#endif

enum rbcolour { red,black};

struct __rbnode {
        enum rbcolour colour;
        struct __rbnode *left, *right, *parent;
        void *key,*value;
};
typedef struct __rbnode * RBNODE;

struct __rbtree {
        RBNODE root;
        int (*compfun)(const void *, const void *);
        void * (*copykey)(const void *);
        void (*freekey)(void *);
};
typedef struct __rbtree * RBTREE;


RBTREE create_rbtree ( int (*compfun)(const void *, const void *), void * (*copykey)(const void *), void (*freekey)(void *) );
void * getelt_rbtree ( const RBTREE tree, const void * key);
bool member_rbtree ( const RBTREE tree, const void * key);
void * insertelt_rbtree ( RBTREE tree, void * key, void * value);
void * removeelt_rbtree (RBTREE tree, const void * key);

void map_rbtree ( RBTREE tree, void * (*mapfun)(const void *, void *) );
void * minelt_rbtree ( const RBTREE tree);
void * maxelt_rbtree ( const RBTREE tree);

int lexo ( const void * pt1, const void * pt2);
void * strcopykey(const void * key);
void strfreekey (void * key);

#endif
