#ifndef _LIKE_H_
#define _LIKE_H_

#define DELTA   1e-6



int CalcLike_Sub ( NODE * node, NODE * parent, TREE * tree, MODEL * model);
int LikeVector ( TREE * tree, MODEL * model, double p[]);
int LikeVectorSub ( TREE * tree, MODEL * model, double p[]);
double LikeDiff ( double scale1, double scale2, double * like1, double * like2, double * freq, int size);
double Like ( double scale, double like[], double freq[], int usize, double * pi , int nsize, int * index);


double CalcLike ( double pt[]);
void DCalcLike ( double pt[], double grad[]);
#endif

