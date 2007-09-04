#ifndef _NUCMODEL_H_
#define _NUCMODEL_H_

#ifndef _MODEL_H_
#include "model.h"
#endif

MODEL *NewJC69Model_full (int nbr);
MODEL * NewNNNModel_full ( const int * desc, const double * params, const int nparam, const double * pi, const int freq_type, const int nbr, const int alt_scale, const int opt_pi);

const int desc_JC69[] = { -1, -1, -1, -1,
                          -1, -1, -1, -1,
                          -1, -1, -1, -1,
                          -1, -1, -1, -1};
      
const int desc_HKY[]  = { -1, -1,  0, -1,
						  -1, -1, -1,  0,
						   0, -1, -1, -1,
						  -1,  0, -1, -1};
						  
const int desc_REV[]  = { -1, -1,  0, 1,
                          -1, -1,  2, 3,
                           0,  2, -1, 4,
                           1,  3,  4, -1};
#endif
