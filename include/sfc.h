
#ifndef SFC_H
#define SFC_H

#include<iostream>
#include<vector>

#include "Point.h"

#define SWAP(a, b) ((&(a) == &(b)) || \
                    (((a) -= (b)), ((b) += (a)), ((a) = (b) - (a))))

extern int G_MAX_DEPTH;
extern int G_dim;




void rotate(int index,char* current,char * rot_index,int dim,bool cal_true_index=true);
//void rotate_table_based(int index,char* current,char * rot_index,int dim);

#endif

