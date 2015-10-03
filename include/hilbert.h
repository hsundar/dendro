#ifndef HILBERT_H
#define HILBERT_H

#include<iostream>

#include "sfc.h"
#include "../include/Point.h"
#include "binUtils.h"
//#include "rotation.h"

#include "hcurvedata.h"


bool hilbert_order_NCA(const Point& p1,const Point& p2) ;
bool hilbert_order(const Point& p1,const Point& p2) ;


#endif