
#ifndef MORTON_H
#define MORTON_H

#include "sfc.h"
#include <iostream>

#include "../point/Point.h"
#include "../binOps/binUtils.h"

bool morton_order(Point p1,Point p2);
bool morton_order_NCA(const Point& p1,const Point& p2);



#endif