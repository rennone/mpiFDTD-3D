#ifndef _CIRCLE_H
#define _CIRCLE_H

#include "bool.h"
extern double (*circleModel_EPS(void))(double, double, double,int,int,int);

extern void circleModel_init(void);
extern bool circleModel_isFinish(void);
extern void circleModel_moveDirectory(void);
#endif
