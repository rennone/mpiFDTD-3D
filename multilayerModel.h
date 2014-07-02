#ifndef MULTI_LAYER_MODEL_H
#define MULTI_LAYER_MODEL_H

#include "bool.h"
#include <stdio.h>
#include <stdlib.h>
extern double ( *multilayerModel_EPS(void))(double x, double y, double z, int, int, int);

extern void multilayerModel_init();
extern bool multilayerModel_isFinish(void);
extern void multilayerModel_setThickness(int thickness1_nm, int thickness2_nm);
extern void multilayerModel_needSize(int *x_nm, int *y_nm,int *z_nm);
extern void multilayerModel_moveDirectory();
#endif
