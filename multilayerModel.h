#ifndef MULTI_LAYER_MODEL_H
#define MULTI_LAYER_MODEL_H

#include "bool.h"
#include <stdio.h>
#include <stdlib.h>
extern double ( *multilayerModel_EPS(void))(double x, double y, double z, int, int, int);
extern bool multilayerModel_isFinish(void);
//extern void (*multiLayerModel_output(void))(FILE *, double complex*);

#endif
