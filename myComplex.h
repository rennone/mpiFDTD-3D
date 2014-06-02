#ifndef MY_COMPLEX_H
#define MY_COMPLEX_H

#include <complex.h>

typedef double complex dcomplex;
extern double* newDouble(int size);
extern void freeDouble(double *array);

extern dcomplex* newDComplex(int size);
extern void freeDComplex(dcomplex* array);

extern double cnorm(dcomplex c);
extern double complex cbilinear(dcomplex *p, double x, double y, int index, int toNextX, int toNextY);
extern double complex cbilinear3D(dcomplex *p, double x, double y, double z, int index, int toNextX, int toNextY, int toNextZ);
#endif
