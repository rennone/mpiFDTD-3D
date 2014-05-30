#include "myComplex.h"
#include <stdlib.h>
double* newDouble(int size)
{
  double *array = (double*)malloc(sizeof(double)*size);
  memset(array, 0, sizeof(double)*size);
  return array;
}

void freeDouble(double *array)
{
  if(array != NULL)
  {
    free(array); array=NULL;
  }
}

dcomplex* newDComplex(int size)
{
  dcomplex *array = (dcomplex*)malloc(sizeof(dcomplex)*size);
  memset(array, 0, sizeof(dcomplex)*size);
  return array;
}

void freeDComplex(dcomplex *array)
{
  if(array != NULL)
  {
    free(array); array=NULL;
  }
}

//norm of complex
double cnorm(dcomplex c){
  double re = creal(c);
  double im = cimag(c);
  return re*re + im*im;
}

double complex cbilinear(dcomplex *p, double x, double y, int width, int height)
{
  int i = floor(x);
  int j = floor(y);
  double dx = x - i;
  double dy = y - j;
  int index = i*height + j;
  return p[index]*(1.0-dx)*(1.0-dy)
       + p[index+height]*dx*(1.0-dy)
       + p[index+1]*(1.0-dx)*dy
       + p[index+height+1]*dx*dy;
}
