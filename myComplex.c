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

double complex cbilinear(dcomplex *p, double x, double y, int index, int toNextX, int toNextY)
{
  int i = floor(x);
  int j = floor(y);
  double dx = x - i;
  double dy = y - j;
  return p[index]*(1.0-dx)*(1.0-dy)
       + p[index+toNextX]*dx*(1.0-dy)
       + p[index+toNextY]*(1.0-dx)*dy
       + p[index+toNextX+toNextY]*dx*dy;
}

extern double complex cbilinear3D(dcomplex *p, double x, double y, double z, int index, int toNextX, int toNextY, int toNextZ)
{
  int i = floor(x);
  int j = floor(y);
  int k = floor(z);
  double dx = x - i;
  double dy = y - j;
  double dz = z - k;

  return p[index]*(1.0-dx)*(1.0-dy)*(1.0-dz)         //p[i,j,k]
    + p[index+toNextX]*dx*(1.0-dy)*(1.0-dz)          //p[i+1,j,k]
    + p[index+toNextY]*(1.0-dx)*dy*(1.0-dz)          //p[i,j+1,k]
    + p[index+toNextZ]*(1.0-dx)*(1.0-dy)*dz          //p[i,j,k+1]
    + p[index+toNextX+toNextY]*dx*dy*(1.0-dz)        //p[i+1,j+1,k]
    + p[index+toNextX+toNextZ]*dx*(1.0-dy)*dz        //p[i+1,j,k+1]
    + p[index+toNextY+toNextZ]*(1.0-dx)*dy*dz        //p[i,j+1,k+1]
    + p[index+toNextX+toNextY+toNextZ]*dx*dy*dz;     //p[i,j+1,k+1]
}
