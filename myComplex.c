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

extern double complex cbilinear3D(dcomplex *p, double x, double y, double z, int width, int height, int depth)
{
  int i = floor(x);
  int j = floor(y);
  int k = floor(z);
  double dx = x - i;
  double dy = y - j;
  double dz = z - k;

  int xOffset = height*depth;
  int yOffset = depth;
  int index = i*xOffset + j*yOffset + k;
  return p[index]*(1.0-dx)*(1.0-dy)*(1.0-dz)   //p[i,j,k]
    + p[index+xOffset]*dx*(1.0-dy)*(1.0-dz)    //p[i+1,j,k]
    + p[index+yOffset]*(1.0-dx)*dy*(1.0-dz)    //p[i,j+1,k]
    + p[index+1]*(1.0-dx)*(1.0-dy)*dz          //p[i,j,k+1]
    + p[index+xOffset+yOffset]*dx*dy*(1.0-dz)  //p[i+1,j+1,k]
    + p[index+xOffset+1]*dx*(1.0-dy)*dz        //p[i+1,j,k+1]
    + p[index+yOffset+1]*(1.0-dx)*dy*dz        //p[i,j+1,k+1]
    + p[index+xOffset+yOffset+1]*dx*dy*dz;     //p[i,j+1,k+1]
}
