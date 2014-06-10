#ifndef _NTFF_3D_H
#define _NTFF_3D_H

#include <stdio.h>
#include "myComplex.h"

extern void ntff3D_Frequency( dcomplex *Ex, dcomplex *Ey,dcomplex *Ez,
                              dcomplex *Hx, dcomplex *Hy, dcomplex *Hz);

extern void ntff3D_SubTimeCalc(dcomplex *Ex,dcomplex *Ey,dcomplex *Ez,
                               dcomplex *Hx,dcomplex *Hy,dcomplex *Hz,
                               dcomplex *Ux,dcomplex *Uy,dcomplex *Uz,
                               dcomplex *Wx,dcomplex *Wy,dcomplex *Wz);
#endif

