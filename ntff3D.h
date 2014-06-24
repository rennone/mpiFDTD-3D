#ifndef _NTFF_3D_H
#define _NTFF_3D_H

#include <stdio.h>
#include "myComplex.h"

extern void ntff3D_Frequency( dcomplex *Ex, dcomplex *Ey,dcomplex *Ez,
                              dcomplex *Hx, dcomplex *Hy, dcomplex *Hz);

extern void ntff3D_SubTimeCalc(dcomplex *Ex,dcomplex *Ey,dcomplex *Ez,
                               dcomplex *Hx,dcomplex *Hy,dcomplex *Hz);

//時間領域の遠方解 Ux, Uy, WzからEth, Ephを求め, fpRe, fpImに書き出す
extern void ntff3D_TimeOutput();

extern void ntff3D_Init(void);
#endif

