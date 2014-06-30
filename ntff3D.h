#ifndef _NTFF_3D_H
#define _NTFF_3D_H

#include <stdio.h>
#include "myComplex.h"

extern void ntff3D_Frequency( dcomplex *Ex, dcomplex *Ey,dcomplex *Ez,
                              dcomplex *Hx, dcomplex *Hy, dcomplex *Hz);

//周波数ドメインの遠方解もSubフィールドでやらないと時間的に厳しい(360*360*4面はシミュレーションまわすより時間がかかる)
extern void ntff3D_SubFrequency( dcomplex *subEx, dcomplex *subEy,dcomplex *subEz,
                              dcomplex *subHx, dcomplex *subHy, dcomplex *subHz);

extern void ntff3D_SubTimeCalc(dcomplex *subEx,dcomplex *subEy,dcomplex *subEz,
                               dcomplex *subHx,dcomplex *subHy,dcomplex *subHz);

//時間領域の遠方解 Ux, Uy, WzからEth, Ephを求め, fpRe, fpImに書き出す
extern void ntff3D_TimeOutput();

extern void ntff3D_Init(void);
#endif

