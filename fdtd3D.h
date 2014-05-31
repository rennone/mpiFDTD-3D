#ifndef _FDTD_3D_H
#define _FDTD_3D_H
#include "myComplex.h"

extern void (* fdtd3D_getUpdate(void))(void);
extern void (* fdtd3D_getFinish(void))(void);
extern void (* fdtd3D_getReset(void))(void);
extern void (* fdtd3D_getInit(void))(void);

extern dcomplex* fdtd3D_getEx(void);
extern dcomplex* fdtd3D_getEy(void);
extern dcomplex* fdtd3D_getEz(void);
extern dcomplex* fdtd3D_getHx(void);
extern dcomplex* fdtd3D_getHy(void);
extern dcomplex* fdtd3D_getHz(void);
extern double*   fdtd3D_getEps(void);

#endif
