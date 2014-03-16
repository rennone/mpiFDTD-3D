#ifndef _MPI_FDTD_3D_H
#define _MPI_FDTD_3D_H
#include <complex.h>

extern void fdtd3D_update(void);
extern void fdtd3D_init(void);
extern void fdtd3D_finish(void);

extern double complex* fdtd3D_upml_getEx(void);
extern double complex* fdtd3D_upml_getEy(void);
extern double complex* fdtd3D_upml_getEz(void);
extern double complex* fdtd3D_upml_getHx(void);
extern double complex* fdtd3D_upml_getHy(void);
extern double complex* fdtd3D_upml_getHz(void);

extern inline int fdtd3D_getSubNx(void);
extern inline int fdtd3D_getSubNy(void);
extern inline int fdtd3D_getSubNpx(void);
extern inline int fdtd3D_getSubNpy(void);
extern inline int fdtd3D_getSubNcell(void);

#endif
