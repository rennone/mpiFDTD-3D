#ifndef _MPI_FDTD_3D_H
#define _MPI_FDTD_3D_H
#include "myComplex.h"

extern void (* mpi_fdtd3D_upml_getUpdate(void))(void);
extern void (* mpi_fdtd3D_upml_getFinish(void))(void);
extern void (* mpi_fdtd3D_upml_getReset(void))(void);
extern void (* mpi_fdtd3D_upml_getInit(void))(void);

extern dcomplex* mpi_fdtd3D_upml_getEx(void);
extern dcomplex* mpi_fdtd3D_upml_getEy(void);
extern dcomplex* mpi_fdtd3D_upml_getEz(void);
extern dcomplex* mpi_fdtd3D_upml_getHx(void);
extern dcomplex* mpi_fdtd3D_upml_getHy(void);
extern dcomplex* mpi_fdtd3D_upml_getHz(void);
extern double*   mpi_fdtd3D_upml_getEps(void);

#endif
