#ifndef _MPI_FDTD_3D_H
#define _MPI_FDTD_3D_H
#include <complex.h>

extern void (* mpi_fdtd3D_upml_getUpdate(void))(void);
extern void (* mpi_fdtd3D_upml_getFinish(void))(void);
extern void (* mpi_fdtd3D_upml_getReset(void))(void);
extern void (* mpi_fdtd3D_upml_getInit(void))(void);

extern double complex* mpi_fdtd3D_upml_getEx(void);
extern double complex* mpi_fdtd3D_upml_getEy(void);
extern double complex* mpi_fdtd3D_upml_getEz(void);
extern double complex* mpi_fdtd3D_upml_getHx(void);
extern double complex* mpi_fdtd3D_upml_getHy(void);
extern double complex* mpi_fdtd3D_upml_getHz(void);

#endif
