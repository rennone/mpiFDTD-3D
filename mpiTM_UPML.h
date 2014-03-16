#ifndef FDTD_TM_UPML_H
#define FDTD_TM_UPML_H
#include <complex.h>

extern void (* fdtdTM_upml_getUpdate(void))(void);
extern void (* fdtdTM_upml_getFinish(void))(void);
extern void (* fdtdTM_upml_getInit(void))(void);
extern double complex* fdtdTM_upml_getHx(void);
extern double complex* fdtdTM_upml_getHy(void);
extern double complex* fdtdTM_upml_getEz(void);
extern double * fdtdTM_upml_getEps(void);

extern inline int fdtdTM_upml_getSubNx(void);
extern inline int fdtdTM_upml_getSubNy(void);
extern inline int fdtdTM_upml_getSubNpx(void);
extern inline int fdtdTM_upml_getSubNpy(void);
extern inline int fdtdTM_upml_getSubNcell(void);

#endif
