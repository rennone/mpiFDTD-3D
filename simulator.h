#ifndef _SIMULATOR_H
#define _SIMULATOR_H

#include "myComplex.h"
#include "bool.h"
#include "field.h"
#include "models.h"

enum SOLVER{
  FDTD_3D,
  MPI_FDTD_3D
};

extern void simulator_init(FieldInfo fInfo, enum MODEL m, enum SOLVER s);
extern void simulator_calc(void);
extern bool simulator_isFinish(void);
extern void simulator_finish(void);
extern void simulator_reset(void);
extern void simulator_solverInit(void); //epsとかも計算し直してから再スタートする用

//モデル(のパラメータ)を変更するので, カレントディレクトリを一段上に行く.
extern void simulator_changeModelAndRestart(void); 
extern double complex* simulator_getDrawingData();
extern double *simulator_getEps();

#endif
