#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "simulator.h"
#include "field.h"
#include "drawer.h"
#include "models.h"
#include <sys/time.h>


static double complex* (*getDataX)() = NULL;
static double complex* (*getDataY)() = NULL;
static double complex* (*getDataZ)() = NULL;

static void (*updateMethod)() = NULL;
static void (* finishMethod)() = NULL;
static void (* initMethod)() = NULL;
static void (* resetMethod)() = NULL;

static double complex* (*getDrawData)() = NULL;
static double* (* getEpsMethod )() = NULL;

static char folderName[256];
static struct timeval timer1, timer2;


static void setMPI3D()
{
  updateMethod = mpi_fdtd3D_upml_getUpdate();
  initMethod   = mpi_fdtd3D_upml_getInit();
  finishMethod = mpi_fdtd3D_upml_getFinish();
}

static void setSolver(enum SOLVER solver)
{
  switch(solver){
  case MPI_FDTD_3D:
  default:
    setMPI3D();
    break;
  }

  
}

void simulator_calc(){
  (*updateMethod)();
  
  field_nextStep();   //時間を一つ進める  
}

void simulator_init(FieldInfo field_info, enum MODEL model, enum SOLVER solver){
  //横幅(nm), 縦幅(nm), 1セルのサイズ(nm), pmlレイヤの数, 波長(nm), 計算ステップ
  field_init(field_info);

  /*NO_MODEL. MIE_CYLINDER, SHELF(todo), NONSHELF(todo) */
  setModel(model);     //次にこれ,モデル(散乱体)を定義

  setSolver(solver);      //Solverの設定

  (*initMethod)();   //Solverの初期化, EPS, Coeffの設定

   /*POINT_LIGHT_IN_CENTER: 中心に点光源   SCATTER: 散乱波*/
  //field_setDefaultIncidence(SCATTER); //入射波の設定

  gettimeofday(&timer1, NULL);
}

void simulator_finish(){
  printf("finish\n");
  gettimeofday(&timer2,NULL);
  printf("time = %lf \n", timer2.tv_sec-timer1.tv_sec+(timer2.tv_usec-timer1.tv_usec)*1e-6);
  char fileName[256];
  strcpy(fileName, folderName);
  strcat(fileName, "mie.txt");
  //field_outputElliptic(fileName, (*getDataZ)());
  (*finishMethod)(); //メモリの解放等  
}

double complex* simulator_getDrawingData(void){
  return (* getDataZ)();
}

bool simulator_isFinish(void)
{
  return field_isFinish();
}
