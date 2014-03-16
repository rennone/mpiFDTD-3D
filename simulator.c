#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "simulator.h"
#include "field.h"
#include "mpiTM_UPML.h"
#include "mpiTE_UPML.h"
#include "drawer.h"
#include "models.h"
#include <sys/time.h>

static void (*update)() = NULL;
static double complex* (*getDataX)() = NULL;
static double complex* (*getDataY)() = NULL;
static double complex* (*getDataZ)() = NULL;
static void (* finishMethod)() = NULL;
static void (* initMethod)() = NULL;
static char folderName[256];
static struct timeval timer1, timer2;

static void setTMupml(){
  update = fdtdTM_upml_getUpdate();
  initMethod = fdtdTM_upml_getInit();
  finishMethod = fdtdTM_upml_getFinish();
  getDataX = fdtdTM_upml_getHx;
  getDataY = fdtdTM_upml_getHy;
  getDataZ = fdtdTM_upml_getEz;
  strcpy(folderName, "data/TMupml/");
}

static void setTEupml(){
  update = fdtdTE_upml_getUpdate();
  initMethod = fdtdTE_upml_getInit();
  finishMethod = fdtdTE_upml_getFinish();
  getDataX = fdtdTE_upml_getEx;
  getDataY = fdtdTE_upml_getEy;
  getDataZ = fdtdTE_upml_getHz;
  strcpy(folderName, "data/TEupml/");
}

static void setSolver(enum SOLVER solver)
{
  switch(solver){
  case TM_UPML_2D:
    setTMupml();
    break;
  case TE_UPML_2D:
    setTEupml();
    break;
  default:
    break;
  }
}

void simulator_calc(){
  (*update)();
  
  field_nextStep();   //時間を一つ進める  
}

void simulator_init(){
  //横幅(nm), 縦幅(nm), 1セルのサイズ(nm), pmlレイヤの数, 波長(nm), 計算ステップ
  setField(256, 256, 1, 10, 60, 2500); //必ず最初にこれ

  /*NO_MODEL. MIE_CYLINDER, SHELF(todo), NONSHELF(todo) */
  setModel(MIE_CYLINDER);  //次にこれ,  モデル(散乱体)を定義

  setSolver(TE_UPML_2D);//Solverの設定

  (*initMethod)(); //Solverの初期化, EPS, Coeffの設定

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
