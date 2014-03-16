#ifndef _FIELD_H
#define _FIELD_H
#include <stdio.h>
#include <complex.h>
#include "bool.h"

//入射波のモード
enum WAVE_MODE{
  POINT_LIGHT_IN_CENTER,  //中心に点光源
  SCATTER //散乱波
};

typedef struct ntffInfo{
  int top, bottom, left, right;
  int cx,cy;
  double RFperC;
  int step;
} NTFFInfo;

//シミュレーション上の物理定数 
static const double LIGHT_SPEED_S = 0.7;
static const double EPSILON_0_S = 1.0;
static const double MU_0_S = 1.0/0.7/0.7;
static const double Z_0_S = 1.42857142857; //√(1.0/0.7/0.7/1.0) = √(μ/ε);

extern int N_X;
extern int N_Y;
extern int N_Z;
extern int N_PML;
extern int N_PX;
extern int N_PY;
extern int N_PZ;
extern int N_CELL;

extern inline int field_getOffsetX();
extern inline int field_getOffsetY();
extern inline int field_getOffsetZ();

extern inline int field_getSubNx();
extern inline int field_getSubNy();
extern inline int field_getSubNz();

extern inline int field_getSubNpx();
extern inline int field_getSubNpy();
extern inline int field_getSubNpz();
extern inline int field_getSubNcell();

//インデックスを取ってくる 
extern inline int ind(const int, const int, const int);

//フィールドの横,縦の大きさ, 1セルのサイズ, pmlレイヤの数, 波長(nm), 計算ステップ
extern void setField(const int wid, const int hei, const int dep, const double h, const int pml, const double lambda, const double step);

//pml用のσを取ってくる
extern inline double field_sigmaX(double x, double y, double z);
extern inline double field_sigmaY(double x, double y, double z);
extern inline double field_sigmaZ(double x, double y, double z);
extern inline double field_toCellUnit(const double);
extern inline double field_toPhisycalUnit(const double);

//---------------入射波---------------
extern inline double complex field_pointLight(void);

extern inline void field_nextStep(void);
extern inline bool field_isFinish(void);


//:getter
extern inline double field_getK(void);
extern inline double field_getRayCoef(void);
extern inline double field_getOmega(void);
extern inline double field_getLambda(void);
extern inline double field_getWaveAngle(void);
extern inline double field_getTime(void);
extern inline double field_getMaxTime(void);
extern inline NTFFInfo field_getNTFFInfo(void);

//output method
extern void field_outputElliptic(const char *fileName,double complex* data); //
extern void field_outputAllDataComplex(const char *fileName,double complex* data); //
extern void field_outputAllDataDouble(const char *fileName,double* data); //

#endif
