#ifndef _FIELD_H
#define _FIELD_H
#include <stdio.h>
#include "myComplex.h"
#include "bool.h"
#include <mpi.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

//各面の同期をとるための型
extern MPI_Datatype MPI_DCOMPLEX_XY_PLANE, MPI_DCOMPLEX_YZ_PLANE, MPI_DCOMPLEX_XZ_PLANE;

//入射波のモード
enum WAVE_MODE{
  POINT_LIGHT_IN_CENTER,  //中心に点光源
  SCATTER //散乱波
};

typedef struct NTFFInfo{
  int top, bottom, left, right, front, back;
  int cx,cy,cz;
  double RFperC;
  int step;
  int arraySize;
} NTFFInfo;

//計算領域に関する物理パラメータ
typedef struct FieldInfo
{
  int width_nm, height_nm, depth_nm; // 領域のサイズ
  int h_u_nm;              //1セルの大きさ
  int pml;                 //pmlレイヤの大きさ(セル数)
  int lambda_nm;           //波長
  int theta_deg, phi_deg;  //入射角度
  int stepNum;             //計算ステップ
}FieldInfo;

//プログラムで扱う計算領域のパラメータ
typedef struct FieldInfo_S
{
  int N_X, N_Y, N_Z;    //領域のサイズ(セル)
  int N_PX, N_PY, N_PZ; //PMLレイヤを含めた領域サイズ(セル)
  int N_CELL;           //全セル数
  int N_PML;            //PMLレイヤの層の数
  int N_PYZ;            // N_PY*N_PZ(あらかじめ計算しておく)
}FieldInfo_S;

//MPI分割後の小領域のパラメータ
typedef struct SubFieldInfo_S
{
  int COORDINATES[3];
  //小領域のインデックスが全体領域のどのインデックスと対応するかのオフセット量
  int OFFSET_X, OFFSET_Y, OFFSET_Z;  // sub.x + OFFSET_X = field.xになる(-1する必要なし).
  int SUB_N_X, SUB_N_Y, SUB_N_Z;
  int SUB_N_PX, SUB_N_PY, SUB_N_PZ;
  int SUB_N_CELL;
  int SUB_N_PYZ;   // SUB_N_PY*SUB_N_PZ(あらかじめ計算しておく)
  int Nproc; //全ランク数
  int Rank; //自身のランク
  int RtRank, LtRank, TpRank, BmRank, FtRank, BkRank; //周りの領域のプロセスランク
}SubFieldInfo_S;

//入射波のパラメータ
typedef struct WaveInfo_S
{
  double Lambda_s; //波長
  double T_s;      //周期
  double Omega_s;  //角周波数
  double K_s;      //波数
  double Theta_deg, Phi_deg;   //入射角
} WaveInfo_S;

//シミュレーション上の物理定数
// 3Dでは1/√3が収束条件
#define C_0_S 0.5735
static const double EPSILON_0_S = 1.0;
static const double MU_0_S = 1.0/C_0_S/C_0_S;

//Z0はm0,e0を計算して求まる値なのでコンパイル時に計算できない.
extern double field_getZ_0_S();

//インデックスを取ってくる 
extern int field_index(const int, const int, const int);
extern int field_left(int ind);
extern int field_right(int ind);
extern int field_top(int ind);
extern int field_bottom(int ind);
extern int field_front(int ind);
extern int field_back(int ind);

extern int field_subIndex(const int, const int, const int);
extern int field_subLeft(int ind);
extern int field_subRight(int ind);
extern int field_subTop(int ind);
extern int field_subBottom(int ind);
extern int field_subFront(int ind);
extern int field_subBack(int ind);

extern int field_subToOneX(int i);
extern int field_subToOneY(int j);
extern int field_subToOneZ(int z);

//フィールドの横,縦の大きさ, 1セルのサイズ, pmlレイヤの数, 波長(nm), 計算ステップ
//フィールドの横,縦の大きさ, 1セルのサイズ, pmlレイヤの数, 波長(nm), 計算ステップ
extern void field_init(FieldInfo field_info);
extern void field_reset(void);

extern void setField(const int wid, const int hei, const int dep, const double h, const int pml, const double lambda, const double step);

//pml用のσを取ってくる
extern double field_sigmaX(double x, double y, double z);
extern double field_sigmaY(double x, double y, double z);
extern double field_sigmaZ(double x, double y, double z);
extern double field_pmlCoef(double x, double y);
extern double field_pmlCoef2(double x, double y);

extern double field_toCellUnit(const double);
extern double field_toPhisycalUnit(const double);


//---------------入射波---------------
extern double complex field_pointLight(void);

//散乱波
// gapX, gapY : Ex-z, Hx-zは格子点からずれた位置に配置され散る為,格子点からのずれを送る必要がある.
extern void field_scatteredPulse(dcomplex *p, double *eps, double gapX, double gapY, double gapZ);

extern void field_nextStep(void);
extern bool field_isFinish(void);

//:getter
extern double field_getK(void);
extern double field_getRayCoef(void);
extern double field_getOmega(void);
extern double field_getLambda(void);
extern double field_getWaveAngle(void);
extern double field_getTime(void);
extern double field_getMaxTime(void);
extern double field_getTheta(void);
extern double field_getPhi(void);

extern int field_getOffsetX();
extern int field_getOffsetY();
extern int field_getOffsetZ();

extern int field_getSubNx();
extern int field_getSubNy();
extern int field_getSubNz();
extern int field_getSubNpx();
extern int field_getSubNpy();
extern int field_getSubNpz();
extern int field_getSubNcell();

extern NTFFInfo field_getNTFFInfo(void);
extern WaveInfo_S field_getWaveInfo_S(void);
extern SubFieldInfo_S field_getSubFieldInfo_S(void);
extern FieldInfo_S field_getFieldInfo_S(void);
extern FieldInfo field_getFieldInfo(void);

//output method
//plane => 0 : XY平面. 1 : YZ平面. 2 : XZ平面
extern void field_outputElliptic(const char *fileName,double complex* data, int plane);
extern void field_outputAllDataComplex(const char *fileName,double complex* data); //
extern void field_outputAllDataDouble(const char *fileName,double* data); //

#endif
