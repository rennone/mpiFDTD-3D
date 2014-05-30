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
  int OFFSET_X, OFFSET_Y, OFFSET_Z; //左下からのオフセット量(セル)
  int SUB_N_X, SUB_N_Y, SUB_N_Z;
  int SUB_N_PX, SUB_N_PY, SUB_N_PZ;
  int SUB_N_CELL;
  int SUB_N_PYZ;   // SUB_N_PY*SUB_N_PZ(あらかじめ計算しておく)
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
#define C_0_S 0.7071
static const double LIGHT_SPEED_S = 0.7;
static const double EPSILON_0_S = 1.0;
static const double MU_0_S = 1.0/C_0_S/C_0_S;
static const double Z_0_S = 1.41422712488; //√(1.0/0.7/0.7/1.0) = √(μ/ε);

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
extern inline int field_index(const int, const int, const int);
extern inline int field_subIndex(const int, const int, const int);

//フィールドの横,縦の大きさ, 1セルのサイズ, pmlレイヤの数, 波長(nm), 計算ステップ
//フィールドの横,縦の大きさ, 1セルのサイズ, pmlレイヤの数, 波長(nm), 計算ステップ
extern void field_init(FieldInfo field_info);
extern void field_reset(void);

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
extern double field_getK(void);
extern double field_getRayCoef(void);
extern double field_getOmega(void);
extern double field_getLambda(void);
extern double field_getWaveAngle(void);
extern double field_getTime(void);
extern double field_getMaxTime(void);

extern NTFFInfo field_getNTFFInfo(void);
extern WaveInfo_S field_getWaveInfo_S(void);
extern SubFieldInfo_S field_getSubFieldInfo_S(void);
extern FieldInfo_S field_getFieldInfo_S(void);
extern FieldInfo field_getFieldInfo(void);

//output method
extern void field_outputElliptic(const char *fileName,double complex* data); //
extern void field_outputAllDataComplex(const char *fileName,double complex* data); //
extern void field_outputAllDataDouble(const char *fileName,double* data); //

#endif
