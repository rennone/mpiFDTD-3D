#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "field.h"
#include "models.h"
#include "mpiFDTD3D.h"
#include "myComplex.h"
#include "function.h"

/* about MPI  */
/*
static MPI_Datatype X_DIRECTION_DOUBLE_COMPLEX;
static MPI_Datatype DCOMPLEX_XY; //XY平面
static MPI_Datatype DCOMPLEX_YZ; //YZ平面
static MPI_Datatype DCOMPLEX_XZ; //XZ平面
*/
//岡田さんの論文と同じ空間配置にしてみる
//h = 1, Δt = 1で計算
//系は右手系
//x(left-, right+)
//y(bottom-, top+)
//z(back-, front+)

//Ex(i    ,j+0.5,k+0.5,t    ) -> Ex[i,j,k]
//Ey(i+0.5,     ,k+0.5,t    ) -> Ey[i,j,k]
//Ez(i+0.5,j+0.5,k    ,t    ) -> Ez[i,j,k]

//Hx(i+0.5,j    ,k    ,t+0.5) -> Hx[i,j,k]
//Hy(i    ,j+0.5,k    ,t+0.5) -> Hy[i,j,k]
//Hz(i    ,j    ,k+0.5,t+0.5) -> Hz[i,j,k]

static dcomplex *Ex = NULL, *Ey = NULL, *Ez = NULL;
static dcomplex *Hx = NULL, *Hy = NULL, *Hz = NULL;
static double *EPS_EX=NULL, *EPS_EY=NULL, *EPS_EZ=NULL;

static double *C_EX_HZHY = NULL, *C_EY_HXHZ = NULL, *C_EZ_HYHZ = NULL;
static double *C_HX_EYEZ = NULL, *C_HY_EZEX = NULL, *C_HZ_EXEY= NULL;

//PMLで使う
static dcomplex *Exy = NULL, *Exz = NULL;
static dcomplex *Eyx = NULL, *Eyz = NULL;;
static dcomplex *Ezx = NULL, *Ezy = NULL;

static dcomplex *Hxy = NULL, *Hxz = NULL;
static dcomplex *Hyx = NULL, *Hyz = NULL;
static dcomplex *Hzx = NULL, *Hzy = NULL;

static double *C_EXY = NULL, *C_EXZ = NULL, *C_EXY_HZ = NULL, *C_EXZ_HY = NULL;
static double *C_EYX = NULL, *C_EYZ = NULL, *C_EYX_HZ = NULL, *C_EYZ_HX = NULL;
static double *C_EZX = NULL, *C_EZY = NULL, *C_EZX_HY = NULL, *C_EZY_HX = NULL;

static double *C_HXY = NULL, *C_HXZ = NULL, *C_HXY_EZ = NULL, *C_HXZ_EY = NULL;;
static double *C_HYX = NULL, *C_HYZ = NULL, *C_HYX_EZ = NULL, *C_HYZ_EX = NULL;
static double *C_HZX = NULL, *C_HZY = NULL, *C_HZX_EY = NULL, *C_HZY_EX = NULL;

static void update(void);
static void finish(void);
static void init(void);
static void reset(void);

static void calcE(void);
static void calcH(void);
static void calcE_pml(void);
static void calcH_pml(void);
static void allocateMemories(void);
static void setCoefficient(void);
static void freeMemories(void);

static void initializeElectroMagneticField(void);

dcomplex* fdtd3D_getEx(void){  return Ex;}
dcomplex* fdtd3D_getEy(void){  return Ey;}
dcomplex* fdtd3D_getEz(void){  return Ez;}
dcomplex* fdtd3D_getHx(void){  return Hx;}
dcomplex* fdtd3D_getHy(void){  return Hy;}
dcomplex* fdtd3D_getHz(void){  return Hz;}

double* fdtd3D_getEps(void){
  return EPS_EZ;
}
void (* fdtd3D_getUpdate(void))(void){
  return update;
}
void (* fdtd3D_getFinish(void))(void){
  return finish;
}
void (* fdtd3D_getReset(void))(void){
  return reset;
}
void (* fdtd3D_getInit(void))(void){
  return init;
}

//Public:
//Update
static void update(void)
{  
//  calcH();    
//  calcE();
  calcH_pml();
  calcE_pml();
  
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();
  int i = field_subIndex(sInfo.SUB_N_PX/2, sInfo.SUB_N_PY/2, sInfo.SUB_N_PZ/2);
//  if(field_getTime() == 3)
    Ez[i] += field_pointLight();

}

static void init(){
  allocateMemories();
  setCoefficient();
}

static void reset()
{
}

static void allocateMemories()
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();
  
  Ex = newDComplex(sInfo.SUB_N_CELL);
  Ey = newDComplex(sInfo.SUB_N_CELL);
  Ez = newDComplex(sInfo.SUB_N_CELL);
  
  Hx = newDComplex(sInfo.SUB_N_CELL);
  Hy = newDComplex(sInfo.SUB_N_CELL);
  Hz = newDComplex(sInfo.SUB_N_CELL);

  EPS_EX = newDouble(sInfo.SUB_N_CELL);
  EPS_EY = newDouble(sInfo.SUB_N_CELL);
  EPS_EZ = newDouble(sInfo.SUB_N_CELL);

  /*
  C_EX_HZHY = newDouble(sInfo.SUB_N_CELL);
  C_EY_HXHZ = newDouble(sInfo.SUB_N_CELL);
  C_EZ_HYHZ = newDouble(sInfo.SUB_N_CELL);
  
  C_HX_EYEZ = newDouble(sInfo.SUB_N_CELL);
  C_HY_EZEX = newDouble(sInfo.SUB_N_CELL);
  C_HZ_EXEY = newDouble(sInfo.SUB_N_CELL);
  */

  //pml用
  Exy = newDComplex(sInfo.SUB_N_CELL);
  Exz = newDComplex(sInfo.SUB_N_CELL);
  Eyx = newDComplex(sInfo.SUB_N_CELL);
  Eyz = newDComplex(sInfo.SUB_N_CELL);
  Ezx = newDComplex(sInfo.SUB_N_CELL);
  Ezy = newDComplex(sInfo.SUB_N_CELL);

  Hxy = newDComplex(sInfo.SUB_N_CELL);
  Hxz = newDComplex(sInfo.SUB_N_CELL);
  Hyx = newDComplex(sInfo.SUB_N_CELL);
  Hyz = newDComplex(sInfo.SUB_N_CELL);
  Hzx = newDComplex(sInfo.SUB_N_CELL);
  Hzy = newDComplex(sInfo.SUB_N_CELL);

  C_EXY    = newDouble(sInfo.SUB_N_CELL);
  C_EXZ    = newDouble(sInfo.SUB_N_CELL);
  C_EXY_HZ = newDouble(sInfo.SUB_N_CELL);
  C_EXZ_HY = newDouble(sInfo.SUB_N_CELL);
 
  C_EYX = newDouble(sInfo.SUB_N_CELL);
  C_EYZ = newDouble(sInfo.SUB_N_CELL);
  C_EYX_HZ = newDouble(sInfo.SUB_N_CELL);
  C_EYZ_HX = newDouble(sInfo.SUB_N_CELL);
  
  C_EZX = newDouble(sInfo.SUB_N_CELL);
  C_EZY =newDouble(sInfo.SUB_N_CELL);
  C_EZX_HY = newDouble(sInfo.SUB_N_CELL);
  C_EZY_HX = newDouble(sInfo.SUB_N_CELL);

  C_HXY    = newDouble(sInfo.SUB_N_CELL);
  C_HXZ    = newDouble(sInfo.SUB_N_CELL);
  C_HXY_EZ = newDouble(sInfo.SUB_N_CELL);
  C_HXZ_EY = newDouble(sInfo.SUB_N_CELL);
 
  C_HYX = newDouble(sInfo.SUB_N_CELL);
  C_HYZ = newDouble(sInfo.SUB_N_CELL);
  C_HYX_EZ = newDouble(sInfo.SUB_N_CELL);
  C_HYZ_EX = newDouble(sInfo.SUB_N_CELL);
  
  C_HZX = newDouble(sInfo.SUB_N_CELL);
  C_HZY =newDouble(sInfo.SUB_N_CELL);
  C_HZX_EY = newDouble(sInfo.SUB_N_CELL);
  C_HZY_EX = newDouble(sInfo.SUB_N_CELL);
}

static void finish()
{
  
}

static void initializeElectroMagneticField()
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();
  memset(Ex, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
  memset(Ey, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
  memset(Ez, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);

  memset(Hx, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
  memset(Hy, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
  memset(Hz, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
}

static void setCoefficient()
{
  
  //Hz, Ex, Eyそれぞれでσx, σx*, σy, σy*が違う(場所が違うから)
  double sig_ex_x, sig_ex_y, sig_ex_z;
  double sig_ey_x, sig_ey_y, sig_ey_z;
  double sig_ez_x, sig_ez_y, sig_ez_z;
  double sig_hx_x, sig_hx_y, sig_hx_z;
  double sig_hy_x, sig_hy_y, sig_hy_z;
  double sig_hz_x, sig_hz_y, sig_hz_z;

  
  double R = 1.0e-8;
  double M = 2.0;
  
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  const double sig_max = -(M+1.0)*EPSILON_0_S*C_0_S/2.0/fInfo_s.N_PML*log(R);

  for(int i=1; i<sInfo.SUB_N_PX-1; i++)
    for(int j=1; j<sInfo.SUB_N_PY-1; j++)
      for(int k=1; k<sInfo.SUB_N_PZ-1; k++)
      {
        int w = field_subIndex(i,j,k);
        int x = i+sInfo.OFFSET_X;
        int y = j+sInfo.OFFSET_Y;
        int z = k+sInfo.OFFSET_Z;
        EPS_EX[w] = models_eps(x    ,y+0.5,z+0.5,D_Y); //todo 
        EPS_EY[w] = models_eps(x+0.5,y    ,z+0.5,D_X);
        EPS_EZ[w] = models_eps(x+0.5,y+0.5,z    ,D_XY);

        /*
        C_EX_HZHY[w] = 1.0/EPS_EX[w];
        C_EY_HXHZ[w] = 1.0/EPS_EY[w];
        C_EZ_HYHZ[w] = 1.0/EPS_EZ[w];
        C_HX_EYEZ[w] = 1.0/MU_0_S;
        C_HY_EZEX[w] = 1.0/MU_0_S;
        C_HZ_EXEY[w] = 1.0/MU_0_S;
        */
//PML用
        sig_ex_x = sig_max*field_sigmaX(x    ,y+0.5,z+0.5);
        sig_ex_y = sig_max*field_sigmaY(x    ,y+0.5,z+0.5);
        sig_ex_z = sig_max*field_sigmaZ(x    ,y+0.5,z+0.5);        
        sig_ey_x = sig_max*field_sigmaX(x+0.5,y    ,z+0.5);
        sig_ey_y = sig_max*field_sigmaY(x+0.5,y    ,z+0.5);
        sig_ey_z = sig_max*field_sigmaZ(x+0.5,y    ,z+0.5);        
        sig_ez_x = sig_max*field_sigmaX(x+0.5,y+0.5,z    );
        sig_ez_y = sig_max*field_sigmaY(x+0.5,y+0.5,z    );
        sig_ez_z = sig_max*field_sigmaZ(x+0.5,y+0.5,z    );

        C_EXY[w]    = field_pmlCoef(EPS_EX[w], sig_ex_y);
        C_EXZ[w]    = field_pmlCoef(EPS_EX[w], sig_ex_z);
        C_EYX[w]    = field_pmlCoef(EPS_EY[w], sig_ey_x);
        C_EYZ[w]    = field_pmlCoef(EPS_EY[w], sig_ey_z);
        C_EZX[w]    = field_pmlCoef(EPS_EZ[w], sig_ez_x);
        C_EZY[w]    = field_pmlCoef(EPS_EZ[w], sig_ez_y);

        C_EXY_HZ[w] = field_pmlCoef2(EPS_EX[w], sig_ex_y);
        C_EXZ_HY[w] = field_pmlCoef2(EPS_EX[w], sig_ex_z);
        C_EYX_HZ[w] = field_pmlCoef2(EPS_EY[w], sig_ey_x);
        C_EYZ_HX[w] = field_pmlCoef2(EPS_EY[w], sig_ey_z);
        C_EZX_HY[w] = field_pmlCoef2(EPS_EZ[w], sig_ez_x);
        C_EZY_HX[w] = field_pmlCoef2(EPS_EZ[w], sig_ez_y);

        //σ*
        sig_hx_x = sig_max*field_sigmaX(x+0.5,y    ,z    ) * MU_0_S/EPSILON_0_S;
        sig_hx_y = sig_max*field_sigmaY(x+0.5,y    ,z    ) * MU_0_S/EPSILON_0_S;
        sig_hx_z = sig_max*field_sigmaZ(x+0.5,y    ,z    ) * MU_0_S/EPSILON_0_S;
        sig_hy_x = sig_max*field_sigmaX(x    ,y+0.5,z    ) * MU_0_S/EPSILON_0_S;
        sig_hy_y = sig_max*field_sigmaY(x    ,y+0.5,z    ) * MU_0_S/EPSILON_0_S;
        sig_hy_z = sig_max*field_sigmaZ(x    ,y+0.5,z    ) * MU_0_S/EPSILON_0_S;
        sig_hz_x = sig_max*field_sigmaX(x    ,y    ,z+0.5) * MU_0_S/EPSILON_0_S;
        sig_hz_y = sig_max*field_sigmaY(x    ,y    ,z+0.5) * MU_0_S/EPSILON_0_S;
        sig_hz_z = sig_max*field_sigmaZ(x    ,y    ,z+0.5) * MU_0_S/EPSILON_0_S;

        C_HXY[w]    = field_pmlCoef(MU_0_S, sig_hx_y);
        C_HXZ[w]    = field_pmlCoef(MU_0_S, sig_hx_z);
        C_HYX[w]    = field_pmlCoef(MU_0_S, sig_hy_x);
        C_HYZ[w]    = field_pmlCoef(MU_0_S, sig_hy_z);
        C_HZX[w]    = field_pmlCoef(MU_0_S, sig_hz_x);
        C_HZY[w]    = field_pmlCoef(MU_0_S, sig_hz_y);

        C_HXY_EZ[w] = field_pmlCoef2(MU_0_S, sig_hx_y);
        C_HXZ_EY[w] = field_pmlCoef2(MU_0_S, sig_hx_z);
        C_HYX_EZ[w] = field_pmlCoef2(MU_0_S, sig_hy_x);
        C_HYZ_EX[w] = field_pmlCoef2(MU_0_S, sig_hy_z);
        C_HZX_EY[w] = field_pmlCoef2(MU_0_S, sig_hz_x);
        C_HZY_EX[w] = field_pmlCoef2(MU_0_S, sig_hz_y);
    }
}

//Private : 
static void calcE()
{
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  for(int i=1; i<subInfo_s.SUB_N_PX-1; i++){
    for(int j=1; j<subInfo_s.SUB_N_PY-1; j++){
      for(int k=1; k<subInfo_s.SUB_N_PZ-1; k++){
        int w = field_subIndex(i,j,k);
        int w_top = field_subTop(w);
        int w_frt = field_subFront(w);
        int w_rht = field_subRight(w);

        Ex[w] = Ex[w] + C_EX_HZHY[w] * ( Hz[w_top] - Hz[w] - Hy[w_frt] + Hy[w]);
        Ey[w] = Ey[w] + C_EY_HXHZ[w] * ( Hx[w_frt] - Hx[w] - Hz[w_rht] + Hz[w]);
        Ez[w] = Ez[w] + C_EZ_HYHZ[w] * ( Hy[w_rht] - Hy[w] - Hx[w_top] + Hx[w]);
      }
    }
  }  
}

static void calcH()
{
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  for(int i=1; i<subInfo_s.SUB_N_PX-1; i++){
    for(int j=1; j<subInfo_s.SUB_N_PY-1; j++){
      for(int k=1; k<subInfo_s.SUB_N_PZ-1; k++){
        int w = field_subIndex(i,j,k);
        int w_bck = field_subBack(w);
        int w_lft = field_subLeft(w);
        int w_btm = field_subBottom(w);

        Hx[w] = Hx[w] + C_HX_EYEZ[w]* ( Ey[w]-Ey[w_bck] -Ez[w]+Ez[w_btm]);
        Hy[w] = Hy[w] + C_HY_EZEX[w]* ( Ez[w]-Ez[w_lft] -Ex[w]+Ex[w_bck]);

        Hz[w] = Hz[w] + C_HZ_EXEY[w]* ( Ex[w]-Ex[w_btm] -Ey[w]*Ey[w_lft]);
      }
    }
  }
}

static void calcE_pml()
{
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  int w;
  int nextX;

  //配列アクセスの関係か,3つに分けた方が13%程度速くなった
  //Ex
  FAST_3FOR(w, subInfo_s, nextX){
    int w_frt = field_subFront(w);
    int w_top = field_subTop(w);
    Exz[w] = C_EXZ[w]*Exz[w] + C_EXZ_HY[w]*( -Hy[w_frt] + Hy[w] );
    Exy[w] = C_EXY[w]*Exy[w] + C_EXY_HZ[w]*(  Hz[w_top] - Hz[w] );
    Ex[w] = Exy[w] + Exz[w];
  }

  //Ey
  FAST_3FOR(w, subInfo_s, nextX){
    int w_rht = field_subRight(w);
    int w_frt = field_subFront(w);
    Eyx[w] = C_EYX[w]*Eyx[w] + C_EYX_HZ[w]*( -Hz[w_rht] + Hz[w] );
    Eyz[w] = C_EYZ[w]*Eyz[w] + C_EYZ_HX[w]*(  Hx[w_frt] - Hx[w] );
    Ey[w] = Eyx[w] + Eyz[w];
  }
  
  //Ez
  FAST_3FOR(w, subInfo_s, nextX){
    int w_top = field_subTop(w);   
    int w_rht = field_subRight(w);
    Ezy[w] = C_EZY[w]*Ezy[w] + C_EZY_HX[w]*( -Hx[w_top] + Hx[w] );
    Ezx[w] = C_EZX[w]*Ezx[w] + C_EZX_HY[w]*(  Hy[w_rht] - Hy[w] );
    Ez[w] = Ezx[w] + Ezy[w];    
  }
}

static void calcH_pml()
{
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  int w;
  int nextX;

  //配列アクセスの関係か, 3つに分けた方が13%程度速くなった
  //Hx
  FAST_3FOR(w, subInfo_s, nextX)
  {
    int w_bck = field_subBack(w);  
    int w_btm = field_subBottom(w);
  
    Hxy[w] = C_HXY[w]*Hxy[w] + C_HXY_EZ[w]*( -Ez[w]+Ez[w_btm] );        
    Hxz[w] = C_HXZ[w]*Hxz[w] + C_HXZ_EY[w]*(  Ey[w]-Ey[w_bck] );
    Hx[w] = Hxy[w] + Hxz[w];
  }

  //Hy
  FAST_3FOR(w, subInfo_s, nextX)
  {
    int w_bck = field_subBack(w);
    int w_lft = field_subLeft(w);  
    Hyz[w] = C_HYZ[w]*Hyz[w] + C_HYZ_EX[w]*( -Ex[w]+Ex[w_bck] );
    Hyx[w] = C_HYX[w]*Hyx[w] + C_HYX_EZ[w]*(  Ez[w]-Ez[w_lft] );
    Hy[w] = Hyx[w] + Hyz[w];
  }

  //Hz
  FAST_3FOR(w, subInfo_s, nextX)
  {
    int w_lft = field_subLeft(w);
    int w_btm = field_subBottom(w);
    Hzx[w] = C_HZX[w]*Hzx[w] + C_HZX_EY[w]*( -Ey[w]*Ey[w_lft] );        
    Hzy[w] = C_HZY[w]*Hzy[w] + C_HZY_EX[w]*(  Ex[w]-Ex[w_btm] );
    Hz[w] = Hzx[w] + Hzy[w];
  }  
}
