#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "field.h"
#include "models.h"
#include "mpiFDTD3D_upml.h"
#include "myComplex.h"
#include "function.h"

/* about MPI  */
static int rank;      //MPIのランク
static int nproc;     //全プロセス数
static MPI_Datatype X_DIRECTION_DOUBLE_COMPLEX;
static MPI_Datatype DCOMPLEX_XY; //XY平面
static MPI_Datatype DCOMPLEX_YZ; //YZ平面
static MPI_Datatype DCOMPLEX_XZ; //XZ平面

//岡田さんの論文と同じ空間配置にしてみる=>遠方解の時の補完が楽になりそう
//h = 1, Δt = 1で計算
//系は右手系
//x(left-, right+)
//y(bottom-, top+)
//z(back-, front+)

//Hx(i,j+0.5,k)         -> Hx[i,j,k]
//Hy(i+0.5,j,k)         -> Hy[i,j,k]
//Hz(i+0.5,j+0.5,k+0.5) -> Hz[i,j,k]

//Ex(i+0.5,j,k+0.5)     -> Ex[i,j,k]
//Ey(i,j+0.5,k+0.5)     -> Ey[i,j,k]
//Ez(i,j,k)             -> Ez[i,j,k]

static dcomplex *Ex = NULL;
static dcomplex *Jx = NULL;
static dcomplex *Dx = NULL;

static dcomplex *Ey = NULL;
static dcomplex *Jy = NULL;
static dcomplex *Dy = NULL;

static dcomplex *Ez = NULL;
static dcomplex *Jz = NULL;
static dcomplex *Dz = NULL;

static dcomplex *Hx = NULL;
static dcomplex *Mx = NULL;
static dcomplex *Bx = NULL;

static dcomplex *Hy = NULL;
static dcomplex *My = NULL;
static dcomplex *By = NULL;

static dcomplex *Hz = NULL;
static dcomplex *Mz = NULL;
static dcomplex *Bz = NULL;

static double *C_JX = NULL, *C_JY = NULL, *C_JZ = NULL;
static double *C_MX = NULL, *C_MY = NULL, *C_MZ= NULL;

static double *C_DX = NULL, *C_DY = NULL, *C_DZ = NULL;
static double *C_BX = NULL, *C_BY = NULL, *C_BZ= NULL;

static double *C_JZHXHY = NULL, *C_JXHYHZ = NULL, *C_JYHXHZ = NULL;
static double *C_MXEYEZ = NULL, *C_MYEXEZ = NULL, *C_MZEXEY = NULL;

static double *C_DXJX0=NULL, *C_DXJX1=NULL;
static double *C_DYJY0=NULL, *C_DYJY1=NULL;
static double *C_DZJZ0=NULL, *C_DZJZ1=NULL;

static double *C_BXMX0=NULL, *C_BXMX1=NULL;
static double *C_BYMY0=NULL, *C_BYMY1=NULL;
static double *C_BZMZ0=NULL, *C_BZMZ1=NULL;

static double *EPS_EX=NULL, *EPS_EY=NULL, *EPS_EZ=NULL;

static void update(void);
static void finish(void);
static void init(void);
static void reset(void);

static void calcJDE(void);
static void calcMBH(void);
static void allocateMemories(void);
static void setCoefficient(void);
static void freeMemories(void);

static void Connection_SendRecvE(void);
static void Connection_SendRecvH(void);
static void scatteredWave(dcomplex *p, double *eps, double gapX, double gapY, double gapZ);

static void miePrint(void);

static void initializeElectroMagneticField(void);

dcomplex* mpi_fdtd3D_upml_getEx(void){  return Ex;}
dcomplex* mpi_fdtd3D_upml_getEy(void){  return Ey;}
dcomplex* mpi_fdtd3D_upml_getEz(void){  return Ez;}
dcomplex* mpi_fdtd3D_upml_getHx(void){  return Hx;}
dcomplex* mpi_fdtd3D_upml_getHy(void){  return Hy;}
dcomplex* mpi_fdtd3D_upml_getHz(void){  return Hz;}

double*   mpi_fdtd3D_upml_getEps(void){  return EPS_EZ;}

void (* mpi_fdtd3D_upml_getUpdate(void))(void){
  return update;
}
void (* mpi_fdtd3D_upml_getFinish(void))(void){
  return finish;
}
void (* mpi_fdtd3D_upml_getReset(void))(void){
  return reset;
}
void (* mpi_fdtd3D_upml_getInit(void))(void){
  return init;
}

static void init(){
  allocateMemories();
  setCoefficient();
}

static void finish(){
  //output();
  miePrint();
  freeMemories();
}

static void reset()
{
}

//Update
static void update(void)
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();
  calcMBH();
  Connection_SendRecvH();
  calcJDE();
  scatteredWave(Ez, EPS_EZ, 0.5, 0.5, 0.0);
  Connection_SendRecvE();  
}

//単一波長の散乱波
// gapX, gapY : Ex-z, Hx-zは格子点からずれた位置に配置され散る為,格子点からのずれを送る必要がある.
// UPML専用
static void scatteredWave(dcomplex *p, double *eps, double gapX, double gapY, double gapZ){
  double ray_coef = field_getRayCoef();  
  double k_s = field_getK();
  double theta_rad = field_getTheta()*M_PI/180.0;
  double phi_rad   = field_getPhi()*M_PI/180.0;

  double ks_cos_cos = cos(phi_rad)*cos(theta_rad)*k_s;
  double ks_cos_sin = cos(phi_rad)*sin(theta_rad)*k_s;
  double ks_sin     = sin(phi_rad)*k_s;
  double w_s_time = field_getOmega() * field_getTime();
  
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  int nextX   = 2*subInfo_s.SUB_N_PZ;
  int endX    = subInfo_s.SUB_N_X;
  int endY    = subInfo_s.SUB_N_Y;
  int endZ    = subInfo_s.SUB_N_Z;
  int offsetX = subInfo_s.OFFSET_X + gapX;
  int offsetY = subInfo_s.OFFSET_Y + gapY;
  int offsetZ = subInfo_s.OFFSET_Z + gapZ;
  int w = field_subIndex(1,1,1);

  double ray_coef_EPS_0 = ray_coef*EPSILON_0_S;
  for(int i=0; i<endX; i++, w+=nextX) 
    for(int j=0; j<endY; j++, w+=2) 
      for(int k=0; k<endZ; k++, w+=1)
      {
        // 空気中は追加の散乱波は0なので無視する.
        if(eps[w] == EPSILON_0_S)
          continue;
        double x = i+offsetX;
        double y = j+offsetY;
        double z = k+offsetZ;

        double kr = x*ks_cos_cos + y*ks_sin + z*ks_cos_sin;
        //p[k] -= かも(岡田さんのメール参照)
        p[w] += (ray_coef_EPS_0/eps[w] - ray_coef)*cexp( I*(kr-w_s_time) );
//        p[w] += ray_coef*(EPSILON_0_S/eps[w] - 1.0)*cexp( I*(kr-w_s_time) );
      }
}

static void Connection_SendRecvE(void)
{
  // (Hの計算に)必要な物は
  // Eのright, top, backなので
  // right, top, back を受け取り
  // left , bottom, frontを送る
  MPI_Status status;
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  
  //左右のランクとの同期
  //this needs only Hy[i-1, j] so send to right and recieve from left
  int rtRecv = field_subIndex(subInfo_s.SUB_N_PX-1, 1, 1);  //最右に格納する
  int ltSend = field_subIndex(1, 1, 1);                     //最左の一つ右を送る
  MPI_Sendrecv(&Ez[ltSend], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.LtRank, 1,
               &Ez[rtRecv], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.RtRank, 1, MPI_COMM_WORLD, &status);
  MPI_Sendrecv(&Ey[ltSend], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.LtRank, 1,
               &Ey[rtRecv], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.RtRank, 1, MPI_COMM_WORLD, &status);

  //上下のランクとの同期
  int bmSend = field_subIndex(1, 1, 1);                      //最下の一つ上を送る
  int tpRecv = field_subIndex(1, subInfo_s.SUB_N_PY-1, 1);  //最上に格納する
  MPI_Sendrecv(&Ez[bmSend], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.BmRank, 1,
               &Ez[tpRecv], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.TpRank, 1, MPI_COMM_WORLD, &status);
  MPI_Sendrecv(&Ex[bmSend], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.BmRank, 1,
               &Ex[tpRecv], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.TpRank, 1, MPI_COMM_WORLD, &status);

  //前後のランクとの同期
  int frSend = field_subIndex(1, 1, 1);                     //最前の一つ上を送る
  int bkRecv = field_subIndex(1, 1, subInfo_s.SUB_N_PZ-1);  //最後に格納する
  MPI_Sendrecv(&Ex[frSend], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.FtRank, 1,
               &Ex[bkRecv], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.BkRank, 1, MPI_COMM_WORLD, &status);
  MPI_Sendrecv(&Ey[frSend], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.FtRank, 1,
               &Ey[bkRecv], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.BkRank, 1, MPI_COMM_WORLD, &status);

}

static void Connection_SendRecvH(void)
{  
  // (Eの計算に)必要な物は
  // Hのleft, bottom, frontなので
  // left, bottom, front を受け取り
  // rihgt, top, backを送る
  
  MPI_Status status;

  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  //左右のランクとの同期
  //this needs only Hy[i-1, j] so send to right and recieve from left 
  int ltRecv = field_subIndex(0                 , 1, 1); //最左の面に格納する為, xのインデックスは0  
  int rtSend = field_subIndex(subInfo_s.SUB_N_PX-2, 1, 1);  //最右-1の値を送るため(最右には何も入っていないから), xのインデックスはSUB_N_PX-2

  //Hz と Hyが左右を使う.
  MPI_Sendrecv(&Hy[rtSend], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.RtRank, 1,
               &Hy[ltRecv], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.LtRank, 1,MPI_COMM_WORLD, &status);
  MPI_Sendrecv(&Hz[rtSend], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.RtRank, 1,
               &Hz[ltRecv], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.LtRank, 1,MPI_COMM_WORLD, &status);

  //上下のランクとの同期
  //this needs only Hy[i,j-1] so send to top and recieve from bottom   
  int bmRecv = field_subIndex(1,0                 , 1);
  int tpSend = field_subIndex(1,subInfo_s.SUB_N_PY-2, 1);
  MPI_Sendrecv(&Hx[tpSend], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.TpRank, 1,
               &Hx[bmRecv], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.BmRank, 1, MPI_COMM_WORLD, &status);
  MPI_Sendrecv(&Hz[tpSend], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.TpRank, 1,
               &Hz[bmRecv], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.BmRank, 1, MPI_COMM_WORLD, &status);

  //前後のランクとの同期
  int ftRecv = field_subIndex(1,1, 0);
  int bkSend = field_subIndex(1,1,subInfo_s.SUB_N_PZ-2);

  MPI_Sendrecv(&Hx[bkSend], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.BkRank, 1,
               &Hx[ftRecv], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.FtRank, 1, MPI_COMM_WORLD, &status);
  MPI_Sendrecv(&Hz[bkSend], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.BkRank, 1,
               &Hz[ftRecv], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.FtRank, 1, MPI_COMM_WORLD, &status);
}

//calculate J and D
static  void calcJDE()
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();

  int w, toNextX;
  //X
  FAST_3FOR(w, sInfo, toNextX)
  {
    const int w_btm = field_subBottom(w);
    const int w_frt = field_subFront(w);

    const dcomplex nowJx = Jx[w];
    Jx[w] = C_JX[w]*Jx[w] + C_JXHYHZ[w]*(Hz[w]-Hz[w_btm] -Hy[w_frt]+Hy[w]);
    Dx[w] = C_DX[w]*Dx[w] + C_DXJX1[w]*Jx[w] - C_DXJX0[w]*nowJx;
    Ex[w] = Dx[w]/EPS_EX[w];
  }

  //Y
  FAST_3FOR(w, sInfo, toNextX)
  {
    const int w_lft = field_subLeft(w);   //一つ左
    const int w_frt = field_subFront(w);  //todo
    const dcomplex nowJy = Jy[w];
    Jy[w] = C_JY[w]*Jy[w] + C_JYHXHZ[w]*( Hx[w_frt]-Hx[w] -Hz[w]+Hz[w_lft] );
    Dy[w] = C_DY[w]*Dy[w] + C_DYJY1[w]*Jy[w] - C_DYJY0[w]*nowJy;
    Ey[w] = Dy[w]/EPS_EY[w];
  }

  //Z
  FAST_3FOR(w, sInfo, toNextX)
  {
    const int w_lft = field_subLeft(w);   //一つ左
    const int w_btm = field_subBottom(w); //一つ下
    const dcomplex nowJz = Jz[w];
    Jz[w] = C_JZ[w]*Jz[w] + C_JZHXHY[w]*(+Hy[w] - Hy[w_lft] -Hx[w]+Hx[w_btm]);
    Dz[w] = C_DZ[w]*Dz[w] + C_DZJZ1[w]*Jz[w] - C_DZJZ0[w]*nowJz;
    Ez[w] = Dz[w]/EPS_EZ[w];
  }
    
}

//calculate M and B
static void calcMBH()
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();
  int w, toNextX;
  
  //X
  FAST_3FOR(w, sInfo, toNextX)
  {
    const int w_top = field_subTop(w);   //一つ上
    const int w_bck = field_subBack(w);  //1つ前        
    const dcomplex nowMx = Mx[w];    
    Mx[w] = C_MX[w]*Mx[w] - C_MXEYEZ[w]*(Ez[w_top]-Ez[w] -Ey[w]+Ey[w_bck]); //原因
    Bx[w] = C_BX[w]*Bx[w] + C_BXMX1[w]*Mx[w] - C_BXMX0[w]*nowMx;
    Hx[w] = Bx[w]/MU_0_S;
  }

  //Y
  FAST_3FOR(w, sInfo, toNextX)
  {
    const int w_rht = field_subRight(w);
    const int w_bck = field_subBack(w);  //1つ前
    const dcomplex nowMy = My[w];
    My[w] = C_MY[w]*My[w] - C_MYEXEZ[w]*( Ex[w]-Ex[w_bck] -Ez[w_rht]+Ez[w]); //原因
    By[w] = C_BY[w]*By[w] + C_BYMY1[w]*My[w] - C_BYMY0[w]*nowMy;
    Hy[w] = By[w]/MU_0_S;
  }

  //Z
  FAST_3FOR(w, sInfo, toNextX)
  {
    const int w_rht = field_subRight(w);
    const int w_top = field_subTop(w);   //一つ上
    const dcomplex nowMz = Mz[w];
    Mz[w] = C_MZ[w]*Mz[w] - C_MZEXEY[w]*( Ey[w_rht]-Ey[w] -Ex[w_top]+Ex[w] );
    Bz[w] = C_BZ[w]*Bz[w] + C_BZMZ1[w]*Mz[w] - C_BZMZ0[w]*nowMz;
    Hz[w] = Bz[w]/MU_0_S;
  }
  
}

//-----------------memory allocate-------------//

static void allocateMemories()
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();
  
  Ex = newDComplex(sInfo.SUB_N_CELL);   //Ex(i+0.5,j    ,k+0.5) -> Ex[i,j,k]
  Ey = newDComplex(sInfo.SUB_N_CELL);   //Ey(i    ,j+0.5,k+0.5) -> Ey[i,j,k]
  Ez = newDComplex(sInfo.SUB_N_CELL);   //Ez(i    ,j    ,k    ) -> Ez[i,j,k]
  
  Jx = newDComplex(sInfo.SUB_N_CELL);
  Jy = newDComplex(sInfo.SUB_N_CELL);
  Jz = newDComplex(sInfo.SUB_N_CELL);
  
  Dx = newDComplex(sInfo.SUB_N_CELL);
  Dy = newDComplex(sInfo.SUB_N_CELL);
  Dz = newDComplex(sInfo.SUB_N_CELL);
  
  Hx = newDComplex(sInfo.SUB_N_CELL);  //Hx(i    , j+0.5, k    )->Hx[i,j,k]
  Hy = newDComplex(sInfo.SUB_N_CELL);  //Hy(i+0.5, j    , k    )->Hy[i,j,k]
  Hz = newDComplex(sInfo.SUB_N_CELL);  //Hz(i+0.5, j+0.5, k+0.5)->Hz[i,j,k]
  
  Mx = newDComplex(sInfo.SUB_N_CELL);
  My = newDComplex(sInfo.SUB_N_CELL);
  Mz = newDComplex(sInfo.SUB_N_CELL);
  
  Bx = newDComplex(sInfo.SUB_N_CELL);
  By = newDComplex(sInfo.SUB_N_CELL);
  Bz = newDComplex(sInfo.SUB_N_CELL);

  C_JX = newDouble(sInfo.SUB_N_CELL);
  C_JY = newDouble(sInfo.SUB_N_CELL);
  C_JZ = newDouble(sInfo.SUB_N_CELL);
  
  C_MX = newDouble(sInfo.SUB_N_CELL);
  C_MY = newDouble(sInfo.SUB_N_CELL);
  C_MZ = newDouble(sInfo.SUB_N_CELL);

  C_DX = newDouble(sInfo.SUB_N_CELL);
  C_DY = newDouble(sInfo.SUB_N_CELL);
  C_DZ = newDouble(sInfo.SUB_N_CELL);
  
  C_BX = newDouble(sInfo.SUB_N_CELL);
  C_BY = newDouble(sInfo.SUB_N_CELL);
  C_BZ = newDouble(sInfo.SUB_N_CELL);

  C_JXHYHZ = newDouble(sInfo.SUB_N_CELL);
  C_JYHXHZ = newDouble(sInfo.SUB_N_CELL);
  C_JZHXHY = newDouble(sInfo.SUB_N_CELL);
  
  C_MXEYEZ = newDouble(sInfo.SUB_N_CELL);
  C_MYEXEZ = newDouble(sInfo.SUB_N_CELL);
  C_MZEXEY = newDouble(sInfo.SUB_N_CELL);
  
  C_DXJX0 = newDouble(sInfo.SUB_N_CELL);
  C_DXJX1 = newDouble(sInfo.SUB_N_CELL);
  C_DYJY0 = newDouble(sInfo.SUB_N_CELL);
  C_DYJY1 = newDouble(sInfo.SUB_N_CELL);
  C_DZJZ0 = newDouble(sInfo.SUB_N_CELL);
  C_DZJZ1 = newDouble(sInfo.SUB_N_CELL);

  C_BXMX1 = newDouble(sInfo.SUB_N_CELL);
  C_BXMX0 = newDouble(sInfo.SUB_N_CELL);
  C_BYMY1 = newDouble(sInfo.SUB_N_CELL);
  C_BYMY0 = newDouble(sInfo.SUB_N_CELL);
  C_BZMZ0 = newDouble(sInfo.SUB_N_CELL);
  C_BZMZ1 = newDouble(sInfo.SUB_N_CELL);

  EPS_EX = newDouble(sInfo.SUB_N_CELL);
  EPS_EY = newDouble(sInfo.SUB_N_CELL);
  EPS_EZ = newDouble(sInfo.SUB_N_CELL);
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

  memset(Jx, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
  memset(Jy, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
  memset(Jz, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
  
  memset(Mx, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
  memset(My, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
  memset(Mz, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);

  memset(Dz, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
  memset(Dx, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
  memset(Dy, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);

  memset(Bx, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
  memset(By, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);
  memset(Bz, 0, sizeof(dcomplex)*sInfo.SUB_N_CELL);    
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
        int x = i-1+sInfo.OFFSET_X;
        int y = j-1+sInfo.OFFSET_Y;
        int z = k-1+sInfo.OFFSET_Z;
        EPS_EX[w] = models_eps(x+0.5,y    ,z+0.5,D_Y); //todo 
        EPS_EY[w] = models_eps(x    ,y+0.5,z+0.5,D_X);
        EPS_EZ[w] = models_eps(x    ,y    ,z    ,D_XY);

        sig_ex_x = sig_max*field_sigmaX(x+0.5,y    ,z+0.5);
        sig_ex_y = sig_max*field_sigmaY(x+0.5,y    ,z+0.5);
        sig_ex_z = sig_max*field_sigmaZ(x+0.5,y    ,z+0.5);
        
        sig_ey_x = sig_max*field_sigmaX(x    ,y+0.5,z+0.5);
        sig_ey_y = sig_max*field_sigmaY(x    ,y+0.5,z+0.5);
        sig_ey_z = sig_max*field_sigmaZ(x    ,y+0.5,z+0.5);
        
        sig_ez_x = sig_max*field_sigmaX(x    ,y    ,z    );
        sig_ez_y = sig_max*field_sigmaY(x    ,y    ,z    );
        sig_ez_z = sig_max*field_sigmaZ(x    ,y    ,z    );
        
        sig_hx_x = sig_max*field_sigmaX(x    ,y+0.5,z    );
        sig_hx_y = sig_max*field_sigmaY(x    ,y+0.5,z    );
        sig_hx_z = sig_max*field_sigmaZ(x    ,y+0.5,z    );
        
        sig_hy_x = sig_max*field_sigmaX(x+0.5,y    ,z    );
        sig_hy_y = sig_max*field_sigmaY(x+0.5,y    ,z    );
        sig_hy_z = sig_max*field_sigmaZ(x+0.5,y    ,z    );
        
        sig_hz_x = sig_max*field_sigmaX(x+0.5,y+0.5,z+0.5);
        sig_hz_y = sig_max*field_sigmaY(x+0.5,y+0.5,z+0.5);
        sig_hz_z = sig_max*field_sigmaZ(x+0.5,y+0.5,z+0.5);

        //Δt = 1 , Κ_i = 1, h = 1
        double eps = EPSILON_0_S;        
        C_JX[w]    = (2*eps - sig_ex_y) / (2*eps + sig_ex_y);
        C_JXHYHZ[w]= (2*eps)            / (2*eps + sig_ex_y);
        C_DX[w]    = (2*eps - sig_ex_z) / (2*eps + sig_ex_z);
        C_DXJX1[w] = (2*eps + sig_ex_x) / (2*eps + sig_ex_z);
        C_DXJX0[w] = (2*eps - sig_ex_x) / (2*eps + sig_ex_z);

        C_JY[w]    = (2*eps - sig_ey_z) / (2*eps + sig_ey_z);        
        C_JYHXHZ[w]= (2*eps)            / (2*eps + sig_ey_z);
        C_DY[w]    = (2*eps - sig_ey_x) / (2*eps + sig_ey_x);
        C_DYJY1[w] = (2*eps + sig_ey_y) / (2*eps + sig_ey_x);
        C_DYJY0[w] = (2*eps - sig_ey_y) / (2*eps + sig_ey_x);

        C_JZ[w]    = (2*eps - sig_ez_x) / (2*eps + sig_ez_x);
        C_JZHXHY[w]= (2*eps )           / (2*eps + sig_ez_x);
        C_DZ[w]    = (2*eps - sig_ez_y) / (2*eps + sig_ez_y);      
        C_DZJZ1[w] = (2*eps + sig_ez_z) / (2*eps + sig_ez_y);
        C_DZJZ0[w] = (2*eps - sig_ez_z) / (2*eps + sig_ez_y);
        
        C_MX[w]    = (2*eps - sig_hx_y) / (2*eps + sig_hx_y);
        C_MXEYEZ[w]= (2*eps)            / (2*eps + sig_hx_y);
        C_BX[w]    = (2*eps - sig_hx_z) / (2*eps + sig_hx_z);
        C_BXMX1[w] = (2*eps + sig_hx_x) / (2*eps + sig_hx_z);
        C_BXMX0[w] = (2*eps - sig_hx_x) / (2*eps + sig_hx_z);

        C_MY[w]    = (2*eps - sig_hy_z) / (2*eps + sig_hy_z);
        C_MYEXEZ[w]= (2*eps)            / (2*eps + sig_hy_z);      
        C_BY[w]    = (2*eps - sig_hy_x) / (2*eps + sig_hy_x);
        C_BYMY1[w] = (2*eps + sig_hy_y) / (2*eps + sig_hy_x);
        C_BYMY0[w] = (2*eps - sig_hy_y) / (2*eps + sig_hy_x);
        
        C_MZ[w]    = (2*eps - sig_hz_x) / (2*eps + sig_hz_x);
        C_MZEXEY[w]= (2*eps)            / (2*eps + sig_hz_x);
        C_BZ[w]    = (2*eps - sig_hz_y) / (2*eps + sig_hz_y);
        C_BZMZ1[w] = (2*eps + sig_hz_z) / (2*eps + sig_hz_y);
        C_BZMZ0[w] = (2*eps - sig_hz_z) / (2*eps + sig_hz_y);
    }

}

static void cpy(dcomplex *entire, dcomplex *region, int dx, int dy, int dz)
{
  //sub領域のregionをentireにコピー
  //SUB_N_PXとかは, 全プロセスで共通(なはず)なので, subInfo_sの値をそのまま使う
  //オフセットはプロセスごとに違うので外部から与える.
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();  
  for(int i=1, x=dx; i<subInfo_s.SUB_N_PX-1; i++, x++)
    for(int j=1, y=dy; j<subInfo_s.SUB_N_PY-1; j++, y++)
      for(int k=1, z=dz; k<subInfo_s.SUB_N_PZ-1; k++, z++)
      {
        entire[field_index(x,y,z)] = region[field_subIndex(i,j,k)];
      }
}

//分割された領域をまとめる.
static dcomplex* unifyToRank0(dcomplex *phi)
{
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  //マスターにすべて集める
  if(subInfo_s.Rank == 0)
  {
    MPI_Status status;
    dcomplex *entire = newDComplex(fInfo_s.N_CELL);
    cpy(entire, Ez, subInfo_s.OFFSET_X, subInfo_s.OFFSET_Y, subInfo_s.OFFSET_Z);

    dcomplex *tmp = newDComplex(subInfo_s.SUB_N_CELL);
    int offset[3];
    for(int i=1; i<subInfo_s.Nproc; i++)
    {
      MPI_Recv(offset, 3, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(tmp, subInfo_s.SUB_N_CELL, MPI_C_DOUBLE_COMPLEX, i, 1, MPI_COMM_WORLD, &status);
      cpy(entire, tmp, offset[0], offset[1], offset[2]);      
    }
    free(tmp);
    return entire;
  } else {
    int offset[3];
    offset[0] = subInfo_s.OFFSET_X;
    offset[1] = subInfo_s.OFFSET_Y;
    offset[2] = subInfo_s.OFFSET_Z;
    MPI_Send(offset, 3, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(Ez, subInfo_s.SUB_N_CELL, MPI_C_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);
    
    return NULL; //マスター以外はNULLを返す.
  }
}

//---------------------メモリの解放--------------------//
static void miePrint()
{
  dcomplex *entire = unifyToRank0(Ez);
  if( entire != NULL){
    field_outputElliptic("Ez.txt",entire);
    free(entire);
  }
  /*
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  //マスターにすべて集める
  if(subInfo_s.Rank == 0)
  {
    MPI_Status status;
    dcomplex *entire = newDComplex(fInfo_s.N_CELL);

    cpy(entire, Ez, subInfo_s.OFFSET_X, subInfo_s.OFFSET_Y, subInfo_s.OFFSET_Z);

    int offset[3];
    for(int i=1; i<subInfo_s.Nproc; i++)
    {
      MPI_Recv(offset, 3, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(Ez, subInfo_s.SUB_N_CELL, MPI_C_DOUBLE_COMPLEX, i, 1, MPI_COMM_WORLD, &status);
      cpy(entire, Ez, offset[0], offset[1], offset[2]);      
    }
    field_outputElliptic("Ez.txt",entire);
    free(entire);
  } else {
    int offset[3];
    offset[0] = subInfo_s.OFFSET_X;
    offset[1] = subInfo_s.OFFSET_Y;
    offset[2] = subInfo_s.OFFSET_Z;
    MPI_Send(offset, 3, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(Ez, subInfo_s.SUB_N_CELL, MPI_C_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);
    }  */
}

static void freeMemories()
{
  if(Ex != NULL){   free(Ex); Ex = NULL;}  
  if(Ey != NULL){   free(Ey); Ey = NULL;}
  if(Ez != NULL){   free(Ez); Ez = NULL;}
  
  if(Jx != NULL){   free(Jx); Jx = NULL;}  
  if(Jy != NULL){   free(Jy); Jy = NULL;}
  if(Jz != NULL){   free(Jz); Jz = NULL;}

  if(Dx != NULL){   free(Dx); Dx = NULL;}  
  if(Dy != NULL){   free(Dy); Dy = NULL;}
  if(Dz != NULL){   free(Dz); Dz = NULL;}

  if(Hx != NULL){   free(Hx); Hx = NULL;}
  if(Hy != NULL){   free(Hy); Hy = NULL;}
  if(Hz != NULL){   free(Hz); Hz = NULL;}

  if(Mx != NULL){   free(Mx); Mx = NULL;}
  if(My != NULL){   free(My); My = NULL;}
  if(Mz != NULL){   free(Mz); Mz = NULL;}

  if(Bx != NULL){   free(Bx); Bx = NULL;}
  if(By != NULL){   free(By); By = NULL;}
  if(Bz != NULL){   free(Bz); Bz = NULL;}

  if(C_JX!= NULL){    free(C_JX);  C_JX = NULL;}
  if(C_JY!= NULL){    free(C_JY);  C_JY = NULL;}
  if(C_JZ!= NULL){    free(C_JZ);  C_JZ = NULL;}  

  if(C_MX!= NULL){    free(C_MX);  C_MX = NULL;}
  if(C_MY!= NULL){    free(C_MY);  C_MY = NULL;}
  if(C_MZ!= NULL){    free(C_MZ);  C_MZ = NULL;}  

  if(C_DX!= NULL){    free(C_DX);  C_DX = NULL;}
  if(C_DY!= NULL){    free(C_DY);  C_DY = NULL;}
  if(C_DZ!= NULL){    free(C_DZ);  C_DZ = NULL;}  

  if(C_BX!= NULL){    free(C_BX);  C_BX = NULL;}
  if(C_BY!= NULL){    free(C_BY);  C_BY = NULL;}
  if(C_BZ!= NULL){    free(C_BZ);  C_BZ = NULL;}  

  
  if(C_JXHYHZ!= NULL){ free(C_JXHYHZ); C_JXHYHZ = NULL;}
  if(C_JYHXHZ!= NULL){ free(C_JYHXHZ); C_JYHXHZ = NULL;}
  if(C_JZHXHY!= NULL){ free(C_JZHXHY); C_JZHXHY = NULL;}

  if(C_MXEYEZ!= NULL){ free(C_MXEYEZ); C_MXEYEZ = NULL;}
  if(C_MYEXEZ!= NULL){ free(C_MYEXEZ); C_MYEXEZ = NULL;}
  if(C_MZEXEY!= NULL){ free(C_MZEXEY); C_MZEXEY = NULL;}
  

  if(C_DXJX0 != NULL){ free(C_DXJX0); C_DXJX0 = NULL;}
  if(C_DXJX1 != NULL){ free(C_DXJX1); C_DXJX1 = NULL;}
  if(C_DYJY0 != NULL){ free(C_DYJY0); C_DYJY0 = NULL;}
  if(C_DYJY1 != NULL){ free(C_DYJY1); C_DYJY1 = NULL;}
  if(C_DZJZ1 != NULL){ free(C_DZJZ1); C_DZJZ1 = NULL;}
  if(C_DZJZ0 != NULL){ free(C_DZJZ0); C_DZJZ0 = NULL;}

  
  if(C_BXMX1 != NULL){ free(C_BXMX1); C_BXMX1 = NULL;}
  if(C_BXMX0 != NULL){ free(C_BXMX0); C_BXMX0 = NULL;}
  if(C_BYMY1 != NULL){ free(C_BYMY1); C_BYMY1 = NULL;}
  if(C_BYMY0 != NULL){ free(C_BYMY0); C_BYMY0 = NULL;}
  if(C_BZMZ0 != NULL){ free(C_BZMZ0); C_BZMZ0 = NULL;}
  if(C_BZMZ1 != NULL){ free(C_BZMZ1); C_BZMZ1 = NULL;}
  
  if(EPS_EX != NULL)   free(EPS_EX);
  if(EPS_EY != NULL)   free(EPS_EY);
  if(EPS_EZ != NULL)   free(EPS_EZ);
}




//============================== Debug ==============================//
static bool debugCheck(dcomplex *p)
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();
  for(int i=1; i<sInfo.SUB_N_PX-1; i++)
    for(int j=1; j<sInfo.SUB_N_PY-1; j++)
      for(int k=1; k<sInfo.SUB_N_PZ-1; k++)
      {
        int w = field_subIndex(i,j,k);
        if( fabs( creal(p[w]) ) > 100 )
        {
          printf("%lf, %lf, %d, %d, %d\n", creal(p[w]), cimag(p[w]), i, j, k);
          return true;
        }
      }
  return false;
}

static void debugPrint()
{  
  int w = field_subIndex(31, 33, 26);
  const int w_lft = field_subLeft(w);   //一つ左
  const int w_btm = field_subBottom(w); //一つ下
  const int w_frt = field_subFront(w);  //todo
  const int w_bck = field_subBack(w);  //todo
  const int w_rht = field_subRight(w);
  const int w_top = field_subTop(w);   //一つ上

  dcomplex dJx1 = (+Hz[w] - Hz[w_btm]);
  dcomplex dJx2 = -Hy[w_frt] + Hy[w];
  dcomplex dJy1 = (+Hx[w_frt] - Hx[w]);
  dcomplex dJy2 = -Hz[w] + Hz[w_lft];
  dcomplex dJz1 = (+Hy[w] - Hy[w_lft]);
  dcomplex dJz2 = -Hx[w] + Hx[w_btm];
  
  dcomplex dMx1 = (Ez[w_top] - Ez[w]);
  dcomplex dMx2 = -Ey[w_frt]+Ey[w];
  dcomplex dMy1 = (Ex[w_frt] - Ex[w]);
  dcomplex dMy2 = -Ez[w_rht]+Ez[w];
  dcomplex dMz1 = (Ey[w_rht] - Ey[w]);
  dcomplex dMz2 = -Ex[w_top]+Ex[w];
  
  printf("dJ1( %lf , %lf,  %lf) , dJ2(%lf , %lf,  %lf ) dM1( %lf, %lf, %lf)  dM2(%lf, %lf, %lf) \n ",
         creal(dJx1), creal(dJy1), creal(dJz1), creal(dJx2), creal(dJy2), creal(dJz2), creal(dMx1), creal(dMy1), creal(dMz1), creal(dMx2), creal(dMy2), creal(dMz2) );

}