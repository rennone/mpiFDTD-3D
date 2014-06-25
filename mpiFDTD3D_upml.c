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
#include "ntff3D.h"
//岡田さんの論文と同じ空間配置にしてみる=>遠方解の時の補完が楽になりそう
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

#ifdef DEBUG
static void debugOutput(void);
#endif

static void Connection_SendRecvE(void);
static void Connection_SendRecvH(void);
static void scatteredWave(dcomplex *p, double *eps, double gapX, double gapY, double gapZ);
static void scatteredPulse(dcomplex *p, double *eps, double gapX, double gapY, double gapZ);
static void pointLightInCenter(dcomplex *p);

static void output(void);

static void initializeElectroMagneticField(void);

dcomplex* mpi_fdtd3D_upml_getEx(void){  return Ex;}
dcomplex* mpi_fdtd3D_upml_getEy(void){  return Ey;}
dcomplex* mpi_fdtd3D_upml_getEz(void){  return Ez;}
dcomplex* mpi_fdtd3D_upml_getHx(void){  return Hx;}
dcomplex* mpi_fdtd3D_upml_getHy(void){  return Hy;}
dcomplex* mpi_fdtd3D_upml_getHz(void){  return Hz;}

double*   mpi_fdtd3D_upml_getEps(void) {  return EPS_EY;}
dcomplex* mpi_fdtd3D_upml_getData(void){  return Ey;}

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

static void readData(dcomplex *ex,dcomplex *ey,dcomplex *ez,dcomplex *hx,dcomplex *hy,dcomplex *hz)
{
  FILE * fpEx= fopen("mieScattering_Ey_h_5nm/Ex.txt", "r");
  FILE * fpEy= fopen("mieScattering_Ey_h_5nm/Ey.txt", "r");
  FILE * fpEz= fopen("mieScattering_Ey_h_5nm/Ez.txt", "r");
  FILE * fpHx= fopen("mieScattering_Ey_h_5nm/Hx.txt", "r");
  FILE * fpHy= fopen("mieScattering_Ey_h_5nm/Hy.txt", "r");
  FILE * fpHz= fopen("mieScattering_Ey_h_5nm/Hz.txt", "r");

  if( (fpEx==NULL) || (fpEy==NULL) || (fpEz==NULL) ||
      (fpHx==NULL) || (fpHy==NULL) || (fpHz==NULL) )
  {
    printf("cannot read data files \n");
    exit(2);
  }
  
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  double re, im;
  for(int i=0; i<fInfo_s.N_CELL; i++)
  {  
    fscanf(fpEx, "%lf, %lf\n", &re, &im);
    ex[i] = re + I*im;
  }

  for(int i=0; i<fInfo_s.N_CELL; i++)
  {    
    fscanf(fpEy, "%lf, %lf\n", &re, &im);
    ey[i] = re + I*im;
  }
  
  for(int i=0; i<fInfo_s.N_CELL; i++)
  {    
    fscanf(fpEz, "%lf, %lf\n", &re, &im);
    ez[i] = re + I*im;
  }
  
  for(int i=0; i<fInfo_s.N_CELL; i++)
  {    
    fscanf(fpHx, "%lf, %lf\n", &re, &im);
    hx[i] = re + I*im;
  }
  
  for(int i=0; i<fInfo_s.N_CELL; i++)
  {    
    fscanf(fpHy, "%lf, %lf\n", &re, &im);
    hy[i] = re + I*im;
  }
  
  for(int i=0; i<fInfo_s.N_CELL; i++)
  {    
    fscanf(fpHz, "%lf, %lf\n", &re, &im);
    hz[i] = re + I*im;
  }
}

static void readDataAndFinish()
{
  printf("read Data at Init \n");
  printf("value of electro magnetic field are set by .txt file, \n");
  printf("but time step is zero yet, so you should set the time by maxStep if you use a time \n");

  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  dcomplex *ex = newDComplex(fInfo_s.N_CELL);
  dcomplex *ey = newDComplex(fInfo_s.N_CELL);
  dcomplex *ez = newDComplex(fInfo_s.N_CELL);
  dcomplex *hx = newDComplex(fInfo_s.N_CELL);
  dcomplex *hy = newDComplex(fInfo_s.N_CELL);
  dcomplex *hz = newDComplex(fInfo_s.N_CELL);
  readData(ex,ey,ez,hx,hy,hz);

  ntff3D_Frequency(ex, ey, ez, hx, hy, hz);
  exit(0);
}

static void init(){
  allocateMemories();
  setCoefficient();

#ifdef USE_FILE_DATA
//  readDataAndFinish();
#endif
  ntff3D_Init();
}

static void finish(){
  output();
//  ntff3D_TimeOutput();  
  freeMemories();
}

static void reset()
{
}

//Update
static void update(void)
{

  calcMBH();
  Connection_SendRecvH();
//  MPI_Barrier(MPI_COMM_WORLD);
  calcJDE();

  scatteredWave(Ey, EPS_EY, 0.5, 0.0, 0.5);
//  scatteredPulse(Ey, EPS_EY, 0.5, 0.0, 0.5);
//  MPI_Barrier(MPI_COMM_WORLD);
  Connection_SendRecvE();

//  ntff3D_SubTimeCalc(Ex, Ey, Ez, Hx, Hy, Hz);
}

static void pointLightInCenter(dcomplex *p)
{
  FieldInfo_S fInfo = field_getFieldInfo_S();
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();
  bool XX = (sInfo.OFFSET_X < fInfo.N_PX/2) && ( fInfo.N_PX/2 <= sInfo.OFFSET_X + sInfo.SUB_N_X);
  bool YY = (sInfo.OFFSET_Y < fInfo.N_PY/2) && ( fInfo.N_PY/2 <= sInfo.OFFSET_Y + sInfo.SUB_N_Y);
  bool ZZ = (sInfo.OFFSET_Z < fInfo.N_PZ/2) && ( fInfo.N_PZ/2 <= sInfo.OFFSET_Z + sInfo.SUB_N_Z);

  if(XX && YY && ZZ){
    int w = field_subIndex(fInfo.N_PX/2 - sInfo.OFFSET_X, fInfo.N_PY/2 - sInfo.OFFSET_Y, fInfo.N_PZ/2 - sInfo.OFFSET_Z);
    p[w] += field_pointLight();
  }
}

//単一波長の散乱波
// gapX, gapY : Ex-z, Hx-zは格子点からずれた位置に配置され散る為,格子点からのずれを送る必要がある.
// UPML専用
static void scatteredWave(dcomplex *p, double *eps, double gapX, double gapY, double gapZ){
  double ray_coef = field_getRayCoef();  
  double k_s = field_getK();
  double theta_rad = field_getTheta()*M_PI/180.0;
  double phi_rad   = field_getPhi()*M_PI/180.0;

  double ks_sin_cos = sin(theta_rad)*cos(phi_rad)*k_s;
  double ks_sin_sin = sin(theta_rad)*sin(phi_rad)*k_s;
  double ks_cos     = cos(theta_rad)*k_s;
  double w_s_time = field_getOmega() * field_getTime();
  
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  int nextX   = 2*subInfo_s.SUB_N_PZ; //最後ののりしろの分１行多くずれる
  int endX    = subInfo_s.SUB_N_PX-1;
  int endY    = subInfo_s.SUB_N_PY-1;
  int endZ    = subInfo_s.SUB_N_PZ-1;
  int offsetX = subInfo_s.OFFSET_X + gapX;
  int offsetY = subInfo_s.OFFSET_Y + gapY;
  int offsetZ = subInfo_s.OFFSET_Z + gapZ;
  int w = field_subIndex(1,1,1);
  double ray_coef_EPS_0 = ray_coef*EPSILON_0_S;
  for(int i=1; i<endX; i++, w+=nextX) 
    for(int j=1; j<endY; j++, w+=2) 
      for(int k=1; k<endZ; k++, w+=1)
      {
        // 空気中は追加の散乱波は0なので無視する.
        if(eps[w] == EPSILON_0_S)
          continue;
        double x = i+offsetX;
        double y = j+offsetY;
        double z = k+offsetZ;

        double kr = x*ks_sin_cos + y*ks_sin_cos + z*ks_cos;
        //p[k] -= かも(岡田さんのメール参照)
        p[w] += (ray_coef_EPS_0/eps[w] - ray_coef)*cexp( I*(kr-w_s_time) );
//        p[w] += ray_coef*(EPSILON_0_S/eps[w] - 1.0)*cexp( I*(kr-w_s_time) );
      }
}

//パルス波の散乱波
// gapX, gapY : Ex-z, Hx-zは格子点からずれた位置に配置され散る為,格子点からのずれを送る必要がある.
// UPML専用
static void scatteredPulse(dcomplex *p, double *eps, double gapX, double gapY, double gapZ)
{
  double theta_rad = field_getTheta()*M_PI/180.0;
  double phi_rad   = field_getPhi()  *M_PI/180.0;
  double sin_cos_per_c = sin(theta_rad)*cos(phi_rad)/C_0_S;
  double sin_sin_per_c = sin(theta_rad)*sin(phi_rad)/C_0_S;
  double cos_per_c     = cos(theta_rad)/C_0_S;
  
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  int nextX   = 2*subInfo_s.SUB_N_PZ; //最後ののりしろの分１行多くずれる
  int endX    = subInfo_s.SUB_N_PX-1;
  int endY    = subInfo_s.SUB_N_PY-1;
  int endZ    = subInfo_s.SUB_N_PZ-1;
  int offsetX = subInfo_s.OFFSET_X + gapX;
  int offsetY = subInfo_s.OFFSET_Y + gapY;
  int offsetZ = subInfo_s.OFFSET_Z + gapZ;
  
  double w_s  = field_getOmega();
  const double beam_width = 50; //パルスの幅
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  
  //waveAngleにより, t0の値を変えないとちょうどいいところにピークが来なため,それを計算.
  const double center_peak = (fInfo_s.N_PX/2.0+gapX)*sin_cos_per_c + (fInfo_s.N_PY/2+gapY)*sin_sin_per_c + (fInfo_s.N_PZ/2+gapZ)*cos_per_c; //中心にピークがくる時間
  const double t_minus_t0 = field_getTime()-center_peak + 100; // t-t0. 常に100ステップの時に,領域の中心にピークが来るようにする.
  int w = field_subIndex(1,1,1);
  for(int i=1; i<endX; i++, w+=nextX) 
    for(int j=1; j<endY; j++, w+=2) 
      for(int k=1; k<endZ; k++, w+=1)
      {
        if(eps[w] == EPSILON_0_S)        // 空気中は追加の散乱波は0なので無視する.
          continue;
        double x = i+offsetX; double y = j+offsetY;  double z = k+offsetZ;

        const double r = x*sin_cos_per_c + y*sin_sin_per_c + z*cos_per_c - (t_minus_t0); // (x*cos+y*sin)/C - (time-t0)
        const double gaussian_coef = exp( -pow(r/beam_width, 2 ) );
        p[w] += gaussian_coef*(EPSILON_0_S/eps[w] - 1)*cexp(I*r*w_s);     //p[k] -= かも(岡田さんのメール参照)
      }
}

static void Connection_SendRecvE(void)
{
  // (Hの計算に)必要な物は
  // Ey,Ezのleft. Ex,Ezのbottm. Ex,Eyのbackなので
  // left, bottom, back を受け取り
  // right , top, frontを送る
  MPI_Status status;
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  
  //左右のランクとの同期
  int ltRecv = field_subIndex(0, 1, 1);                    //最左に格納する
  int rtSend = field_subIndex(subInfo_s.SUB_N_PX-2, 1, 1); //最右の一つ左を送る
  MPI_Sendrecv(&Ez[rtSend], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.RtRank, 1,
               &Ez[ltRecv], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.LtRank, 1, MPI_COMM_WORLD, &status);
  MPI_Sendrecv(&Ey[rtSend], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.RtRank, 1,
               &Ey[ltRecv], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.LtRank, 1, MPI_COMM_WORLD, &status);

  //上下のランクとの同期
  int tpSend = field_subIndex(1, subInfo_s.SUB_N_PY-2, 1); //最上の一つ下を送る
  int bmRecv = field_subIndex(1, 0, 1);                    //最下に格納する
  MPI_Sendrecv(&Ez[tpSend], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.TpRank, 1,
               &Ez[bmRecv], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.BmRank, 1, MPI_COMM_WORLD, &status);
  MPI_Sendrecv(&Ex[tpSend], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.TpRank, 1,
               &Ex[bmRecv], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.BmRank, 1, MPI_COMM_WORLD, &status);

  //前後のランクとの同期
  int frSend = field_subIndex(0, 0, subInfo_s.SUB_N_PZ-2);//field_subIndex(1, 1, subInfo_s.SUB_N_PZ-2);  //最前面の一つ手前を送る(右手系なので手前が大きい)
  int bkRecv = field_subIndex(0, 0, 0);//field_subIndex(1, 1, 0);                     //最背面に格納する
  MPI_Sendrecv(&Ex[frSend], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.FtRank, 1,
               &Ex[bkRecv], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.BkRank, 1, MPI_COMM_WORLD, &status);
  MPI_Sendrecv(&Ey[frSend], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.FtRank, 1,
               &Ey[bkRecv], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.BkRank, 1, MPI_COMM_WORLD, &status);

}

static void Connection_SendRecvH(void)
{  
  // (Eの計算に)必要な物は
  // Hy,Hzのright. Hx,Hzのtop. Hx,Hyのfrontなので
  // right, top, front を受け取り
  // left, bottom, backを送る
  
  MPI_Status status;
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  
  //左右のランクとの同期
  int rtRecv = field_subIndex(subInfo_s.SUB_N_PX-1 , 1, 1); //最右の面に格納する為, xのインデックスはN_PX-1
  int ltSend = field_subIndex(1, 1, 1);           //最右+1の値を送るため(最右には何も入っていないから), xのインデックスはSUB_N_PX-2

  //Hz と Hyが左右を使う.
  MPI_Sendrecv(&Hy[ltSend], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.LtRank, 1,
               &Hy[rtRecv], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.RtRank, 1,MPI_COMM_WORLD, &status);
  MPI_Sendrecv(&Hz[ltSend], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.LtRank, 1,
               &Hz[rtRecv], 1, MPI_DCOMPLEX_YZ_PLANE, subInfo_s.RtRank, 1,MPI_COMM_WORLD, &status);

  //上下のランクとの同期
  //this needs only Hy[i,j-1] so send to top and recieve from bottom   
  int tpRecv = field_subIndex(1,subInfo_s.SUB_N_PY-1, 1);
  int bmSend = field_subIndex(1, 1, 1);
  MPI_Sendrecv(&Hx[bmSend], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.BmRank, 1,
               &Hx[tpRecv], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.TpRank, 1, MPI_COMM_WORLD, &status);
  MPI_Sendrecv(&Hz[bmSend], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.BmRank, 1,
               &Hz[tpRecv], 1, MPI_DCOMPLEX_XZ_PLANE, subInfo_s.TpRank, 1, MPI_COMM_WORLD, &status);

  //前後のランクとの同期=>これは, 配列が2段階に隙間があるから.のりしろを含めた平面全部を同期しないと行けない
  int ftRecv = field_subIndex(0, 0, subInfo_s.SUB_N_PZ-1);//field_subIndex(1, 1, subInfo_s.SUB_N_PZ-1);
  int bkSend = field_subIndex(0, 0, 1);

  MPI_Sendrecv(&Hx[bkSend], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.BkRank, 1,
               &Hx[ftRecv], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.FtRank, 1, MPI_COMM_WORLD, &status);
  MPI_Sendrecv(&Hy[bkSend], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.BkRank, 1,
               &Hy[ftRecv], 1, MPI_DCOMPLEX_XY_PLANE, subInfo_s.FtRank, 1, MPI_COMM_WORLD, &status);
}

//calculate J and D
static  void calcJDE()
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();

  int w, toNextX;
  //X
  FAST_3FOR(w, sInfo, toNextX)
  {
    const int w_top = field_subTop(w);
    const int w_frt = field_subFront(w);

    const dcomplex nowJx = Jx[w];
    Jx[w] = C_JX[w]*Jx[w] + C_JXHYHZ[w]*(Hz[w_top]-Hz[w] -Hy[w_frt]+Hy[w]);
    Dx[w] = C_DX[w]*Dx[w] + C_DXJX1[w]*Jx[w] - C_DXJX0[w]*nowJx;
    Ex[w] = Dx[w]/EPS_EX[w];
  }

  //Y
  FAST_3FOR(w, sInfo, toNextX)
  {
    const int w_rht = field_subRight(w);   //一つ左
    const int w_frt = field_subFront(w);  //todo
    const dcomplex nowJy = Jy[w];
    Jy[w] = C_JY[w]*Jy[w] + C_JYHXHZ[w]*( Hx[w_frt]-Hx[w] -Hz[w_rht]+Hz[w] );
    Dy[w] = C_DY[w]*Dy[w] + C_DYJY1[w]*Jy[w] - C_DYJY0[w]*nowJy;
    Ey[w] = Dy[w]/EPS_EY[w];
  }

  //Z
  FAST_3FOR(w, sInfo, toNextX)
  {
    const int w_rht = field_subRight(w);   //一つ左
    const int w_top = field_subTop(w); //一つ下
    const dcomplex nowJz = Jz[w];
    Jz[w] = C_JZ[w]*Jz[w] + C_JZHXHY[w]*(+Hy[w_rht] - Hy[w] -Hx[w_top]+Hx[w]);
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
    const int w_btm = field_subBottom(w);   //一つ上
    const int w_bck = field_subBack(w);  //1つ前        
    const dcomplex nowMx = Mx[w];    
    Mx[w] = C_MX[w]*Mx[w] - C_MXEYEZ[w]*(Ez[w]-Ez[w_btm] -Ey[w]+Ey[w_bck]); //原因
    Bx[w] = C_BX[w]*Bx[w] + C_BXMX1[w]*Mx[w] - C_BXMX0[w]*nowMx;
    Hx[w] = Bx[w]/MU_0_S;
  }

  //Y
  FAST_3FOR(w, sInfo, toNextX)
  {
    const int w_lft = field_subLeft(w);
    const int w_bck = field_subBack(w);  //1つ前
    const dcomplex nowMy = My[w];
    My[w] = C_MY[w]*My[w] - C_MYEXEZ[w]*( Ex[w]-Ex[w_bck] -Ez[w]+Ez[w_lft]); //原因
    By[w] = C_BY[w]*By[w] + C_BYMY1[w]*My[w] - C_BYMY0[w]*nowMy;
    Hy[w] = By[w]/MU_0_S;
  }

  //Z
  FAST_3FOR(w, sInfo, toNextX)
  {
    const int w_lft = field_subLeft(w);
    const int w_btm = field_subBottom(w);   //一つ上
    const dcomplex nowMz = Mz[w];
    Mz[w] = C_MZ[w]*Mz[w] - C_MZEXEY[w]*( Ey[w]-Ey[w_lft] -Ex[w]+Ex[w_btm] );
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
        int x = i+sInfo.OFFSET_X; //-1する必要はない.
        int y = j+sInfo.OFFSET_Y;
        int z = k+sInfo.OFFSET_Z;
        
        EPS_EX[w] = models_eps(x    ,y+0.5,z+0.5,D_Y); //todo 
        EPS_EY[w] = models_eps(x+0.5,y    ,z+0.5,D_X);
        EPS_EZ[w] = models_eps(x+0.5,y+0.5,z    ,D_XY);

        sig_ex_x = sig_max*field_sigmaX(x    ,y+0.5,z+0.5);
        sig_ex_y = sig_max*field_sigmaY(x    ,y+0.5,z+0.5);
        sig_ex_z = sig_max*field_sigmaZ(x    ,y+0.5,z+0.5);        
        sig_ey_x = sig_max*field_sigmaX(x+0.5,y    ,z+0.5);
        sig_ey_y = sig_max*field_sigmaY(x+0.5,y    ,z+0.5);
        sig_ey_z = sig_max*field_sigmaZ(x+0.5,y    ,z+0.5);        
        sig_ez_x = sig_max*field_sigmaX(x+0.5,y+0.5,z    );
        sig_ez_y = sig_max*field_sigmaY(x+0.5,y+0.5,z    );
        sig_ez_z = sig_max*field_sigmaZ(x+0.5,y+0.5,z    );

        sig_hx_x = sig_max*field_sigmaX(x+0.5,y    ,z    );
        sig_hx_y = sig_max*field_sigmaY(x+0.5,y    ,z    );
        sig_hx_z = sig_max*field_sigmaZ(x+0.5,y    ,z    );
        sig_hy_x = sig_max*field_sigmaX(x    ,y+0.5,z    );
        sig_hy_y = sig_max*field_sigmaY(x    ,y+0.5,z    );
        sig_hy_z = sig_max*field_sigmaZ(x    ,y+0.5,z    );
        sig_hz_x = sig_max*field_sigmaX(x    ,y    ,z+0.5);
        sig_hz_y = sig_max*field_sigmaY(x    ,y    ,z+0.5);
        sig_hz_z = sig_max*field_sigmaZ(x    ,y    ,z+0.5);

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
  for(int i=1; i<subInfo_s.SUB_N_PX-1; i++)  
    for(int j=1; j<subInfo_s.SUB_N_PY-1; j++)
      for(int k=1; k<subInfo_s.SUB_N_PZ-1; k++)
      {
        int w = field_index(i+dx,j+dy,k+dz);
        entire[w] = region[field_subIndex(i,j,k)];
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
    cpy(entire, phi, subInfo_s.OFFSET_X, subInfo_s.OFFSET_Y, subInfo_s.OFFSET_Z);

    dcomplex *tmp = newDComplex(subInfo_s.SUB_N_CELL);
    int offset[3];
    for(int i=1; i<subInfo_s.Nproc; i++)
    {
      MPI_Recv(offset, 3, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(tmp, subInfo_s.SUB_N_CELL, MPI_C_DOUBLE_COMPLEX, i, 0, MPI_COMM_WORLD, &status);

      cpy(entire, tmp, offset[0], offset[1], offset[2]);      
    }
    free(tmp);
    return entire;
  }
  else {
    int offset[3];
    offset[0] = subInfo_s.OFFSET_X;
    offset[1] = subInfo_s.OFFSET_Y;
    offset[2] = subInfo_s.OFFSET_Z;
    MPI_Send(offset, 3, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(phi, subInfo_s.SUB_N_CELL, MPI_C_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD);
    
    return NULL; //マスター以外はNULLを返す.
  }
}

//---------------------メモリの解放--------------------//

static void output()
{  
  dcomplex *entireEx = unifyToRank0(Ex);
  dcomplex *entireEy = unifyToRank0(Ey);
  dcomplex *entireEz = unifyToRank0(Ez);
  dcomplex *entireHx = unifyToRank0(Hx);
  dcomplex *entireHy = unifyToRank0(Hy);
  dcomplex *entireHz = unifyToRank0(Hz);

  MPI_Barrier(MPI_COMM_WORLD);
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  if(subInfo_s.Rank == 0)
  {
    field_outputAllDataComplex("Ex.txt", entireEx);
    field_outputAllDataComplex("Ey.txt", entireEy);    
    field_outputAllDataComplex("Ez.txt", entireEz);
    field_outputAllDataComplex("Hx.txt", entireHx);
    field_outputAllDataComplex("Hy.txt", entireHy);
    field_outputAllDataComplex("Hz.txt", entireHz);
    
    field_outputElliptic("Ex_xy.txt", entireEx, 0);
    field_outputElliptic("Ex_zy.txt", entireEx, 1);
    field_outputElliptic("Ex_xz.txt", entireEx, 2);    
    field_outputElliptic("Ey_xy.txt", entireEy, 0);
    field_outputElliptic("Ey_zy.txt", entireEy, 1);
    field_outputElliptic("Ey_xz.txt", entireEy, 2);
    field_outputElliptic("Ez_xy.txt", entireEz, 0);
    field_outputElliptic("Ez_zy.txt", entireEz, 1);
    field_outputElliptic("Ez_xz.txt", entireEz, 2);    
  
    ntff3D_Frequency(entireEx,entireEy,entireEz,entireHx,entireHy,entireHz);
    free(entireEx);
    free(entireEy);
    free(entireEz);
    free(entireHx);
    free(entireHy);
    free(entireHz);
  }  
  //勝手にfreeしないように吐き出しが終わるまでは,待つ.
  MPI_Barrier(MPI_COMM_WORLD);
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
#ifdef DEBUG
static void outputAllDataComplex(const char *fileName, dcomplex* data)
{
  FILE *fp = openFile(fileName);
  SubFieldInfo_S sInfo_s = field_getSubFieldInfo_S();
  for(int i=0; i<sInfo_s.SUB_N_CELL; i++)
  {
    fprintf(fp,  "%.18lf, %.18lf \n", creal(data[i]), cimag(data[i]));
  }
  fclose(fp);
}

static void debugOutput()
{
}

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
#endif
