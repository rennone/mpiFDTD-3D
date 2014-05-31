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
static int rank;      //MPIのランク
static int nproc;     //全プロセス数
static MPI_Datatype X_DIRECTION_DOUBLE_COMPLEX;
static MPI_Datatype DCOMPLEX_XY; //XY平面
static MPI_Datatype DCOMPLEX_YZ; //YZ平面
static MPI_Datatype DCOMPLEX_XZ; //XZ平面

//系は右手系
//x(left-, right+)
//y(bottom-, top+)
//z(back-, front+)
//Hx(i,j+0.5,k)         -> Hx[i,j,k]
//Hy(i+0.5,j,k)         -> Hy[i,j,k]
//Ez(i,j,k)             -> Ez[i,j,k]

//Ex(i+0.5,j,k)         -> Ex[i,j,k]
//Ey(i,j+0.5,k)         -> Ey[i,j,k]
//Hz(i+0.5,j+0.5,k+0.5) -> Hz[i,j,k]

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
static double *EPS_HX=NULL, *EPS_HY=NULL, *EPS_HZ=NULL;

static void update(void);
static void finish(void);
static void init(void);
static void reset(void);

static void calcJD(void);
static void calcE(void);
static void calcMB(void);
static void calcH(void);
static void allocateMemories(void);
static void setCoefficient(void);
static void freeMemories(void);

static void initializeElectroMagneticField(void);

//static void initMpi(void);

dcomplex* mpi_fdtd3D_upml_getEx(void){  return Ex;}
dcomplex* mpi_fdtd3D_upml_getEy(void){  return Ey;}
dcomplex* mpi_fdtd3D_upml_getEz(void){  return Ez;}
dcomplex* mpi_fdtd3D_upml_getHx(void){  return Hx;}
dcomplex* mpi_fdtd3D_upml_getHy(void){  return Hy;}
dcomplex* mpi_fdtd3D_upml_getHz(void){  return Hz;}

double*   mpi_fdtd3D_upml_getEps(void){  return EPS_EZ;}

void (* mpi_fdtd3D_upml_getUpdate(void))(void)
{
  return update;
}
void (* mpi_fdtd3D_upml_getFinish(void))(void)
{
  return finish;
}
void (* mpi_fdtd3D_upml_getReset(void))(void)
{
  return reset;
}
void (* mpi_fdtd3D_upml_getInit(void))(void)
{
  return init;
}

//-------------------- index method --------------------//
/*
static  int subInd(const int i, const int j, const int k)
{
  return i*SUB_N_PY*SUB_N_PZ + j*SUB_N_PZ + k;
}
static  int subIndLeft(const int i)
{
  return i - SUB_N_PY*SUB_N_PZ;//左
}
static  int subIndRight(const int i)
{
  return i + SUB_N_PY*SUB_N_PZ;//右
}
static  int subIndTop(const int i)
{
  return i + SUB_N_PZ;//上
}
//
static  int subIndBottom(const int i)
{
  return i - SUB_N_PZ;//下
}

static  int subIndBack(const int i)
{
  return i + 1;//奥
}

static  int subIndFront(const int i)
{
  NOT_DONE("subIndFront => this is left hand\n");
  return i - 1;//手前
}
*/
#define SubIndexLeft(sInfo, ind)   ind - sInfo.SUB_N_PYZ
#define SubIndexRight(sInfo, ind)  ind + sInfo.SUB_N_PYZ
#define SubIndexTop(sInfo, ind)    ind + sInfo.SUB_N_PZ
#define SubIndexBottom(sInfo, ind) ind - sInfo.SUB_N_PZ
#define SubIndexFront(sInfo, ind)  ind + 1
#define SubIndexBack(sInfo, ind)   ind - 1

//Standard Scattered Wave
static void scatteredWave(dcomplex *p, double *eps){
  double time     = field_getTime();
  double w_s      = field_getOmega();
  double ray_coef = field_getRayCoef();
  double k_s      = field_getK();  
  double rad      = field_getWaveAngle()*M_PI/180;	//ラジアン変換
  double ks_cos = cos(rad)*k_s, ks_sin = sin(rad)*k_s;	//毎回計算すると時間かかりそうだから,代入しておく

  SubFieldInfo_S sInfo_s = field_getSubFieldInfo_S();
  for(int i=1; i<sInfo_s.SUB_N_PX-1; i++)
    for(int j=1; j<sInfo_s.SUB_N_PY-1; j++)
      for(int k=1; k<sInfo_s.SUB_N_PZ; k++)
      {
        int l = field_subIndex(i,j,k);
        int x = i-1+sInfo_s.OFFSET_X;
        int y = j-1+sInfo_s.OFFSET_Y;
        int z = k-1+sInfo_s.OFFSET_Z;

        NOT_DONE("check scattered wave");
        double ikx = x*ks_cos + y*ks_sin; //k_s*(i*cos + j*sin)
        p[l] += ray_coef*(EPSILON_0_S/eps[l] - 1)*(cos(ikx-w_s*time) + I*sin(ikx-w_s*time));
      }   
  
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
          //printf("%lf, %lf, %d, %d, %d\n", creal(p[w]), cimag(p[w]), i, j, k);
          return false;
//return true;
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
//  const int w_bck = field_subBack(w);  //1つ前
//  const int w_frt = field_subFront(w);  //1つ前

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
/*
  printf("Ex[frt-w] %lf %lf  Ez[rht-w] %lf %lf \n",
         creal(Ex[field_subFront(i)]), creal(Ex[i]),
         creal(Ez[field_subRight(i)]), creal(Ez[i])),
  
  printf("E(%lf,%lf,%lf) H(%lf,%lf,%lf)  ",
         creal(Ex[i]), creal(Ey[i]), creal(Ez[i]),
    creal(Hx[i]), creal(Hy[i]), creal(Hz[i]));

  printf("J(%lf,%lf,%lf) M(%lf,%lf,%lf)%\n ",
         creal(Jx[i]), creal(Jy[i]), creal(Jz[i]),
         creal(Bx[i]), creal(By[i]), creal(Bz[i]));*/
}

//Update
static void update(void)
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();
/*
  //軸を描く為の処理
  for(int i=1; i<sInfo.SUB_N_PY-1; i++)
  {
    int indexX = field_subIndex(i, sInfo.SUB_N_PY/2, sInfo.SUB_N_PZ/2);
    int indexY = field_subIndex(sInfo.SUB_N_PX/2, i,  sInfo.SUB_N_PZ/2);
    int indexZ = field_subIndex(sInfo.SUB_N_PX/2, sInfo.SUB_N_PY/2, i);
    Ez[indexZ] = 1.0;
  }
  return;
*/
  
  calcJD();
  if( debugCheck(Jz) )
  {
    STOP("Jz\n");
  }
  if( debugCheck(Jx) )
  {
    STOP("Jx\n");
  }
  if( debugCheck(Jy) )
  {
    STOP("Jy\n");
  }
  calcE();
  int i = field_subIndex(sInfo.SUB_N_PX/2, sInfo.SUB_N_PY/2, sInfo.SUB_N_PZ/2);
  dcomplex wave = 10.0*field_getRayCoef()*cexp(-I*field_getOmega()*field_getTime());
  Ez[i] += wave;

  calcMB();
  if( debugCheck(My) )
  {
    STOP("My\n");
  }
  if( debugCheck(Mx) )
  {
    STOP("Mx\n");
  }  
  if( debugCheck(Mz) )
  {
    STOP("Mz\n");
  }

  calcH();
  // Hz[i] += wave;

  debugPrint();
  //scatteredWave(Ez, EPS_EZ);
}

//calculate J and D
static  void calcJD()
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();

  for(int i=1; i<sInfo.SUB_N_PX-1; i++)
    for(int j=1; j<sInfo.SUB_N_PY-1; j++)
      for(int k=1; k<sInfo.SUB_N_PZ-1; k++)        
      {
        const int w = field_subIndex(i,j,k);
        const int w_lft = field_subLeft(w);   //一つ左
        const int w_btm = field_subBottom(w); //一つ下
        const int w_frt = field_subFront(w);  //todo
        const int w_bck = field_subBack(w);  //todo
        const dcomplex nowJx = Jx[w];
        const dcomplex nowJy = Jy[w];
        const dcomplex nowJz = Jz[w];

        Jx[w] = C_JX[w]*Jx[w] + C_JXHYHZ[w]*(Hz[w]-Hz[w_btm] -Hy[w_frt]+Hy[w]);
        Dx[w] = C_DX[w]*Dx[w] + C_DXJX1[w]*Jx[w] - C_DXJX0[w]*nowJx;
        
        Jy[w] = C_JY[w]*Jy[w] + C_JYHXHZ[w]*( Hx[w_frt]-Hx[w] -Hz[w]+Hz[w_lft] );
        Dy[w] = C_DY[w]*Dy[w] + C_DYJY1[w]*Jy[w] - C_DYJY0[w]*nowJy;

        Jz[w] = C_JZ[w]*Jz[w] + C_JZHXHY[w]*(+Hy[w] - Hy[w_lft] -Hx[w]+Hx[w_btm]);
        Dz[w] = C_DZ[w]*Dz[w] + C_DZJZ1[w]*Jz[w] - C_DZJZ0[w]*nowJz;
      } 
}

//calculate E 
static void calcE()
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();
  for(int i=1; i<sInfo.SUB_N_PX-1; i++)
    for(int j=1; j<sInfo.SUB_N_PY-1; j++)
      for(int k=1; k<sInfo.SUB_N_PZ-1; k++)
      {
        const int w = field_subIndex(i,j,k);
        Ex[w] = Dx[w]/EPS_EX[w];
        Ey[w] = Dy[w]/EPS_EY[w];
        Ez[w] = Dz[w]/EPS_EZ[w];
      }
}

//calculate M and B
static void calcMB()
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();

  for(int i=1; i<sInfo.SUB_N_PX-1; i++)
    for(int j=1; j<sInfo.SUB_N_PY-1; j++)
      for(int k=1; k<sInfo.SUB_N_PZ-1; k++)
      {
        const int w     = field_subIndex(i,j,k);
        const int w_rht = field_subRight(w);
        const int w_top = field_subTop(w);   //一つ上
        const int w_bck = field_subBack(w);  //1つ前
        const int w_frt = field_subFront(w);  //1つ前
        
        const dcomplex nowMx = Mx[w];
        const dcomplex nowMy = My[w];
        const dcomplex nowMz = Mz[w];
    
        Mx[w] = C_MX[w]*Mx[w] - C_MXEYEZ[w]*(Ez[w_top]-Ez[w] -Ey[w_frt]+Ey[w]); //原因
        Bx[w] = C_BX[w]*Bx[w] + C_BXMX1[w]*Mx[w] - C_BXMX0[w]*nowMx;
      
        My[w] = C_MY[w]*My[w] - C_MYEXEZ[w]*( Ex[w_frt]-Ex[w] -Ez[w_rht]+Ez[w]); //原因
        By[w] = C_BY[w]*By[w] + C_BYMY1[w]*My[w] - C_BYMY0[w]*nowMy;
        
        Mz[w] = C_MZ[w]*Mz[w] - C_MZEXEY[w]*( Ey[w_rht]-Ey[w] -Ex[w_top]+Ex[w] );
        Bz[w] = C_BZ[w]*Bz[w] + C_BZMZ1[w]*Mz[w] - C_BZMZ0[w]*nowMz;
      }
}

//calculate H
static void calcH()
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();
  for(int i=1; i<sInfo.SUB_N_PX-1; i++)
    for(int j=1; j<sInfo.SUB_N_PY-1; j++)
      for(int k=1; k<sInfo.SUB_N_PZ-1; k++)
      {
        const int w = field_subIndex(i,j,k);
        Hx[w] = Bx[w]/MU_0_S;
        Hy[w] = By[w]/MU_0_S;
        Hz[w] = Bz[w]/MU_0_S;
      }
}

//-----------------memory allocate-------------//
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

  EPS_HY = newDouble(sInfo.SUB_N_CELL);
  EPS_HX = newDouble(sInfo.SUB_N_CELL);
  EPS_HZ = newDouble(sInfo.SUB_N_CELL);

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

        EPS_HX[w] = models_eps(x    ,y+0.5,z    ,D_Y);
        EPS_HY[w] = models_eps(x+0.5,y    ,z    ,D_X);
        EPS_HZ[w] = 0.5*(models_eps(x+0.5,y+0.5,z+0.5, D_X) + models_eps(x+0.5,y+0.5,z+0.5, D_Y));

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

//---------------------メモリの解放--------------------//
static void finish(){
  //output();
  freeMemories();
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
  if(EPS_EZ != NULL)   free(EPS_HZ);  

  if(EPS_HX != NULL)   free(EPS_EX);
  if(EPS_HY != NULL)   free(EPS_EY);
  if(EPS_HZ != NULL)   free(EPS_HZ);  
}
