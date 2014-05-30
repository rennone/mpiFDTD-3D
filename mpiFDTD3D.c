#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "field.h"
#include "models.h"
#include "mpiFDTD3D.h"
#include "myComplex.h"

/* about MPI  */
static int rank;      //MPIのランク
static int nproc;     //全プロセス数
static int offsetX, offsetY, offsetZ;  //計算領域における, このミニフィールドのオフセット量
static int ltRank, rtRank, tpRank, bmRank; //左右上下のランク
static int SUB_N_X, SUB_N_Y, SUB_N_Z;
static int SUB_N_PX, SUB_N_PY, SUB_N_PZ;
static int SUB_N_CELL;
static MPI_Datatype X_DIRECTION_DOUBLE_COMPLEX;


//Ex(i+0.5,j) -> Ex[i,j]
//Ey(i,j+0.5) -> Ey[i,j]
//Hz(i+0.5,j+0.5) -> Hz[i,j]
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

static void calcJD(void);
static void calcE(void);
static void calcMB(void);
static void calcH(void);
static void allocateMemories(void);
static void setCoefficient(void);
static void freeMemories(void);

static void initializeElectroMagneticField(void);

static void initMpi(void);

dcomplex* fdtd3D_getEx(void)
{
  return Ex;
}
dcomplex* fdtd3D_getEy(void)
{
  return Ey;
}
dcomplex* fdtd3D_getEz(void)
{
  return Ez;
}
dcomplex* fdtd3D_getHz(void)
{
  return Hz;
}
dcomplex* fdtd3D_getHx(void)
{
  return Hx;
}
dcomplex* fdtd3D_getHy(void)
{
  return Hy;
}

//-------------------- index method --------------------//
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
  return i - 1;//手前
}

//Standard Scattered Wave
static void scatteredWave(dcomplex *p, double *eps){
  double time     = field_getTime();
  double w_s      = field_getOmega();
  double ray_coef = field_getRayCoef();
  double k_s      = field_getK();  
  double rad      = field_getWaveAngle()*M_PI/180;	//ラジアン変換
  double ks_cos = cos(rad)*k_s, ks_sin = sin(rad)*k_s;	//毎回計算すると時間かかりそうだから,代入しておく

  SubFieldInfo_S sInfo_s = field_getSubFieldInfo_S();
  int i,j,k;
  for(i=1; i<SUB_N_PX-1; i++)
    for(j=1; j<SUB_N_PY-1; j++)
      for(k=1; k<SUB_N_PZ; k++)
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

//Update
static void update(void)
{
  calcJD();
  calcE();
  scatteredWave(Ez, EPS_EZ);
//  MPI_Barrier(MPI_COMM_WORLD);  
  //Connection_SendRecvE();
//  Connection_ISend_IRecvE();
  calcMB();
  calcH();
//  Connection_ISend_IRecvH();
  //Connection_SendRecvH();
}

//calculate J and D
static  void calcJD()
{
  int i,j,k;
  for(i=1; i<SUB_N_PX-1; i++)
    for(j=1; j<SUB_N_PY-1; j++)
      for(k=1; j<SUB_N_PZ-1; k++)
      {
        const int w = subInd(i,j,k);
        const int w_lft = subIndLeft(w);    //一つ左
        const int w_btm = subIndBottom(w);  //一つ下
        const int w_bck = subIndBack(w);
        const dcomplex nowJx = Jx[w];
        const dcomplex nowJy = Jy[k];
        const dcomplex nowJz = Jz[w];

        Jx[w] = C_JX[w]*Jx[w] + C_JXHYHZ[w]*(Hz[w]-Hz[w_btm] -Hy[w_bck]+Hy[w]);
        Dx[w] = C_DX[w]*Dx[w] + C_DXJX1[w]*Jx[w] - C_DXJX0[w]*nowJx;
        
        Jy[w] = C_JY[w]*Jy[w] + C_JYHXHZ[w]*(Hx[w_bck]-Hx[w] -Hz[w]+Hz[w_lft] );
        Dy[w] = C_DY[w]*Dy[w] + C_DYJY1[w]*Jy[w] - C_DYJY0[w]*nowJy;

        Jz[w] = C_JZ[w]*Jz[w] + C_JZHXHY[w]*(+Hy[w] - Hy[w_lft] -Hx[w]+Hx[w_btm]);
        Dz[w] = C_DZ[w]*Dz[w] + C_DZJZ1[w]*Jz[w] - C_DZJZ0[w]*nowJz;
      } 
}

//calculate E 
static void calcE()
{
  int i,j,k;
  for(i=1; i<SUB_N_PX-1; i++)
    for(j=1; j<SUB_N_PY-1; j++)
      for(k=1; k<SUB_N_PZ-1; k++)
      {
        const int w = subInd(i,j,k);
        Ex[w] = Dx[w]/EPS_EX[w];
        Ey[w] = Dy[w]/EPS_EY[w];
        Ez[w] = Dz[w]/EPS_EZ[w];
      } 
}

//calculate M and B
static void calcMB()
{
  int i,j,k;
  for(i=1; i<SUB_N_PX-1; i++)
    for(j=1; j<SUB_N_PY-1; j++)
      for(k=1; j<SUB_N_PZ-1; k++)
      {
        const int w = subInd(i,j,w);
        const int w_rht = subIndRight(w);
        const int w_top = subIndTop(w); //一つ上
        const int w_frt = subIndFront(w);
        const dcomplex nowMx = Mx[w];
        const dcomplex nowMy = My[w];
        const dcomplex nowMz = Mz[w];
      
        Mx[w] = C_MX[w]*Mx[w] - C_MXEYEZ[w]*(Ez[w_top]-Ez[w] -Ey[w]+Ey[w_frt]);
        Bx[w] = C_BX[w]*Bx[w] + C_BXMX1[w]*Mx[w] - C_BXMX0[w]*nowMx;
      
        My[w] = C_MY[w]*My[w] - C_MYEXEZ[w]*(Ex[w]-Ex[w_frt] -Ez[w_rht]+Ez[w]);
        By[w] = C_BY[w]*By[w] + C_BYMY1[w]*My[w] - C_BYMY0[w]*nowMy;

        Mz[w] = C_MZ[w]*Mz[w] - C_MZEXEY[w]*(Ey[w_rht]-Ey[w] -Ex[w_top]+Ex[w]);
        Bz[w] = C_BZ[w]*Bz[w] + C_BZMZ1[w]*Mz[w] - C_BZMZ0[w]*nowMz;
      }
}

//calculate H
static void calcH()
{
  int i,j,k;
  for(i=1; i<SUB_N_PX-1; i++)
    for(j=1; j<SUB_N_PY-1; j++)
      for(k=1; j<SUB_N_PZ-1; k++)
      {
        const int w = subInd(i,j,k);
        Hx[w] = Bx[w]/MU_0_S;
        Hy[w] = By[w]/MU_0_S;
        Hz[w] = Bz[w]/MU_0_S;
      }
}

//-----------------memory allocate-------------//
static void init(){
  initMpi();
  allocateMemories();
  setCoefficient();
}

static void initMpi()
{
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  int dim = 3;          //number of dimension is 2
  int procs[3] = {0,0,0};       //[0]: x方向の分割数, [1]:y方向の分割数, [2]:z方向の分割数
  int period[3] = {0,0,0};//境界条件, 固定境界
  MPI_Comm grid_comm;
  int reorder = 1;   //re-distribute rank flag

  MPI_Dims_create(nproc, dim, procs);
  MPI_Cart_create(MPI_COMM_WORLD, dim, procs, period, reorder, &grid_comm);
  MPI_Cart_shift(grid_comm, 0, 1, &ltRank, &rtRank);
  MPI_Cart_shift(grid_comm, 1, 1, &bmRank, &tpRank);

  //プロセス座標において自分がどの位置に居るのか求める(何行何列に居るか)
  int coordinates[2];
  MPI_Comm_rank(grid_comm, &rank);
  MPI_Cart_coords(grid_comm, rank, 2, coordinates);
  
  SUB_N_X = N_PX / procs[0];
  SUB_N_Y = N_PY / procs[1];
  SUB_N_Z = N_PZ / procs[1];
  
  SUB_N_PX = SUB_N_X + 2; //のりしろの分2大きい
  SUB_N_PY = SUB_N_Y + 2; //のりしろの分2大きい
  SUB_N_PZ = SUB_N_Z + 2; //のりしろの分2大きい
  SUB_N_CELL = SUB_N_PX*SUB_N_PY*SUB_N_PZ;
  
  offsetX = coordinates[0] * SUB_N_X; //ランクのインデックスではなく, セル単位のオフセットなのでSUB_N_Xずれる
  offsetY = coordinates[1] * SUB_N_Y;
  offsetY = coordinates[1] * SUB_N_Z;
  
/*これだと, 1個のデータをSUB_N_PY跳び(次のデータまでSUB_N_PY-1個隙間がある),SUB_N_X行ぶん取ってくる事になる */
  //todo 
  MPI_Type_vector(SUB_N_X, 1, SUB_N_PY, MPI_C_DOUBLE_COMPLEX, &X_DIRECTION_DOUBLE_COMPLEX); 
  MPI_Type_commit(&X_DIRECTION_DOUBLE_COMPLEX);
}

static void allocateMemories()
{
  Ex = newDComplex(SUB_N_CELL);   //Ex(i+0.5,j    ,k+0.5) -> Ex[i,j,k]
  Ey = newDComplex(SUB_N_CELL);   //Ey(i    ,j+0.5,k+0.5) -> Ey[i,j,k]
  Ez = newDComplex(SUB_N_CELL);   //Ez(i    ,j    ,k    ) -> Ez[i,j,k]
  
  Jx = newDComplex(SUB_N_CELL);
  Jy = newDComplex(SUB_N_CELL);
  Jz = newDComplex(SUB_N_CELL);
  
  Dx = newDComplex(SUB_N_CELL);
  Dy = newDComplex(SUB_N_CELL);
  Dz = newDComplex(SUB_N_CELL);
  
  Hx = newDComplex(SUB_N_CELL);  //Hx(i    , j+0.5, k    )->Hx[i,j,k]
  Hy = newDComplex(SUB_N_CELL);  //Hy(i+0.5, j    , k    )->Hy[i,j,k]
  Hz = newDComplex(SUB_N_CELL);  //Hz(i+0.5, j+0.5, k+0.5)->Hz[i,j,k]
  
  Mx = newDComplex(SUB_N_CELL);
  My = newDComplex(SUB_N_CELL);
  Mz = newDComplex(SUB_N_CELL);
  
  Bx = newDComplex(SUB_N_CELL);
  By = newDComplex(SUB_N_CELL);
  Bz = newDComplex(SUB_N_CELL);

  C_JX = newDouble(SUB_N_CELL);
  C_JY = newDouble(SUB_N_CELL);
  C_JZ = newDouble(SUB_N_CELL);
  
  C_MX = newDouble(SUB_N_CELL);
  C_MY = newDouble(SUB_N_CELL);
  C_MZ = newDouble(SUB_N_CELL);

  C_DX = newDouble(SUB_N_CELL);
  C_DY = newDouble(SUB_N_CELL);
  C_DZ = newDouble(SUB_N_CELL);
  
  C_BX = newDouble(SUB_N_CELL);
  C_BY = newDouble(SUB_N_CELL);
  C_BZ = newDouble(SUB_N_CELL);

  C_JXHYHZ = newDouble(SUB_N_CELL);
  C_JYHXHZ = newDouble(SUB_N_CELL);
  C_JZHXHY = newDouble(SUB_N_CELL);
  
  C_MXEYEZ = newDouble(SUB_N_CELL);
  C_MYEXEZ = newDouble(SUB_N_CELL);
  C_MZEXEY = newDouble(SUB_N_CELL);
  
  C_DXJX0 = newDouble(SUB_N_CELL);
  C_DXJX1 = newDouble(SUB_N_CELL);
  C_DYJY0 = newDouble(SUB_N_CELL);
  C_DYJY1 = newDouble(SUB_N_CELL);
  C_DZJZ0 = newDouble(SUB_N_CELL);
  C_DZJZ1 = newDouble(SUB_N_CELL);

  C_BXMX1 = newDouble(SUB_N_CELL);
  C_BXMX0 = newDouble(SUB_N_CELL);
  C_BYMY1 = newDouble(SUB_N_CELL);
  C_BYMY0 = newDouble(SUB_N_CELL);
  C_BZMZ0 = newDouble(SUB_N_CELL);
  C_BZMZ1 = newDouble(SUB_N_CELL);

  EPS_EX = newDouble(SUB_N_CELL);
  EPS_EY = newDouble(SUB_N_CELL);
  EPS_EZ = newDouble(SUB_N_CELL);

  EPS_HY = newDouble(SUB_N_CELL);
  EPS_HX = newDouble(SUB_N_CELL);
  EPS_HZ = newDouble(SUB_N_CELL);

}

static void initializeElectroMagneticField()
{
  memset(Ex, 0, sizeof(dcomplex)*SUB_N_CELL);
  memset(Ey, 0, sizeof(dcomplex)*SUB_N_CELL);
  memset(Ez, 0, sizeof(dcomplex)*SUB_N_CELL);

  memset(Hx, 0, sizeof(dcomplex)*SUB_N_CELL);
  memset(Hy, 0, sizeof(dcomplex)*SUB_N_CELL);
  memset(Hz, 0, sizeof(dcomplex)*SUB_N_CELL);

  memset(Jx, 0, sizeof(dcomplex)*SUB_N_CELL);
  memset(Jy, 0, sizeof(dcomplex)*SUB_N_CELL);
  memset(Jz, 0, sizeof(dcomplex)*SUB_N_CELL);
  
  memset(Mx, 0, sizeof(dcomplex)*SUB_N_CELL);
  memset(My, 0, sizeof(dcomplex)*SUB_N_CELL);
  memset(Mz, 0, sizeof(dcomplex)*SUB_N_CELL);

  memset(Dz, 0, sizeof(dcomplex)*SUB_N_CELL);
  memset(Dx, 0, sizeof(dcomplex)*SUB_N_CELL);
  memset(Dy, 0, sizeof(dcomplex)*SUB_N_CELL);

  memset(Bx, 0, sizeof(dcomplex)*SUB_N_CELL);
  memset(By, 0, sizeof(dcomplex)*SUB_N_CELL);
  memset(Bz, 0, sizeof(dcomplex)*SUB_N_CELL);    
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
  
  const double sig_max = -(M+1.0)*EPSILON_0_S*LIGHT_SPEED_S/2.0/N_PML*log(R);

  int i,j,k;
  for(i=1; i<SUB_N_PX-1; i++)
    for(j=1; j<SUB_N_PY-1; j++)
      for(k=1; k<SUB_N_PZ; j++)
      {
        int w = subInd(i,j,k);
        int x = i-1+offsetX;
        int y = j-1+offsetY;
        int z = z-1+offsetZ;
        EPS_EX[k] = models_eps(x+0.5,y    ,z+0.5,D_Y);
        EPS_EY[k] = models_eps(x    ,y+0.5,z+0.5,D_X);
        EPS_EZ[k] = models_eps(x    ,y    ,z    ,D_XY);

        EPS_HX[k] = models_eps(x    ,y+0.5,z    ,D_Y);
        EPS_HY[k] = models_eps(x+0.5,y    ,z    ,D_X);
        EPS_HZ[k] = 0.5*(models_eps(x+0.5,y+0.5,z+0.5, D_X) + models_eps(x+0.5,y+0.5,z+0.5, D_Y));

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
        C_MXEYEZ[w]  = (2*eps)            / (2*eps + sig_hx_y);
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
