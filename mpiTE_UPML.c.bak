#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "mpiTE_UPML.h"
#include "field.h"
#include "models.h"

/* about MPI  */
static int rank;      //MPIのランク
static int nproc;     //全プロセス数
static int offsetX, offsetY;  //計算領域における, このミニフィールドのオフセット量
static int ltRank, rtRank, tpRank, bmRank; //左右上下のランク
static int SUB_N_X, SUB_N_Y;
static int SUB_N_PX, SUB_N_PY;
static int SUB_N_CELL;
static MPI_Datatype X_DIRECTION_DOUBLE_COMPLEX;

//Ex(i+0.5,j) -> Ex[i,j]
//Ey(i,j+0.5) -> Ey[i,j]
//Hz(i+0.5,j+0.5) -> Hz[i,j]
static double complex *Ex = NULL;
static double complex *Jx = NULL;
static double complex *Dx = NULL;

static double complex *Ey = NULL;
static double complex *Jy = NULL;
static double complex *Dy = NULL;

static double complex *Hz = NULL;
static double complex *Mz = NULL;
static double complex *Bz = NULL;

static double *C_JX = NULL, *C_JY = NULL, *C_MZ= NULL;
static double *C_JXHZ = NULL, *C_JYHZ = NULL, *C_MZEXEY= NULL;
static double *C_DX=NULL, *C_DY=NULL, *C_BZ=NULL;
static double *C_DXJX0=NULL, *C_DXJX1=NULL;
static double *C_DYJY0=NULL, *C_DYJY1=NULL;
static double *C_BZMZ0=NULL, *C_BZMZ1=NULL;

static double *EPS_EX=NULL, *EPS_EY=NULL, *EPS_HZ=NULL;

static const int timeStep = 500;
static double complex Uz[360][timeStep],Wx[360][timeStep], Wy[360][timeStep];

//------prototype--------//
static void update(void);
static void finish(void);
static void init(void);
static void initMpi();
static void allocateMemories();
static void freeMemories();
static void setCoefficient();
static void output();

static inline void calcE(void);
static inline void calcJD(void);
static inline void calcH(void);
static inline void calcMB(void);

//:public-------------------------------//
void (* fdtdTE_upml_getUpdate(void))(void)
{
  return update;
}

void (* fdtdTE_upml_getFinish(void))(void)
{
  return finish;
}

void (* fdtdTE_upml_getInit(void))(void)
{
  return init;
}

double complex* fdtdTE_upml_getEx(void){
  return Ex;
}

double complex* fdtdTE_upml_getEy(void){
  return Ey;
}

double complex* fdtdTE_upml_getHz(void){
  return Hz;
}

double* fdtdTE_upml_getEps(void){
  return EPS_EY;
}

inline int fdtdTE_upml_getSubNx(void){
  return SUB_N_X;
}

inline int fdtdTE_upml_getSubNy(void){
  return SUB_N_Y;
}

inline int fdtdTE_upml_getSubNpx(void){
  return SUB_N_PX;
}

inline int fdtdTE_upml_getSubNpy(void){
  return SUB_N_PY;
}

inline int fdtdTE_upml_getSubNcell(void){
  return SUB_N_CELL;
}
//:private ----------------------------//
static inline int subInd(const int i, const int j)
{
  return (i)*SUB_N_PY + (j);
}
static inline int subIndLeft(const int i){
  return i - SUB_N_PY;
}
static inline int subIndRight(const int i){
  return i + SUB_N_PY;
}
static inline int subIndTop(const int i){
  return i+1;
}
static inline int subIndBottom(const int i){
  return i-1;
}

//--------------------connection--------------------
static inline void Connection_ISend_IRecvE(void)
{
  MPI_Status status;
  MPI_Request req1, req2, req3, req4;
  //left & right
  //this needs only Ey[i-1, j] so send to right and recieve from left 
  int ltRecv = subInd(0,1);
  int rtSend = subInd(SUB_N_PX-2, 1);  
  MPI_Isend(&Ey[rtSend], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, rtRank, 1, MPI_COMM_WORLD, &req1);
  MPI_Irecv(&Ey[ltRecv], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, ltRank, 1, MPI_COMM_WORLD, &req2);
  
  //this needs only Ey[i,j-1] so send to top and recieve from bottom   
  int bmRecv = subInd(1,0);
  int tpSend = subInd(1,SUB_N_PY-2);
  MPI_Isend(&Ex[tpSend], 1, X_DIRECTION_DOUBLE_COMPLEX, tpRank, 1, MPI_COMM_WORLD, &req3);
  MPI_Irecv(&Ex[bmRecv], 1, X_DIRECTION_DOUBLE_COMPLEX, bmRank, 1, MPI_COMM_WORLD, &req4);

  MPI_Wait(&req1, &status);
  MPI_Wait(&req2, &status);
  MPI_Wait(&req3, &status);
  MPI_Wait(&req4, &status); 
}

static inline void Connection_ISend_IRecvH(void)
{
  MPI_Status status;
  MPI_Request req1, req2, req3, req4;

  //left & right
  //this needs only Hz[i+1, j] so send to left and recieve from right 
  int rtRecv = subInd(SUB_N_PX-1,1);
  int ltSend = subInd(         1,1);
  MPI_Isend(&Hz[ltSend], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, ltRank, 1, MPI_COMM_WORLD, &req1);
  MPI_Irecv(&Hz[rtRecv], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, rtRank, 1, MPI_COMM_WORLD, &req2);
  
  //bottom & top
//this needs only Hz[i, j+1] so send to bottom and recieve from top
  int bmSend = subInd(1,1);
  int tpRecv = subInd(1,SUB_N_PY-1);
  MPI_Isend(&Hz[bmSend], 1, X_DIRECTION_DOUBLE_COMPLEX, bmRank, 1, MPI_COMM_WORLD, &req3);
  MPI_Irecv(&Hz[tpRecv], 1, X_DIRECTION_DOUBLE_COMPLEX, tpRank, 1, MPI_COMM_WORLD, &req4);  
  MPI_Wait(&req1, &status);
  MPI_Wait(&req2, &status);
  MPI_Wait(&req3, &status);
  MPI_Wait(&req4, &status);
}

static inline void Connection_SendRecvE(void)
{
  MPI_Status status;  
  //left & right
  //this needs only Ez[i+1, j] so send to left and recieve from right 
  int rtRecv = subInd(SUB_N_PX-1,1);
  int ltSend = subInd(         1,1);
  MPI_Sendrecv(&Ey[ltSend], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, ltRank, 1,
               &Ey[rtRecv], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, rtRank, 1, MPI_COMM_WORLD, &status);
  
  //bottom & top
//this needs only Ez[i, j+1] so send to bottom and recieve from top
  int bmSend = subInd(1,1);
  int tpRecv = subInd(1,SUB_N_PY-1);
  MPI_Sendrecv(&Ex[bmSend], 1, X_DIRECTION_DOUBLE_COMPLEX, bmRank, 1,
               &Ex[tpRecv], 1, X_DIRECTION_DOUBLE_COMPLEX, tpRank, 1, MPI_COMM_WORLD, &status);
}

static inline void Connection_SendRecvH(void)
{
  MPI_Status status;
  //left & right
  //this needs only Hy[i-1, j] so send to right and recieve from left 
  int ltRecv = subInd(0,1);
  int rtSend = subInd(SUB_N_PX-2, 1);
  MPI_Sendrecv(&Hz[rtSend], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, rtRank, 1,
               &Hz[ltRecv], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, ltRank, 1,MPI_COMM_WORLD, &status);
  
  //this needs only Hy[i,j-1] so send to top and recieve from bottom   
  int bmRecv = subInd(1,0);
  int tpSend = subInd(1,SUB_N_PY-2);
  MPI_Sendrecv(&Hz[tpSend], 1, X_DIRECTION_DOUBLE_COMPLEX, tpRank, 1,
               &Hz[bmRecv], 1, X_DIRECTION_DOUBLE_COMPLEX, bmRank, 1, MPI_COMM_WORLD, &status); 
}

//Standard Scattered Wave
static inline void scatteredWave(double complex *p, double *eps){
  double time = field_getTime();
  double w_s  = field_getOmega();
  double ray_coef = field_getRayCoef();
  double k_s = field_getK();  
  double rad = field_getWaveAngle()*M_PI/180;	//ラジアン変換
  double ks_cos = cos(rad)*k_s, ks_sin = sin(rad)*k_s;	//毎回計算すると時間かかりそうだから,代入しておく
  int i,j;
  for(i=1; i<SUB_N_PX-1; i++){
    for(j=1; j<SUB_N_PY-1; j++){
      int k = subInd(i,j); 
      int x = i-1+offsetX;
      int y = j-1+offsetY;
      double ikx = x*ks_cos + y*ks_sin; //k_s*(i*cos + j*sin)
      p[k] += ray_coef*(EPSILON_0_S/eps[k] - 1)*(cos(ikx-w_s*time) + I*sin(ikx-w_s*time));
    }
  }
}

static inline void update(void)
{
  calcJD();
  calcE();
  scatteredWave(Ey, EPS_EY);
  Connection_SendRecvE();
  calcMB();
  calcH();
  Connection_SendRecvH();
}

static inline void calcJD(void)
{
  for(int i=1; i<SUB_N_PX-1; i++){
    for(int j=1; j<SUB_N_PY-1; j++){
      int k = subInd(i,j);
      int k_btm = subIndBottom(k);
      double complex nowJx = Jx[k];
      Jx[k] = C_JX[k]*Jx[k] + C_JXHZ[k]*(Hz[k] - Hz[k_btm]);
      Dx[k] = C_DX[k]*Dx[k] + C_DXJX1[k]*Jx[k] - C_DXJX0[k]*nowJx;
    }
  }
  
  for(int i=1; i<SUB_N_PX-1; i++){
    for(int j=1; j<SUB_N_PY-1; j++){
      int k = subInd(i,j);
      int k_lft = subIndLeft(k);
      double complex nowJy = Jy[k];
      Jy[k] = C_JY[k]*Jy[k] + C_JYHZ[k]*(-Hz[k] + Hz[k_lft]);
      Dy[k] = C_DY[k]*Dy[k] + C_DYJY1[k]*Jy[k] - C_DYJY0[k]*nowJy;
    }
  }
}

static inline void calcE(void)
{
  for(int i=1; i<SUB_N_PX-1; i++){
    for(int j=1; j<SUB_N_PY-1; j++){
      const int k = subInd(i,j);
      Ex[k] = Dx[k]/EPS_EX[k];
    }
  }
  
  for(int i=1; i<SUB_N_PX-1; i++){
    for(int j=1; j<SUB_N_PY-1; j++){
      const int k = subInd(i,j);
      Ey[k] = Dy[k]/EPS_EY[k];
    }
  }
}

static inline void calcMB(void)
{
  for(int i=1; i<SUB_N_PX-1; i++){
    for(int j=1; j<SUB_N_PY-1; j++){
      int k = subInd(i,j);
      int k_rht = subIndRight(k);
      int k_top = subIndTop(k);
      double complex nowMz = Mz[k];
      Mz[k] = C_MZ[k]*Mz[k] - C_MZEXEY[k]*(Ey[k_rht] - Ey[k] - Ex[k_top] + Ex[k]);
      Bz[k] = C_BZ[k]*Bz[k] + C_BZMZ1[k]*Mz[k] - C_BZMZ0[k]*nowMz; 
    }
  }
}

static inline void calcH(void)
{
  for(int i=1; i<SUB_N_PX-1; i++){
    for(int j=1; j<SUB_N_PY-1; j++){
      const int k = subInd(i,j);
      Hz[k] = Bz[k]/MU_0_S;
    }
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
  int dim = 2;          //number of dimension is 2
  int procs[2] = {0,0};         //[0]: x方向の分割数, [1]:y方向の分割数
  int period[2] = {0,0};//境界条件, 固定境界
  MPI_Comm grid_comm;
  int reorder = 1;   //re-distribute rank flag

  MPI_Dims_create(nproc, dim, procs);
  MPI_Cart_create(MPI_COMM_WORLD, 2, procs, period, reorder, &grid_comm);
  MPI_Cart_shift(grid_comm, 0, 1, &ltRank, &rtRank);
  MPI_Cart_shift(grid_comm, 1, 1, &bmRank, &tpRank);

  //プロセス座標において自分がどの位置に居るのか求める(何行何列に居るか)
  int coordinates[2];
  MPI_Comm_rank(grid_comm, &rank);
  MPI_Cart_coords(grid_comm, rank, 2, coordinates);
  
  SUB_N_X = N_PX / procs[0];
  SUB_N_Y = N_PY / procs[1];
  SUB_N_PX = SUB_N_X + 2; //のりしろの分2大きい
  SUB_N_PY = SUB_N_Y + 2; //のりしろの分2大きい
  SUB_N_CELL = SUB_N_PX*SUB_N_PY;  
  offsetX = coordinates[0] * SUB_N_X; //ランクのインデックスではなく, セル単位のオフセットなのでSUB_N_Xずれる
  offsetY = coordinates[1] * SUB_N_Y;

/*これだと, 1個のデータをSUB_N_PY跳び(次のデータまでSUB_N_PY-1個隙間がある),SUB_N_X行ぶん取ってくる事になる */
  MPI_Type_vector(SUB_N_X, 1, SUB_N_PY, MPI_C_DOUBLE_COMPLEX, &X_DIRECTION_DOUBLE_COMPLEX); 
  MPI_Type_commit(&X_DIRECTION_DOUBLE_COMPLEX);
}

static void allocateMemories()
{
  Ex = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);   //Ex(i+0.5,j) -> Ex[i,j]
  Ey = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);   //Ey(i,j+0.5) -> Ey[i,j]
  Hz = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);  //Hz(i+0.5, j+0.5)->Hz[i,j]

  Jx = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);   //Ex(i+0.5,j) -> Ex[i,j]
  Jy = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);   //Ey(i,j+0.5) -> Ey[i,j]
  Mz = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);  //Hz(i+0.5, j+0.5)->Hz[i,j]

  Dx = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);   //Ex(i+0.5,j) -> Ex[i,j]
  Dy = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);   //Ey(i,j+0.5) -> Ey[i,j]
  Bz = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);  //Hz(i+0.5, j+0.5)->Hz[i,j]

  C_JX = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_JY = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_MZ = (double *)malloc(sizeof(double)*SUB_N_CELL);  
  C_JXHZ = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_JYHZ = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_MZEXEY = (double *)malloc(sizeof(double)*SUB_N_CELL);

  C_DX = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_DY = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_BZ = (double *)malloc(sizeof(double)*SUB_N_CELL);  
  C_DXJX0 = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_DXJX1 = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_DYJY0 = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_DYJY1 = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_BZMZ0 = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_BZMZ1 = (double *)malloc(sizeof(double)*SUB_N_CELL);
  
  EPS_EX = (double *)malloc(sizeof(double)*SUB_N_CELL);
  EPS_EY = (double *)malloc(sizeof(double)*SUB_N_CELL);
  EPS_HZ = (double *)malloc(sizeof(double)*SUB_N_CELL);

  memset(Ex, 0, sizeof(double complex)*SUB_N_CELL);
  memset(Ey, 0, sizeof(double complex)*SUB_N_CELL);
  memset(Hz,0, sizeof(double complex)*SUB_N_CELL);
  memset(Jx, 0, sizeof(double complex)*SUB_N_CELL);
  memset(Jy, 0, sizeof(double complex)*SUB_N_CELL);
  memset(Mz,0, sizeof(double complex)*SUB_N_CELL);
  memset(Dx, 0, sizeof(double complex)*SUB_N_CELL);
  memset(Dy, 0, sizeof(double complex)*SUB_N_CELL);
  memset(Bz,0, sizeof(double complex)*SUB_N_CELL);
}

static void setCoefficient()
{
  //Hz, Ex, Eyそれぞれでσx, σx*, σy, σy*が違う(場所が違うから)
  double sig_ex_x, sig_ex_y;
  double sig_ey_x, sig_ey_y;
  double sig_hz_x, sig_hz_y;
  double R = 1.0e-8;
  double M = 2.0;
  const double sig_max = -(M+1.0)*EPSILON_0_S*LIGHT_SPEED_S/2.0/N_PML*log(R);
  for(int i=1; i<SUB_N_PX-1; i++){
    for(int j=1; j<SUB_N_PY-1; j++){      
      int k = subInd(i,j);
      int x = i-1+offsetX;
      int y = j-1+offsetY;
      EPS_EX[k] = models_eps(x+0.5,y, D_Y);
      EPS_EY[k] = models_eps(x,y+0.5, D_X);
      EPS_HZ[k] = 0.5*(models_eps(x+0.5,y+0.5, D_X) + models_eps(x+0.5,y+0.5, D_Y));

      sig_ex_x = sig_max*field_sigmaX(x+0.5,y);
      sig_ex_y = sig_max*field_sigmaY(x+0.5,y);
      sig_ey_x = sig_max*field_sigmaX(x,y+0.5);
      sig_ey_y = sig_max*field_sigmaY(x,y+0.5);
      sig_hz_x = sig_max*field_sigmaX(x+0.5,y+0.5);
      sig_hz_y = sig_max*field_sigmaY(x+0.5,y+0.5);

      double eps = EPSILON_0_S;
      double sig_z = 0;

      C_JX[k] = (2*eps - sig_ex_y)/(2*eps + sig_ex_y);
      C_JXHZ[k] = (2*eps)/(2*eps + sig_ex_y);
      C_DX[k] = (2*eps - sig_z) / (2*eps + sig_z);
      C_DXJX1[k] = (2*eps + sig_ex_x) / (2*eps + sig_z);
      C_DXJX0[k] = (2*eps - sig_ex_x) / (2*eps + sig_z);

      C_JY[k] = (2*eps - sig_z) / (2*eps + sig_z);
      C_JYHZ[k] = (2*eps)/(2*eps + sig_z);
      C_DY[k] = (2*eps - sig_ey_x) / (2*eps + sig_ey_x);
      C_DYJY1[k] = (2*eps + sig_ey_y) / (2*eps + sig_ey_x);
      C_DYJY0[k] = (2*eps - sig_ey_y) / (2*eps + sig_ey_x);

      C_MZ[k] = (2*eps - sig_hz_x) / (2*eps + sig_hz_x );
      C_MZEXEY[k] = (2*eps) / (2*eps + sig_hz_x);
      C_BZ[k] = (2*eps - sig_hz_y) / (2*eps + sig_hz_y);
      C_BZMZ1[k] = (2*eps + sig_z) / (2*eps + sig_hz_y);
      C_BZMZ0[k] = (2*eps - sig_z) / (2*eps + sig_hz_y);
    }
  }  
}

//---------------------メモリの解放--------------------//
static void finish(){
  output();
  freeMemories();
}

static void freeMemories()
{
  if(Ex != NULL){    free(Ex); Ex = NULL;}  
  if(Ey != NULL){    free(Ey); Ey = NULL;}  
  if(Hz != NULL){    free(Hz); Hz = NULL;}

  if(C_JX!= NULL){    free(C_JX);  C_JX = NULL;}
  if(C_JXHZ!= NULL){   free(C_JXHZ); C_JXHZ = NULL;}
  if(C_DX!= NULL){   free(C_DX); C_DX = NULL;}
  if(C_DXJX0 != NULL){   free(C_DXJX0); C_DXJX0 = NULL;}
  if(C_DXJX1 != NULL){   free(C_DXJX1); C_DXJX1 = NULL;}
  
  if(C_JY!= NULL){    free(C_JY);  C_JY = NULL;}
  if(C_JYHZ!= NULL){   free(C_JYHZ); C_JYHZ = NULL;}
  if(C_DY!= NULL){   free(C_DY); C_DY = NULL;}
  if(C_DYJY0 != NULL){   free(C_DYJY0); C_DYJY0 = NULL;}
  if(C_DYJY1 != NULL){   free(C_DYJY1); C_DYJY1 = NULL;}
  
  if(C_MZ!= NULL){    free(C_MZ);  C_MZ = NULL;}
  if(C_MZEXEY!= NULL){   free(C_MZEXEY); C_MZEXEY = NULL;}
  if(C_BZ != NULL){   free(C_BZ); C_BZ = NULL;}
  if(C_BZMZ0 != NULL){   free(C_BZMZ0); C_BZMZ0 = NULL;}
  if(C_BZMZ1 != NULL){   free(C_BZMZ1); C_BZMZ1 = NULL;}
  
  if(EPS_EX != NULL)   free(EPS_EX);
  if(EPS_EY != NULL)   free(EPS_EY);
  if(EPS_HZ != NULL)   free(EPS_HZ);
}

static void output()
{
  if(rank == 0){
    int i,j;
    int x,y;
    int p;
    MPI_Status status;
    double complex *entire = (double complex*)malloc(sizeof(double complex)*N_CELL);
    memset(entire, 0, sizeof(double complex)*N_CELL);
    for(i=1, x=offsetX; i<SUB_N_PX-1; i++, x++)
      for(j=1, y=offsetY; j<SUB_N_PY-1; j++, y++)
        entire[ind(x,y)] = Ey[subInd(i, j)];
    int offsets[2];
    for(p=1; p<nproc; p++){     
      MPI_Recv(offsets, 2, MPI_INT, p, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(Ey, SUB_N_CELL, MPI_C_DOUBLE_COMPLEX, p, 1, MPI_COMM_WORLD, &status);      
      int offX = offsets[0];
      int offY = offsets[1];
      for(i=1, x=offsets[0]; i<SUB_N_PX-1; i++, x++)
        for(j=1, y=offsets[1]; j<SUB_N_PY-1; j++, y++)
          entire[ind(x,y)] = Ey[subInd(i, j)];      
    }    
    field_outputElliptic("mpi_mieTE.txt", entire);
    free(entire);
  }
  else{
    int offsets[2];
    offsets[0] = offsetX; offsets[1] = offsetY;
    MPI_Send(offsets, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);
    MPI_Send(Ey, SUB_N_CELL, MPI_C_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

//---------------------- ntff--------------------//
/*
static inline void ntffCoef(double time, double timeShift, int *m, double *a, double *b, double *ab)
{  
  double t = time + timeShift;
  *m = floor(t + 0.5);
  *a = (0.5 + t - *m);
  *b = 1.0-*a;
  *ab = *a-*b;
}

static void ntff(void)
{
  const double C = LIGHT_SPEED_S;  
  const int cx = SUB_N_PX/2; //a center of field is origin
  const int cy = SUB_N_PY/2;  
  const int offset = 5;		// closed line offset
  const int tp = N_PY - N_PML - offset;	//上から-5
  const int bm = N_PML + offset;	//下から5
  const int rt = SUB_N_PX - N_PML - offset;	//右から-5
  const int lt = N_PML + offset;	//左から5
  const int len = tp-bm+1; //length
  const double RFperC = len*2;  //C = 1/√2,Rf=len*√2 => Rf/C = len*2  
  const double coef = 1.0/(4*M_PI*C);
  
  int m_e, m_h;
  double timeE = field_getTime() - 1;  //t - Δt  todo
  double timeH = field_getTime() - 0.5;  //t - Δt/2  todo
  double a_e, b_e, ab_e, a_h, b_h, ab_h;
  
  for(int ang=0; ang<360; ang++){
    double rad = ang*M_PI/180.0;
    double r1x = cos(rad), r1y = sin(rad);    
    //bottom side
    //normal vector n is (0,-1)
    //Js = n × H = (-Hz, 0, 0)  Ms = E × n = (0, 0, -Ex)
    for(int i=lt; i<rt; i++){
      double r2x = i-cx+0.5, r2y = bm-cy;            //distance between origin to cell
      double timeShift = -(r1x*r2x + r1y*r2y)/C + RFperC;
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_e, &a_e, &b_e, &ab_e);
      
      double complex ex = -EX(i,bm);
      double complex hz = -0.5*( HZ(i,bm) + HZ(i,bm-1) );
          
      Uz[ang][m_e-1] += ex*b_e*coef;
      Wx[ang][m_h-1] += hz*b_h*coef;
      Uz[ang][m_e]   += ex*ab_e*coef;
      Wx[ang][m_h]   += hz*ab_h*coef;
      Uz[ang][m_e+1] -= ex*a_e*coef;
      Wx[ang][m_h+1] -= hz*a_h*coef;      
    }

    //top side
    //normal vector n is (0,1)
    //Js = n × H = (Hz, 0, 0)  Ms = E × n = (0, 0, Ex)
    for(int i=lt; i<rt; i++){
      double r2x = i-cx+0.5, r2y = tp-cy;
      double timeShift = -(r1x*r2x + r1y*r2y)/C + RFperC;      
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

      double complex ex = EX(i,bm);
      double complex hz = 0.5*( HZ(i,tp) + HZ(i,tp-1) );

      Uz[ang][m_e-1] += ex*b_e*coef;
      Wx[ang][m_h-1] += hz*b_h*coef;
      Uz[ang][m_e]   += ex*ab_e*coef;
      Wx[ang][m_h]   += hz*ab_h*coef;
      Uz[ang][m_e+1] -= ex*a_e*coef;
      Wx[ang][m_h+1] -= hz*a_h*coef;            
    }
    
    //left side
    //normal vector n is (-1,0)
    //Js = n × H = (0,Hz,0)    Ms = E × n = (0,0,Ey)
    for(int j=bm; j<tp; j++){
      double r2x = lt-cx, r2y = j-cy+0.5;
      double timeShift = -(r1x*r2x + r1y*r2y)/C + RFperC;
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);
      
      double complex ey = EY(lt,j);
      double complex hz = 0.5*( HZ(lt,j) + HZ(lt-1,j) );
      
      Uz[ang][m_e-1] += ey*b_e*coef;
      Wy[ang][m_h-1] += hz*b_h*coef;
      Uz[ang][m_e]   += ey*ab_e*coef;
      Wy[ang][m_h]   += hz*ab_h*coef;
      Uz[ang][m_e+1] -= ey*a_e*coef;
      Wy[ang][m_h+1] -= hz*a_h*coef;
    }

    //right side
    //normal vector n is (1,0)
    //Js = n × H = (0,-Hz,0)    Ms = E × n = (0,0,-Ey)    
    for(int j=bm; j<tp; j++){
      double r2x = rt-cx, r2y = j-cy+0.5;
      double timeShift = -(r1x*r2x + r1y*r2y)/C + RFperC;
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

      double complex ey = -EY(lt,j);
      double complex hz = -0.5*( HZ(lt,j) + HZ(lt-1,j) );
      Uz[ang][m_e-1] += ey*b_e*coef;
      Wy[ang][m_h-1] += hz*b_h*coef;
      Uz[ang][m_e]   += ey*ab_e*coef;
      Wy[ang][m_h]   += hz*ab_h*coef;
      Uz[ang][m_e+1] -= ey*a_e*coef;
      Wy[ang][m_h+1] -= hz*a_h*coef;      
    }

  }
}
*/
