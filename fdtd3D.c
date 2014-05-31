#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "field.h"
#include "models.h"
#include "fdtd3D.h"
#include "myComplex.h"
#include "function.h"

/* about MPI  */
static int rank;      //MPIのランク
static int nproc;     //全プロセス数
static MPI_Datatype X_DIRECTION_DOUBLE_COMPLEX;
static MPI_Datatype DCOMPLEX_XY; //XY平面
static MPI_Datatype DCOMPLEX_YZ; //YZ平面
static MPI_Datatype DCOMPLEX_XZ; //XZ平面

//岡田さんの論文と同じ空間配置にしてみる
//h = 1, Δt = 1で計算
//系は右手系
//x(left-, right+)
//y(bottom-, top+)
//z(back-, front+)

//Hx(i+0.5,j    ,k    ,t+0.5) -> Hx[i,j,k]
//Hy(i    ,j+0.5,k    ,t+0.5) -> Hy[i,j,k]
//Hz(i    ,j    ,k+0.5,t+0.5) -> Hz[i,j,k]

//Ex(i    ,j+0.5,k+0.5,t    ) -> Ex[i,j,k]
//Ey(i+0.5,     ,k+0.5,t    ) -> Ey[i,j,k]
//Ez(i+0.5,j+0.5,k    ,t    ) -> Ez[i,j,k]


static dcomplex *Ex = NULL;
static dcomplex *Ey = NULL;
static dcomplex *Ez = NULL;
static dcomplex *Hx = NULL;
static dcomplex *Hy = NULL;
static dcomplex *Hz = NULL;

static double *C_EX_HZHY = NULL, *C_EY_HXHZ = NULL, *C_EZ_HYHZ = NULL;
static double *C_HX_EYEZ = NULL, *C_HY_EZEX = NULL, *C_HZ_EXEY= NULL;

static double *EPS_EX=NULL, *EPS_EY=NULL, *EPS_EZ=NULL;
static double *EPS_HX=NULL, *EPS_HY=NULL, *EPS_HZ=NULL;

static void update(void);
static void finish(void);
static void init(void);
static void reset(void);

static void calcE(void);
static void calcH(void);
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
  calcH();    
  calcE();
  
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

  C_EX_HZHY = newDouble(sInfo.SUB_N_CELL);
  C_EY_HXHZ = newDouble(sInfo.SUB_N_CELL);
  C_EZ_HYHZ = newDouble(sInfo.SUB_N_CELL);
  
  C_HX_EYEZ = newDouble(sInfo.SUB_N_CELL);
  C_HY_EZEX = newDouble(sInfo.SUB_N_CELL);
  C_HZ_EXEY = newDouble(sInfo.SUB_N_CELL);

  EPS_EX = newDouble(sInfo.SUB_N_CELL);
  EPS_EY = newDouble(sInfo.SUB_N_CELL);
  EPS_EZ = newDouble(sInfo.SUB_N_CELL);
  EPS_HY = newDouble(sInfo.SUB_N_CELL);
  EPS_HX = newDouble(sInfo.SUB_N_CELL);
  EPS_HZ = newDouble(sInfo.SUB_N_CELL);

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
  /*
  //Hz, Ex, Eyそれぞれでσx, σx*, σy, σy*が違う(場所が違うから)
  double sig_ex_x, sig_ex_y, sig_ex_z;
  double sig_ey_x, sig_ey_y, sig_ey_z;
  double sig_ez_x, sig_ez_y, sig_ez_z;
  double sig_hx_x, sig_hx_y, sig_hx_z;
  double sig_hy_x, sig_hy_y, sig_hy_z;
  double sig_hz_x, sig_hz_y, sig_hz_z;
  */
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
        EPS_EX[w] = models_eps(x    ,y+0.5,z+0.5,D_Y); //todo 
        EPS_EY[w] = models_eps(x    ,y+0.5,z+0.5,D_X);
        EPS_EZ[w] = models_eps(x+0.5,y+0.5,z    ,D_XY);

        EPS_HX[w] = models_eps(x+0.5,y    ,z    ,D_Y);
        EPS_HY[w] = models_eps(x    ,y+0.5,z    ,D_X);
        EPS_HZ[w] =(models_eps(x,y,z+0.5, D_X)+models_eps(x,y,z+0.5, D_Y))/2.0;

        C_EX_HZHY[w] = 1.0/EPS_EX[w];
        C_EY_HXHZ[w] = 1.0/EPS_EY[w];
        C_EZ_HYHZ[w] = 1.0/EPS_EZ[w];

        C_HX_EYEZ[w] = 1.0/MU_0_S;
        C_HY_EZEX[w] = 1.0/MU_0_S;
        C_HZ_EXEY[w] = 1.0/MU_0_S;
    }
}

#define IndexCheck(w, sInfo, msg) if(w>=sInfo.SUB_N_CELL) { printf(msg); exit(2);}
//Private : 
static void calcE()
{
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();

  int c = field_subIndex(subInfo_s.SUB_N_PX/2,subInfo_s.SUB_N_PY/2,subInfo_s.SUB_N_PZ/2);
  for(int i=1; i<subInfo_s.SUB_N_PX-1; i++){
    for(int j=1; j<subInfo_s.SUB_N_PY-1; j++){
      for(int k=1; k<subInfo_s.SUB_N_PZ-1; k++){
        int w = field_subIndex(i,j,k);
        int w_top = field_subIndex(i  ,j+1, k  );
        int w_frt = field_subIndex(i  ,j  , k+1);
        int w_rht = field_subIndex(i+1,j  , k  );

        Ex[w] = Ex[w] + /*C_EX_HZHY[w]* */( Hz[w_top] - Hz[w] - Hy[w_frt] + Hy[w]);
        if( w==c)
        {
          //printf("(%lf, %lf)\n", creal(-Hy[w_frt]), creal(+Hy[w]));          
        }
        Ey[w] = Ey[w] + /*C_EY_HXHZ[w]* */( Hx[w_frt] - Hx[w] - Hz[w_rht] + Hz[w]);
        Ez[w] = Ez[w] + /*C_EZ_HYHZ[w]* */( Hy[w_rht] - Hy[w] - Hx[w_top] + Hx[w]);
      }
    }
  }  
}

static void calcH()
{
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  int c = field_subIndex(subInfo_s.SUB_N_PX/2,subInfo_s.SUB_N_PY/2,subInfo_s.SUB_N_PZ/2);
  for(int i=1; i<subInfo_s.SUB_N_PX-1; i++){
    for(int j=1; j<subInfo_s.SUB_N_PY-1; j++){
      for(int k=1; k<subInfo_s.SUB_N_PZ-1; k++){
        int w = field_subIndex(i,j,k);
        int w_bck = field_subIndex(i  ,j  ,k-1);
        int w_lft = field_subIndex(i-1,j  ,k  );
        int w_btm = field_subIndex(i  ,j-1,k  );

        Hx[w] = Hx[w] + C_HX_EYEZ[w]* ( Ey[w]-Ey[w_bck] -Ez[w]+Ez[w_btm]);
        Hy[w] = Hy[w] + C_HY_EZEX[w]* ( Ez[w]-Ez[w_lft] -Ex[w]+Ex[w_bck]);

        if( w==c )
        {
          printf("(%lf, %lf)\n", creal(Ez[w]), creal(+Ez[w_bck]));          
        }
        Hz[w] = Hz[w] + C_HZ_EXEY[w]* ( Ex[w]-Ex[w_btm] -Ey[w]*Ey[w_lft]);          
      }
    }
  }
}

