#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "field.h"
#include "models.h"
#include "function.h"

FieldInfo fieldInfo;
FieldInfo_S fieldInfo_s;
SubFieldInfo_S subFieldInfo_s;
WaveInfo_S waveInfo_s;

//MPI分割したときのフィールドパラメータ
static int SUB_N_X;
static int SUB_N_Y;
static int SUB_N_Z;

static int SUB_N_PX;
static int SUB_N_PY;
static int SUB_N_PZ;
static int SUB_N_CELL;
static int RANK, NPROC;
static int LTRANK, RTRANK, TPRANK,BMRANK;
static int OFFSET_X, OFFSET_Y, OFFSET_Z;
static MPI_Datatype MPI_DOUBLE_COMPLEX_COLUM;

//_u : 物理量変換単位, _s:シミュレーション単位
static const int H_s = 1;
static double H_u;         //1セルの大きさ(nm)
static double time;     //ステップ数
static double ray_coef; //波をゆっくり入れる為の係数;
static double waveAngle;//入射角
static double Lambda_s; //波長 
static double k_s;      //波数 
static double w_s;      //角周波数
static double T_s;      //周期 
static void (*defaultWave)(double complex* p, double* eps);
static double maxTime;

static NTFFInfo ntff_info;

static void field_setNTFF(int);
static void mpiSplit(void);

//:public------------------------------------//
 double field_toCellUnit(const double phisycalUnit){
  return phisycalUnit/fieldInfo.h_u_nm;   //セル単位に変換 
}

 double field_toPhisycalUnit(const double cellUnit){
  return cellUnit*fieldInfo.h_u_nm;    //物理単位(nm)に変換
}

int field_getOffsetX(){  return subFieldInfo_s.OFFSET_X;}
int field_getOffsetY(){  return subFieldInfo_s.OFFSET_Y;}
int field_getOffsetZ(){  return subFieldInfo_s.OFFSET_Z;}
int field_getSubNx(){  return subFieldInfo_s.SUB_N_X;}
int field_getSubNy(){  return subFieldInfo_s.SUB_N_Y;}
int field_getSubNz(){  return subFieldInfo_s.SUB_N_Z;}
int field_getSubNpx(){  return subFieldInfo_s.SUB_N_PX;}
int field_getSubNpy(){  return subFieldInfo_s.SUB_N_PY;}
int field_getSubNpz(){  return subFieldInfo_s.SUB_N_PZ;}
int field_getSubNcell(){  return subFieldInfo_s.SUB_N_CELL;}
int field_getSubNpyz(){  return subFieldInfo_s.SUB_N_PYZ;}

WaveInfo_S field_getWaveInfo_S()         { return waveInfo_s;}
SubFieldInfo_S field_getSubFieldInfo_S() { return subFieldInfo_s;}
FieldInfo_S field_getFieldInfo_S()       { return fieldInfo_s;}
FieldInfo field_getFieldInfo()   { return fieldInfo;}

int field_subIndex(int i, int j, int k){
    return i*subFieldInfo_s.SUB_N_PYZ + j*subFieldInfo_s.SUB_N_PZ + k;
}

int field_index(int i, int j, int k){
    return i*fieldInfo_s.N_PYZ + j*fieldInfo_s.N_PZ+k;
}

int field_subLeft(int ind)
{
  return ind - subFieldInfo_s.SUB_N_PYZ;
}
int field_subRight(int ind)
{
  return ind + subFieldInfo_s.SUB_N_PYZ;
}
int field_subTop(int ind)
{
  return ind + subFieldInfo_s.SUB_N_PZ;
}
int field_subBottom(int ind)
{
  return ind - subFieldInfo_s.SUB_N_PZ;
}
int field_subFront(int ind)
{
  return ind+1;
}
int field_subBack(int ind)
{
  return ind-1;
}

void field_init(FieldInfo field_info)
{
  //フィールド情報の保存(最初にしないとtoCellUnit, PhisicalUnitが使えない.
  fieldInfo = field_info;
  
  //領域のシミュレータパラメータを計算
  fieldInfo_s.N_PX  = field_toCellUnit(fieldInfo.width_nm);
  fieldInfo_s.N_PY  = field_toCellUnit(fieldInfo.height_nm);
  fieldInfo_s.N_PZ  = field_toCellUnit(fieldInfo.depth_nm);
  
  fieldInfo_s.N_PML = fieldInfo.pml;
  fieldInfo_s.N_X   = fieldInfo_s.N_PX - 2*fieldInfo_s.N_PML;
  fieldInfo_s.N_Y   = fieldInfo_s.N_PY - 2*fieldInfo_s.N_PML;
  fieldInfo_s.N_Z   = fieldInfo_s.N_PZ - 2*fieldInfo_s.N_PML;
  fieldInfo_s.N_CELL= fieldInfo_s.N_PY*fieldInfo_s.N_PX*fieldInfo_s.N_PZ;
  fieldInfo_s.N_PYZ = fieldInfo_s.N_Y*fieldInfo_s.N_Z;

  //入射波パラメータの計算
  waveInfo_s.Lambda_s  = field_toCellUnit(fieldInfo.lambda_nm);
  waveInfo_s.T_s       = waveInfo_s.Lambda_s/C_0_S;
  waveInfo_s.K_s       = 2*M_PI/waveInfo_s.Lambda_s;
  waveInfo_s.Omega_s   = C_0_S*waveInfo_s.K_s;
  waveInfo_s.Theta_deg = fieldInfo.theta_deg;
  waveInfo_s.Phi_deg   = fieldInfo.phi_deg;
  
  mpiSplit();  //小領域のパラメータを計算

  time = 0;
  maxTime = fieldInfo.stepNum;
  ray_coef  = 0;
  
  /* NTFF設定 */
  ntff_info.cx     = fieldInfo_s.N_PX/2;
  ntff_info.cy     = fieldInfo_s.N_PY/2;
  ntff_info.cz     = fieldInfo_s.N_PY/2;
  ntff_info.top    = fieldInfo_s.N_PY - fieldInfo_s.N_PML - 5;
  ntff_info.bottom = fieldInfo_s.N_PML + 5;
  ntff_info.left   = fieldInfo_s.N_PML + 5;
  ntff_info.right  = fieldInfo_s.N_PX - fieldInfo_s.N_PML - 5;
  ntff_info.front  = fieldInfo_s.N_PML + 5;
  ntff_info.back   = fieldInfo_s.N_PZ - fieldInfo_s.N_PML - 5;

  // todo
//  NOT_DONE("you have to check RFperC in 3D\n");
  double len = (ntff_info.top - ntff_info.bottom)/2;
  ntff_info.RFperC = len*2;
  ntff_info.arraySize = maxTime + 2*ntff_info.RFperC;
}

//-------------------getter-------------------//
double  field_getK()
{
  return waveInfo_s.K_s;
}

double  field_getRayCoef()
{
  return ray_coef;
}

double  field_getOmega()
{
  return waveInfo_s.Omega_s;
}

double  field_getLambda()
{
  return waveInfo_s.Lambda_s;
}

double  field_getWaveAngle()
{
  return waveAngle;
}

double  field_getTime()
{
  return time;
}

double  field_getMaxTime()
{
  return maxTime;
}

NTFFInfo  field_getNTFFInfo()
{
  return ntff_info;
}

//----------------------------------------//
 double field_sigmaX(const double x, const double __y, const double __z)
{
  const int M = 2;
  if(x<fieldInfo_s.N_PML)
    return pow(1.0*(fieldInfo_s.N_PML-x)/fieldInfo_s.N_PML, M);
  
  else if(fieldInfo_s.N_PML <= x && x < (fieldInfo_s.N_X+fieldInfo_s.N_PML))    
    return 0;
  
  else
    return pow(1.0*(x - (fieldInfo_s.N_PX-fieldInfo_s.N_PML-1))/fieldInfo_s.N_PML, M);
}

 double field_sigmaY(const double __x, const double y, const double __z)
{
  const int M = 2;
  if(y<fieldInfo_s.N_PML)
    return pow(1.0*(fieldInfo_s.N_PML - y)/fieldInfo_s.N_PML,M);
  
  else if(y>=fieldInfo_s.N_PML && y<(fieldInfo_s.N_Y+fieldInfo_s.N_PML))
    return 0.0;

  else
    return pow(1.0*(y - (fieldInfo_s.N_PY-fieldInfo_s.N_PML-1))/fieldInfo_s.N_PML,M);
}

double field_sigmaZ(const double __x, const double __y, const double z)
{
  const int M = 2;
  if(z<fieldInfo_s.N_PML)
    return pow(1.0*(fieldInfo_s.N_PML-z)/fieldInfo_s.N_PML, M);
  
  else if( z>=fieldInfo_s.N_PML && z<(fieldInfo_s.N_Z+fieldInfo_s.N_PML))
    return 0.0;
  
  else
    return pow(1.0*(z - (fieldInfo_s.N_PZ-fieldInfo_s.N_PML-1))/fieldInfo_s.N_PML,M);
}

//------------------getter-------------------------//

//------------------light method----------------------//
//点光源を返す
 double complex field_pointLight(void)
{
  return ray_coef * (cos(w_s*time) + sin(w_s*time)*I);
}

//------------------light method----------------------//
 void field_nextStep(void){
  time+=1.0;
  ray_coef = 1.0*(1.0 - exp(-0.0001*time*time));
}

bool field_isFinish(void){
  return time >= maxTime;
}

static void mpiSplit(void)
{
  MPI_Comm_size(MPI_COMM_WORLD, &NPROC);
  int dim = 3;          //number of dimension
  int procs[3] = {0,0,0};         //[0]: x方向の分割数, [1]:y方向の分割数
  int period[3] = {0,0,0};//境界条件, 固定境界
  MPI_Comm grid_comm;
  int reorder = 1;   //re-distribute rank flag

  MPI_Dims_create(NPROC, dim, procs);
  MPI_Cart_create(MPI_COMM_WORLD, dim, procs, period, reorder, &grid_comm);
  MPI_Cart_shift(grid_comm, 0, 1, &subFieldInfo_s.LtRank, &subFieldInfo_s.RtRank); //x方向の隣
  MPI_Cart_shift(grid_comm, 1, 1, &subFieldInfo_s.BmRank, &subFieldInfo_s.TpRank); //y方向

//  NOT_DONE("field.c check the order of FtRank and BkRank ");
  MPI_Cart_shift(grid_comm, 2, 1, &subFieldInfo_s.BkRank, &subFieldInfo_s.FtRank); //z方向 todo 逆かもしれない

  //プロセス座標において自分がどの位置に居るのか求める(何行何列に居るか)
  int coordinates[3];
  MPI_Comm_rank(grid_comm, &RANK);
  MPI_Cart_coords(grid_comm, RANK, dim, coordinates);
  
/*これだと, 1個のデータをSUB_N_PY跳び(次のデータまでSUB_N_PY-1個隙間がある),SUB_N_X行ぶん取ってくる事になる */
  MPI_Type_vector(SUB_N_X, 1, SUB_N_PY, MPI_C_DOUBLE_COMPLEX, &MPI_DOUBLE_COMPLEX_COLUM); 
  MPI_Type_commit(&MPI_DOUBLE_COMPLEX_COLUM);


  //小領域のパラメータ
  subFieldInfo_s.SUB_N_X    = fieldInfo_s.N_PX / procs[0];
  subFieldInfo_s.SUB_N_Y    = fieldInfo_s.N_PY / procs[1];
  subFieldInfo_s.SUB_N_Z    = fieldInfo_s.N_PZ / procs[2];
  
  subFieldInfo_s.SUB_N_PX   = subFieldInfo_s.SUB_N_X + 2; //のりしろの分2大きい
  subFieldInfo_s.SUB_N_PY   = subFieldInfo_s.SUB_N_Y + 2; //のりしろの分2大きい
  subFieldInfo_s.SUB_N_PZ   = subFieldInfo_s.SUB_N_Z + 2; //のりしろの分2大きい
  
  subFieldInfo_s.SUB_N_CELL = subFieldInfo_s.SUB_N_PX*subFieldInfo_s.SUB_N_PY*subFieldInfo_s.SUB_N_PZ;
  subFieldInfo_s.SUB_N_PYZ  = subFieldInfo_s.SUB_N_PY*subFieldInfo_s.SUB_N_PZ;

  //ランクのインデックスではなく, セル単位のオフセットなのでSUB_N_Xずれる
  subFieldInfo_s.OFFSET_X  = coordinates[0] * subFieldInfo_s.SUB_N_X;
  subFieldInfo_s.OFFSET_Y  = coordinates[1] * subFieldInfo_s.SUB_N_Y;
  subFieldInfo_s.OFFSET_Z  = coordinates[2] * subFieldInfo_s.SUB_N_Z;
}
