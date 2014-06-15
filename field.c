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

static double Z_0_S;

//MPI分割したときのフィールドパラメータ
static MPI_Datatype MPI_DCOMPLEX_YZ_COL; //YZ平面の列に対応する(Zが列, Yが行), これをsize(X)繰り返せばXY平面全部とって来れるようにする.

MPI_Datatype MPI_DCOMPLEX_YZ_PLANE; //yz平面全部をまとめた型
MPI_Datatype MPI_DCOMPLEX_XZ_PLANE; //xz平面全部をまとめた型
MPI_Datatype MPI_DCOMPLEX_XY_PLANE; //xy平面全部をまとめた型


static int RANK, NPROC;

//_u : 物理量変換単位, _s:シミュレーション単位
//static const int H_s = 1;
static double time;     //ステップ数
static double ray_coef; //波をゆっくり入れる為の係数;
static double maxTime;
static NTFFInfo ntff_info;

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
NTFFInfo  field_getNTFFInfo(){  return ntff_info;}

double field_getZ_0_S(){ return Z_0_S;}
double  field_getK(){  return waveInfo_s.K_s;}
double  field_getRayCoef(){  return ray_coef;}
double  field_getOmega(){  return waveInfo_s.Omega_s;}
double  field_getLambda(){  return waveInfo_s.Lambda_s;}
double field_getTheta(){  return waveInfo_s.Theta_deg;}
double field_getPhi(){  return waveInfo_s.Phi_deg;}
double  field_getTime(){  return time;}
double  field_getMaxTime(){  return maxTime;}

int field_index(int i, int j, int k){
    return i*fieldInfo_s.N_PYZ + j*fieldInfo_s.N_PZ+k;
}
int field_left(int ind){
  return ind - fieldInfo_s.N_PYZ;
}
int field_right(int ind){
  return ind + fieldInfo_s.N_PYZ;
}
int field_top(int ind){
  return ind + fieldInfo_s.N_PZ;
}
int field_bottom(int ind){
  return ind - fieldInfo_s.N_PZ;
}
int field_front(int ind){
  return ind+1;
}
int field_back(int ind){
  return ind-1;
}

int field_subIndex(int i, int j, int k){
    return i*subFieldInfo_s.SUB_N_PYZ + j*subFieldInfo_s.SUB_N_PZ + k;
}
int field_subLeft(int ind){
  return ind - subFieldInfo_s.SUB_N_PYZ;
}
int field_subRight(int ind){
  return ind + subFieldInfo_s.SUB_N_PYZ;
}
int field_subTop(int ind){
  return ind + subFieldInfo_s.SUB_N_PZ;
}
int field_subBottom(int ind){
  return ind - subFieldInfo_s.SUB_N_PZ;
}
int field_subFront(int ind){
  return ind+1;
}
int field_subBack(int ind){
  return ind-1;
}
int field_subToOneX(int i){
  return i-1+subFieldInfo_s.OFFSET_X;
}
int field_subToOneY(int j){
  return j-1+subFieldInfo_s.OFFSET_Y;
}
int field_subToOneZ(int k){
  return k-1+subFieldInfo_s.OFFSET_Z;
}

void field_init(FieldInfo field_info)
{
  Z_0_S = sqrt(MU_0_S/EPSILON_0_S);
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
  fieldInfo_s.N_PYZ = fieldInfo_s.N_PY*fieldInfo_s.N_PZ;

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
  ntff_info.back   = fieldInfo_s.N_PML + 5;
  ntff_info.front  = fieldInfo_s.N_PZ - fieldInfo_s.N_PML - 5;


  // todo
//  NOT_DONE("you have to check RFperC in 3D\n");
  printf("you have to check RFperC in 3D\n");
  double len = (ntff_info.top - ntff_info.bottom)/2;
  ntff_info.RFperC = len*2;
  ntff_info.arraySize = maxTime + 2*ntff_info.RFperC;
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

//pml用の係数のひな形 Δt = 1
//ep_mu εかμ(Eの係数->ε, Hの係数-> μ
//sig  σ
double field_pmlCoef(double ep_mu, double sig)
{
  return (1.0 - sig/ep_mu)/(1.0+sig/ep_mu);
}
double field_pmlCoef2(double ep_mu, double sig)
{
  return 1.0/(ep_mu + sig); // 1.0/{ep_mu(1.0 + sig/ep_mu)}と同じ
}


//------------------light method----------------------//
//点光源を返す
 double complex field_pointLight(void)
{
  return ray_coef * (cos(waveInfo_s.Omega_s*time) + sin(waveInfo_s.Omega_s*time)*I);
}

//ガウシアンパルス
// gapX, gapY : Ex-z, Hx-zは格子点からずれた位置に配置され散る為,格子点からのずれを送る必要がある.
// 分割領域では使えない->gapをさらに領域のオフセットだけずらせば使えそう
// UPML専用
void field_scatteredPulse(dcomplex *p, double *eps, double gapX, double gapY, double gapZ)
{
/*
  double time = field_getTime();
  double w_s  = field_getOmega();
  double rad = field_getWaveAngle()*M_PI/180.0;	//ラジアン変換  

  double cos_per_c = cos(rad)/C_0_S, sin_per_c = sin(rad)/C_0_S;
  const double beam_width = 50; //パルスの幅

  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  
  //waveAngleにより, t0の値を変えないとちょうどいいところにピークが来なため,それを計算.
  const double center_peak = (fInfo_s.N_PX/2.0+gapX)*cos_per_c+(fInfo_s.N_PY/2+gapY)*sin_per_c; //中心にピークがくる時間
  const double t0 = -center_peak + 100; //常に100ステップの時に,領域の中心にピークが来るようにする.
*/
  /*
  for(int i=1; i<fInfo_s.N_PX-1; i++) {
    for(int j=1; j<fInfo_s.N_PY-1; j++) {
      int k = field_index(i,j);
      const double r = (i+gapX)*cos_per_c+(j+gapY)*sin_per_c-(time-t0); // (x*cos+y*sin)/C - (time-t0)
      const double gaussian_coef = exp( -pow(r/beam_width, 2 ) );
      p[k] += gaussian_coef*(EPSILON_0_S/eps[k] - 1)*cexp(I*r*w_s);     //p[k] -= かも(岡田さんのメール参照)
    }
    }*/
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

  //z方向 todo 逆かもしれない=>たぶん合ってる
  MPI_Cart_shift(grid_comm, 2, 1, &subFieldInfo_s.BkRank, &subFieldInfo_s.FtRank); 

  //プロセス座標において自分がどの位置に居るのか求める(何行何列に居るか)
  int coordinates[3];
  MPI_Comm_rank(grid_comm, &RANK);
  MPI_Cart_coords(grid_comm, RANK, dim, coordinates);

  
  //ワールドタブにおける自分のランクを求める.
  MPI_Comm_rank(MPI_COMM_WORLD, &subFieldInfo_s.Rank);
  MPI_Comm_size(MPI_COMM_WORLD, &subFieldInfo_s.Nproc);
  
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
  // -1はのりしろのがあるため, sub領域の0番目は,隣の領域(左,下,後ろ)の領域をさすため
  subFieldInfo_s.OFFSET_X  = coordinates[0] * subFieldInfo_s.SUB_N_X - 1;
  subFieldInfo_s.OFFSET_Y  = coordinates[1] * subFieldInfo_s.SUB_N_Y - 1;
  subFieldInfo_s.OFFSET_Z  = coordinates[2] * subFieldInfo_s.SUB_N_Z - 1;

  //YZ平面の同期をとるための型を定義
  //SUB_N_Z個の連続したデータ, SUB_N_PZ跳び(次のデータまでpz-p = 2個の隙間がある)に, SUB_N_Y行 取ってくる事になる
  MPI_Type_vector(subFieldInfo_s.SUB_N_Y, subFieldInfo_s.SUB_N_Z,
                  subFieldInfo_s.SUB_N_PZ, MPI_C_DOUBLE_COMPLEX, &MPI_DCOMPLEX_YZ_PLANE);
  MPI_Type_commit(&MPI_DCOMPLEX_YZ_PLANE);

  //XZ平面の同期をとるための型を定義
  MPI_Type_vector(subFieldInfo_s.SUB_N_X, subFieldInfo_s.SUB_N_Z,
                  subFieldInfo_s.SUB_N_PYZ, MPI_C_DOUBLE_COMPLEX, &MPI_DCOMPLEX_XZ_PLANE );
  MPI_Type_commit(&MPI_DCOMPLEX_XZ_PLANE);


  //XY平面は連続する領域が無い(隙間が2種類ある)ので, のりしろも含めた全領域を同期する必要がある.
  MPI_Type_vector(subFieldInfo_s.SUB_N_PX*subFieldInfo_s.SUB_N_PY, 1, subFieldInfo_s.SUB_N_PZ, MPI_C_DOUBLE_COMPLEX, &MPI_DCOMPLEX_XY_PLANE);
  MPI_Type_commit(&MPI_DCOMPLEX_XY_PLANE);
  
  /*
  //XYは隙間が2段階あるので,一気に登録は出来ない.  
  MPI_Type_vector(subFieldInfo_s.SUB_N_Y, 1, subFieldInfo_s.SUB_N_PZ, MPI_C_DOUBLE_COMPLEX, &MPI_DCOMPLEX_YZ_COL);
  MPI_Type_commit(&MPI_DCOMPLEX_YZ_COL);

  MPI_Type_vector(subFieldInfo_s.SUB_N_X, 1, subFieldInfo_s.SUB_N_PY, MPI_DCOMPLEX_YZ_COL, &MPI_DCOMPLEX_XY_PLANE);
  MPI_Type_commit(&MPI_DCOMPLEX_XY_PLANE);  */
  
  printf("field.c rank=%d, offset(%d, %d, %d)\n", subFieldInfo_s.Rank, subFieldInfo_s.OFFSET_X, subFieldInfo_s.OFFSET_Y, subFieldInfo_s.OFFSET_Z);

}

//plane => 0 : XY平面(z=dep/2)
//1 : YZ平面(x=wid/2)
//2 : XZ平面(y=hei/2)
void field_outputElliptic(const char *fileName, dcomplex* data, int plane)
{
  printf("output start\n");
  //file open
  FILE *fp = openFile(fileName);
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  WaveInfo_S wInfo_s = field_getWaveInfo_S();


  double r = 1.5*wInfo_s.Lambda_s;
  double ToRad = M_PI/180.0;

  int ox = fInfo_s.N_PX/2;
  int oy = fInfo_s.N_PY/2;
  int oz = fInfo_s.N_PZ/2;
  for(int ang=360; ang; ang--)
  {
    double rad = ang*ToRad;
    dcomplex phi;
    if(plane == 0)
    {
      double x =  ox + r*cos(rad);
      double y =  oy + r*sin(rad);
      int index = field_index((int)x, (int)y, oz);
      phi = cbilinear(data, x, y, index, fInfo_s.N_PYZ, fInfo_s.N_PZ);
    } else if( plane == 1) {
      double z =  oz + r*cos(rad);
      double y =  oy + r*sin(rad);
      int index = field_index( ox, (int)y, (int)z);
      phi = cbilinear(data, z, y, index, 1, fInfo_s.N_PZ);
    }else {
      double x =  ox + r*cos(rad);
      double z =  oz + r*sin(rad);
      int index = field_index( (int)x, oy, (int)z);
      phi = cbilinear(data, x, z, index, fInfo_s.N_PYZ, 1);
    }
    fprintf(fp, "%d, %.18lf \n", 360-ang, cnorm(phi));
  }
  
  fclose(fp);
  printf("output to %s end\n", fileName);
}

void field_outputAllDataComplex(const char *fileName,dcomplex* data)
{
  FILE *fp = openFile(fileName);
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  for(int i=0; i<fInfo_s.N_CELL; i++)
  {
    fprintf(fp,  "%.18lf, %.18lf \n", creal(data[i]), cimag(data[i]));
  }
  fclose(fp);
}
