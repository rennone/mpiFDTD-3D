#include <math.h>
#include "circleModel.h"
#include "field.h"
#include "function.h"

#define ST_RADIUS_NM 100
#define EN_RADIUS_NM 100
#define DELTA_RADIUS_NM 50

static int radius_nm = ST_RADIUS_NM;
static double radius;
static double epsilon;

static double posx;
static double posy;
static double posz;

static double eps(double, double, double, int, int, int);

double (*circleModel_EPS())(double, double, double, int , int, int)
{
  return eps;
}

void circleModel_init()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  posx = fInfo_s.N_PX/2;
  posy = fInfo_s.N_PY/2;
  posz = fInfo_s.N_PZ/2;

  radius = field_toCellUnit(radius_nm);

  double n = 1.6;
  epsilon = n*n*EPSILON_0_S;
}

//col : D_Xモード, row : D_Yモード
static double eps(double x, double y, double z, int col, int row, int dep)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  //PML領域には物質を置かない
  if(x < fInfo_s.N_PML || x > fInfo_s.N_X + fInfo_s.N_PML ||
     y < fInfo_s.N_PML || y > fInfo_s.N_Y + fInfo_s.N_PML ||
     z < fInfo_s.N_PML || z > fInfo_s.N_Z + fInfo_s.N_PML )
    return EPSILON_0_S;

  double dx = x-posx;
  double dy = y-posy;
  double dz = z-posz;
  
  //2乗距離
  double lenSqr = dx*dx+dy*dy+dz*dz;

  //中心との距離がr+1.2セル以上なら,そのセルは完全に媒質の外
  //0.5*√3 = 0.8 以上なら完全に媒質の外
  if(lenSqr >= (radius+1.0)*(radius+1.0))
    return EPSILON_0_S;

  //中心との距離がr-1セル以下なら,そのセルは完全に媒質の外 
  if(lenSqr <= (radius-1.0)*(radius-1.0))
    return epsilon;

  
  //さらに分割した領域で媒質内と媒質外の数を求めepsilonを細かに決定する
  double split = 10; //分割数
  double sum=0;
//  NOT_DONE("i dont need its loop if flag(col row dep) is false\n");
  for(int i=-split/2+0.5; i<split/2; i+=1)
    for(int j=-split/2+0.5; j<split/2; j+=1)
      for(int k=-split/2+0.5; k<split/2; k+=1)
      {
        if(pow(dx+col*i/split, 2.0) + pow(dy+row*j/split, 2.0) + pow(dz+dep*k/split, 2.0) <= radius*radius)
	sum+=1;
      }
  
  sum /= split*split*split;
  return epsilon*sum + EPSILON_0_S*(1-sum);
}

bool circleModel_isFinish()
{
  //50nm増やす
  radius += field_toCellUnit(50);
  
  return radius > field_toCellUnit(EN_RADIUS_NM);
}

void circleModel_moveDirectory(void)
{
  char buf[512];
  sprintf(buf, "%d_nm", (int)field_toPhysicalUnit(radius));
  makeDirectory(buf);
  moveDirectory(buf); 
}
