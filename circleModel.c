#include <math.h>
#include "circleModel.h"
#include "field.h"
static double radius;
static double epsilon;
static double posx;
static double posy;
static double posz;

static double eps(double, double, int, int);

double (*circleModel_EPS(double x, double y, double z, double r))(double, double, double, int , int)
{
  radius = r;
  posx = x;
  posy = y;
  posz = z;
  epsilon = 1.6*1.6*EPSILON_0_S;
  return eps;
}

//col : D_Xモード, row : D_Yモード
static double eps(double x, double y, double z, int col, int row, int dep)
{
  if(x < N_PML || y < N_PML || x > N_X+N_PML || y > N_Y + N_PML)
    return EPSILON_0_S;

  double dx = x-posx;
  double dy = y-posy;
  double dz = z-posz;
  //2乗距離
  double len = dx*dx+dy*dy+dz*dz;

  //中心との距離がr+1セル以上なら,そのセルは完全に媒質の外
  //0.5*√3 = 0.8 以上なら完全に媒質の外
  if(len >= (radius+1)*(radius+1))
    return EPSILON_0_S;

  //中心との距離がr-1セル以下なら,そのセルは完全に媒質の外 
  if(len <= (radius-1)*(radius-1))
    return epsilon;

  //さらに32*32分割し媒質内と媒質外の数を求めepsilonを決定する
  double sum=0;
  double i, j, k;
  for(i=-16+0.5; i<16; i+=1)
    for(j=-16+0.5; j<16; j+=1)
      for(k=-16+0.5; k<16; k+=1)
      {
        if(pow(dx+col*i/32.0, 2.0) + pow(dy+row*j/32.0, 2.0) + pow(dz+dep*k/32.0, 2.0) <= radius*radius)
	sum+=1;
      }  
  
  sum /= 32.0*32.0*32.0;
  return epsilon*sum + EPSILON_0_S*(1-sum);
}
