#include "multilayerModel.h"
#include "field.h"
#include "function.h"
#include <math.h>

double width_s[2];     //幅
double thickness_s[2]; //厚さ
double n[2];            //屈折率
double ep[2];           //誘電率 = n*n*ep0
int layerNum;          //枚数
bool asymmetry;      //左右比対称

//col : D_Xモード row : D_Yモード
//x,yを中心に, 計算領域のセルと同じ大きさの領域を調べる
static double eps(double x, double y, double z, int col, int row, int dep)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();

  double width  = max(width_s[0], width_s[1]);
  double thick  = thickness_s[0] + thickness_s[1];
  double height = thick*layerNum;

  //領域の中心から, 下にheight/2ずれた位置がレイヤの下部
  int oy = fInfo_s.N_PY/2 - height/2;
  int ox = fInfo_s.N_PX/2;
  int oz = fInfo_s.N_PZ/2;
  
  double _x = x-ox;	//ox,oyを座標の原点に
  double _y = y-oy;
  double _z = z-oz;

  //上下左右に飛び出ていないか確認(細分化したセルがあるため, 0.5の余白をとっている)
  if( fabs(_x) > (width/2+0.5) ||
      fabs(_z) > (width/2+0.5) ||
      _y < -0.5 || _y > height+0.5  )  
    return EPSILON_0_S;

  double s[2]={0,0}; //n1,n2それぞれの分割セルの数が入る
  double split = 10; //分割数
  for(int i=-split/2+0.5; i<split/2; i+=1)
    for(int j=-split/2+0.5; j<split/2; j+=1)
      for(int k=-split/2+0.5; k<split/2; k+=1)
      {
        double sx = _x + col*i/split; //細分化したセルの位置
        double sy = _y + row*j/split;
        double sz = _z + dep*k/split;

        //上下に飛び出ていないか確認
        if(sy < 0 || sy >= height)
          continue;
      
        //thickで割ったあまり(double型なのでこんなやり方をしている)
        double modY = sy - floor(sy/thick)*thick;

        //境界上のときは両方の平均になる(普通は無い).
        if(modY == thickness_s[0]) {
          s[0] += 0.5*(fabs(sx) < width_s[0]/2);
          s[1] += 0.5*(fabs(sx) < width_s[1]/2);
          continue;
        }

        int k = (modY > thickness_s[0]); //どっちの屈折率にいるか調べる

        if (sx < 0 && asymmetry)
          k = 1-k;		//左右で反転, 互い違いでなかったら反転しない

        //正方形の内側ならその媒質内にある
        if( abs(sx) < width_s[k]/2 && abs(sz) < width_s[k] )
          s[k] +=1;
      }

  s[0] /= split*split*split;
  s[1] /= split*split*split;
  return EPSILON_0_S*(1-s[0]-s[1]) + ep[0]*s[0] + ep[1]*s[1];
}

double ( *multilayerModel_EPS(void))(double, double, double, int, int, int)
{
  width_s[0]     = field_toCellUnit(300);
  width_s[1]     = field_toCellUnit(300);
  thickness_s[0] = field_toCellUnit(90);
  thickness_s[1] = field_toCellUnit(90);
  layerNum = 8;
  n[0] = 1.56;
  n[1] = 1.0;
  ep[0] = n[0]*n[0]*EPSILON_0_S;
  ep[1] = n[1]*n[1]*EPSILON_0_S;

  asymmetry = true;

  return eps;
}

bool multilayerModel_isFinish(void)
{
  return true;
}
