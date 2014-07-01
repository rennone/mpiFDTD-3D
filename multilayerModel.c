#include "multilayerModel.h"
#include "field.h"
#include "function.h"
#include <stdlib.h>
#include <math.h>

static double depth_s[2];     //奥行き z
static double width_s[2];     //幅 x
static double thickness_s[2]; //厚さ y
static double n[2];            //屈折率
static double ep[2];           //誘電率 = n*n*ep0
static int layerNum;          //枚数
static bool asymmetry;      //左右比対称

static int start_depth_nm =  300;
static int end_depth_nm   = 1000;
static int delta_depth_nm = 100;

static int start_thickness_nm = 80;
static int end_thickness_nm   = 120;
static int delta_thickness_nm = 10;

//col : D_Xモード row : D_Yモード
//x,yを中心に, 計算領域のセルと同じ大きさの領域を調べる
static double eps(double x, double y, double z, int col, int row, int dep)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();

  double depth  = max(depth_s[0], depth_s[1]);
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
      fabs(_z) > (depth/2+0.5) ||
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
          s[0] += 0.5*(fabs(sx) < width_s[0]/2 && fabs(sz) < depth_s[0]/2);
          s[1] += 0.5*(fabs(sx) < width_s[1]/2 && fabs(sz) < depth_s[1]/2);
          continue;
        }

        int k = (modY > thickness_s[0]); //どっちの屈折率にいるか調べる

        if (sx < 0 && asymmetry)
          k = 1-k;		//左右で反転, 互い違いでなかったら反転しない

        //正方形の内側ならその媒質内にある
        if( abs(sx) < width_s[k]/2 && abs(sz) < depth_s[k] )
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
  depth_s[0]     = field_toCellUnit(start_depth_nm);
  depth_s[1]     = field_toCellUnit(start_depth_nm);
  thickness_s[0] = field_toCellUnit(start_thickness_nm);
  thickness_s[1] = field_toCellUnit(start_thickness_nm);
  layerNum = 8;
  n[0] = 1.56;
  n[1] = 1.0;
  ep[0] = n[0]*n[0]*EPSILON_0_S;
  ep[1] = n[1]*n[1]*EPSILON_0_S;

  asymmetry = false;

  return eps;
}

void multilayerModel_setThickness(int thickness1_nm, int thickness2_nm)
{
  thickness_s[0] = field_toCellUnit(thickness1_nm);
  thickness_s[1] = field_toCellUnit(thickness2_nm);
}

// 構造が終わったか確認
bool multilayerModel_isFinish()
{
  depth_s[0] += field_toCellUnit(delta_depth_nm);
  depth_s[1] += field_toCellUnit(delta_depth_nm);

  //奥行きが大きくなったら,厚みを増やす
  if(depth_s[0] > field_toCellUnit(end_depth_nm))
  {
    thickness_s[0] += field_toCellUnit(delta_thickness_nm);
    thickness_s[1] += field_toCellUnit(delta_thickness_nm);
    
    depth_s[0] = field_toCellUnit(start_depth_nm);
    depth_s[1] = field_toCellUnit(start_depth_nm);

    //厚みがend_thicknessに達したら終了
    if(thickness_s[0] > field_toCellUnit(end_thickness_nm))
      return true;    
  }

  return false;
}

void multilayerModel_needSize(int *x_nm, int *y_nm,int *z_nm)
{
  double width = max(width_s[0], width_s[1]);
  double depth = max(depth_s[0], depth_s[1]);
  double thick = thickness_s[0] + thickness_s[1];
  (*x_nm) = field_toPhysicalUnit(width);
  (*z_nm) = field_toPhysicalUnit(depth);
  (*y_nm) = field_toPhysicalUnit(thick*layerNum);
/*
  printf("w=%d,d=%d,t=%d\n",(int)field_toPhysicalUnit(width),
         (int)field_toPhysicalUnit(depth),
         (int)field_toPhysicalUnit(thick*layerNum));
*/
}

void multilayerModel_moveDirectory()
{
  char buf[512];

  sprintf(buf, "thick_%dnm_depth_%dnm",(int)field_toPhysicalUnit(thickness_s[0]),
    (int)field_toPhysicalUnit(depth_s[0]));
  makeDirectory(buf);
  moveDirectory(buf);
}
