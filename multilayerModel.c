#include "multilayerModel.h"
#include "field.h"
#include "function.h"
#include <stdlib.h>
#include <math.h>

//円形の多層膜
#define CIRCLE_LAYER true

//幅(x)
#define ST_WIDTH_NM 300
#define EN_WIDTH_NM 300
#define DELTA_WIDTH_NM 100

//奥行き(z)
#define ST_DEPTH_NM 300
#define EN_DEPTH_NM 300
#define DELTA_DEPTH_NM 100

//ラメラの厚さ(y)
#define ST_THICK_NM 80
#define EN_THICK_NM 120
#define DELTA_THICK_NM 10

//ラメラの枚数
#define LAYER_NUM 8

//互い違い
#define ASYMMETRY false

//中心に以下の幅で軸となる枝を入れる => 軸の屈折率はN_1になる
#define BRANCH_NM 50

//屈折率
#define N_0 1.0
#define N_1 1.56

//先端における横幅の割合
#define ST_EDGE_RATE 1.0
#define EN_EDGE_RATE 1.0
#define DELTA_EDGE_RATE 0.1

//ラメラの先端を丸める曲率 (1で四角形のまま, 0.0で最もカーブする)
#define CURVE 0.8

static int depth_nm[2]     = {ST_DEPTH_NM, ST_DEPTH_NM};
static int width_nm[2]     = {ST_WIDTH_NM, ST_WIDTH_NM};
static int thickness_nm[2] = {ST_THICK_NM, ST_THICK_NM};
static int layerNum        = LAYER_NUM;
static int branch_width_nm = BRANCH_NM;

static double depth_s[2];     //奥行き z
static double width_s[2];     //幅 x
static double thickness_s[2]; //厚さ y

static double edge_width_rate = ST_EDGE_RATE;

static double n[2];            //屈折率
static double ep[2];           //誘電率 = n*n*ep0

//原点位置(中心の下)
static double posx, posy, posz;

//2次関数の比例定数
static double c0_x, c1_x, c0_z, c1_z;

//楕円形の多層膜
static double eps_circle(double x, double y, double z, int col, int row, int dep)
{
  double depth  = max(depth_s[0], depth_s[1]);
  double width  = max(width_s[0], width_s[1]);
  double thick  = thickness_s[0] + thickness_s[1];
  double height = thick*LAYER_NUM;

  double _x = x-posx;	//pox,poyを座標の原点に
  double _y = y-posy;
  double _z = z-posz;

  //上下左右に飛び出ていないか確認(細分化したセルがあるため, 0.5の余白をとっている)
  if( pow(_x*depth/2.0, 2) + pow(_z*width/2.0, 2) > pow(width*depth/4.0,2)+1.0 ||
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
          s[0] += 0.5*(pow(sx*depth_s[0]/2.0, 2) + pow(sz*width_s[0]/2.0, 2) <= pow(width_s[0]*depth_s[0]/4.0,2));
          s[1] += 0.5*(pow(sx*depth_s[1]/2.0, 2) + pow(sz*width_s[1]/2.0, 2) <= pow(width_s[1]*depth_s[1]/4.0,2));
          continue;
          }

        int k = (modY > thickness_s[0]); //どっちの屈折率にいるか調べる

        if (sx < 0 && ASYMMETRY)
          k = 1-k;		//左右で反転, 互い違いでなかったら反転しない

        //楕円系の内側ならその媒質内にある
        if( pow(sx*depth_s[k]/2.0, 2) + pow(sz*width_s[k]/2.0, 2) <= pow(width_s[k]*depth_s[k]/4.0,2))
          s[k] += 1;
      }

  s[0] /= split*split*split;
  s[1] /= split*split*split;
  return EPSILON_0_S*(1-s[0]-s[1]) + ep[0]*s[0] + ep[1]*s[1];
}

//col : D_Xモード row : D_Yモード
//x,yを中心に, 計算領域のセルと同じ大きさの領域を調べる
static double eps(double x, double y, double z, int col, int row, int dep)
{
  double depth  = max(depth_s[0], depth_s[1]);
  double width  = max(width_s[0], width_s[1]);
  double thick  = thickness_s[0] + thickness_s[1];
  double height = thick*LAYER_NUM;

  double _x = x-posx;	//ox,oyを座標の原点に
  double _y = y-posy;
  double _z = z-posz;

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

        if (sx < 0 && ASYMMETRY)
          k = 1-k;		//左右で反転, 互い違いでなかったら反転しない

        //正方形の内側ならその媒質内にある
        if( abs(sx) < width_s[k]/2 && abs(sz) < depth_s[k] )
          s[k] +=1;
      }

  s[0] /= split*split*split;
  s[1] /= split*split*split;
  return EPSILON_0_S*(1-s[0]-s[1]) + ep[0]*s[0] + ep[1]*s[1];
}

double ( *multilayerModel_EPS(void))(double x, double y, double z, int, int, int)
{
  if(CIRCLE_LAYER)
    return eps_circle;
  else
    return eps;
}

void multilayerModel_init()
{
  width_s[0]     = field_toCellUnit(width_nm[0]);
  width_s[1]     = field_toCellUnit(width_nm[0]);
  depth_s[0]     = field_toCellUnit(depth_nm[0]);
  depth_s[1]     = field_toCellUnit(depth_nm[0]);
  thickness_s[0] = field_toCellUnit(thickness_nm[0]);
  thickness_s[1] = field_toCellUnit(thickness_nm[0]);
  n[0] = 1.56;
  n[1] = 1.0;
  ep[0] = n[0]*n[0]*EPSILON_0_S;
  ep[1] = n[1]*n[1]*EPSILON_0_S;

  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  //領域の中心から, 下にheight/2ずれた位置がレイヤの下部
  double height = (thickness_s[0] + thickness_s[1])*LAYER_NUM;
  posy = fInfo_s.N_PY/2 - height/2;
  posx = fInfo_s.N_PX/2;
  posz = fInfo_s.N_PZ/2;  
}


// 構造が終わったか確認
bool multilayerModel_isFinish()
{
  depth_nm[0] += DELTA_DEPTH_NM;
  depth_nm[1] += DELTA_DEPTH_NM;

  //奥行きが大きくなったら,厚みを増やす
  if(depth_nm[0] > EN_DEPTH_NM)
  {
    thickness_nm[0] += DELTA_THICK_NM;
    thickness_nm[1] += DELTA_THICK_NM;

    depth_nm[0] = ST_DEPTH_NM;
    depth_nm[1] = ST_DEPTH_NM;
    
    //厚みがend_thicknessに達したら終了
    if(thickness_nm[0] > EN_THICK_NM)
      return true;    
  }

  return false;
}

void multilayerModel_needSize(int *x_nm, int *y_nm,int *z_nm)
{ 
  (*x_nm) = max(width_nm[0], width_nm[1]);
  (*z_nm) = max(depth_nm[0], depth_nm[1]);
  (*y_nm) = (thickness_nm[0] + thickness_nm[1])*LAYER_NUM;

  printf("w=%d,d=%d,t=%d\n",*x_nm, *y_nm, *z_nm);
}

void multilayerModel_moveDirectory()
{
  char buf[512];

  sprintf(buf, "thick_%dnm_depth_%dnm",thickness_nm[0],depth_nm[0]);
  makeDirectory(buf);
  moveDirectory(buf);
}
