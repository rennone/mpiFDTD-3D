#ifdef USE_OPENGL
#include "drawer.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <GL/glew.h>
#include <GLUT/glut.h>
#include "function.h"
#include "myComplex.h"
#include "field.h"

typedef struct {
 GLfloat r,g,b;
}colorf;

#define TEX_SIZE 512
#define TEX_NX 512
#define TEX_NY 512

static const int vertexNum = 4; //頂点数
static colorf texColor[TEX_NX][TEX_NY]={};
static GLuint ver_buf, tex_buf;
static GLuint texId;
//static GLuint texIds[3];

static double (*colorMode)( dcomplex );
static void colorTransform(double p, colorf *c);


static GLfloat vertices[] =
  {-1.0f, -1.0f, 0.0f,
   +1.0f, -1.0f, 0.0f, 
   +1.0f, +1.0f, 0.0f, 
   -1.0f, +1.0f, 0.0f};

static GLfloat texCoords[] =
  { 0.0f, 0.0f,
    0.0f, 1.0f,
    1.0f, 1.0f,
    1.0f, 0.0f };

//--------------------prototype--------------------//
void drawer_paintImage(int l, int b,int r, int t, int wid,int hei, dcomplex*);
void drawer_paintModel(int l, int b,int r, int t, int wid,int hei, double *);
void drawer_draw();
//--------------------------------------//


//--------------public Method-----------------//
void (*drawer_getDraw(void))(void)
{
  return drawer_draw;
}
//--------------------------------------//

void drawer_init(enum COLOR_MODE cm)
{
  if(cm == CREAL)
    colorMode = creal;
  else
    colorMode = cnorm;//cabs;

  glGenBuffers(1, &ver_buf);

  glBindBuffer(GL_ARRAY_BUFFER, ver_buf);
  glBufferData(GL_ARRAY_BUFFER, 3*vertexNum*sizeof(GLfloat), vertices, GL_STATIC_DRAW);

  glGenBuffers(1, &tex_buf);
  glBindBuffer(GL_ARRAY_BUFFER, tex_buf);
  glBufferData(GL_ARRAY_BUFFER, 2*vertexNum*sizeof(GLfloat), texCoords, GL_STATIC_DRAW);

  //
  glEnable( GL_TEXTURE_2D );
  glGenTextures( 1, &texId );

  glActiveTexture( GL_TEXTURE0 );

  glBindTexture( GL_TEXTURE_2D, texId );
  glTexImage2D( GL_TEXTURE_2D, 0, 3, TEX_NX, TEX_NY, 0, GL_RGB, GL_FLOAT, texColor);

    //min, maxフィルタ
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );    //min, maxフィルター

}

void drawer_draw()
{  
  glBindBuffer( GL_ARRAY_BUFFER, ver_buf);
  glVertexPointer( 3, GL_FLOAT, 0, 0);

  glBindBuffer( GL_ARRAY_BUFFER, tex_buf);
  glTexCoordPointer(2, GL_FLOAT, 0, 0);
  
  glBindTexture( GL_TEXTURE_2D, texId );

  glTexImage2D( GL_TEXTURE_2D, 0, 3, TEX_NX, TEX_NY, 0, GL_RGB, GL_FLOAT, texColor);  

  glDrawArrays( GL_POLYGON, 0, vertexNum);  
}

static double _dbilinear(double *p, double x, double y, int index, int toNextX, int toNextY)
{
  int i = floor(x);
  int j = floor(y);
  double dx = x - i;
  double dy = y - j;
  return p[index]*(1.0-dx)*(1.0-dy)
    + p[index+toNextX]*dx*(1.0-dy)
    + p[index+toNextY]*(1.0-dx)*dy
    + p[index+toNextX+toNextY]*dx*dy;
}

static dcomplex _cbilinear(dcomplex *p, double x, double y, int index, int toNextX, int toNextY)
{
  int i = floor(x);
  int j = floor(y);
  double dx = x - i;
  double dy = y - j;
  return p[index]*(1.0-dx)*(1.0-dy)
    + p[index+toNextX]*dx*(1.0-dy)
    + p[index+toNextY]*(1.0-dx)*dy
    + p[index+toNextX+toNextY]*dx*dy;
}

void drawer_paintImage(int left, int bottom, int right, int top, int width, int height, dcomplex *phis)
{
  colorf c;
  dcomplex cphi;
  double ux = 1.0*(right-left)/TEX_NX;
  double uy = 1.0*(top-bottom)/TEX_NY;
  double u = max(ux,uy);
  int i,j;
  double x,y;
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();

  for(i=0,x=left; i<TEX_NX && x<right; i++, x+=u){
    for(j=0,y=bottom; j<TEX_NY && y<top; j++, y+=u){
      int index = field_subIndex( (int)x, (int)y ,sInfo.SUB_N_PZ/2);
      cphi = _cbilinear(phis,x, y, index, sInfo.SUB_N_PYZ, sInfo.SUB_N_PZ);
//      cphi = cbilinear(phis,x,y,width,height);
      colorTransform(colorMode(cphi), &c);
      texColor[i][j] = c;
    }
  }
}

/*
  3次元空間を3面図で描画する. [0]は横, [1]は縦に描画する
  phis : 3次元空間の複素数データ.
  startIndex       : テクスチャの左下に対応する3D空間のインデックス
  length[2]        : 平面の長さ
  offsetIndex[2]   : 各次元における.1次元配列で保存している為, となりのインデックスまでのオフセット量( 基本的にxはheight*depth, yはdepth, zは1 )
  fixedIndex       : 各面のx,y,zの固定点( yz平面を描画する場合はx座標のインデックスを渡す)
  fixedIndexOffset : 固定点に対する隣のインデックスへのオフセット量
  quadrant         : 描画位置(象限)
 */

static void paint(dcomplex *phis, int startIndex, int length[2], int offsetIndex[2], int fixedIndex, int fixedIndexOffset, int quadrant)
{
  //3D空間とテクスチャのサイズ比
  double ux = 1.0*(length[0])/TEX_SIZE;
  double uy = 1.0*(length[1])/TEX_SIZE;
  double u = max(ux, uy); //一番大きなスケールに合わせる(はみ出ら無いように).
  
  colorf c;
  dcomplex cphi;
  int dx = (quadrant&1)*TEX_SIZE;       //横の象限
  int dy = ((quadrant>>1)&1)*TEX_SIZE;  //縦
  
  int i,j;
  double x,y;
  //第一象限にxy平面を描画 xが横軸
  for(i=0, x=0; i<TEX_SIZE && x<=length[0]; i++, x+=u){
    for(j=0, y=0; j<TEX_SIZE && y<length[1]; j++, y+=u){
      int index = (int)x*offsetIndex[0] + (int)y*offsetIndex[1] + startIndex;
      cphi = _cbilinear(phis, x, y, index, offsetIndex[0], offsetIndex[1]);
      colorTransform(colorMode(cphi), &c);      
      texColor[i+dx][j+dy] = c;
    }
  }
}

void drawer_paintImage3(dcomplex *phis)
{
  FieldInfo_S sInfo = field_getFieldInfo_S();
  double ux = 1.0*sInfo.N_X/TEX_SIZE;
  double uy = 1.0*sInfo.N_Y/TEX_SIZE;
  double uz = 1.0*sInfo.N_Z/TEX_SIZE;
  double u = max(max(ux,uy),uz);
  
  colorf c;
  dcomplex cphi;

  int i,j;
  double x,y,z;

  //第一象限にxy平面を描画 xが横軸
  for(i=0,x=sInfo.N_PML; i<TEX_SIZE && x<sInfo.N_PX; i++, x+=u){
    for(j=0,y=sInfo.N_PML; j<TEX_SIZE && y<sInfo.N_PY; j++, y+=u){
      int index = field_index( (int)x, (int)y ,sInfo.N_PZ/2);
      cphi = _cbilinear(phis, x, y, index, sInfo.N_PYZ, sInfo.N_PZ);
      colorTransform(colorMode(cphi), &c);

      texColor[i+TEX_SIZE][j+TEX_SIZE] = c;
    }
  }  

  //第二象限にxz平面を描画  xが横軸
  for(i=0,x=sInfo.N_PML; i<TEX_SIZE && x<sInfo.N_PX; i++, x+=u){
    for(j=0,z=sInfo.N_PML; j<TEX_SIZE && z<sInfo.N_PZ; j++, z+=u){
      int index = field_index( (int)x, sInfo.N_PY/2, (int)z);
      cphi = _cbilinear(phis, x, z, index, sInfo.N_PYZ, 1);
      colorTransform(colorMode(cphi), &c);
      texColor[i][j+TEX_SIZE] = c;
    }
  }

  //第三象限にzy平面を描画  (zが横軸)
  for(i=0,y=sInfo.N_PML; i<TEX_SIZE && y<sInfo.N_PY; i++, y+=u){
    for(j=0,z=sInfo.N_PML; j<TEX_SIZE && z<sInfo.N_PZ; j++, z+=u){
      int index = field_subIndex( sInfo.N_PX/2, (int)y, (int)z);
      cphi = _cbilinear(phis, z, y, index, 1, sInfo.N_PZ);
      colorTransform(colorMode(cphi), &c);
      texColor[j][i] = c;
    }
  }
}

void drawer_subFieldPaintImage3(dcomplex *phis)
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();

  //第一象限に
  //線の長さは,格子点の数よりも1小さい.
  double ux = (sInfo.SUB_N_X-1.0)/TEX_SIZE;
  double uy = (sInfo.SUB_N_Y-1.0)/TEX_SIZE;
  double uz = (sInfo.SUB_N_Z-1.0)/TEX_SIZE;
  double u = max(max(ux,uy),uz);
  
  colorf c;
  dcomplex cphi;

  int i,j;
  double x,y,z;

  double amp = 10; //描画するときは見やすいように増幅する
  int offsetT, offsetS;

  offsetT = 0;//TEX_SIZE;
  offsetS = 0;//TEX_SIZE;
  //第一象限にxy平面を描画 xが横軸.線形補完で一つ後ろの値も使う為,終了条件はPX-2になる.
  for(i=0,x=1; i<TEX_SIZE && x<sInfo.SUB_N_PX-2; i++, x+=u){
    for(j=0,y=1; j<TEX_SIZE && y<sInfo.SUB_N_PY-2; j++, y+=u){
      int index = field_subIndex( (int)x, (int)y ,sInfo.SUB_N_PZ/2);
      cphi = _cbilinear(phis,x, y, index, sInfo.SUB_N_PYZ, sInfo.SUB_N_PZ);
      colorTransform(colorMode(cphi) * amp, &c);    

      texColor[i+offsetT][j+offsetS] = c;
    }
  }  
/*
  
  //第二象限にxz平面を描画  xが横軸

  offsetT = 0;
  offsetS = TEX_SIZE;
  for(i=0,x=1; i<TEX_SIZE && x<sInfo.SUB_N_PX-2; i++, x+=u){
    for(j=0,z=1; j<TEX_SIZE && z<sInfo.SUB_N_PZ-2; j++, z+=u){
      int index = field_subIndex( (int)x, sInfo.SUB_N_PY/2, (int)z);
      cphi = _cbilinear(phis, x, z, index, sInfo.SUB_N_PYZ, 1);
      colorTransform(colorMode(cphi) * amp, &c);
      texColor[i+offsetT][j+offsetS] = c;
    }
  }

  //第三象限にzy平面を描画  (zが横軸)
  offsetT = 0;
  offsetS = 0;
  for(i=0,z=1; i<TEX_SIZE && z<sInfo.SUB_N_PZ-2; i++, z+=u){
    for(j=0,y=1; j<TEX_SIZE && y<sInfo.SUB_N_PY-2; j++, y+=u){
      int index = field_subIndex( sInfo.SUB_N_PX/2, (int)y, (int)z);
      cphi = _cbilinear(phis, z, y, index, 1, sInfo.SUB_N_PZ);
      colorTransform(colorMode(cphi) * amp, &c);
      texColor[i+offsetT][j+offsetS] = c;
    }
  }*/
}

void drawer_subFieldPaintModel3(double *phis)
{
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();
  
  //第一象限に
  double ux = (sInfo.SUB_N_X-1.0)/TEX_SIZE;
  double uy = (sInfo.SUB_N_Y-1.0)/TEX_SIZE;
  double uz = (sInfo.SUB_N_Z-1.0)/TEX_SIZE;
  double u = max(max(ux,uy),uz);
  
  colorf c;
  double dphi;

  int i,j;
  double x,y,z;

  //第一象限にxy平面を描画 xが横軸
  //モデルは端っこの
  for(i=0,x=1; i<TEX_SIZE && x<sInfo.SUB_N_PX-2; i++, x+=u){
    for(j=0,y=1; j<TEX_SIZE && y<sInfo.SUB_N_PY-2; j++, y+=u){
      int index = field_subIndex( (int)x, (int)y ,sInfo.SUB_N_PZ/2);
      dphi = _dbilinear(phis,x, y, index, sInfo.SUB_N_PYZ, sInfo.SUB_N_PZ);

      double n = 1-1.0/dphi;

      //中心には線を引く
      if(i==TEX_SIZE/2 || j==TEX_SIZE/2) {
        texColor[i][j].r = 1;
        texColor[i][j].g = 1;
        texColor[i][j].b = 1;
      } else {
        texColor[i][j].r -= n;
        texColor[i][j].g -= n;
        texColor[i][j].b -= n;
      }
    }
  }
  
/*
  //第二象限にxz平面を描画  xが横軸
  for(i=0,x=1; i<TEX_SIZE && x<sInfo.SUB_N_PX-1; i++, x+=u){
    for(j=0,z=1; j<TEX_SIZE && z<sInfo.SUB_N_PZ-1; j++, z+=u){
      int index = field_subIndex( (int)x, sInfo.SUB_N_PY/2, (int)z);
      dphi = _dbilinear(phis,x, y, index, sInfo.SUB_N_PYZ, sInfo.SUB_N_PZ);
      
      double n = 1-1.0/dphi;
      texColor[i][j+TEX_SIZE].r -= n;
      texColor[i][j+TEX_SIZE].g -= n;
      texColor[i][j+TEX_SIZE].b -= n;
    }
  }

  //第三象限にzy平面を描画  (zが横軸)
  for(i=0,y=1; i<TEX_SIZE && y<sInfo.SUB_N_PY-1; i++, y+=u){
    for(j=0,z=1; j<TEX_SIZE && z<sInfo.SUB_N_PZ-1; j++, z+=u){
      int index = field_subIndex( sInfo.SUB_N_PX/2, (int)y, (int)z);
      dphi = _dbilinear(phis,x, y, index, sInfo.SUB_N_PYZ, sInfo.SUB_N_PZ);
      double n = 1-1.0/dphi;
      texColor[i][j].r -= n;
      texColor[i][j].g -= n;
      texColor[i][j].b -= n; 
    }
    }*/
}

void drawer_paintModel(int left, int bottom, int right, int top, int width, int height, double *phis)
{
  double dphi;
  double ux = 1.0*(right-left)/TEX_NX;
  double uy = 1.0*(top-bottom)/TEX_NY;
  double u = max(ux,uy);
  int i,j;
  double x,y;
  
  for(i=0,x=left; i<TEX_NX && x<right; i++, x+=u)
  {
    for(j=0,y=bottom; j<TEX_NY && y<top; j++, y+=u)
    {      
      dphi = dbilinear(phis,x,y,width,height);
      double n = 1-1.0/dphi;
      texColor[i][j].r -= n;
      texColor[i][j].g -= n;
      texColor[i][j].b -= n;
    }
  }
}

void drawer_paintTest(void){
  colorf c;
  for(int i=0; i<TEX_NX; i++){
    for(int j=0; j<TEX_NY; j++){
      colorTransform(1.0*i/TEX_NX, &c);
      texColor[i][j] = c;
    }
  }
}

void drawer_finish()
{
  printf("drawer_finish not implemented\n");
}

//--------------------public Method--------------------//

//--------------------Color Trancform---------------------//
static void colorTransform(double phi, colorf *col)
{
  double range = 1.0; //波の振幅  
  double ab_phi = phi < 0 ? -phi : phi;
  double a = ab_phi < range ? (ab_phi <  range/3.0 ? 3.0/range*ab_phi : (-3.0/4.0/range*ab_phi+1.25) ) : 0.5;
  
  col->r = phi > 0 ? a:0;
  col->b = phi < 0 ? a:0;
  col->g = min(1.0, max(0.0, -3*ab_phi+2));
}

#endif



