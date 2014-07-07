#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "myComplex.h"
#include "simulator.h"
#include "field.h"
#include "multilayerModel.h"
#include "models.h"
#include "function.h"
// 以下 OPEN_GLの関数
#ifdef USE_OPENGL

#include "drawer.h"
#define WINDOW_WIDTH 300
#define WINDOW_HEIGHT 300

#include <GL/glew.h>

//Macの場合
#ifdef MAC_OS
#include <GLUT/glut.h>
#endif

//Mac以外
#ifndef MAC_OS
#include <GL/glut.h>
#endif

//プロトタイプ宣言
static void drawField();
static void drawSubField();
static void display(void);
static void idle(void);

#endif //USE_OPENGL


char root[512];

#define ST_LAMBDA_NM 380
#define EN_LAMBDA_NM 380
#define DELTA_LAMBDA_NM 10
static int start_lambda_nm = 380;
static int end_lambda_nm   = 640;
static int lambda_nm;

static void cpy(double *entire, double *region, int dx, int dy, int dz)
{
  //sub領域のregionをentireにコピー
  //SUB_N_PXとかは, 全プロセスで共通(なはず)なので, subInfo_sの値をそのまま使う
  //オフセットはプロセスごとに違うので外部から与える.  
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();  
  for(int i=1; i<subInfo_s.SUB_N_PX-1; i++)  
    for(int j=1; j<subInfo_s.SUB_N_PY-1; j++)
      for(int k=1; k<subInfo_s.SUB_N_PZ-1; k++)
      {
        int w = field_index(i+dx,j+dy,k+dz);
        entire[w] = region[field_subIndex(i,j,k)];
      }  
}

//分割された領域をまとめる.
static double* unifyToRank0(double *phi)
{
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  
  //マスターにすべて集める
  if(subInfo_s.Rank == 0)
  {
    MPI_Status status;
    double *entire = newDouble(fInfo_s.N_CELL);
    cpy(entire, phi, subInfo_s.OFFSET_X, subInfo_s.OFFSET_Y, subInfo_s.OFFSET_Z);

    double *tmp = newDouble(subInfo_s.SUB_N_CELL);
    int offset[3];
    for(int i=1; i<subInfo_s.Nproc; i++)
    {
      MPI_Recv(offset, 3, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(tmp, subInfo_s.SUB_N_CELL, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

      cpy(entire, tmp, offset[0], offset[1], offset[2]);      
    }
    free(tmp);
    return entire;
  }
  else {
    int offset[3];
    offset[0] = subInfo_s.OFFSET_X;
    offset[1] = subInfo_s.OFFSET_Y;
    offset[2] = subInfo_s.OFFSET_Z;
    MPI_Send(offset, 3, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(phi, subInfo_s.SUB_N_CELL, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    
    return NULL; //マスター以外はNULLを返す.
  }
}


static void outputImage()
{
  double *entire = unifyToRank0(simulator_getEps());
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  if(entire != NULL)
    drawer_outputImage("image.bmp", entire, fInfo_s.N_PX, fInfo_s.N_PY, fInfo_s.N_PZ, field_index);

  MPI_Barrier(MPI_COMM_WORLD);
}

/*
static void move(enum MODEL modelType)
{
  switch(modelType)
  {
  case LAYER:
    makeDirectory("Multilayer");
    moveDirectory("Multilayer");
    break;
  case MIE_SPHERE:
    makeDirectory("Mie");
    moveDirectory("Mie");
    break;
  default :
    printf("set Model Layer or Sphere");
    exit(2);
  }
  }*/

//モデルとは別に必要なフィールドの大きさ
static void setFieldSize(FieldInfo *fInfo, int x_nm, int y_nm, int z_nm)
{
  //pmlレイヤ + 遠方解の積分路(端から5ずつ) + 余白
  int extends = fInfo->h_u_nm*(fInfo->pml + 5)*2 + 200;

  fInfo->width_nm  = x_nm + extends;
  fInfo->height_nm = y_nm + extends;
  fInfo->depth_nm  = z_nm + extends;
}

int main( int argc, char *argv[] )
{
  getcwd(root, 512); //カレントディレクトリを保存
  
  enum MODEL modelType   = LAYER;//MIE_SPHERE;//NO_MODEL;
  models_setModel(modelType);
  
  enum SOLVER solberType = MPI_FDTD_3D;
  simulator_setSolver(solberType);

  //モデルに必要な大きさを求める
  int x_nm, y_nm, z_nm;
  models_needSize(&x_nm, &y_nm, &z_nm);

  FieldInfo fInfo;

  fInfo.h_u_nm    = 10;
  fInfo.pml       = 10;
  fInfo.lambda_nm = start_lambda_nm;
  fInfo.stepNum   = 10;
  fInfo.theta_deg = 0;
  fInfo.phi_deg   = 90;
  
  setFieldSize(&fInfo, x_nm, y_nm, z_nm);

  printf("%d, %d, %d\n",fInfo.width_nm, fInfo.height_nm, fInfo.depth_nm);
  
//  move(modelType);
  models_moveDirectory();
  
  MPI_Init( 0, 0 );
  simulator_init(fInfo);
  outputImage();
  
#ifndef USE_OPENGL    //only calculate mode
  while(1)
  {  
    lambda_nm = field_toPhysicalUnit(field_getLambda_S());
    while(lambda_nm <= end_lambda_nm)
    {
      while(!simulator_isFinish())
      {
        simulator_calc();
      }
      MPI_Barrier(MPI_COMM_WORLD);
      lambda_nm += 10;
      
      if(lambda_nm > end_lambda_nm)
        break;
      
      simulator_reset();
      field_setLambda(lambda_nm);
      SubFieldInfo_S sInfo_s = field_getSubFieldInfo_S();
      if(sInfo_s.Rank == 0)
        printf("next Simulation. lambda = %d\n", lambda_nm);
    }
    simulator_finish(); //構造を変えるのでシミュレーションを終わらせる.
    
    //構造が終了か調べる
    if(models_isFinish())
      break;

    //移動し直し
    moveDirectory(root);
    models_moveDirectory();
    
    int x_nm, y_nm, z_nm;
    //モデルに必要な大きさを求める
    models_needSize(&x_nm, &y_nm, &z_nm);
    setFieldSize(&fInfo, x_nm, y_nm, z_nm);
    fInfo.lambda_nm = ST_LAMBDA_NM;
    simulator_init(fInfo);

    // outputImage();
    SubFieldInfo_S sInfo_s = field_getSubFieldInfo_S();
    if(sInfo_s.Rank == 0)
      printf("next Simulation. lambda = %d, size(%d, %d, %d)\n", lambda_nm, fInfo.width_nm, fInfo.height_nm, fInfo.depth_nm);    
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Finalize();
#endif

#ifdef USE_OPENGL
  SubFieldInfo_S subInfo = field_getSubFieldInfo_S(); 
  int windowX = 1.0*subInfo.OFFSET_X / subInfo.SUB_N_X * WINDOW_WIDTH;
  int windowY = 800-1.0*subInfo.OFFSET_Y/subInfo.SUB_N_Y * WINDOW_HEIGHT - WINDOW_HEIGHT;
  enum COLOR_MODE colorMode = CREAL;

  glutInit(&argc, argv);
  glutInitWindowPosition(windowX,windowY);
  glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
  glutCreateWindow("FDTD Simulator");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  glewInit();
  drawer_init(colorMode);
  glutMainLoop();
#endif

  return 0;
}

#ifdef USE_OPENGL
static void drawField()
{
  FieldInfo_S fInfo = field_getFieldInfo_S();
  dcomplex *data3D = simulator_getDrawingData();
  double *eps3D = simulator_getEps();
}

static void drawSubField()
{

  dcomplex *data3D = simulator_getDrawingData();
  double *eps3D    = simulator_getEps();

  FieldInfo_S fInfo = field_getFieldInfo_S();
  SubFieldInfo_S sInfo = field_getSubFieldInfo_S();

  drawer_subFieldPaintImage3(data3D, eps3D, ALL_PLANE);
}

void display(void)
{
  glEnableClientState( GL_VERTEX_ARRAY );
  glEnableClientState( GL_TEXTURE_COORD_ARRAY );

  drawSubField();
  drawer_draw();
    
  glDisableClientState( GL_VERTEX_ARRAY );
  glDisableClientState( GL_TEXTURE_COORD_ARRAY );
  glutSwapBuffers();
}


void idle(void)
{
  simulator_calc();
  
  if( !simulator_isFinish() )
    return;
  
  //終了したときの処理

  //波長を変える
  lambda_nm += 10;

  //すべての波長を計算したら->構造を変える
  if(lambda_nm > EN_LAMBDA_NM)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    simulator_finish(); //シミュレーションを終了させる(メモリを解放する)

//    moveDirectory("../"); //一つ上のディレクトリに戻る
    
    //構造が終了か調べる
    if( !models_isFinish())
    {
      FieldInfo fInfo = field_getFieldInfo();

      int x_nm, y_nm, z_nm;
      models_needSize(&x_nm, &y_nm, &z_nm);      //フィールドの大きさを構造に合わせて変化させる
      setFieldSize(&fInfo, x_nm, y_nm, z_nm);
      
      fInfo.lambda_nm = ST_LAMBDA_NM;      //lambdaも元に戻してる.
      simulator_init(fInfo);

      outputImage();
      
      SubFieldInfo_S sInfo_s = field_getSubFieldInfo_S();
      if(sInfo_s.Rank == 0)
        printf("next Simulation. lambda = %d, size(%d, %d, %d)\n",
               lambda_nm, fInfo.width_nm, fInfo.height_nm, fInfo.depth_nm);

      models_moveDirectory();
    }
  } else {
    field_setLambda(lambda_nm);
    simulator_reset(); //電磁波の状態とステップ数だけリセット
    return;
  }  

  MPI_Finalize();
  exit(0);
  
  glutPostRedisplay();  //再描画
}

#endif
