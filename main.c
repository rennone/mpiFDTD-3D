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

#endif

static int rank;
static int numProc;

static int start_lambda_nm = 500;
static int end_lambda_nm   = 500;
void move(enum MODEL modelType)
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
  }
}

int main( int argc, char *argv[] )
{
  FieldInfo fInfo;
  fInfo.width_nm  = 500;
  fInfo.height_nm = 500;
  fInfo.depth_nm  = 500;
  fInfo.h_u_nm    = 5;
  fInfo.pml       = 10;
  fInfo.lambda_nm = start_lambda_nm;
  fInfo.stepNum   = 1500;
  fInfo.theta_deg = 0;
  fInfo.phi_deg   = 90;
  enum MODEL modelType   = MIE_SPHERE;//NO_MODEL;LAYER;//
  enum SOLVER solberType = MPI_FDTD_3D;
  move(modelType);
  MPI_Init( 0, 0 );  
  simulator_init(fInfo, modelType, solberType);

#ifndef USE_OPENGL    //only calculate mode
  while(1)
  {
    models_moveDirectory();
    int lambda_nm = field_toPhysicalUnit(field_getLambda_S());
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
    moveDirectory("../");
    //構造が終了か調べる
    if(models_isFinish())
      break;
    
    // 10nm * 2 * 8レイヤ増える todo ハードコーディングはやめる
    fInfo.height_nm += 160;
    //lambdaも元に戻してる.
    simulator_resetField(fInfo);
    SubFieldInfo_S sInfo_s = field_getSubFieldInfo_S();
    if(sInfo_s.Rank == 0)
      printf("next Simulation. lambda = %d, size(%d, %d, %d)\n",
           lambda_nm, fInfo.width_nm, fInfo.height_nm, fInfo.depth_nm);    
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

  if( simulator_isFinish() ){
    MPI_Barrier(MPI_COMM_WORLD);
    simulator_finish();
    MPI_Finalize();
    exit(0);
  }
  glutPostRedisplay();  //再描画
  MPI_Barrier(MPI_COMM_WORLD);
}

#endif
