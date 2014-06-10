#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "myComplex.h"
#include "simulator.h"
#include "field.h"
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

int rank;
int numProc;

int main( int argc, char *argv[] )
{
  FieldInfo fInfo;
  fInfo.width_nm  = 640;
  fInfo.height_nm = 640;
  fInfo.depth_nm  = 640;
  fInfo.h_u_nm    = 10;
  fInfo.pml       = 10;
  fInfo.lambda_nm = 80;
  fInfo.stepNum   = 1200;
  fInfo.theta_deg = 0;
  fInfo.phi_deg   = 0;
  enum MODEL modelType   = MIE_SPHERE;//NO_MODEL;
  enum SOLVER solberType = MPI_FDTD_3D;
  MPI_Init( 0, 0 );
  simulator_init(fInfo, modelType, solberType);

#ifndef USE_OPENGL    //only calculate mode
    while(!simulator_isFinish())
    {
       simulator_calc();       
       //   MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    simulator_finish();
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

    return 1;
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

  drawer_subFieldPaintImage3(data3D, eps3D);
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
