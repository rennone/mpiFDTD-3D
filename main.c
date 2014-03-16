#include "drawer.h" //ここに_USE_OPENGLを定義
#include <stdio.h>
#include <stdlib.h>
#ifdef _USE_OPENGL
  #include <GL/glew.h>
  #include <GLUT/glut.h>
#endif
#include <mpi.h>
#include <complex.h>
#include "simulator.h"
#include "mpiTM_UPML.h"
#include "mpiTE_UPML.h"
int windowX = 100;
int windowY = 100;
int windowWidth = 300;
int windowHeight=300;

void display(void)
{
#ifdef _USE_OPENGL
  glEnableClientState( GL_VERTEX_ARRAY );
  glEnableClientState( GL_TEXTURE_COORD_ARRAY );

  drawer_paintImage(1,1, fdtdTE_upml_getSubNx(), fdtdTE_upml_getSubNy(),
                         fdtdTE_upml_getSubNpx(), fdtdTE_upml_getSubNpy(), simulator_getDrawingData());
  //drawer_paintImage2(1,1, fdtdTE_upml_getSubNx(), fdtdTE_upml_getSubNy(), fdtdTE_upml_getSubNpx(), fdtdTE_upml_getSubNpy(), fdtdTE_upml_getEps());
  drawer_draw();
    
  glDisableClientState( GL_VERTEX_ARRAY );
  glDisableClientState( GL_TEXTURE_COORD_ARRAY );
  glutSwapBuffers();
#endif
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
#ifdef _USE_OPENGL
  glutPostRedisplay();  //再描画
#endif  
  MPI_Barrier(MPI_COMM_WORLD);
}

int main( int argc, char *argv[] )
{
    MPI_Init( 0, 0 ); 
    simulator_init();
#ifdef _USE_OPENGL
    glutInit(&argc, argv);
    glutInitWindowPosition(windowX,windowY);
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutCreateWindow("FDTD Simulator");
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glewInit();
    drawer_init(CABS);
    glutMainLoop();
    MPI_Finalize();
#endif

#ifndef _USE_OPENGL    //only calculate mode

    while(!simulator_isFinish()){
       simulator_calc();
       //   MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    simulator_finish();
    MPI_Finalize();
#endif

    return 1;
}
