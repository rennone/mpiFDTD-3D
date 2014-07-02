#include "field.h"
#include "models.h"
#include "noModel.h"
#include "circleModel.h"
#include "multilayerModel.h"
//#include "shelf.h"
//#include "nonshelf.h"

static void (*initMethod)();
static bool (*isFinishMethod)(void);
static bool (*moveDirectoryMethod)(void);
static double (*epsMethod)(double, double, double, int, int, int);
static void (*needSizeMethod)(int*, int*, int*);

static void noModel(void)
{
  epsMethod = noModel_EPS();
}

static void circleModel(void)
{ 
  epsMethod = circleModel_EPS();
  initMethod = multilayerModel_init;
  isFinishMethod = circleModel_isFinish;
  moveDirectoryMethod = circleModel_moveDirectory;
}

static void multilayerModel()
{
  epsMethod      = multilayerModel_EPS();
  initMethod     = multilayerModel_init;
  isFinishMethod = multilayerModel_isFinish;
  needSizeMethod = multilayerModel_needSize;
  moveDirectoryMethod = multilayerModel_moveDirectory;
}

void models_setModel(enum MODEL model)
{
  switch(model){
  case NO_MODEL:
    noModel();
    break;
  case MIE_SPHERE:
    circleModel();
    break;
  case LAYER:
    multilayerModel();
    break;
  }
}

double models_eps(double x, double y, double z, enum MODE mode){
  double epsilon = EPSILON_0_S;
  switch(mode){
  case D_X :
    epsilon = (*epsMethod)(x, y, z, 1, 0, 0);
  case D_Y :
    epsilon = (*epsMethod)(x, y, z, 0, 1, 0);
  case D_Z :
    epsilon = (*epsMethod)(x, y, z, 0, 0, 1);
  case D_XY :
    epsilon = (*epsMethod)(x, y, z, 1, 1, 0);
  case D_XZ :
    epsilon = (*epsMethod)(x, y, z, 1, 0, 1);
  case D_YZ :
    epsilon = (*epsMethod)(x, y, z, 0, 1, 1);
  case D_XYZ :
    epsilon = (*epsMethod)(x, y, z, 1, 1, 1);
  }  
  return epsilon;
}

void models_init()
{
  (*initMethod)();
}

bool models_isFinish()
{
  return (*isFinishMethod)();
}

//必要な計算領域のサイズをとってくる
void models_needSize(int *x_nm, int *y_nm,int *z_nm)
{
  needSizeMethod(x_nm, y_nm, z_nm);
}

void models_moveDirectory(void)
{
  (*moveDirectoryMethod)();
}
