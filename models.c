#include "field.h"
#include "models.h"
#include "noModel.h"
#include "circleModel.h"
//#include "shelf.h"
//#include "nonshelf.h"


static double (*epsMethod)(double, double, double, int, int, int);

static void noModel(void)
{
  //no material
  epsMethod = noModel_EPS();
}

static void circleModel(void)
{
  //cylinder material whitch radius = lambda, origin = center of field
  epsMethod = circleModel_EPS();
}

void setModel(enum MODEL model)
{
  switch(model){
  case NO_MODEL:
    noModel();
    break;
  case MIE_CYLINDER:
    circleModel();
    break;
  case SHELF :
  case NONSHELF:
  case LAYER:
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
  
  if(epsilon != EPSILON_0_S)
  {
    printf("%lf %lf, %lf, %lf", epsilon, x, y, z);
  }
  return epsilon;
}
