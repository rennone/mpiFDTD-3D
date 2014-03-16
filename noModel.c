#include "noModel.h"
#include "field.h"

static double eps(double x, double y, double z, int col, int row, int dep){
  return EPSILON_0_S;
}

double (*noModel_EPS(void))(double, double, double, int, int, int){
  return eps;
}

