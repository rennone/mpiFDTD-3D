#ifndef _MODELS_H
#define _MODELS_H

enum MODEL{
  NO_MODEL,
  MIE_SPHERE,
  LAYER
};

enum MODE{
  D_X, //x方向に線積分
  D_Y, //y方向に線積分
  D_Z, //z方向に線積分
  D_XY, //xy面積分
  D_XZ, //xy面積分
  D_YZ, //xy面積分
  D_XYZ, //xy面積分 
};

extern void models_setModel(enum MODEL model);
extern void models_init();
extern double models_eps(double x, double y, double z, enum MODE mode);
extern bool models_isFinish();
extern void models_needSize(int *x_nm, int *y_nm,int *z_nm);
extern void models_moveDirectory(void);
#endif
