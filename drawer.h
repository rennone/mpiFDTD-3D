#ifdef USE_OPENGL

#ifndef _DRAWER_H
#define _DRAWER_H
#include "myComplex.h"

enum COLOR_MODE{
  CREAL,
  CABS
};

enum DISPLAY_PLANE{
  XY_PLANE,
  XZ_PLANE,
  ZY_PLANE,
  ALL_PLANE
};

extern void (*drawer_getDraw(void))(void);
extern void drawer_paintImage(int l, int b, int r, int t,int wid, int hei, double complex*);
extern void drawer_paintImage3(dcomplex *phis);
extern void drawer_subFieldPaintImage3(dcomplex *phis, double *eps, enum DISPLAY_PLANE);
extern void drawer_paintModel(int l, int b, int r, int t,int wid, int hei, double *);
extern void drawer_paintTest(void);
extern void drawer_init(enum COLOR_MODE);
extern void drawer_finish(void);
extern void drawer_draw(void);
#endif

extern void drawer_outputImage(char *fileName, double *model, int width, int height, int depth, int (*indexMethod)(int, int, int));
extern void drawer_screenshot(const char* filename);
#endif
