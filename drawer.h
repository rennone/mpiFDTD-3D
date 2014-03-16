#define _USE_OPENGL

#ifdef _USE_OPENGL
#ifndef _DRAWER_H
#define _DRAWER_H
#include <complex.h>

enum COLOR_MODE{
  CREAL,
  CABS
};

extern void (*drawer_getDraw(void))(void);
extern void drawer_paintImage(int l, int b, int r, int t,int wid, int hei, double complex*, ...);
extern void drawer_paintImage2(int l, int b, int r, int t,int wid, int hei, double *, ...);
extern void drawer_paintTest(void);
extern void drawer_init(enum COLOR_MODE);
extern void drawer_finish(void);
extern void drawer_draw(void);
#endif
#endif
