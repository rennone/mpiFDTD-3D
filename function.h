#ifndef _FUNCTION_H
#define _FUNCTION_H
#include "myComplex.h"
#include <stdio.h>
#include "bool.h"
#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

extern double dbilinear(double *p, double x, double y, int width, int height);
extern FILE* openFile(const char* file_name);
extern void makeDirectory(const char*);
extern void moveDirectory(const char*);

// w : int (配列のインデックスが入る), sInfo_s : subFieldIndo_S型
#define FAST_3FOR(w, sInfo_s, nextX)  w = field_subIndex(1,1,1); \
  nextX = 2*sInfo_s.SUB_N_PZ;                                    \
  for(int i=sInfo_s.SUB_N_X; i--;  w+=nextX)                     \
    for(int j=sInfo_s.SUB_N_Y; j--;  w+=2)                       \
      for(int k=sInfo_s.SUB_N_Z; k--;  w+=1)                     \

//まだ実装していない部分を,やり忘れないようにexitさせたいときに使う.
#define NOT_DONE(msg) printf(msg); exit(2)

#define STOP(msg) printf(msg); exit(2)
#endif
