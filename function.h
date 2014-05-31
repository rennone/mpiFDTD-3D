#ifndef _FUNCTION_H
#define _FUNCTION_H
#include "myComplex.h"
#include <stdio.h>
#include "bool.h"
#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

extern double dbilinear(double *p, double x, double y, int width, int height);
extern FILE* openFile(const char* file_name);
extern bool makeDirectory(const char*);
extern void moveDirectory(const char*);


//まだ実装していない部分を,やり忘れないようにexitさせたいときに使う.
#define NOT_DONE(msg) printf(msg); exit(2)

#define STOP(msg) printf(msg); exit(2)
#endif
