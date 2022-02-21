#ifndef GLOBAL
#define GLOBAL

#include "phg.h"
#include <math.h>

#define R (Pow(x*x + y*y + z*z, 0.5))
#define Power(x,y) (Pow((FLOAT) (x), (FLOAT) (y)))

extern INT SPACETIME;
extern FLOAT M;
extern FLOAT A;
extern INT gamma0;
extern INT gamma1;
extern INT gamma2;
extern INT NVAR;

#endif


