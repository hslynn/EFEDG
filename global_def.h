#include "phg.h"
#include <math.h>

#define GLOBAL
#define R (Pow(x*x + y*y + z*z, 0.5))
#define Power(x,y) (Pow((FLOAT) (x), (FLOAT) (y)))
#define spacetime 1

extern FLOAT M;
extern INT gamma0;
extern INT gamma1;
extern INT gamma2;
extern INT NVAR;
