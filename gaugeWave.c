#include "phg.h"
#include "global_def.h"
#define PI 3.1415926535897932384

void
func_Psi00(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -1. - AA*Sin((x-TIME)*(2*PI));
}

void
func_Psi01(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Psi02(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Psi03(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Psi11(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 1. + AA*Sin((x-TIME)*(2*PI));
}

void
func_Psi12(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Psi13(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}


void
func_Psi22(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 1;
}

void
func_Psi23(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Psi33(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{

    *value = 1;
}

//nonzero funcs of Pi
void
func_Pi00(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -AA*Cos((x-TIME)*(2*PI))*(2*PI)/Pow(1 + AA*Sin((x-TIME)*(2*PI)), 0.5);
}

void
func_Pi01(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Pi02(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Pi03(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Pi11(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = AA*Cos((x-TIME)*(2*PI))*(2*PI)/Pow(1 + AA*Sin((x-TIME)*(2*PI)), 0.5);
}

void
func_Pi12(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Pi13(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Pi22(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Pi23(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Pi33(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

//nonzero funcs of Phi
void
func_Phi100(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -AA*(2*PI)*Cos((2*PI)*(x-TIME));
}

void
func_Phi101(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0; 
}

void
func_Phi102(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0; 
}

void
func_Phi103(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0; 
}

void
func_Phi111(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = AA*(2*PI)*Cos((2*PI)*(x-TIME));
} 

void
func_Phi112(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0; 
}

void
func_Phi113(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi122(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi123(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi133(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi200(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi201(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi202(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi203(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi211(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0; 
}

void
func_Phi212(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi213(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi222(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi223(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi233(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi300(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi301(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi302(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi303(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi311(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi312(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi313(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0; 
}

void
func_Phi322(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi323(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_Phi333(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

/*functions of H_a*/
void
func_H0(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0; 
}

void
func_H1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0; 
}

void
func_H2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0; 
}

void
func_H3(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0; 
}

/*functions of deriH*/
void
func_deriH00(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0.; 
}

void
func_deriH01(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0.;
}

void
func_deriH02(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0.;
}

void
func_deriH03(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0.;
}

void
func_deriH10(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_deriH11(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_deriH12(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_deriH13(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_deriH20(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_deriH21(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_deriH22(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_deriH23(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_deriH30(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_deriH31(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_deriH32(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void
func_deriH33(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}


