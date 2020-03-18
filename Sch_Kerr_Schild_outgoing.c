#include "phg.h"
#include "global_def.h"

void
func_Psi00(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -1. + (2*M)/R;
}

void
func_Psi01(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -(2*M*x)/Pow(R, 2.);
}

void
func_Psi02(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -(2*M*y)/Pow(R, 2.);
}

void
func_Psi03(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -(2*M*z)/Pow(R, 2.);
}

void
func_Psi11(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 1 + (2*M*Pow(x,2))/Pow(R, 3.);
}

void
func_Psi12(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (2*M*x*y)/Pow(R, 3.);
}

void
func_Psi13(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (2*M*x*z)/Pow(R, 3.);
}


void
func_Psi22(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 1 + (2*M*Pow(y,2))/Pow(R, 3.);
}

void
func_Psi23(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (2*M*y*z)/Pow(R, 3.);
}

void
func_Psi33(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{

    *value = 1 + (2*M*Pow(z,2))/Pow(R, 3.);
}

//nonzero funcs of Pi
void
func_Pi00(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*(2.*Power(x,2) + 2.*Power(y,2) + 2.*Power(z,2)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),1.)*
     (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
       0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. - (1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.))/
        (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
          0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)),0.5));
}

void
func_Pi01(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*x*((-4.*Power(x,2) - 4.*Power(y,2) - 4.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
       2.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)*
     (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
       0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. - (1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.))/
        (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
          0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)),0.5));
}

void
func_Pi02(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*y*((-4.*Power(x,2) - 4.*Power(y,2) - 4.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
       2.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)*
     (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
       0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. - (1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.))/
        (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
          0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)),0.5));
}

void
func_Pi03(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*z*((-4.*Power(x,2) - 4.*Power(y,2) - 4.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
       2.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)*
     (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
       0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. - (1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.))/
        (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
          0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)),0.5));
}

void
func_Pi11(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*Power(x,2)*(-4.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
       Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
        (6.*Power(x,2) + 6.*Power(y,2) + 6.*Power(z,2))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)*
     (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
       0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. - (1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.))/
        (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
          0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)),0.5));
}

void
func_Pi12(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*x*y*(-4.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
       Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
        (6.*Power(x,2) + 6.*Power(y,2) + 6.*Power(z,2))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)*
     (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
       0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. - (1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.))/
        (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
          0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)),0.5));
}

void
func_Pi13(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*x*z*(-4.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
       Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
        (6.*Power(x,2) + 6.*Power(y,2) + 6.*Power(z,2))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)*
     (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
       0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. - (1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.))/
        (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
          0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)),0.5));
}

void
func_Pi22(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*Power(y,2)*(-4.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
       Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
        (6.*Power(x,2) + 6.*Power(y,2) + 6.*Power(z,2))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)*
     (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
       0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. - (1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.))/
        (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
          0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)),0.5));
}

void
func_Pi23(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*y*z*(-4.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
       Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
        (6.*Power(x,2) + 6.*Power(y,2) + 6.*Power(z,2))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)*
     (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
       0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. - (1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.))/
        (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
          0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)),0.5));
}

void
func_Pi33(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*Power(z,2)*(-4.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
       Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
        (6.*Power(x,2) + 6.*Power(y,2) + 6.*Power(z,2))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)*
     (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
       0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. - (1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.))/
        (1.*M*Power(x,2) + 1.*M*Power(y,2) + 1.*M*Power(z,2) + 
          0.5*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)),0.5));
}

//nonzero funcs of Phi
void
func_Phi100(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2.*M*x)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5);
}

void
func_Phi101(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -M*((-4.*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
     2/Power(Power(x,2) + Power(y,2) + Power(z,2),1.)); 
}

void
func_Phi102(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -(-4.*M*x*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.); 
}

void
func_Phi103(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -(-4.*M*x*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.); 
}

void
func_Phi111(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = M*((-6.*Power(x,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
     (4*x)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5));
} 

void
func_Phi112(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
   *value = M*y*((-6.*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
     2/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)); 
}

void
func_Phi113(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = M*z*((-6.*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
     2/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5));
}

void
func_Phi122(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-6.*M*x*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_Phi123(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-6.*M*x*y*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_Phi133(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-6.*M*x*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_Phi200(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2.*M*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5);
}

void
func_Phi201(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -(-4.*M*x*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
}

void
func_Phi202(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -M*((-4.*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
     2/Power(Power(x,2) + Power(y,2) + Power(z,2),1.));
}

void
func_Phi203(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -(-4.*M*y*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
}

void
func_Phi211(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-6.*M*Power(x,2)*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5); 
}

void
func_Phi212(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = M*x*((-6.*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
     2/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5));
}

void
func_Phi213(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-6.*M*x*y*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_Phi222(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = M*((-6.*Power(y,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
     (4*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5));
}

void
func_Phi223(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = M*z*((-6.*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
     2/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5));
}

void
func_Phi233(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-6.*M*y*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_Phi300(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2.*M*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5);
}

void
func_Phi301(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -(-4.*M*x*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
}

void
func_Phi302(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -(-4.*M*y*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
}

void
func_Phi303(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -M*((-4.*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
     2/Power(Power(x,2) + Power(y,2) + Power(z,2),1.));
}

void
func_Phi311(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-6.*M*Power(x,2)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_Phi312(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-6.*M*x*y*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_Phi313(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = M*x*((-6.*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
     2/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)); 
}

void
func_Phi322(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-6.*M*Power(y,2)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_Phi323(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = M*y*((-6.*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
     2/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5));
}

void
func_Phi333(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = M*((-6.*Power(z,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
     (4*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5));
}

/*functions of H_a*/
void
func_H0(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (2*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.);
}

void
func_H1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2*M*x)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5); 
}

void
func_H2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2*M*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5); 
}

void
func_H3(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2*M*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5); 
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
    *value = (-4.*M*x)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
}

void
func_deriH11(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (6.*M*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
   (2*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5);
}

void
func_deriH12(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (6.*M*x*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_deriH13(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (6.*M*x*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_deriH20(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4.*M*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
}

void
func_deriH21(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (6.*M*x*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_deriH22(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (6.*M*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
   (2*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5);
}

void
func_deriH23(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (6.*M*y*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_deriH30(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4.*M*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
}

void
func_deriH31(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (6.*M*x*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_deriH32(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (6.*M*y*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_deriH33(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (6.*M*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
   (2*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5);
}

