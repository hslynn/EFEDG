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
    *value = (2*M*x)/Pow(R, 2.);
}

void
func_Psi02(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (2*M*y)/Pow(R, 2.);
}

void
func_Psi03(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (2*M*z)/Pow(R, 2.);
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
    *value = (Power(M,2)*(-2.*Power(x,2) - 2.*Power(y,2) - 2.*Power(z,2)))/
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
    *value = (Power(M,2)*Power(x,2)*((-6.*Power(x,2) - 6.*Power(y,2) - 6.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
       4.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)))/
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
    *value = (Power(M,2)*x*y*((-6.*Power(x,2) - 6.*Power(y,2) - 6.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
       4.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)))/
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
    *value = (Power(M,2)*x*z*((-6.*Power(x,2) - 6.*Power(y,2) - 6.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
       4.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)))/
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
    *value = (Power(M,2)*Power(y,2)*((-6.*Power(x,2) - 6.*Power(y,2) - 6.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
       4.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)))/
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
    *value = (Power(M,2)*y*z*((-6.*Power(x,2) - 6.*Power(y,2) - 6.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
       4.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)))/
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
    *value = (Power(M,2)*Power(z,2)*((-6.*Power(x,2) - 6.*Power(y,2) - 6.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
       4.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)))/
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
    *value = M*((-4.*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
     2/Power(Power(x,2) + Power(y,2) + Power(z,2),1.)); 
}

void
func_Phi102(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4.*M*x*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.); 
}

void
func_Phi103(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4.*M*x*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.); 
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
    *value = (-4.*M*x*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
}

void
func_Phi202(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = M*((-4.*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
     2/Power(Power(x,2) + Power(y,2) + Power(z,2),1.));
}

void
func_Phi203(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4.*M*y*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
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
    *value = (-4.*M*x*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
}

void
func_Phi302(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4.*M*y*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
}

void
func_Phi303(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = M*((-4.*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
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
    *value = (-4*Power(M,2)*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
   (4*Power(M,2)*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
   (4*Power(M,2)*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
   (2*M*(-1. + (2*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),1.); 
}

void
func_H1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4*Power(M,2)*x*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) - 
   (4*Power(M,2)*x*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
   (4*Power(M,2)*x)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
   (2*M*x*(1 + (2*M*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),1.5); 
}

void
func_H2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4*Power(M,2)*Power(x,2)*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) - 
   (4*Power(M,2)*y*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
   (4*Power(M,2)*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
   (2*M*y*(1 + (2*M*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),1.5); 
}

void
func_H3(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4*Power(M,2)*Power(x,2)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) - 
   (4*Power(M,2)*Power(y,2)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
   (4*Power(M,2)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
   (2*M*z*(1 + (2*M*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),1.5); 
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
    *value = (20.*Power(M,2)*Power(x,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.5) + 
   (20.*Power(M,2)*x*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.5) + 
   (20.*Power(M,2)*x*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.5) - 
   (12.*Power(M,2)*x)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
   (4.*M*x*(-1. + (2*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
}

void
func_deriH11(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (24.*Power(M,2)*Power(x,2)*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (24.*Power(M,2)*Power(x,2)*Power(z,2))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),4.) - 
   (16.*Power(M,2)*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) - 
   (4*Power(M,2)*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) - 
   (4*Power(M,2)*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
   (4*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
   (2*M*x*((-6.*M*Power(x,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
        (4*M*x)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
   (6.*M*Power(x,2)*(1 + (2*M*Power(x,2))/
         Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
   (2*M*(1 + (2*M*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),1.5);
}

void
func_deriH12(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (24.*Power(M,2)*Power(x,3)*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (12.*Power(M,2)*x*Power(y,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (24.*Power(M,2)*x*y*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) - 
   (24.*Power(M,2)*x*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
   (6.*M*x*y*(1 + (2*M*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_deriH13(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (24.*Power(M,2)*Power(x,3)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (24.*Power(M,2)*x*Power(y,2)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (12.*Power(M,2)*x*Power(z,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) - 
   (24.*Power(M,2)*x*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
   (6.*M*x*z*(1 + (2*M*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_deriH20(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (20.*Power(M,2)*Power(x,2)*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),3.5) + 
   (20.*Power(M,2)*Power(y,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.5) + 
   (20.*Power(M,2)*y*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.5) - 
   (12.*Power(M,2)*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
   (4.*M*y*(-1. + (2*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
}

void
func_deriH21(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (12.*Power(M,2)*Power(x,3)*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (24.*Power(M,2)*x*Power(y,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (24.*Power(M,2)*x*y*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) - 
   (24.*Power(M,2)*x*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
   (6.*M*x*y*(1 + (2*M*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_deriH22(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (24.*Power(M,2)*Power(x,2)*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (24.*Power(M,2)*Power(y,2)*Power(z,2))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),4.) - 
   (4*Power(M,2)*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) - 
   (16.*Power(M,2)*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) - 
   (4*Power(M,2)*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
   (4*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
   (2*M*y*((-6.*M*Power(y,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
        (4*M*y)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
   (6.*M*Power(y,2)*(1 + (2*M*Power(y,2))/
         Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
   (2*M*(1 + (2*M*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),1.5);
}

void
func_deriH23(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (24.*Power(M,2)*Power(x,2)*y*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (24.*Power(M,2)*Power(y,3)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (12.*Power(M,2)*y*Power(z,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) - 
   (24.*Power(M,2)*y*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
   (6.*M*y*z*(1 + (2*M*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_deriH30(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (20.*Power(M,2)*Power(x,2)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),3.5) + 
   (20.*Power(M,2)*Power(y,2)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),3.5) + 
   (20.*Power(M,2)*Power(z,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.5) - 
   (12.*Power(M,2)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
   (4.*M*z*(-1. + (2*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),2.);
}

void
func_deriH31(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (12.*Power(M,2)*Power(x,3)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (24.*Power(M,2)*x*Power(y,2)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (24.*Power(M,2)*x*Power(z,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) - 
   (24.*Power(M,2)*x*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
   (6.*M*x*z*(1 + (2*M*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_deriH32(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (24.*Power(M,2)*Power(x,2)*y*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (12.*Power(M,2)*Power(y,3)*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (24.*Power(M,2)*y*Power(z,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) - 
   (24.*Power(M,2)*y*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
   (6.*M*y*z*(1 + (2*M*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),2.5);
}

void
func_deriH33(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (24.*Power(M,2)*Power(x,2)*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),4.) + 
   (24.*Power(M,2)*Power(y,2)*Power(z,2))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),4.) - 
   (4*Power(M,2)*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) - 
   (4*Power(M,2)*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) - 
   (16.*Power(M,2)*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
   (4*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
   (2*M*z*((-6.*M*Power(z,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
        (4*M*z)/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
   (6.*M*Power(z,2)*(1 + (2*M*Power(z,2))/
         Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
   (2*M*(1 + (2*M*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
    Power(Power(x,2) + Power(y,2) + Power(z,2),1.5);
}

/*set dof data using above functions*/
void
set_data_var(DOF **dofs_var)
{
    int i;
    DOF_USER_FUNC funcs_Psi[10] = {func_Psi00, func_Psi01, func_Psi02, func_Psi03,
                                                    func_Psi11, func_Psi12, func_Psi13,
                                                                func_Psi22, func_Psi23,
                                                                            func_Psi33};
    DOF_USER_FUNC funcs_Pi[10] = {func_Pi00, func_Pi01, func_Pi02, func_Pi03,
                                             func_Pi11, func_Pi12, func_Pi13,
                                                        func_Pi22, func_Pi23,
                                                                   func_Pi33};
    DOF_USER_FUNC funcs_Phi[30] = {func_Phi100, func_Phi101, func_Phi102, func_Phi103,        
                                               func_Phi111, func_Phi112, func_Phi113,
                                                          func_Phi122, func_Phi123,
                                                                      func_Phi133,
                                   func_Phi200, func_Phi201, func_Phi202, func_Phi203,        
                                               func_Phi211, func_Phi212, func_Phi213,
                                                          func_Phi222, func_Phi223,
                                                                      func_Phi233,
                                   func_Phi300, func_Phi301, func_Phi302, func_Phi303,        
                                               func_Phi311, func_Phi312, func_Phi313,
                                                          func_Phi322, func_Phi323,
                                                                      func_Phi333};
    for(i=0; i<10; i++){
        phgDofSetDataByFunction(dofs_var[i], funcs_Psi[i]);
        phgDofSetDataByFunction(dofs_var[10 + i], funcs_Pi[i]);
        
        phgDofSetDataByFunction(dofs_var[20 + i], funcs_Phi[i]);
        phgDofSetDataByFunction(dofs_var[30 + i], funcs_Phi[10 + i]);
        phgDofSetDataByFunction(dofs_var[40 + i], funcs_Phi[20 + i]);
    }
}

void
set_data_H(DOF **dofs_H)
{
    phgPrintf("M:%lf\n", M);
    int i;
    DOF_USER_FUNC funcs_H[4] = {func_H0, func_H1, func_H2, func_H3};
    for(i=0;i<4;i++){
        phgDofSetDataByFunction(dofs_H[i], funcs_H[i]);
    }
}

void
set_data_deriH(DOF **dofs_deriH)
{
    int i;
    DOF_USER_FUNC funcs_deriH[16] = {func_deriH00, func_deriH01, func_deriH02, func_deriH03,
                                     func_deriH10, func_deriH11, func_deriH12, func_deriH13,
                                     func_deriH20, func_deriH21, func_deriH22, func_deriH23,
                                     func_deriH30, func_deriH31, func_deriH32, func_deriH33};
    for(i=0;i<16;i++){
        phgDofSetDataByFunction(dofs_deriH[i], funcs_deriH[i]);
    }
}
