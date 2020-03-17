#include "phg.h"
#include "global_def.h"

static void
func_Psi00(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -(1.0-M/R)/(1.0+M/R);
}

static void
func_Psi01(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 4*Power(M, 2.0)*x/Power(R+M, 2.0)/R;
}

static void
func_Psi02(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 4*Power(M, 2.0)*y/Power(R+M, 2.0)/R;
}

static void
func_Psi03(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 4*Power(M, 2.0)*z/Power(R+M, 2.0)/R;
}

static void
func_Psi11(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = ((15*Power(M,3) + 11*Power(M,2)*R + 5*M*Power(R,2) + Power(R,3))*Power(x,2))/
    (Power(R,2)*Power(M + R,3)) + (Power(M + R,2)*(Power(y,2) + Power(z,2)))/Power(R,4);
}

static void
func_Psi12(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -((Power(M + R,2)*x*y)/Power(R,4)) + 
   ((15*Power(M,3) + 11*Power(M,2)*R + 5*M*Power(R,2) + Power(R,3))*x*y)/
    (Power(R,2)*Power(M + R,3)); 
}

static void
func_Psi13(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -((Power(M + R,2)*x*z)/Power(R,4)) + 
   ((15*Power(M,3) + 11*Power(M,2)*R + 5*M*Power(R,2) + Power(R,3))*x*z)/
    (Power(R,2)*Power(M + R,3));  
}


static void
func_Psi22(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = ((15*Power(M,3) + 11*Power(M,2)*R + 5*M*Power(R,2) + Power(R,3))*Power(y,2))/
    (Power(R,2)*Power(M + R,3)) + (Power(M + R,2)*(Power(x,2) + Power(z,2)))/Power(R,4);
}

static void
func_Psi23(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -((Power(M + R,2)*y*z)/Power(R,4)) + 
   ((15*Power(M,3) + 11*Power(M,2)*R + 5*M*Power(R,2) + Power(R,3))*y*z)/
    (Power(R,2)*Power(M + R,3));
}

static void
func_Psi33(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M + R,2)*(Power(x,2) + Power(y,2)))/Power(R,4) + 
   ((15*Power(M,3) + 11*Power(M,2)*R + 5*M*Power(R,2) + Power(R,3))*Power(z,2))/
    (Power(R,2)*Power(M + R,3));
}

//nonzero funcs of Pi
static void
func_Pi00(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-0.5333333333333334*Power(M,3)*(1.*Power(x,2) + 1.*Power(y,2) + 1.*Power(z,2))*
     (M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.)*
     Power(1.*M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)*
     (1.*Power(M,3) + 0.7333333333333334*Power(M,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.33333333333333337*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.06666666666666668*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. + (M*(-2.333333333333333*Power(M,3) - 
            3.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) - 
            1.6666666666666665*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
            0.3333333333333333*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
        (2.5*Power(M,4) + 4.333333333333333*Power(M,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
          2.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          0.16666666666666666*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi01(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*x*(M*(-4.2666666666666675*Power(x,2) - 4.2666666666666675*Power(y,2) - 
          4.2666666666666675*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
       (-3.2*Power(x,2) - 3.2*Power(y,2) - 3.2*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
       2.1333333333333337*M*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
       1.0666666666666669*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
       Power(M,2)*((-1.0666666666666669*Power(x,2) - 1.0666666666666669*Power(y,2) - 
             1.0666666666666669*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)\
           + 1.0666666666666669*Power(Power(x,2) + Power(y,2) + Power(z,2),6.))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),7.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3)*
     (1.*Power(M,3) + 0.7333333333333334*Power(M,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.33333333333333337*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.06666666666666668*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. + (M*(-2.333333333333333*Power(M,3) - 
            3.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) - 
            1.6666666666666665*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
            0.3333333333333333*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
        (2.5*Power(M,4) + 4.333333333333333*Power(M,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
          2.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          0.16666666666666666*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi02(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*y*(M*(-4.2666666666666675*Power(x,2) - 4.2666666666666675*Power(y,2) - 
          4.2666666666666675*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
       (-3.2*Power(x,2) - 3.2*Power(y,2) - 3.2*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
       2.1333333333333337*M*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
       1.0666666666666669*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
       Power(M,2)*((-1.0666666666666669*Power(x,2) - 1.0666666666666669*Power(y,2) - 
             1.0666666666666669*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)\
           + 1.0666666666666669*Power(Power(x,2) + Power(y,2) + Power(z,2),6.))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),7.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3)*
     (1.*Power(M,3) + 0.7333333333333334*Power(M,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.33333333333333337*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.06666666666666668*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. + (M*(-2.333333333333333*Power(M,3) - 
            3.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) - 
            1.6666666666666665*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
            0.3333333333333333*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
        (2.5*Power(M,4) + 4.333333333333333*Power(M,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
          2.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          0.16666666666666666*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi03(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*z*(M*(-4.2666666666666675*Power(x,2) - 4.2666666666666675*Power(y,2) - 
          4.2666666666666675*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
       (-3.2*Power(x,2) - 3.2*Power(y,2) - 3.2*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
       2.1333333333333337*M*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
       1.0666666666666669*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
       Power(M,2)*((-1.0666666666666669*Power(x,2) - 1.0666666666666669*Power(y,2) - 
             1.0666666666666669*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)\
           + 1.0666666666666669*Power(Power(x,2) + Power(y,2) + Power(z,2),6.))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),7.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3)*
     (1.*Power(M,3) + 0.7333333333333334*Power(M,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.33333333333333337*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.06666666666666668*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. + (M*(-2.333333333333333*Power(M,3) - 
            3.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) - 
            1.6666666666666665*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
            0.3333333333333333*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
        (2.5*Power(M,4) + 4.333333333333333*Power(M,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
          2.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          0.16666666666666666*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi11(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*((0.5333333333333334*Power(x,2) + 0.5333333333333334*Power(y,2) + 
          0.5333333333333334*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       M*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)*
        (3.7333333333333334*Power(x,2) + 3.7333333333333334*Power(y,2) + 
          3.7333333333333334*Power(z,2)) + 
       M*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
        (-4.2666666666666675*Power(x,4) - 4.2666666666666675*Power(y,4) - 
          8.533333333333335*Power(y,2)*Power(z,2) - 4.2666666666666675*Power(z,4) + 
          Power(x,2)*(-8.533333333333335*Power(y,2) - 8.533333333333335*Power(z,2))) + 
       Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
        (-0.5333333333333334*Power(x,4) - 0.5333333333333334*Power(y,4) - 
          1.0666666666666669*Power(y,2)*Power(z,2) - 0.5333333333333334*Power(z,4) + 
          Power(x,2)*(-1.0666666666666669*Power(y,2) - 1.0666666666666669*Power(z,2))) + 
       Power(M,3)*(Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
           (22.400000000000002*Power(x,2) + 18.666666666666668*Power(y,2) + 
             18.666666666666668*Power(z,2)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*
           (-34.66666666666667*Power(x,4) - 26.66666666666667*Power(y,4) - 
             53.33333333333334*Power(y,2)*Power(z,2) - 26.66666666666667*Power(z,4) + 
             Power(x,2)*(-61.33333333333334*Power(y,2) - 61.33333333333334*Power(z,2)))) + 
       Power(M,4)*(Power(Power(x,2) + Power(y,2) + Power(z,2),10.)*
           (21.866666666666667*Power(x,2) + 18.666666666666668*Power(y,2) + 
             18.666666666666668*Power(z,2)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (-30.933333333333334*Power(x,4) - 29.333333333333336*Power(y,4) - 
             58.66666666666667*Power(y,2)*Power(z,2) - 29.333333333333336*Power(z,4) + 
             Power(x,2)*(-60.26666666666667*Power(y,2) - 60.26666666666667*Power(z,2)))) + 
       Power(M,2)*(Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
           (11.733333333333334*Power(x,2) + 11.200000000000001*Power(y,2) + 
             11.200000000000001*Power(z,2)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),10.)*
           (-15.466666666666667*Power(x,4) - 14.400000000000002*Power(y,4) - 
             28.800000000000004*Power(y,2)*Power(z,2) - 14.400000000000002*Power(z,4) + 
             Power(x,2)*(-29.866666666666667*Power(y,2) - 29.866666666666667*Power(z,2))))\
        + Power(M,5)*(Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*
           (8.*Power(x,2) + 11.200000000000001*Power(y,2) + 11.200000000000001*Power(z,2))\
           + Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*
           (-8.*Power(x,4) - 19.200000000000003*Power(y,4) - 
             38.400000000000006*Power(y,2)*Power(z,2) - 19.200000000000003*Power(z,4) + 
             Power(x,2)*(-27.200000000000003*Power(y,2) - 27.200000000000003*Power(z,2))))\
        + Power(M,6)*(Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (3.7333333333333334*Power(y,2) + 3.7333333333333334*Power(z,2)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),8.)*
           (-6.9333333333333345*Power(y,4) - 13.866666666666669*Power(y,2)*Power(z,2) - 
             6.9333333333333345*Power(z,4) + 
             Power(x,2)*(-6.9333333333333345*Power(y,2) - 6.9333333333333345*Power(z,2))))\
        + Power(M,7)*((0.5333333333333334*Power(y,2) + 0.5333333333333334*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)*
           (-1.0666666666666669*Power(y,4) - 2.1333333333333337*Power(y,2)*Power(z,2) - 
             1.0666666666666669*Power(z,4) + 
             Power(x,2)*(-1.0666666666666669*Power(y,2) - 1.0666666666666669*Power(z,2))))))
    /(Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)*
     (1.*Power(M,3) + 0.7333333333333334*Power(M,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.33333333333333337*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.06666666666666668*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. + (M*(-2.333333333333333*Power(M,3) - 
            3.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) - 
            1.6666666666666665*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
            0.3333333333333333*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
        (2.5*Power(M,4) + 4.333333333333333*Power(M,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
          2.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          0.16666666666666666*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi12(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*x*y*(M*(-8.000000000000002*Power(x,2) - 8.000000000000002*Power(y,2) - 
          8.000000000000002*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
       (-1.0666666666666669*Power(x,2) - 1.0666666666666669*Power(y,2) - 
          1.0666666666666669*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
       3.7333333333333343*M*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
       0.5333333333333334*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       Power(M,2)*((-1.6*Power(x,2) - 1.6*Power(y,2) - 1.6*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          3.2*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
       Power(M,5)*(-0.5333333333333334*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*
           (1.0666666666666669*Power(x,2) + 1.0666666666666669*Power(y,2) + 
             1.0666666666666669*Power(z,2))) + 
       Power(M,4)*(-3.7333333333333343*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (6.9333333333333345*Power(x,2) + 6.9333333333333345*Power(y,2) + 
             6.9333333333333345*Power(z,2))) + 
       Power(M,3)*(-3.2*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*
           (11.200000000000003*Power(x,2) + 11.200000000000003*Power(y,2) + 
             11.200000000000003*Power(z,2)))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),12.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)*
     (1.*Power(M,3) + 0.7333333333333334*Power(M,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.33333333333333337*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.06666666666666668*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. + (M*(-2.333333333333333*Power(M,3) - 
            3.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) - 
            1.6666666666666665*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
            0.3333333333333333*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
        (2.5*Power(M,4) + 4.333333333333333*Power(M,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
          2.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          0.16666666666666666*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi13(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*x*z*(M*(-8.000000000000002*Power(x,2) - 8.000000000000002*Power(y,2) - 
          8.000000000000002*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
       (-1.0666666666666669*Power(x,2) - 1.0666666666666669*Power(y,2) - 
          1.0666666666666669*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
       3.7333333333333343*M*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
       0.5333333333333334*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       Power(M,2)*((-1.6*Power(x,2) - 1.6*Power(y,2) - 1.6*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          3.2*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
       Power(M,5)*(-0.5333333333333334*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*
           (1.0666666666666669*Power(x,2) + 1.0666666666666669*Power(y,2) + 
             1.0666666666666669*Power(z,2))) + 
       Power(M,4)*(-3.7333333333333343*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (6.9333333333333345*Power(x,2) + 6.9333333333333345*Power(y,2) + 
             6.9333333333333345*Power(z,2))) + 
       Power(M,3)*(-3.2*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*
           (11.200000000000003*Power(x,2) + 11.200000000000003*Power(y,2) + 
             11.200000000000003*Power(z,2)))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),12.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)*
     (1.*Power(M,3) + 0.7333333333333334*Power(M,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.33333333333333337*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.06666666666666668*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. + (M*(-2.333333333333333*Power(M,3) - 
            3.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) - 
            1.6666666666666665*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
            0.3333333333333333*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
        (2.5*Power(M,4) + 4.333333333333333*Power(M,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
          2.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          0.16666666666666666*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi22(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*((0.5333333333333333*Power(x,2) + 0.5333333333333333*Power(y,2) + 
          0.5333333333333333*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       M*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)*
        (3.7333333333333334*Power(x,2) + 3.7333333333333334*Power(y,2) + 
          3.7333333333333334*Power(z,2)) + 
       M*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
        (-4.266666666666667*Power(x,4) - 4.266666666666667*Power(y,4) - 
          8.533333333333333*Power(y,2)*Power(z,2) - 4.266666666666667*Power(z,4) + 
          Power(x,2)*(-8.533333333333333*Power(y,2) - 8.533333333333333*Power(z,2))) + 
       Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
        (-0.5333333333333333*Power(x,4) - 0.5333333333333333*Power(y,4) - 
          1.0666666666666667*Power(y,2)*Power(z,2) - 0.5333333333333333*Power(z,4) + 
          Power(x,2)*(-1.0666666666666667*Power(y,2) - 1.0666666666666667*Power(z,2))) + 
       Power(M,4)*(Power(Power(x,2) + Power(y,2) + Power(z,2),10.)*
           (18.666666666666668*Power(x,2) + 21.86666666666667*Power(y,2) + 
             18.666666666666668*Power(z,2)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (-29.333333333333336*Power(x,4) - 30.93333333333334*Power(y,4) - 
             60.26666666666667*Power(y,2)*Power(z,2) - 29.333333333333336*Power(z,4) + 
             Power(x,2)*(-60.26666666666667*Power(y,2) - 58.66666666666667*Power(z,2)))) + 
       Power(M,3)*(Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
           (18.666666666666668*Power(x,2) + 22.400000000000002*Power(y,2) + 
             18.666666666666668*Power(z,2)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*
           (-26.666666666666668*Power(x,4) - 34.66666666666667*Power(y,4) - 
             61.33333333333334*Power(y,2)*Power(z,2) - 26.666666666666668*Power(z,4) + 
             Power(x,2)*(-61.33333333333334*Power(y,2) - 53.333333333333336*Power(z,2)))) + 
       Power(M,5)*(Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*
           (11.200000000000001*Power(x,2) + 8.000000000000002*Power(y,2) + 
             11.200000000000001*Power(z,2)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*
           (-19.2*Power(x,4) - 8.000000000000002*Power(y,4) - 27.2*Power(y,2)*Power(z,2) - 
             19.2*Power(z,4) + Power(x,2)*(-27.2*Power(y,2) - 38.4*Power(z,2)))) + 
       Power(M,2)*(Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
           (11.200000000000001*Power(x,2) + 11.733333333333334*Power(y,2) + 
             11.200000000000001*Power(z,2)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),10.)*
           (-14.400000000000002*Power(x,4) - 15.466666666666667*Power(y,4) - 
             29.86666666666667*Power(y,2)*Power(z,2) - 14.400000000000002*Power(z,4) + 
             Power(x,2)*(-29.86666666666667*Power(y,2) - 28.800000000000004*Power(z,2)))) + 
       Power(M,6)*(Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (3.7333333333333334*Power(x,2) + 3.7333333333333334*Power(z,2)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),8.)*
           (-6.9333333333333345*Power(x,4) - 6.9333333333333345*Power(y,2)*Power(z,2) - 
             6.9333333333333345*Power(z,4) + 
             Power(x,2)*(-6.9333333333333345*Power(y,2) - 13.866666666666669*Power(z,2))))\
        + Power(M,7)*((0.5333333333333333*Power(x,2) + 0.5333333333333333*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)*
           (-1.0666666666666667*Power(x,4) - 1.0666666666666667*Power(y,2)*Power(z,2) - 
             1.0666666666666667*Power(z,4) + 
             Power(x,2)*(-1.0666666666666667*Power(y,2) - 2.1333333333333333*Power(z,2))))))
    /(Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)*
     (1.*Power(M,3) + 0.7333333333333334*Power(M,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.33333333333333337*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.06666666666666668*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. + (M*(-2.333333333333333*Power(M,3) - 
            3.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) - 
            1.6666666666666665*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
            0.3333333333333333*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
        (2.5*Power(M,4) + 4.333333333333333*Power(M,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
          2.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          0.16666666666666666*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi23(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*y*z*(M*(-8.000000000000002*Power(x,2) - 8.000000000000002*Power(y,2) - 
          8.000000000000002*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
       (-1.0666666666666669*Power(x,2) - 1.0666666666666669*Power(y,2) - 
          1.0666666666666669*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
       3.7333333333333343*M*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
       0.5333333333333334*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       Power(M,2)*((-1.6*Power(x,2) - 1.6*Power(y,2) - 1.6*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          3.2*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
       Power(M,5)*(-0.5333333333333334*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*
           (1.0666666666666669*Power(x,2) + 1.0666666666666669*Power(y,2) + 
             1.0666666666666669*Power(z,2))) + 
       Power(M,4)*(-3.7333333333333343*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (6.9333333333333345*Power(x,2) + 6.9333333333333345*Power(y,2) + 
             6.9333333333333345*Power(z,2))) + 
       Power(M,3)*(-3.2*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*
           (11.200000000000003*Power(x,2) + 11.200000000000003*Power(y,2) + 
             11.200000000000003*Power(z,2)))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),12.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)*
     (1.*Power(M,3) + 0.7333333333333334*Power(M,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.33333333333333337*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.06666666666666668*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. + (M*(-2.333333333333333*Power(M,3) - 
            3.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) - 
            1.6666666666666665*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
            0.3333333333333333*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
        (2.5*Power(M,4) + 4.333333333333333*Power(M,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
          2.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          0.16666666666666666*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi33(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*((0.5333333333333334*Power(x,2) + 0.5333333333333334*Power(y,2) + 
          0.5333333333333334*Power(z,2))*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       M*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)*
        (3.7333333333333343*Power(x,2) + 3.7333333333333343*Power(y,2) + 
          3.7333333333333343*Power(z,2)) + 
       M*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
        (-4.2666666666666675*Power(x,4) - 4.2666666666666675*Power(y,4) - 
          8.533333333333335*Power(y,2)*Power(z,2) - 4.2666666666666675*Power(z,4) + 
          Power(x,2)*(-8.533333333333335*Power(y,2) - 8.533333333333335*Power(z,2))) + 
       Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
        (-0.5333333333333334*Power(x,4) - 0.5333333333333334*Power(y,4) - 
          1.0666666666666669*Power(y,2)*Power(z,2) - 0.5333333333333334*Power(z,4) + 
          Power(x,2)*(-1.0666666666666669*Power(y,2) - 1.0666666666666669*Power(z,2))) + 
       Power(M,3)*(Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
           (18.66666666666667*Power(x,2) + 18.66666666666667*Power(y,2) + 
             22.400000000000002*Power(z,2)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*
           (-26.666666666666675*Power(x,4) - 26.666666666666675*Power(y,4) - 
             61.33333333333335*Power(y,2)*Power(z,2) - 34.66666666666667*Power(z,4) + 
             Power(x,2)*(-53.33333333333335*Power(y,2) - 61.33333333333335*Power(z,2)))) + 
       Power(M,4)*(Power(Power(x,2) + Power(y,2) + Power(z,2),10.)*
           (18.66666666666667*Power(x,2) + 18.66666666666667*Power(y,2) + 
             21.86666666666667*Power(z,2)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (-29.333333333333343*Power(x,4) - 29.333333333333343*Power(y,4) - 
             60.26666666666668*Power(y,2)*Power(z,2) - 30.93333333333334*Power(z,4) + 
             Power(x,2)*(-58.666666666666686*Power(y,2) - 60.26666666666668*Power(z,2)))) + 
       Power(M,2)*(Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
           (11.200000000000001*Power(x,2) + 11.200000000000001*Power(y,2) + 
             11.733333333333334*Power(z,2)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),10.)*
           (-14.4*Power(x,4) - 14.4*Power(y,4) - 29.86666666666667*Power(y,2)*Power(z,2) - 
             15.466666666666669*Power(z,4) + 
             Power(x,2)*(-28.8*Power(y,2) - 29.86666666666667*Power(z,2)))) + 
       Power(M,5)*(Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*
           (11.200000000000001*Power(x,2) + 11.200000000000001*Power(y,2) + 
             8.000000000000002*Power(z,2)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*
           (-19.200000000000003*Power(x,4) - 19.200000000000003*Power(y,4) - 
             27.200000000000006*Power(y,2)*Power(z,2) - 8.000000000000002*Power(z,4) + 
             Power(x,2)*(-38.400000000000006*Power(y,2) - 27.200000000000006*Power(z,2))))\
        + Power(M,6)*((3.7333333333333343*Power(x,2) + 3.7333333333333343*Power(y,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),8.)*
           (-6.9333333333333345*Power(x,4) - 6.9333333333333345*Power(y,4) - 
             6.9333333333333345*Power(y,2)*Power(z,2) + 
             Power(x,2)*(-13.866666666666669*Power(y,2) - 6.9333333333333345*Power(z,2))))\
        + Power(M,7)*((0.5333333333333334*Power(x,2) + 0.5333333333333334*Power(y,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)*
           (-1.0666666666666669*Power(x,4) - 1.0666666666666669*Power(y,4) - 
             1.0666666666666669*Power(y,2)*Power(z,2) + 
             Power(x,2)*(-2.1333333333333337*Power(y,2) - 1.0666666666666669*Power(z,2))))))
    /(Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)*
     (1.*Power(M,3) + 0.7333333333333334*Power(M,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.33333333333333337*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.06666666666666668*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power(1. + (M*(-2.333333333333333*Power(M,3) - 
            3.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) - 
            1.6666666666666665*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
            0.3333333333333333*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
        (2.5*Power(M,4) + 4.333333333333333*Power(M,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
          2.6666666666666665*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          0.16666666666666666*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

//nonzero funcs of Phi
static void
func_Phi100(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2.*M*x)/(Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)*
     Power(M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi101(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*(-4.*M*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
       12.*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
       4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
       4.*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3)); 
}

static void
func_Phi102(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*x*y*(-4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
       12.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3)); 
}

static void
func_Phi103(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*x*z*(-4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
       12.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3)); 
}

static void
func_Phi111(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (x*(Power(M,6)*(-4.*Power(y,2) - 4.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
       Power(M,5)*(-22.*Power(y,2) - 22.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       M*(-14.*Power(x,2) - 14.*Power(y,2) - 14.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
       (-2.*Power(x,2) - 2.*Power(y,2) - 2.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
       12.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
       2.*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
       Power(M,4)*((-30.*Power(x,2) - 50.*Power(y,2) - 50.*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          30.*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
       Power(M,3)*((-86.*Power(x,2) - 60.*Power(y,2) - 60.*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          52.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
       Power(M,2)*((-44.*Power(x,2) - 40.*Power(y,2) - 40.*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          32.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4));
} 

static void
func_Phi112(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
   *value = (Power(M,2)*y*(20.*Power(M,2)*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
       4.*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
       1.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
       Power(M,4)*(4.*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
          1.*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)) + 
       Power(M,3)*(22.*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
          6.*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)) + 
       M*(-26.*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          6.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)); 
}

static void
func_Phi113(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*z*(20.*Power(M,2)*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
       4.*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
       1.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
       Power(M,4)*(4.*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
          1.*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)) + 
       Power(M,3)*(22.*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
          6.*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)) + 
       M*(-26.*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          6.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4));
}

static void
func_Phi122(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (x*(Power(y,2)*(-3. - (45.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (33.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
          (15.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       Power(y,2)*(-2. - (30.*Power(M,4))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
          (52.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (32.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
          (12.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       Power(y,2)*(3. + (11.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          (21.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          (13.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (2.*(Power(x,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),5))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
       (4.*(Power(x,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
       (2*Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.)))/
   Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4);
}

static void
func_Phi123(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*x*y*z*(4.*Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
       22.*Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       20.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
       26.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
       4.*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4));
}

static void
func_Phi133(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (x*(Power(z,2)*(-3. - (45.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (33.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
          (15.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       Power(z,2)*(-2. - (30.*Power(M,4))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
          (52.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (32.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
          (12.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       Power(z,2)*(3. + (11.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          (21.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          (13.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (2.*(Power(x,2) + Power(y,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),5))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
       (4.*(Power(x,2) + Power(y,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
       (2*Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.)))/
   Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4);
}

static void
func_Phi200(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2.*M*y)/(Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)*
     Power(M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi201(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*x*y*(-4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
       12.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3));
}

static void
func_Phi202(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*(-4.*M*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
       12.*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
       4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
       4.*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3));
}

static void
func_Phi203(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*y*z*(-4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
       12.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3));
}

static void
func_Phi211(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (y*(Power(x,2)*(-3. - (45.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (33.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
          (15.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       Power(x,2)*(-2. - (30.*Power(M,4))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
          (52.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (32.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
          (12.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       Power(x,2)*(3. + (11.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          (21.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          (13.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (2.*(Power(y,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),5))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
       (4.*(Power(y,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
       (2*Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.)))/
   Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4); 
}

static void
func_Phi212(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*x*(20.*Power(M,2)*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
       4.*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
       1.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
       Power(M,4)*(4.*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
          1.*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)) + 
       Power(M,3)*(22.*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
          6.*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)) + 
       M*(-26.*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          6.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4));
}

static void
func_Phi213(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*x*y*z*(4.*Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
       22.*Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       20.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
       26.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
       4.*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4));
}

static void
func_Phi222(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (y*(Power(M,6)*(-4.*Power(x,2) - 4.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
       Power(M,5)*(-22.*Power(x,2) - 22.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       M*(-14.*Power(x,2) - 14.*Power(y,2) - 14.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
       (-2.*Power(x,2) - 2.*Power(y,2) - 2.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
       12.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
       2.*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
       Power(M,4)*((-50.*Power(x,2) - 30.*Power(y,2) - 50.*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          30.*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
       Power(M,3)*((-60.*Power(x,2) - 86.*Power(y,2) - 60.*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          52.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
       Power(M,2)*((-40.*Power(x,2) - 44.*Power(y,2) - 40.*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          32.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4));
}

static void
func_Phi223(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*z*(20.*Power(M,2)*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
       4.*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
       1.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
       Power(M,4)*(4.*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
          1.*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)) + 
       Power(M,3)*(22.*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
          6.*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)) + 
       M*(-26.*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          6.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4));
}

static void
func_Phi233(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (y*(Power(z,2)*(-3. - (45.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (33.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
          (15.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       Power(z,2)*(-2. - (30.*Power(M,4))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
          (52.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (32.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
          (12.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       Power(z,2)*(3. + (11.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          (21.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          (13.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (2.*(Power(x,2) + Power(y,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),5))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
       (4.*(Power(x,2) + Power(y,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
       (2*Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.)))/
   Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4);
}

static void
func_Phi300(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2.*M*z)/(Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)*
     Power(M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi301(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*x*z*(-4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
       12.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3));
}

static void
func_Phi302(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*y*z*(-4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
       12.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3));
}

static void
func_Phi303(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*(-4.*M*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
       12.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
       4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
       4.*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3));
}

static void
func_Phi311(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (z*(Power(x,2)*(-3. - (45.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (33.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
          (15.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       Power(x,2)*(-2. - (30.*Power(M,4))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
          (52.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (32.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
          (12.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       Power(x,2)*(3. + (11.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          (21.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          (13.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (2.*(Power(y,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),5))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
       (4.*(Power(y,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
       (2*Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.)))/
   Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4);
}

static void
func_Phi312(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*x*y*z*(4.*Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
       22.*Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       20.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
       26.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
       4.*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4));
}

static void
func_Phi313(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*x*(20.*Power(M,2)*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
       4.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
       1.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
       Power(M,4)*(4.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
          1.*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)) + 
       Power(M,3)*(22.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
          6.*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)) + 
       M*(-26.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          6.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)); 
}

static void
func_Phi322(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (z*(Power(y,2)*(-3. - (45.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (33.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
          (15.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       Power(y,2)*(-2. - (30.*Power(M,4))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
          (52.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (32.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
          (12.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       Power(y,2)*(3. + (11.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          (21.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          (13.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (2.*(Power(x,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),5))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
       (4.*(Power(x,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
       (2*Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.)))/
   Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4);
}

static void
func_Phi323(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*y*(20.*Power(M,2)*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
       4.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
       1.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
       Power(M,4)*(4.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
          1.*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)) + 
       Power(M,3)*(22.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
          6.*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)) + 
       M*(-26.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          6.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4));
}

static void
func_Phi333(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (z*(Power(M,6)*(-4.*Power(x,2) - 4.*Power(y,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
       Power(M,5)*(-22.*Power(x,2) - 22.*Power(y,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       M*(-14.*Power(x,2) - 14.*Power(y,2) - 14.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
       (-2.*Power(x,2) - 2.*Power(y,2) - 2.*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
       12.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
       2.*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
       Power(M,4)*((-50.*Power(x,2) - 50.*Power(y,2) - 30.*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          30.*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
       Power(M,3)*((-60.*Power(x,2) - 60.*Power(y,2) - 86.*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          52.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
       Power(M,2)*((-40.*Power(x,2) - 40.*Power(y,2) - 44.*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          32.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4));
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
    int i;
    for(i=0; i<4; i++){
        phgDofSetDataByValue(dofs_H[i], 0);
    }
}

void
set_data_deriH(DOF **dofs_deriH)
{
    int i;
    for(i=0; i<16; i++){
        phgDofSetDataByValue(dofs_deriH[i], 0);
    }
}
