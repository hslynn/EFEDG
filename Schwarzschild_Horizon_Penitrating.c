#define M 0.2
#define R (Pow(x*x + y*y + z*z, 0.5))

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
    *value = (Power(M,3)*(-0.7272727272727273*Power(x,2) - 0.7272727272727273*Power(y,2) - 
       0.7272727272727273*Power(z,2))*
     (M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.)*
     Power(1.*M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((1.*Power(M,4) + 4.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         6.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         1.*Power(Power(x,2) + Power(y,2) + Power(z,2),2.))/
       (15.*Power(M,4) + 26.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         16.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         6.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi01(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*x*(-4.363636363636363*Power(x,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),6.) - 
       4.363636363636363*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) - 
       4.363636363636363*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
       1.4545454545454546*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
       Power(M,2)*(-1.4545454545454546*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
          1.4545454545454546*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
          1.4545454545454546*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
          1.4545454545454546*Power(Power(x,2) + Power(y,2) + Power(z,2),6.)) + 
       M*(-5.818181818181818*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
          5.818181818181818*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
          5.818181818181818*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
          2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),7.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((1.*Power(M,4) + 4.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         6.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         1.*Power(Power(x,2) + Power(y,2) + Power(z,2),2.))/
       (15.*Power(M,4) + 26.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         16.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         6.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi02(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*y*(-4.363636363636363*Power(x,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),6.) - 
       4.363636363636363*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) - 
       4.363636363636363*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
       1.4545454545454546*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
       Power(M,2)*(-1.4545454545454546*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
          1.4545454545454546*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
          1.4545454545454546*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
          1.4545454545454546*Power(Power(x,2) + Power(y,2) + Power(z,2),6.)) + 
       M*(-5.818181818181818*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
          5.818181818181818*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
          5.818181818181818*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
          2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),7.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((1.*Power(M,4) + 4.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         6.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         1.*Power(Power(x,2) + Power(y,2) + Power(z,2),2.))/
       (15.*Power(M,4) + 26.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         16.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         6.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi03(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*z*(-4.363636363636363*Power(x,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),6.) - 
       4.363636363636363*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) - 
       4.363636363636363*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
       1.4545454545454546*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
       Power(M,2)*(-1.4545454545454546*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
          1.4545454545454546*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
          1.4545454545454546*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
          1.4545454545454546*Power(Power(x,2) + Power(y,2) + Power(z,2),6.)) + 
       M*(-5.818181818181818*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
          5.818181818181818*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
          5.818181818181818*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
          2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),7.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((1.*Power(M,4) + 4.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         6.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         1.*Power(Power(x,2) + Power(y,2) + Power(z,2),2.))/
       (15.*Power(M,4) + 26.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         16.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         6.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi11(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*(-0.7272727272727273*Power(x,4)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       0.7272727272727273*Power(y,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       1.4545454545454546*Power(y,2)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       0.7272727272727273*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
       0.7272727272727273*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       0.7272727272727273*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       Power(x,2)*(-1.4545454545454546*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          1.4545454545454546*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
       Power(M,7)*(-1.4545454545454546*Power(y,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
          1.4545454545454546*Power(z,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          Power(x,2)*(-1.4545454545454546*Power(y,2) - 1.4545454545454546*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          0.7272727272727273*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          Power(y,2)*(-2.909090909090909*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5))) + 
       Power(M,6)*(-9.454545454545455*Power(y,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
          9.454545454545455*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
          Power(x,2)*(-9.454545454545455*Power(y,2) - 9.454545454545455*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
          5.090909090909091*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          Power(y,2)*(-18.90909090909091*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
             5.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),9.))) + 
       Power(M,5)*(-10.909090909090908*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          26.181818181818183*Power(y,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          26.181818181818183*Power(z,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          15.272727272727273*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          Power(x,2)*(-37.09090909090909*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             37.09090909090909*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             10.909090909090908*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
          Power(y,2)*(-52.36363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             15.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5))) + 
       Power(M,4)*(-42.18181818181818*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          40.*Power(y,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          40.*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          25.454545454545453*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          Power(y,2)*(-80.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             25.454545454545453*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
          Power(x,2)*(-82.18181818181819*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             82.18181818181819*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             29.818181818181817*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))) + 
       Power(M,3)*(-47.27272727272727*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          36.36363636363637*Power(y,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          36.36363636363637*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          25.454545454545453*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          Power(y,2)*(-72.72727272727273*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             25.454545454545453*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(x,2)*(-83.63636363636364*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             83.63636363636364*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             30.545454545454547*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5))) + 
       Power(M,2)*(-21.09090909090909*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          19.636363636363637*Power(y,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          19.636363636363637*Power(z,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          15.272727272727273*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          Power(y,2)*(-39.27272727272727*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             15.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(x,2)*(-40.72727272727273*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             40.72727272727273*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             16.*Power(Power(x,2) + Power(y,2) + Power(z,2),11.))) + 
       M*(-5.818181818181818*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          5.818181818181818*Power(y,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          5.818181818181818*Power(z,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          5.090909090909091*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          Power(y,2)*(-11.636363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             5.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(x,2)*(-11.636363636363637*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             11.636363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             5.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((1.*Power(M,4) + 4.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         6.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         1.*Power(Power(x,2) + Power(y,2) + Power(z,2),2.))/
       (15.*Power(M,4) + 26.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         16.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         6.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi12(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*x*y*(-1.4545454545454546*Power(x,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       1.4545454545454546*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       1.4545454545454546*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
       0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       Power(M,5)*(1.4545454545454546*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          1.4545454545454546*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          1.4545454545454546*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
       Power(M,4)*(9.454545454545455*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          9.454545454545455*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          9.454545454545455*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          5.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
       Power(M,3)*(15.272727272727273*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          15.272727272727273*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          15.272727272727273*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          4.363636363636363*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
       Power(M,2)*(-2.1818181818181817*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          2.1818181818181817*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          2.1818181818181817*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          4.363636363636363*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
       M*(-10.90909090909091*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          10.90909090909091*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          10.90909090909091*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          5.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),12.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((1.*Power(M,4) + 4.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         6.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         1.*Power(Power(x,2) + Power(y,2) + Power(z,2),2.))/
       (15.*Power(M,4) + 26.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         16.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         6.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi13(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*x*z*(-1.4545454545454546*Power(x,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       1.4545454545454546*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       1.4545454545454546*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
       0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       Power(M,5)*(1.4545454545454546*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          1.4545454545454546*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          1.4545454545454546*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
       Power(M,4)*(9.454545454545455*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          9.454545454545455*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          9.454545454545455*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          5.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
       Power(M,3)*(15.272727272727273*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          15.272727272727273*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          15.272727272727273*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          4.363636363636363*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
       Power(M,2)*(-2.1818181818181817*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          2.1818181818181817*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          2.1818181818181817*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          4.363636363636363*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
       M*(-10.90909090909091*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          10.90909090909091*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          10.90909090909091*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          5.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),12.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((1.*Power(M,4) + 4.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         6.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         1.*Power(Power(x,2) + Power(y,2) + Power(z,2),2.))/
       (15.*Power(M,4) + 26.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         16.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         6.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi22(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*(-0.7272727272727273*Power(x,4)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       0.7272727272727273*Power(y,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       1.4545454545454546*Power(y,2)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       0.7272727272727273*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
       0.7272727272727273*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       0.7272727272727273*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       Power(x,2)*(-1.4545454545454546*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          1.4545454545454546*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
       Power(M,7)*(-1.4545454545454546*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          Power(x,2)*(-1.4545454545454546*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
             2.909090909090909*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)) + 
          Power(z,2)*(-1.4545454545454546*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
             1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5))) + 
       Power(M,6)*(-9.454545454545455*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
          Power(x,2)*(-9.454545454545455*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
             18.90909090909091*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
             5.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)) + 
          Power(z,2)*(-9.454545454545455*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
             9.454545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
             5.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),9.))) + 
       Power(M,5)*(-26.18181818181818*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          10.90909090909091*Power(y,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          26.18181818181818*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          15.272727272727273*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          Power(y,2)*(-37.090909090909086*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             10.90909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
          Power(x,2)*(-37.090909090909086*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             52.36363636363636*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             15.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5))) + 
       Power(M,4)*(-40.*Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          42.18181818181818*Power(y,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          40.*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          25.454545454545457*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          Power(x,2)*(-82.18181818181819*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             80.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             25.454545454545457*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
          Power(y,2)*(-82.18181818181819*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             29.81818181818182*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))) + 
       Power(M,3)*(-36.36363636363637*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          47.27272727272727*Power(y,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          36.36363636363637*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          25.454545454545457*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          Power(x,2)*(-83.63636363636364*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             72.72727272727273*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             25.454545454545457*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(y,2)*(-83.63636363636364*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             30.545454545454547*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5))) + 
       Power(M,2)*(-19.63636363636364*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          21.09090909090909*Power(y,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          19.63636363636364*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          15.272727272727273*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          Power(x,2)*(-40.72727272727273*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             39.27272727272728*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             15.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(y,2)*(-40.72727272727273*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             16.*Power(Power(x,2) + Power(y,2) + Power(z,2),11.))) + 
       M*(-5.818181818181818*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          5.818181818181818*Power(y,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          5.818181818181818*Power(z,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          5.090909090909091*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          Power(y,2)*(-11.636363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             5.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(x,2)*(-11.636363636363637*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             11.636363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             5.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((1.*Power(M,4) + 4.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         6.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         1.*Power(Power(x,2) + Power(y,2) + Power(z,2),2.))/
       (15.*Power(M,4) + 26.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         16.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         6.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi23(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*y*z*(-1.4545454545454546*Power(x,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       1.4545454545454546*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       1.4545454545454546*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
       0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       Power(M,5)*(1.4545454545454546*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          1.4545454545454546*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          1.4545454545454546*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
       Power(M,4)*(9.454545454545455*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          9.454545454545455*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          9.454545454545455*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          5.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
       Power(M,3)*(15.272727272727273*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          15.272727272727273*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          15.272727272727273*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          4.363636363636363*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
       Power(M,2)*(-2.1818181818181817*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          2.1818181818181817*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          2.1818181818181817*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          4.363636363636363*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
       M*(-10.90909090909091*Power(x,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          10.90909090909091*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          10.90909090909091*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          5.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),12.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((1.*Power(M,4) + 4.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         6.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         1.*Power(Power(x,2) + Power(y,2) + Power(z,2),2.))/
       (15.*Power(M,4) + 26.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         16.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         6.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
}

static void
func_Pi33(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*(-0.7272727272727274*Power(x,4)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       0.7272727272727274*Power(y,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       1.4545454545454548*Power(y,2)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
       0.7272727272727274*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
       0.7272727272727274*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       0.7272727272727274*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
       Power(x,2)*(-1.4545454545454548*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          1.4545454545454548*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          0.7272727272727274*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
       Power(M,7)*(-1.4545454545454548*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          Power(x,2)*(-2.9090909090909096*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
             1.4545454545454548*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
             0.7272727272727274*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)) + 
          Power(y,2)*(-1.4545454545454548*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
             1.4545454545454548*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
             0.7272727272727274*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5))) + 
       Power(M,6)*(-9.454545454545455*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
          Power(x,2)*(-18.90909090909091*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
             9.454545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
             5.090909090909092*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)) + 
          Power(y,2)*(-9.454545454545455*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
             9.454545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
             5.090909090909092*Power(Power(x,2) + Power(y,2) + Power(z,2),9.))) + 
       Power(M,5)*(-26.181818181818183*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          26.181818181818183*Power(y,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          10.90909090909091*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          10.90909090909091*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          Power(y,2)*(-37.09090909090909*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             15.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
          Power(x,2)*(-52.36363636363637*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             37.09090909090909*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             15.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5))) + 
       Power(M,4)*(-40.00000000000001*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          40.00000000000001*Power(y,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          42.18181818181819*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          29.818181818181824*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          Power(y,2)*(-82.18181818181819*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             25.45454545454546*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
          Power(x,2)*(-80.00000000000001*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             82.18181818181819*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             25.45454545454546*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))) + 
       Power(M,3)*(-36.363636363636374*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          36.363636363636374*Power(y,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          47.27272727272727*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          30.545454545454547*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          Power(y,2)*(-83.63636363636364*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             25.45454545454546*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(x,2)*(-72.72727272727275*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             83.63636363636364*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             25.45454545454546*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5))) + 
       Power(M,2)*(-19.636363636363637*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          19.636363636363637*Power(y,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          21.09090909090909*Power(z,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          16.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          Power(y,2)*(-40.727272727272734*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             15.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(x,2)*(-39.27272727272727*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             40.727272727272734*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             15.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),11.))) + 
       M*(-5.818181818181819*Power(x,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          5.818181818181819*Power(y,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          5.818181818181819*Power(z,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          5.090909090909092*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          Power(y,2)*(-11.636363636363638*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             5.090909090909092*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(x,2)*(-11.636363636363638*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             11.636363636363638*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             5.090909090909092*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((1.*Power(M,4) + 4.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         6.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         4.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         1.*Power(Power(x,2) + Power(y,2) + Power(z,2),2.))/
       (15.*Power(M,4) + 26.*Power(M,3)*
          Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
         16.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
         6.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
         Power(Power(x,2) + Power(y,2) + Power(z,2),2.)),0.5));
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
    *value = (-4.*Power(M,2)*x*y*(1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       3.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3)); 
}

static void
func_Phi103(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4.*Power(M,2)*x*z*(1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       3.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3)); 
}

static void
func_Phi111(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (x*(-2.*Power(x,2) - 2.*Power(y,2) - 2.*Power(z,2) + 
       2*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
       (4.*Power(M,6)*(1.*Power(y,2) + 1.*Power(z,2)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) - 
       (22.*Power(M,5)*(1.*Power(y,2) + 1.*Power(z,2)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
       Power(M,4)*((-30.*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
          (50.*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
          (50.*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
          30/Power(Power(x,2) + Power(y,2) + Power(z,2),1.)) + 
       Power(M,3)*((-86.*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (60.*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (60.*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          52/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (12.*M*(-1.1666666666666667*Power(x,2) - 1.1666666666666667*Power(y,2) - 
            1.1666666666666667*Power(z,2) + 
            Power(Power(x,2) + Power(y,2) + Power(z,2),1.)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       (32.*Power(M,2)*(-1.375*Power(x,2) - 1.25*Power(y,2) - 1.25*Power(z,2) + 
            1.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),1.)))/
   Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4);
} 

static void
func_Phi112(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
   *value = (Power(M,2)*y*(20.*Power(M,2)*Power(x,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
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
    *value = (Power(M,2)*z*(20.*Power(M,2)*Power(x,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
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
       Power(y,2)*(3. + (11.*Power(M,3))/
           Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          (21.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          (13.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (2.*(Power(x,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),5))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
       (4.*(Power(x,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
       (2*Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
       (2.*Power(y,2)*(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5))*
          (15*Power(M,3) + 11*Power(M,2)*
             Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
            5*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
            Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
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
       Power(z,2)*(3. + (11.*Power(M,3))/
           Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          (21.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          (13.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (2.*(Power(x,2) + Power(y,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),5))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
       (4.*(Power(x,2) + Power(y,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
       (2*Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
       (2.*Power(z,2)*(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5))*
          (15*Power(M,3) + 11*Power(M,2)*
             Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
            5*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
            Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
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
    *value = (-4.*Power(M,2)*x*y*(1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       3.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
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
    *value = (-4.*Power(M,2)*y*z*(1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       3.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3));
}

static void
func_Phi211(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (y*(Power(x,2)*(-3. - (45.*Power(M,3))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (33.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
          (15.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       Power(x,2)*(3. + (11.*Power(M,3))/
           Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          (21.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          (13.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (2.*(Power(y,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),5))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
       (4.*(Power(y,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
       (2*Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
       (2.*Power(x,2)*(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5))*
          (15*Power(M,3) + 11*Power(M,2)*
             Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
            5*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
            Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.)))/
   Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4); 
}

static void
func_Phi212(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*x*(20.*Power(M,2)*Power(y,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
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
    *value = (y*(-2.*Power(x,2) - 2.*Power(y,2) - 2.*Power(z,2) + 
       2*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) - 
       (4.*Power(M,6)*(1.*Power(x,2) + 1.*Power(z,2)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) - 
       (22.*Power(M,5)*(1.*Power(x,2) + 1.*Power(z,2)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
       Power(M,4)*((-50.*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
          (30.*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
          (50.*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
          30/Power(Power(x,2) + Power(y,2) + Power(z,2),1.)) + 
       Power(M,3)*((-60.*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (86.*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (60.*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          52/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (12.*M*(-1.1666666666666667*Power(x,2) - 1.1666666666666667*Power(y,2) - 
            1.1666666666666667*Power(z,2) + 
            Power(Power(x,2) + Power(y,2) + Power(z,2),1.)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       (32.*Power(M,2)*(-1.25*Power(x,2) - 1.375*Power(y,2) - 1.25*Power(z,2) + 
            1.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),1.)))/
   Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4);
}

static void
func_Phi223(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*z*(20.*Power(M,2)*Power(y,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
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
       Power(z,2)*(3. + (11.*Power(M,3))/
           Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          (21.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          (13.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (2.*(Power(x,2) + Power(y,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),5))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
       (4.*(Power(x,2) + Power(y,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
       (2*Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
       (2.*Power(z,2)*(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5))*
          (15*Power(M,3) + 11*Power(M,2)*
             Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
            5*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
            Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
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
    *value = (-4.*Power(M,2)*x*z*(1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       3.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3));
}

static void
func_Phi302(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4.*Power(M,2)*y*z*(1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       3.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
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
       Power(x,2)*(3. + (11.*Power(M,3))/
           Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          (21.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          (13.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (2.*(Power(y,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),5))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
       (4.*(Power(y,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
       (2*Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
       (2.*Power(x,2)*(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5))*
          (15*Power(M,3) + 11*Power(M,2)*
             Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
            5*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
            Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
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
    *value = (Power(M,2)*x*(20.*Power(M,2)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
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
       Power(y,2)*(3. + (11.*Power(M,3))/
           Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          (21.*Power(M,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
          (13.*M)/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (2.*(Power(x,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),5))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) - 
       (4.*(Power(x,2) + Power(z,2))*
          Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) + 
       (2*Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
       (2.*Power(y,2)*(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5))*
          (15*Power(M,3) + 11*Power(M,2)*
             Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
            5*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
            Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.)))/
   Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4);
}

static void
func_Phi323(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*y*(20.*Power(M,2)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
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
    *value = (z*(-2.*Power(x,2) - 2.*Power(y,2) - 2.*Power(z,2) - 
       (4.*Power(M,6)*(1.*Power(x,2) + 1.*Power(y,2)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),3.) - 
       (22.*Power(M,5)*(1.*Power(x,2) + 1.*Power(y,2)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
       2*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       Power(M,4)*((-50.*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
          (50.*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) - 
          (30.*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
          30/Power(Power(x,2) + Power(y,2) + Power(z,2),1.)) + 
       Power(M,3)*((-60.*Power(x,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (60.*Power(y,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
          (86.*Power(z,2))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
          52/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)) + 
       (12.*M*(-1.1666666666666667*Power(x,2) - 1.1666666666666667*Power(y,2) - 
            1.1666666666666667*Power(z,2) + 
            Power(Power(x,2) + Power(y,2) + Power(z,2),1.)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       (32.*Power(M,2)*(-1.25*Power(x,2) - 1.25*Power(y,2) - 1.375*Power(z,2) + 
            1.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)))/
        Power(Power(x,2) + Power(y,2) + Power(z,2),1.)))/
   Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),4);
}

/*set dof data using above functions*/
static void
set_data_dofs(DOF **dofs_var)
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

