#define M 0.5
#define R (Pow(x*x + y*y + z*z, 0.5))

static void
func_Psi00(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -(1.0-M/R)/(1.0+M/R);
}

static void
func_Psi01(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 4*Power(M, 2.0)*x*x/Power(R+M, 2.0)/R;
}

static void
func_Psi02(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 4*Power(M, 2.0)*y*y/Power(R+M, 2.0)/R;
}

static void
func_Psi03(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 4*Power(M, 2.0)*z*z/Power(R+M, 2.0)/R;
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
    *value = (Power(M,3)*(Power(M,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)*
        (0.7272727272727274*Power(x,5) + 0.7272727272727274*Power(y,5) + 
          0.7272727272727274*Power(y,3)*Power(z,2) + 
          0.7272727272727274*Power(y,2)*Power(z,3) + 0.7272727272727274*Power(z,5) + 
          Power(x,3)*(0.7272727272727274*Power(y,2) + 0.7272727272727274*Power(z,2)) + 
          Power(x,2)*(0.7272727272727274*Power(y,3) + 0.7272727272727274*Power(z,3))) + 
       Power(Power(x,2) + Power(y,2) + Power(z,2),3.5)*
        (0.7272727272727274*Power(x,5) + 0.7272727272727274*Power(y,5) + 
          0.7272727272727274*Power(y,3)*Power(z,2) + 
          0.7272727272727274*Power(y,2)*Power(z,3) + 0.7272727272727274*Power(z,5) + 
          Power(x,3)*(0.7272727272727274*Power(y,2) + 0.7272727272727274*Power(z,2)) + 
          Power(x,2)*(0.7272727272727274*Power(y,3) + 0.7272727272727274*Power(z,3))) + 
       Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)*
        (3.6363636363636367*Power(x,5) + 3.6363636363636367*Power(y,5) + 
          3.6363636363636367*Power(y,3)*Power(z,2) + 
          3.6363636363636367*Power(y,2)*Power(z,3) + 3.6363636363636367*Power(z,5) + 
          Power(x,3)*(3.6363636363636367*Power(y,2) + 3.6363636363636367*Power(z,2)) + 
          Power(x,2)*(3.6363636363636367*Power(y,3) + 3.6363636363636367*Power(z,3))) + 
       M*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
        (3.6363636363636367*Power(x,5) + 3.6363636363636367*Power(y,5) + 
          3.6363636363636367*Power(y,3)*Power(z,2) + 
          3.6363636363636367*Power(y,2)*Power(z,3) + 3.6363636363636367*Power(z,5) + 
          Power(x,3)*(3.6363636363636367*Power(y,2) + 3.6363636363636367*Power(z,2)) + 
          Power(x,2)*(3.6363636363636367*Power(y,3) + 3.6363636363636367*Power(z,3))) + 
       Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)*
        (7.272727272727273*Power(x,5) + 7.272727272727273*Power(y,5) + 
          7.272727272727273*Power(y,3)*Power(z,2) + 
          7.272727272727273*Power(y,2)*Power(z,3) + 7.272727272727273*Power(z,5) + 
          Power(x,3)*(7.272727272727273*Power(y,2) + 7.272727272727273*Power(z,2)) + 
          Power(x,2)*(7.272727272727273*Power(y,3) + 7.272727272727273*Power(z,3))) + 
       Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
        (7.272727272727273*Power(x,5) + 7.272727272727273*Power(y,5) + 
          7.272727272727273*Power(y,3)*Power(z,2) + 
          7.272727272727273*Power(y,2)*Power(z,3) + 7.272727272727273*Power(z,5) + 
          Power(x,3)*(7.272727272727273*Power(y,2) + 7.272727272727273*Power(z,2)) + 
          Power(x,2)*(7.272727272727273*Power(y,3) + 7.272727272727273*Power(z,3)))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
          (0.09090909090909091*Power(x,4) + 0.09090909090909091*Power(y,4) + 
            0.18181818181818182*Power(y,2)*Power(z,2) + 0.09090909090909091*Power(z,4) + 
            Power(x,2)*(0.18181818181818182*Power(y,2) + 0.18181818181818182*Power(z,2)))
           + M*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)*
          (0.8181818181818181*Power(x,4) + 0.8181818181818181*Power(y,4) + 
            1.6363636363636362*Power(y,2)*Power(z,2) + 0.8181818181818181*Power(z,4) + 
            Power(x,2)*(1.6363636363636362*Power(y,2) + 1.6363636363636362*Power(z,2)))\
          + Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.5)*
          (3.2727272727272725*Power(x,4) + 3.2727272727272725*Power(y,4) + 
            6.545454545454545*Power(y,2)*Power(z,2) + 3.2727272727272725*Power(z,4) + 
            Power(x,2)*(6.545454545454545*Power(y,2) + 6.545454545454545*Power(z,2))) + 
         Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.)*
          (7.636363636363637*Power(x,4) + 7.636363636363637*Power(y,4) + 
            15.272727272727273*Power(y,2)*Power(z,2) + 7.636363636363637*Power(z,4) + 
            Power(x,2)*(15.272727272727273*Power(y,2) + 15.272727272727273*Power(z,2)))\
          + Power(M,9)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)*
          (-1.3636363636363635*Power(x,4) + 1.4545454545454546*Power(x,6) - 
            1.3636363636363635*Power(y,4) + 1.4545454545454546*Power(y,6) - 
            2.727272727272727*Power(y,2)*Power(z,2) + 
            2.909090909090909*Power(y,3)*Power(z,3) - 1.3636363636363635*Power(z,4) + 
            1.4545454545454546*Power(z,6) + 
            Power(x,2)*(-2.727272727272727*Power(y,2) - 2.727272727272727*Power(z,2)) + 
            Power(x,3)*(2.909090909090909*Power(y,3) + 2.909090909090909*Power(z,3))) + 
         Power(M,8)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)*
          (-6.454545454545454*Power(x,4) + 7.2727272727272725*Power(x,6) - 
            6.454545454545454*Power(y,4) + 7.2727272727272725*Power(y,6) - 
            12.909090909090908*Power(y,2)*Power(z,2) + 
            14.545454545454545*Power(y,3)*Power(z,3) - 6.454545454545454*Power(z,4) + 
            7.2727272727272725*Power(z,6) + 
            Power(x,2)*(-12.909090909090908*Power(y,2) - 
               12.909090909090908*Power(z,2)) + 
            Power(x,3)*(14.545454545454545*Power(y,3) + 14.545454545454545*Power(z,3)))\
          + Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.5)*
          (1.4545454545454546*Power(x,6) + 1.4545454545454546*Power(y,6) + 
            10.*Power(z,4) + 1.4545454545454546*Power(z,6) + 
            Power(y,4)*(10. + 1.4545454545454546*Power(z,2)) + 
            Power(x,4)*(10. + 1.4545454545454546*Power(y,2) + 
               1.4545454545454546*Power(z,2)) + 
            Power(y,2)*(20.*Power(z,2) + 1.4545454545454546*Power(z,4)) + 
            Power(x,2)*(20.*Power(y,2) + 1.4545454545454546*Power(y,4) + 
               20.*Power(z,2) + 1.4545454545454546*Power(z,4))) + 
         Power(M,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
          (7.2727272727272725*Power(x,6) + 7.2727272727272725*Power(y,6) + 
            4.181818181818182*Power(z,4) + 7.2727272727272725*Power(z,6) + 
            Power(y,4)*(4.181818181818182 + 7.2727272727272725*Power(z,2)) + 
            Power(x,4)*(4.181818181818182 + 7.2727272727272725*Power(y,2) + 
               7.2727272727272725*Power(z,2)) + 
            Power(y,2)*(8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4)) + 
            Power(x,2)*(8.363636363636363*Power(y,2) + 7.2727272727272725*Power(y,4) + 
               8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4))) + 
         Power(M,6)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            2.90909090909091*Power(y,3)*Power(z,3) - 6.909090909090909*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-6.909090909090909 + 16.*Power(z,2)) + 
            Power(x,4)*(-6.909090909090909 + 16.*Power(y,2) + 16.*Power(z,2)) + 
            Power(x,3)*(-2.90909090909091*Power(y,3) - 2.90909090909091*Power(z,3)) + 
            Power(y,2)*(-13.818181818181818*Power(z,2) + 16.*Power(z,4)) + 
            Power(x,2)*(-13.818181818181818*Power(y,2) + 16.*Power(y,4) - 
               13.818181818181818*Power(z,2) + 16.*Power(z,4))) + 
         Power(M,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            14.545454545454543*Power(y,3)*Power(z,3) - 11.272727272727273*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-11.272727272727273 + 21.818181818181817*Power(z,2)) + 
            Power(x,4)*(-11.272727272727273 + 21.818181818181817*Power(y,2) + 
               21.818181818181817*Power(z,2)) + 
            Power(x,3)*(-14.545454545454543*Power(y,3) - 
               14.545454545454543*Power(z,3)) + 
            Power(y,2)*(-22.545454545454547*Power(z,2) + 
               21.818181818181817*Power(z,4)) + 
            Power(x,2)*(-22.545454545454547*Power(y,2) + 21.818181818181817*Power(y,4) - 
               22.545454545454547*Power(z,2) + 21.818181818181817*Power(z,4))))/
       (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
         Power(1.*M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6)*
         (1.3636363636363635*Power(M,3) + 
           1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
           0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
           0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))),0.5));
}

static void
func_Pi01(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*Power(x,2)*(4.363636363636364*Power(x,5)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       x*(-2.909090909090909*Power(y,2) - 2.909090909090909*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
       Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.)*
        (4.363636363636364*Power(y,3) + 4.363636363636364*Power(z,3)) + 
       Power(Power(x,2) + Power(y,2) + Power(z,2),8.)*
        (4.363636363636364*Power(y,5) + 4.363636363636364*Power(y,3)*Power(z,2) + 
          4.363636363636364*Power(y,2)*Power(z,3) + 4.363636363636364*Power(z,5)) + 
       Power(x,3)*(4.363636363636364*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
          4.363636363636364*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
          2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)) + 
       Power(M,6)*(1.4545454545454546*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
          1.4545454545454546*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
          1.4545454545454546*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
          1.4545454545454546*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
          2.909090909090909*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)*
           (1.4545454545454546*Power(y,3) + 1.4545454545454546*Power(z,3)) + 
          Power(y,3)*(1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
             2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),6.)) + 
          Power(x,3)*(1.4545454545454546*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
             1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
             2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),6.))) + 
       Power(M,5)*(11.636363636363637*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
          11.636363636363637*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
          11.636363636363637*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
          11.636363636363637*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
          17.454545454545457*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
           (11.636363636363637*Power(y,3) + 11.636363636363637*Power(z,3)) + 
          Power(y,3)*(11.636363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
             17.454545454545457*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)) + 
          Power(x,3)*(11.636363636363637*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
             11.636363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
             17.454545454545457*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5))) + 
       Power(M,4)*(36.36363636363637*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
          x*(-43.63636363636364*Power(y,2) - 43.63636363636364*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.)*
           (36.36363636363637*Power(y,3) + 36.36363636363637*Power(z,3)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),6.)*
           (36.36363636363637*Power(y,5) + 36.36363636363637*Power(y,3)*Power(z,2) + 
             36.36363636363637*Power(y,2)*Power(z,3) + 36.36363636363637*Power(z,5)) + 
          Power(x,3)*(36.36363636363637*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
             36.36363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.) - 
             43.63636363636364*Power(Power(x,2) + Power(y,2) + Power(z,2),7.))) + 
       Power(M,3)*(58.18181818181819*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
          58.18181818181819*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
          58.18181818181819*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
          58.18181818181819*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
          17.454545454545457*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          x*(-75.63636363636364*Power(y,2) - 75.63636363636364*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)*
           (58.18181818181819*Power(y,3) + 58.18181818181819*Power(z,3)) + 
          Power(x,3)*(58.18181818181819*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
             58.18181818181819*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) - 
             58.18181818181819*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)) + 
          Power(y,3)*(58.18181818181819*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
             17.454545454545457*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5))) + 
       Power(M,2)*(50.909090909090914*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
          50.909090909090914*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
          50.909090909090914*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
          50.909090909090914*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
          2.909090909090909*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
          x*(-46.54545454545455*Power(y,2) - 46.54545454545455*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.)*
           (50.909090909090914*Power(y,3) + 50.909090909090914*Power(z,3)) + 
          Power(x,3)*(50.909090909090914*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
             50.909090909090914*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) - 
             43.63636363636364*Power(Power(x,2) + Power(y,2) + Power(z,2),8.)) + 
          Power(y,3)*(50.909090909090914*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
             2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),8.))) + 
       M*(23.272727272727273*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          x*(-17.454545454545457*Power(y,2) - 17.454545454545457*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)*
           (23.272727272727273*Power(y,3) + 23.272727272727273*Power(z,3)) + 
          Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)*
           (23.272727272727273*Power(y,5) + 23.272727272727273*Power(y,3)*Power(z,2) + 
             23.272727272727273*Power(y,2)*Power(z,3) + 23.272727272727273*Power(z,5)) + 
          Power(x,3)*(23.272727272727273*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
             23.272727272727273*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
             17.454545454545457*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),8.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),7)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
          (0.09090909090909091*Power(x,4) + 0.09090909090909091*Power(y,4) + 
            0.18181818181818182*Power(y,2)*Power(z,2) + 0.09090909090909091*Power(z,4) + 
            Power(x,2)*(0.18181818181818182*Power(y,2) + 0.18181818181818182*Power(z,2)))
           + M*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)*
          (0.8181818181818181*Power(x,4) + 0.8181818181818181*Power(y,4) + 
            1.6363636363636362*Power(y,2)*Power(z,2) + 0.8181818181818181*Power(z,4) + 
            Power(x,2)*(1.6363636363636362*Power(y,2) + 1.6363636363636362*Power(z,2)))\
          + Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.5)*
          (3.2727272727272725*Power(x,4) + 3.2727272727272725*Power(y,4) + 
            6.545454545454545*Power(y,2)*Power(z,2) + 3.2727272727272725*Power(z,4) + 
            Power(x,2)*(6.545454545454545*Power(y,2) + 6.545454545454545*Power(z,2))) + 
         Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.)*
          (7.636363636363637*Power(x,4) + 7.636363636363637*Power(y,4) + 
            15.272727272727273*Power(y,2)*Power(z,2) + 7.636363636363637*Power(z,4) + 
            Power(x,2)*(15.272727272727273*Power(y,2) + 15.272727272727273*Power(z,2)))\
          + Power(M,9)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)*
          (-1.3636363636363635*Power(x,4) + 1.4545454545454546*Power(x,6) - 
            1.3636363636363635*Power(y,4) + 1.4545454545454546*Power(y,6) - 
            2.727272727272727*Power(y,2)*Power(z,2) + 
            2.909090909090909*Power(y,3)*Power(z,3) - 1.3636363636363635*Power(z,4) + 
            1.4545454545454546*Power(z,6) + 
            Power(x,2)*(-2.727272727272727*Power(y,2) - 2.727272727272727*Power(z,2)) + 
            Power(x,3)*(2.909090909090909*Power(y,3) + 2.909090909090909*Power(z,3))) + 
         Power(M,8)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)*
          (-6.454545454545454*Power(x,4) + 7.2727272727272725*Power(x,6) - 
            6.454545454545454*Power(y,4) + 7.2727272727272725*Power(y,6) - 
            12.909090909090908*Power(y,2)*Power(z,2) + 
            14.545454545454545*Power(y,3)*Power(z,3) - 6.454545454545454*Power(z,4) + 
            7.2727272727272725*Power(z,6) + 
            Power(x,2)*(-12.909090909090908*Power(y,2) - 
               12.909090909090908*Power(z,2)) + 
            Power(x,3)*(14.545454545454545*Power(y,3) + 14.545454545454545*Power(z,3)))\
          + Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.5)*
          (1.4545454545454546*Power(x,6) + 1.4545454545454546*Power(y,6) + 
            10.*Power(z,4) + 1.4545454545454546*Power(z,6) + 
            Power(y,4)*(10. + 1.4545454545454546*Power(z,2)) + 
            Power(x,4)*(10. + 1.4545454545454546*Power(y,2) + 
               1.4545454545454546*Power(z,2)) + 
            Power(y,2)*(20.*Power(z,2) + 1.4545454545454546*Power(z,4)) + 
            Power(x,2)*(20.*Power(y,2) + 1.4545454545454546*Power(y,4) + 
               20.*Power(z,2) + 1.4545454545454546*Power(z,4))) + 
         Power(M,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
          (7.2727272727272725*Power(x,6) + 7.2727272727272725*Power(y,6) + 
            4.181818181818182*Power(z,4) + 7.2727272727272725*Power(z,6) + 
            Power(y,4)*(4.181818181818182 + 7.2727272727272725*Power(z,2)) + 
            Power(x,4)*(4.181818181818182 + 7.2727272727272725*Power(y,2) + 
               7.2727272727272725*Power(z,2)) + 
            Power(y,2)*(8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4)) + 
            Power(x,2)*(8.363636363636363*Power(y,2) + 7.2727272727272725*Power(y,4) + 
               8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4))) + 
         Power(M,6)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            2.90909090909091*Power(y,3)*Power(z,3) - 6.909090909090909*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-6.909090909090909 + 16.*Power(z,2)) + 
            Power(x,4)*(-6.909090909090909 + 16.*Power(y,2) + 16.*Power(z,2)) + 
            Power(x,3)*(-2.90909090909091*Power(y,3) - 2.90909090909091*Power(z,3)) + 
            Power(y,2)*(-13.818181818181818*Power(z,2) + 16.*Power(z,4)) + 
            Power(x,2)*(-13.818181818181818*Power(y,2) + 16.*Power(y,4) - 
               13.818181818181818*Power(z,2) + 16.*Power(z,4))) + 
         Power(M,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            14.545454545454543*Power(y,3)*Power(z,3) - 11.272727272727273*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-11.272727272727273 + 21.818181818181817*Power(z,2)) + 
            Power(x,4)*(-11.272727272727273 + 21.818181818181817*Power(y,2) + 
               21.818181818181817*Power(z,2)) + 
            Power(x,3)*(-14.545454545454543*Power(y,3) - 
               14.545454545454543*Power(z,3)) + 
            Power(y,2)*(-22.545454545454547*Power(z,2) + 
               21.818181818181817*Power(z,4)) + 
            Power(x,2)*(-22.545454545454547*Power(y,2) + 21.818181818181817*Power(y,4) - 
               22.545454545454547*Power(z,2) + 21.818181818181817*Power(z,4))))/
       (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
         Power(1.*M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6)*
         (1.3636363636363635*Power(M,3) + 
           1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
           0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
           0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))),0.5));
}

static void
func_Pi02(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*Power(y,2)*(4.363636363636363*Power(x,5)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       4.363636363636363*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       4.363636363636363*Power(y,3)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       4.363636363636363*Power(y,2)*Power(z,3)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       4.363636363636363*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
       2.909090909090909*Power(y,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
       2.909090909090909*y*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
       Power(x,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.)*
        (4.363636363636363*Power(y,2) + 4.363636363636363*Power(z,2)) + 
       Power(x,2)*(4.363636363636363*Power(y,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
          4.363636363636363*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
          2.909090909090909*y*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)) + 
       Power(M,6)*(1.4545454545454546*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
          1.4545454545454546*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
          1.4545454545454546*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
          1.4545454545454546*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
          2.909090909090909*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)*
           (1.4545454545454546*Power(y,3) + 1.4545454545454546*Power(z,3)) + 
          Power(y,3)*(1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
             2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),6.)) + 
          Power(x,3)*(1.4545454545454546*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
             1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
             2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),6.))) + 
       Power(M,5)*(11.636363636363637*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
          11.636363636363637*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
          11.636363636363637*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
          11.636363636363637*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
          17.454545454545453*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
           (11.636363636363637*Power(y,3) + 11.636363636363637*Power(z,3)) + 
          Power(y,3)*(11.636363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
             17.454545454545453*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)) + 
          Power(x,3)*(11.636363636363637*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
             11.636363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
             17.454545454545453*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5))) + 
       Power(M,4)*(36.36363636363637*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
          36.36363636363637*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
          36.36363636363637*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
          36.36363636363637*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) - 
          43.63636363636364*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
          Power(x,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.)*
           (36.36363636363637*Power(y,2) + 36.36363636363637*Power(z,2)) + 
          Power(y,3)*(36.36363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.) - 
             43.63636363636364*Power(Power(x,2) + Power(y,2) + Power(z,2),7.)) + 
          Power(x,2)*(36.36363636363637*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
             36.36363636363637*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.) - 
             43.63636363636364*y*Power(Power(x,2) + Power(y,2) + Power(z,2),7.))) + 
       Power(M,3)*(58.18181818181819*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
          58.18181818181819*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
          58.18181818181819*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
          58.18181818181819*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) - 
          75.63636363636364*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          17.454545454545453*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          Power(y,3)*(58.18181818181819*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) - 
             58.18181818181819*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)) + 
          Power(x,3)*(58.18181818181819*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
             58.18181818181819*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
             17.454545454545453*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)) + 
          Power(x,2)*(58.18181818181819*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
             58.18181818181819*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) - 
             75.63636363636364*y*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5))) + 
       Power(M,2)*(50.909090909090914*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
          50.909090909090914*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
          50.909090909090914*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
          50.909090909090914*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) - 
          46.54545454545455*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
          2.909090909090909*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
          Power(y,3)*(50.909090909090914*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) - 
             43.63636363636364*Power(Power(x,2) + Power(y,2) + Power(z,2),8.)) + 
          Power(x,3)*(50.909090909090914*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
             50.909090909090914*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
             2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),8.)) + 
          Power(x,2)*(50.909090909090914*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
             50.909090909090914*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) - 
             46.54545454545455*y*Power(Power(x,2) + Power(y,2) + Power(z,2),8.))) + 
       M*(23.272727272727273*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          23.272727272727273*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          23.272727272727273*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          23.272727272727273*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
          17.454545454545453*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          Power(x,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)*
           (23.272727272727273*Power(y,2) + 23.272727272727273*Power(z,2)) + 
          Power(y,3)*(23.272727272727273*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
             17.454545454545453*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)) + 
          Power(x,2)*(23.272727272727273*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
             23.272727272727273*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
             17.454545454545453*y*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),8.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),7)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
          (0.09090909090909091*Power(x,4) + 0.09090909090909091*Power(y,4) + 
            0.18181818181818182*Power(y,2)*Power(z,2) + 0.09090909090909091*Power(z,4) + 
            Power(x,2)*(0.18181818181818182*Power(y,2) + 0.18181818181818182*Power(z,2)))
           + M*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)*
          (0.8181818181818181*Power(x,4) + 0.8181818181818181*Power(y,4) + 
            1.6363636363636362*Power(y,2)*Power(z,2) + 0.8181818181818181*Power(z,4) + 
            Power(x,2)*(1.6363636363636362*Power(y,2) + 1.6363636363636362*Power(z,2)))\
          + Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.5)*
          (3.2727272727272725*Power(x,4) + 3.2727272727272725*Power(y,4) + 
            6.545454545454545*Power(y,2)*Power(z,2) + 3.2727272727272725*Power(z,4) + 
            Power(x,2)*(6.545454545454545*Power(y,2) + 6.545454545454545*Power(z,2))) + 
         Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.)*
          (7.636363636363637*Power(x,4) + 7.636363636363637*Power(y,4) + 
            15.272727272727273*Power(y,2)*Power(z,2) + 7.636363636363637*Power(z,4) + 
            Power(x,2)*(15.272727272727273*Power(y,2) + 15.272727272727273*Power(z,2)))\
          + Power(M,9)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)*
          (-1.3636363636363635*Power(x,4) + 1.4545454545454546*Power(x,6) - 
            1.3636363636363635*Power(y,4) + 1.4545454545454546*Power(y,6) - 
            2.727272727272727*Power(y,2)*Power(z,2) + 
            2.909090909090909*Power(y,3)*Power(z,3) - 1.3636363636363635*Power(z,4) + 
            1.4545454545454546*Power(z,6) + 
            Power(x,2)*(-2.727272727272727*Power(y,2) - 2.727272727272727*Power(z,2)) + 
            Power(x,3)*(2.909090909090909*Power(y,3) + 2.909090909090909*Power(z,3))) + 
         Power(M,8)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)*
          (-6.454545454545454*Power(x,4) + 7.2727272727272725*Power(x,6) - 
            6.454545454545454*Power(y,4) + 7.2727272727272725*Power(y,6) - 
            12.909090909090908*Power(y,2)*Power(z,2) + 
            14.545454545454545*Power(y,3)*Power(z,3) - 6.454545454545454*Power(z,4) + 
            7.2727272727272725*Power(z,6) + 
            Power(x,2)*(-12.909090909090908*Power(y,2) - 
               12.909090909090908*Power(z,2)) + 
            Power(x,3)*(14.545454545454545*Power(y,3) + 14.545454545454545*Power(z,3)))\
          + Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.5)*
          (1.4545454545454546*Power(x,6) + 1.4545454545454546*Power(y,6) + 
            10.*Power(z,4) + 1.4545454545454546*Power(z,6) + 
            Power(y,4)*(10. + 1.4545454545454546*Power(z,2)) + 
            Power(x,4)*(10. + 1.4545454545454546*Power(y,2) + 
               1.4545454545454546*Power(z,2)) + 
            Power(y,2)*(20.*Power(z,2) + 1.4545454545454546*Power(z,4)) + 
            Power(x,2)*(20.*Power(y,2) + 1.4545454545454546*Power(y,4) + 
               20.*Power(z,2) + 1.4545454545454546*Power(z,4))) + 
         Power(M,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
          (7.2727272727272725*Power(x,6) + 7.2727272727272725*Power(y,6) + 
            4.181818181818182*Power(z,4) + 7.2727272727272725*Power(z,6) + 
            Power(y,4)*(4.181818181818182 + 7.2727272727272725*Power(z,2)) + 
            Power(x,4)*(4.181818181818182 + 7.2727272727272725*Power(y,2) + 
               7.2727272727272725*Power(z,2)) + 
            Power(y,2)*(8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4)) + 
            Power(x,2)*(8.363636363636363*Power(y,2) + 7.2727272727272725*Power(y,4) + 
               8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4))) + 
         Power(M,6)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            2.90909090909091*Power(y,3)*Power(z,3) - 6.909090909090909*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-6.909090909090909 + 16.*Power(z,2)) + 
            Power(x,4)*(-6.909090909090909 + 16.*Power(y,2) + 16.*Power(z,2)) + 
            Power(x,3)*(-2.90909090909091*Power(y,3) - 2.90909090909091*Power(z,3)) + 
            Power(y,2)*(-13.818181818181818*Power(z,2) + 16.*Power(z,4)) + 
            Power(x,2)*(-13.818181818181818*Power(y,2) + 16.*Power(y,4) - 
               13.818181818181818*Power(z,2) + 16.*Power(z,4))) + 
         Power(M,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            14.545454545454543*Power(y,3)*Power(z,3) - 11.272727272727273*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-11.272727272727273 + 21.818181818181817*Power(z,2)) + 
            Power(x,4)*(-11.272727272727273 + 21.818181818181817*Power(y,2) + 
               21.818181818181817*Power(z,2)) + 
            Power(x,3)*(-14.545454545454543*Power(y,3) - 
               14.545454545454543*Power(z,3)) + 
            Power(y,2)*(-22.545454545454547*Power(z,2) + 
               21.818181818181817*Power(z,4)) + 
            Power(x,2)*(-22.545454545454547*Power(y,2) + 21.818181818181817*Power(y,4) - 
               22.545454545454547*Power(z,2) + 21.818181818181817*Power(z,4))))/
       (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
         Power(1.*M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6)*
         (1.3636363636363635*Power(M,3) + 
           1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
           0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
           0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))),0.5));
}

static void
func_Pi03(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*Power(z,2)*(4.363636363636363*Power(x,5)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       4.363636363636363*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       4.363636363636363*Power(y,3)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       4.363636363636363*Power(y,2)*Power(z,3)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
       4.363636363636363*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
       2.909090909090909*Power(y,2)*z*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
       2.909090909090909*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
       Power(x,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.)*
        (4.363636363636363*Power(y,2) + 4.363636363636363*Power(z,2)) + 
       Power(x,2)*(4.363636363636363*Power(y,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
          4.363636363636363*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) - 
          2.909090909090909*z*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)) + 
       Power(M,6)*(1.4545454545454546*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
          1.4545454545454546*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
          1.4545454545454546*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
          1.4545454545454546*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
          2.909090909090909*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)*
           (1.4545454545454546*Power(y,3) + 1.4545454545454546*Power(z,3)) + 
          Power(y,3)*(1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
             2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),6.)) + 
          Power(x,3)*(1.4545454545454546*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.) + 
             1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.) - 
             2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),6.))) + 
       Power(M,5)*(11.636363636363637*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
          11.636363636363637*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
          11.636363636363637*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
          11.636363636363637*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
          17.454545454545453*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
           (11.636363636363637*Power(y,3) + 11.636363636363637*Power(z,3)) + 
          Power(y,3)*(11.636363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
             17.454545454545453*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)) + 
          Power(x,3)*(11.636363636363637*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) + 
             11.636363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),5.5) - 
             17.454545454545453*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5))) + 
       Power(M,4)*(36.36363636363637*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
          36.36363636363637*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
          36.36363636363637*Power(y,3)*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
          36.36363636363637*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) - 
          43.63636363636364*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
          Power(x,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.)*
           (36.36363636363637*Power(y,2) + 36.36363636363637*Power(z,2)) + 
          Power(y,2)*(36.36363636363637*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.) - 
             43.63636363636364*z*Power(Power(x,2) + Power(y,2) + Power(z,2),7.)) + 
          Power(x,2)*(36.36363636363637*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + 
             36.36363636363637*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.) - 
             43.63636363636364*z*Power(Power(x,2) + Power(y,2) + Power(z,2),7.))) + 
       Power(M,3)*(58.18181818181819*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
          58.18181818181819*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
          58.18181818181819*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) - 
          58.18181818181819*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          Power(y,3)*(58.18181818181819*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
             17.454545454545453*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)) + 
          Power(x,3)*(58.18181818181819*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
             58.18181818181819*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
             17.454545454545453*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)) + 
          Power(y,2)*(58.18181818181819*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) - 
             75.63636363636364*z*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)) + 
          Power(x,2)*(58.18181818181819*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 
             58.18181818181819*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) - 
             75.63636363636364*z*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5))) + 
       Power(M,2)*(50.909090909090914*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
          50.909090909090914*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
          50.909090909090914*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) - 
          43.63636363636364*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 
          Power(y,3)*(50.909090909090914*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
             2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),8.)) + 
          Power(x,3)*(50.909090909090914*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
             50.909090909090914*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
             2.909090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),8.)) + 
          Power(y,2)*(50.909090909090914*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) - 
             46.54545454545455*z*Power(Power(x,2) + Power(y,2) + Power(z,2),8.)) + 
          Power(x,2)*(50.909090909090914*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) + 
             50.909090909090914*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.) - 
             46.54545454545455*z*Power(Power(x,2) + Power(y,2) + Power(z,2),8.))) + 
       M*(23.272727272727273*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          23.272727272727273*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          23.272727272727273*Power(y,3)*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
          23.272727272727273*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
          17.454545454545453*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          Power(x,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5)*
           (23.272727272727273*Power(y,2) + 23.272727272727273*Power(z,2)) + 
          Power(y,2)*(23.272727272727273*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
             17.454545454545453*z*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)) + 
          Power(x,2)*(23.272727272727273*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) + 
             23.272727272727273*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 
             17.454545454545453*z*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),8.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),7)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
          (0.09090909090909091*Power(x,4) + 0.09090909090909091*Power(y,4) + 
            0.18181818181818182*Power(y,2)*Power(z,2) + 0.09090909090909091*Power(z,4) + 
            Power(x,2)*(0.18181818181818182*Power(y,2) + 0.18181818181818182*Power(z,2)))
           + M*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)*
          (0.8181818181818181*Power(x,4) + 0.8181818181818181*Power(y,4) + 
            1.6363636363636362*Power(y,2)*Power(z,2) + 0.8181818181818181*Power(z,4) + 
            Power(x,2)*(1.6363636363636362*Power(y,2) + 1.6363636363636362*Power(z,2)))\
          + Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.5)*
          (3.2727272727272725*Power(x,4) + 3.2727272727272725*Power(y,4) + 
            6.545454545454545*Power(y,2)*Power(z,2) + 3.2727272727272725*Power(z,4) + 
            Power(x,2)*(6.545454545454545*Power(y,2) + 6.545454545454545*Power(z,2))) + 
         Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.)*
          (7.636363636363637*Power(x,4) + 7.636363636363637*Power(y,4) + 
            15.272727272727273*Power(y,2)*Power(z,2) + 7.636363636363637*Power(z,4) + 
            Power(x,2)*(15.272727272727273*Power(y,2) + 15.272727272727273*Power(z,2)))\
          + Power(M,9)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)*
          (-1.3636363636363635*Power(x,4) + 1.4545454545454546*Power(x,6) - 
            1.3636363636363635*Power(y,4) + 1.4545454545454546*Power(y,6) - 
            2.727272727272727*Power(y,2)*Power(z,2) + 
            2.909090909090909*Power(y,3)*Power(z,3) - 1.3636363636363635*Power(z,4) + 
            1.4545454545454546*Power(z,6) + 
            Power(x,2)*(-2.727272727272727*Power(y,2) - 2.727272727272727*Power(z,2)) + 
            Power(x,3)*(2.909090909090909*Power(y,3) + 2.909090909090909*Power(z,3))) + 
         Power(M,8)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)*
          (-6.454545454545454*Power(x,4) + 7.2727272727272725*Power(x,6) - 
            6.454545454545454*Power(y,4) + 7.2727272727272725*Power(y,6) - 
            12.909090909090908*Power(y,2)*Power(z,2) + 
            14.545454545454545*Power(y,3)*Power(z,3) - 6.454545454545454*Power(z,4) + 
            7.2727272727272725*Power(z,6) + 
            Power(x,2)*(-12.909090909090908*Power(y,2) - 
               12.909090909090908*Power(z,2)) + 
            Power(x,3)*(14.545454545454545*Power(y,3) + 14.545454545454545*Power(z,3)))\
          + Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.5)*
          (1.4545454545454546*Power(x,6) + 1.4545454545454546*Power(y,6) + 
            10.*Power(z,4) + 1.4545454545454546*Power(z,6) + 
            Power(y,4)*(10. + 1.4545454545454546*Power(z,2)) + 
            Power(x,4)*(10. + 1.4545454545454546*Power(y,2) + 
               1.4545454545454546*Power(z,2)) + 
            Power(y,2)*(20.*Power(z,2) + 1.4545454545454546*Power(z,4)) + 
            Power(x,2)*(20.*Power(y,2) + 1.4545454545454546*Power(y,4) + 
               20.*Power(z,2) + 1.4545454545454546*Power(z,4))) + 
         Power(M,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
          (7.2727272727272725*Power(x,6) + 7.2727272727272725*Power(y,6) + 
            4.181818181818182*Power(z,4) + 7.2727272727272725*Power(z,6) + 
            Power(y,4)*(4.181818181818182 + 7.2727272727272725*Power(z,2)) + 
            Power(x,4)*(4.181818181818182 + 7.2727272727272725*Power(y,2) + 
               7.2727272727272725*Power(z,2)) + 
            Power(y,2)*(8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4)) + 
            Power(x,2)*(8.363636363636363*Power(y,2) + 7.2727272727272725*Power(y,4) + 
               8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4))) + 
         Power(M,6)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            2.90909090909091*Power(y,3)*Power(z,3) - 6.909090909090909*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-6.909090909090909 + 16.*Power(z,2)) + 
            Power(x,4)*(-6.909090909090909 + 16.*Power(y,2) + 16.*Power(z,2)) + 
            Power(x,3)*(-2.90909090909091*Power(y,3) - 2.90909090909091*Power(z,3)) + 
            Power(y,2)*(-13.818181818181818*Power(z,2) + 16.*Power(z,4)) + 
            Power(x,2)*(-13.818181818181818*Power(y,2) + 16.*Power(y,4) - 
               13.818181818181818*Power(z,2) + 16.*Power(z,4))) + 
         Power(M,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            14.545454545454543*Power(y,3)*Power(z,3) - 11.272727272727273*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-11.272727272727273 + 21.818181818181817*Power(z,2)) + 
            Power(x,4)*(-11.272727272727273 + 21.818181818181817*Power(y,2) + 
               21.818181818181817*Power(z,2)) + 
            Power(x,3)*(-14.545454545454543*Power(y,3) - 
               14.545454545454543*Power(z,3)) + 
            Power(y,2)*(-22.545454545454547*Power(z,2) + 
               21.818181818181817*Power(z,4)) + 
            Power(x,2)*(-22.545454545454547*Power(y,2) + 21.818181818181817*Power(y,4) - 
               22.545454545454547*Power(z,2) + 21.818181818181817*Power(z,4))))/
       (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
         Power(1.*M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6)*
         (1.3636363636363635*Power(M,3) + 
           1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
           0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
           0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))),0.5));
}

static void
func_Pi11(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*(0.7272727272727274*Power(x,7)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       0.7272727272727274*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       1.4545454545454548*Power(y,5)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       0.7272727272727274*Power(y,4)*Power(z,3)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       0.7272727272727274*Power(y,3)*Power(z,4)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       1.4545454545454548*Power(y,2)*Power(z,5)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       0.7272727272727274*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
       0.7272727272727274*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),15.) - 
       0.7272727272727274*Power(y,3)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),15.) - 
       0.7272727272727274*Power(y,2)*Power(z,3)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),15.) - 
       0.7272727272727274*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),15.) + 
       Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)*
        (0.7272727272727274*Power(y,3) + 0.7272727272727274*Power(z,3)) + 
       Power(x,5)*(1.4545454545454548*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          1.4545454545454548*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
          0.7272727272727274*Power(Power(x,2) + Power(y,2) + Power(z,2),15.)) + 
       Power(x,3)*(0.7272727272727274*Power(y,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          0.7272727272727274*Power(z,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
          0.7272727272727274*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),15.) + 
          Power(y,2)*(1.4545454545454548*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
             0.7272727272727274*Power(Power(x,2) + Power(y,2) + Power(z,2),15.))) + 
       Power(x,2)*(1.4545454545454548*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          1.4545454545454548*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          1.4545454545454548*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
          0.7272727272727274*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),15.) + 
          Power(y,3)*(1.4545454545454548*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
             0.7272727272727274*Power(Power(x,2) + Power(y,2) + Power(z,2),15.))) + 
       Power(M,11)*(1.4545454545454548*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          1.4545454545454548*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          1.4545454545454548*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          0.7272727272727274*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          Power(x,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*
           (1.4545454545454548*Power(y,2) + 1.4545454545454548*Power(z,2)) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*
           (1.4545454545454548*Power(y,5) + 1.4545454545454548*Power(y,3)*Power(z,2) + 
             1.4545454545454548*Power(y,2)*Power(z,3) + 1.4545454545454548*Power(z,5)) + 
          Power(y,5)*(2.9090909090909096*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             0.7272727272727274*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
          Power(y,3)*(1.4545454545454548*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             0.7272727272727274*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
          Power(y,2)*(2.9090909090909096*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             0.7272727272727274*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
          Power(x,3)*(1.4545454545454548*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             1.4545454545454548*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             0.7272727272727274*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             Power(y,2)*(2.9090909090909096*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
                0.7272727272727274*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)))) + 
       Power(M,10)*(15.272727272727275*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          15.272727272727275*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          15.272727272727275*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          8.000000000000002*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          Power(x,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (15.272727272727275*Power(y,2) + 15.272727272727275*Power(z,2)) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (15.272727272727275*Power(y,5) + 15.272727272727275*Power(y,3)*Power(z,2) + 
             15.272727272727275*Power(y,2)*Power(z,3) + 15.272727272727275*Power(z,5)) + 
          Power(y,5)*(30.54545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             8.000000000000002*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
          Power(y,3)*(15.272727272727275*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             8.000000000000002*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)
             ) + Power(y,2)*(30.54545454545455*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             8.000000000000002*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)
             ) + Power(x,3)*(15.272727272727275*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             15.272727272727275*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             8.000000000000002*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             Power(y,2)*(30.54545454545455*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
                8.000000000000002*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)))) + 
       Power(M,9)*(10.909090909090912*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          72.72727272727273*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          72.72727272727273*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          72.72727272727273*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          40.00000000000001*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*
           (10.909090909090912*Power(y,3) + 10.909090909090912*Power(z,3)) + 
          Power(y,5)*(145.45454545454547*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             40.00000000000001*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(x,5)*(83.63636363636364*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             83.63636363636364*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             10.909090909090912*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(y,3)*(72.72727272727273*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             40.00000000000001*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(y,2)*(145.45454545454547*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             40.00000000000001*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(x,3)*(72.72727272727273*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             72.72727272727273*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             29.090909090909097*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             Power(y,2)*(145.45454545454547*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
                29.090909090909097*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5))) + 
          Power(x,2)*(83.63636363636364*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             83.63636363636364*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             83.63636363636364*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             21.818181818181824*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             Power(y,3)*(83.63636363636364*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
                21.818181818181824*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)))) + 
       Power(M,8)*(85.81818181818183*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          207.2727272727273*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          207.2727272727273*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          207.2727272727273*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          120.00000000000001*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)*
           (85.81818181818183*Power(y,3) + 85.81818181818183*Power(z,3)) + 
          Power(y,5)*(414.5454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             120.00000000000001*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(x,5)*(293.0909090909091*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             293.0909090909091*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             73.45454545454545*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(y,3)*(207.2727272727273*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             120.00000000000001*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(y,2)*(414.5454545454546*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             120.00000000000001*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(x,2)*(293.0909090909091*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             293.0909090909091*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             293.0909090909091*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             146.9090909090909*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             Power(y,3)*(293.0909090909091*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
                146.9090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),11.))) + 
          Power(x,3)*(207.2727272727273*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             207.2727272727273*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             46.545454545454554*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             Power(y,2)*(414.5454545454546*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
                46.545454545454554*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)))) + 
       Power(M,7)*(281.4545454545455*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          392.72727272727275*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          392.72727272727275*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          392.72727272727275*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          240.00000000000003*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
           (281.4545454545455*Power(y,3) + 281.4545454545455*Power(z,3)) + 
          Power(y,5)*(785.4545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             240.00000000000003*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(x,5)*(674.1818181818182*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             674.1818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             215.2727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(y,3)*(392.72727272727275*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             240.00000000000003*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(y,2)*(785.4545454545455*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             240.00000000000003*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(x,2)*(674.1818181818182*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             674.1818181818182*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             674.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             266.90909090909093*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             Power(y,3)*(674.1818181818182*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
                266.90909090909093*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5))) + 
          Power(x,3)*(392.72727272727275*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             392.72727272727275*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             188.36363636363637*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             Power(y,2)*(785.4545454545455*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
                188.36363636363637*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)))) + 
       Power(M,6)*(506.909090909091*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          519.2727272727274*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          519.2727272727274*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          519.2727272727274*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          336.*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
           (506.909090909091*Power(y,3) + 506.909090909091*Power(z,3)) + 
          Power(x,5)*(1026.1818181818182*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             1026.1818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             360.7272727272728*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(y,5)*(1038.5454545454547*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             336.*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(y,3)*(519.2727272727274*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             336.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(y,2)*(1038.5454545454547*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             336.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(x,3)*(519.2727272727274*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             519.2727272727274*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             378.909090909091*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             Power(y,2)*(1038.5454545454547*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
                378.909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),12.))) + 
          Power(x,2)*(1026.1818181818182*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             1026.1818181818182*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             1026.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             317.8181818181818*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             Power(y,3)*(1026.1818181818182*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
                317.8181818181818*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)))) + 
       Power(M,5)*(553.4545454545455*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          488.7272727272728*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          488.7272727272728*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          488.7272727272728*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
          336.*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)*
           (553.4545454545455*Power(y,3) + 553.4545454545455*Power(z,3)) + 
          Power(x,5)*(1042.1818181818182*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             1042.1818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             382.54545454545456*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(y,5)*(977.4545454545456*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             336.*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(y,3)*(488.7272727272728*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             336.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(y,2)*(977.4545454545456*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             336.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(x,3)*(488.7272727272728*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             488.7272727272728*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             390.5454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             Power(y,2)*(977.4545454545456*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
                390.5454545454546*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5))) + 
          Power(x,2)*(1042.1818181818182*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             1042.1818181818182*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             1042.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             328.00000000000006*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             Power(y,3)*(1042.1818181818182*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
                328.00000000000006*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)))) + 
       Power(M,4)*(381.81818181818187*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          327.2727272727273*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          327.2727272727273*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          327.2727272727273*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
          240.00000000000003*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)*
           (381.81818181818187*Power(y,3) + 381.81818181818187*Power(z,3)) + 
          Power(x,5)*(709.0909090909091*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             709.0909090909091*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             269.0909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(y,5)*(654.5454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             240.00000000000003*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(y,3)*(327.2727272727273*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             240.00000000000003*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(y,2)*(654.5454545454546*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             240.00000000000003*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(x,3)*(327.2727272727273*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             327.2727272727273*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             269.81818181818187*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             Power(y,2)*(654.5454545454546*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
                269.81818181818187*Power(Power(x,2) + Power(y,2) + Power(z,2),13.))) + 
          Power(x,2)*(709.0909090909091*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             709.0909090909091*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             709.0909090909091*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             239.2727272727273*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             Power(y,3)*(709.0909090909091*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
                239.2727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)))) + 
       Power(M,3)*(169.45454545454547*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          152.72727272727275*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          152.72727272727275*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          152.72727272727275*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
          120.00000000000001*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)*
           (169.45454545454547*Power(y,3) + 169.45454545454547*Power(z,3)) + 
          Power(x,5)*(322.1818181818182*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             322.1818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             128.00000000000003*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(y,5)*(305.4545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             120.00000000000001*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(y,3)*(152.72727272727275*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             120.00000000000001*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(y,2)*(305.4545454545455*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             120.00000000000001*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(x,3)*(152.72727272727275*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             152.72727272727275*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             128.00000000000003*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             Power(y,2)*(305.4545454545455*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
                128.00000000000003*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5))) + 
          Power(x,2)*(322.1818181818182*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             322.1818181818182*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             322.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             120.00000000000001*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             Power(y,3)*(322.1818181818182*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
                120.00000000000001*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)))) + 
       Power(M,2)*(48.727272727272734*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          47.27272727272728*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          47.27272727272728*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          47.27272727272728*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
          40.00000000000001*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)*
           (48.727272727272734*Power(y,3) + 48.727272727272734*Power(z,3)) + 
          Power(x,5)*(96.00000000000001*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             96.00000000000001*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.727272727272734*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)) + 
          Power(y,5)*(94.54545454545456*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.00000000000001*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)) + 
          Power(y,3)*(47.27272727272728*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.00000000000001*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)
             ) + Power(y,2)*(94.54545454545456*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.00000000000001*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)
             ) + Power(x,3)*(47.27272727272728*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             47.27272727272728*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.727272727272734*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
             Power(y,2)*(94.54545454545456*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
                40.727272727272734*Power(Power(x,2) + Power(y,2) + Power(z,2),14.))) + 
          Power(x,2)*(96.00000000000001*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             96.00000000000001*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             96.00000000000001*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.00000000000001*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
             Power(y,3)*(96.00000000000001*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
                40.00000000000001*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)))) + 
       M*(8.727272727272728*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          8.727272727272728*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          8.727272727272728*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          8.727272727272728*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
          8.000000000000002*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)*
           (8.727272727272728*Power(y,3) + 8.727272727272728*Power(z,3)) + 
          Power(y,5)*(17.454545454545457*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.000000000000002*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5)) + 
          Power(x,5)*(17.454545454545457*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             17.454545454545457*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.000000000000002*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5)) + 
          Power(y,3)*(8.727272727272728*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.000000000000002*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.5)) + 
          Power(y,2)*(17.454545454545457*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.000000000000002*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.5)) + 
          Power(x,3)*(8.727272727272728*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             8.727272727272728*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.000000000000002*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.5) + 
             Power(y,2)*(17.454545454545457*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
                8.000000000000002*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5))) + 
          Power(x,2)*(17.454545454545457*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             17.454545454545457*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             17.454545454545457*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.000000000000002*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.5) + 
             Power(y,3)*(17.454545454545457*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
                8.000000000000002*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5))))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),13.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),8)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
          (0.09090909090909091*Power(x,4) + 0.09090909090909091*Power(y,4) + 
            0.18181818181818182*Power(y,2)*Power(z,2) + 0.09090909090909091*Power(z,4) + 
            Power(x,2)*(0.18181818181818182*Power(y,2) + 0.18181818181818182*Power(z,2)))
           + M*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)*
          (0.8181818181818181*Power(x,4) + 0.8181818181818181*Power(y,4) + 
            1.6363636363636362*Power(y,2)*Power(z,2) + 0.8181818181818181*Power(z,4) + 
            Power(x,2)*(1.6363636363636362*Power(y,2) + 1.6363636363636362*Power(z,2)))\
          + Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.5)*
          (3.2727272727272725*Power(x,4) + 3.2727272727272725*Power(y,4) + 
            6.545454545454545*Power(y,2)*Power(z,2) + 3.2727272727272725*Power(z,4) + 
            Power(x,2)*(6.545454545454545*Power(y,2) + 6.545454545454545*Power(z,2))) + 
         Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.)*
          (7.636363636363637*Power(x,4) + 7.636363636363637*Power(y,4) + 
            15.272727272727273*Power(y,2)*Power(z,2) + 7.636363636363637*Power(z,4) + 
            Power(x,2)*(15.272727272727273*Power(y,2) + 15.272727272727273*Power(z,2)))\
          + Power(M,9)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)*
          (-1.3636363636363635*Power(x,4) + 1.4545454545454546*Power(x,6) - 
            1.3636363636363635*Power(y,4) + 1.4545454545454546*Power(y,6) - 
            2.727272727272727*Power(y,2)*Power(z,2) + 
            2.909090909090909*Power(y,3)*Power(z,3) - 1.3636363636363635*Power(z,4) + 
            1.4545454545454546*Power(z,6) + 
            Power(x,2)*(-2.727272727272727*Power(y,2) - 2.727272727272727*Power(z,2)) + 
            Power(x,3)*(2.909090909090909*Power(y,3) + 2.909090909090909*Power(z,3))) + 
         Power(M,8)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)*
          (-6.454545454545454*Power(x,4) + 7.2727272727272725*Power(x,6) - 
            6.454545454545454*Power(y,4) + 7.2727272727272725*Power(y,6) - 
            12.909090909090908*Power(y,2)*Power(z,2) + 
            14.545454545454545*Power(y,3)*Power(z,3) - 6.454545454545454*Power(z,4) + 
            7.2727272727272725*Power(z,6) + 
            Power(x,2)*(-12.909090909090908*Power(y,2) - 
               12.909090909090908*Power(z,2)) + 
            Power(x,3)*(14.545454545454545*Power(y,3) + 14.545454545454545*Power(z,3)))\
          + Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.5)*
          (1.4545454545454546*Power(x,6) + 1.4545454545454546*Power(y,6) + 
            10.*Power(z,4) + 1.4545454545454546*Power(z,6) + 
            Power(y,4)*(10. + 1.4545454545454546*Power(z,2)) + 
            Power(x,4)*(10. + 1.4545454545454546*Power(y,2) + 
               1.4545454545454546*Power(z,2)) + 
            Power(y,2)*(20.*Power(z,2) + 1.4545454545454546*Power(z,4)) + 
            Power(x,2)*(20.*Power(y,2) + 1.4545454545454546*Power(y,4) + 
               20.*Power(z,2) + 1.4545454545454546*Power(z,4))) + 
         Power(M,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
          (7.2727272727272725*Power(x,6) + 7.2727272727272725*Power(y,6) + 
            4.181818181818182*Power(z,4) + 7.2727272727272725*Power(z,6) + 
            Power(y,4)*(4.181818181818182 + 7.2727272727272725*Power(z,2)) + 
            Power(x,4)*(4.181818181818182 + 7.2727272727272725*Power(y,2) + 
               7.2727272727272725*Power(z,2)) + 
            Power(y,2)*(8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4)) + 
            Power(x,2)*(8.363636363636363*Power(y,2) + 7.2727272727272725*Power(y,4) + 
               8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4))) + 
         Power(M,6)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            2.90909090909091*Power(y,3)*Power(z,3) - 6.909090909090909*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-6.909090909090909 + 16.*Power(z,2)) + 
            Power(x,4)*(-6.909090909090909 + 16.*Power(y,2) + 16.*Power(z,2)) + 
            Power(x,3)*(-2.90909090909091*Power(y,3) - 2.90909090909091*Power(z,3)) + 
            Power(y,2)*(-13.818181818181818*Power(z,2) + 16.*Power(z,4)) + 
            Power(x,2)*(-13.818181818181818*Power(y,2) + 16.*Power(y,4) - 
               13.818181818181818*Power(z,2) + 16.*Power(z,4))) + 
         Power(M,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            14.545454545454543*Power(y,3)*Power(z,3) - 11.272727272727273*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-11.272727272727273 + 21.818181818181817*Power(z,2)) + 
            Power(x,4)*(-11.272727272727273 + 21.818181818181817*Power(y,2) + 
               21.818181818181817*Power(z,2)) + 
            Power(x,3)*(-14.545454545454543*Power(y,3) - 
               14.545454545454543*Power(z,3)) + 
            Power(y,2)*(-22.545454545454547*Power(z,2) + 
               21.818181818181817*Power(z,4)) + 
            Power(x,2)*(-22.545454545454547*Power(y,2) + 21.818181818181817*Power(y,4) - 
               22.545454545454547*Power(z,2) + 21.818181818181817*Power(z,4))))/
       (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
         Power(1.*M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6)*
         (1.3636363636363635*Power(M,3) + 
           1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
           0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
           0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))),0.5));
}

static void
func_Pi12(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*x*y*(1.4545454545454546*Power(x,5)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
       1.4545454545454546*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
       1.4545454545454546*Power(y,3)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
       1.4545454545454546*Power(y,2)*Power(z,3)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
       1.4545454545454546*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
       0.36363636363636365*Power(y,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
       0.36363636363636365*y*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       x*(-0.36363636363636365*Power(y,2) - 0.36363636363636365*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       Power(x,3)*(1.4545454545454546*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          1.4545454545454546*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
          0.36363636363636365*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)) + 
       Power(x,2)*(1.4545454545454546*Power(y,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          1.4545454545454546*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
          0.36363636363636365*y*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)) + 
       Power(M,9)*(-1.4545454545454546*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          1.4545454545454546*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          1.4545454545454546*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          1.4545454545454546*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          0.7272727272727273*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*
           (-1.4545454545454546*Power(y,3) - 1.4545454545454546*Power(z,3)) + 
          Power(y,3)*(-1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
          Power(x,3)*(-1.4545454545454546*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5))) + 
       Power(M,8)*(-15.272727272727272*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          15.272727272727272*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          15.272727272727272*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          15.272727272727272*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          8.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (-15.272727272727272*Power(y,3) - 15.272727272727272*Power(z,3)) + 
          Power(y,3)*(-15.272727272727272*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
          Power(x,3)*(-15.272727272727272*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             15.272727272727272*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))) + 
       Power(M,7)*(-61.81818181818181*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          61.81818181818181*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          61.81818181818181*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          61.81818181818181*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          5.454545454545455*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          18.181818181818183*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          x*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
           (5.454545454545455*Power(y,2) + 5.454545454545455*Power(z,2)) + 
          Power(y,3)*(-61.81818181818181*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             23.636363636363637*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(x,3)*(-61.81818181818181*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             61.81818181818181*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             23.636363636363637*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(x,2)*(-61.81818181818181*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             61.81818181818181*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             5.454545454545455*y*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5))) + 
       Power(M,6)*(-121.45454545454545*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          121.45454545454545*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          121.45454545454545*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          121.45454545454545*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          36.72727272727273*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          26.909090909090907*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          x*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
           (36.72727272727273*Power(y,2) + 36.72727272727273*Power(z,2)) + 
          Power(y,3)*(-121.45454545454545*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             9.818181818181818*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(x,3)*(-121.45454545454545*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             121.45454545454545*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             9.818181818181818*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(x,2)*(-121.45454545454545*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             121.45454545454545*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             36.72727272727273*y*Power(Power(x,2) + Power(y,2) + Power(z,2),11.))) + 
       Power(M,5)*(-111.27272727272727*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          111.27272727272727*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          111.27272727272727*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          111.27272727272727*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          25.818181818181817*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
          26.909090909090907*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          x*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)*
           (25.818181818181817*Power(y,2) + 25.818181818181817*Power(z,2)) + 
          Power(y,3)*(-111.27272727272727*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             1.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(x,3)*(-111.27272727272727*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             111.27272727272727*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             1.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(x,2)*(-111.27272727272727*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             111.27272727272727*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             25.818181818181817*y*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5))) + 
       Power(M,4)*(-12.363636363636363*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          12.363636363636363*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          12.363636363636363*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          12.363636363636363*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          21.454545454545453*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          18.181818181818183*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          x*(-21.454545454545453*Power(y,2) - 21.454545454545453*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          Power(y,3)*(-12.363636363636363*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             3.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(x,3)*(-12.363636363636363*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             12.363636363636363*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             3.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(x,2)*(-12.363636363636363*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             12.363636363636363*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             21.454545454545453*y*Power(Power(x,2) + Power(y,2) + Power(z,2),12.))) + 
       Power(M,3)*(64.72727272727272*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          64.72727272727272*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          64.72727272727272*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          64.72727272727272*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
          27.272727272727273*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          8.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          x*(-27.272727272727273*Power(y,2) - 27.272727272727273*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          Power(y,3)*(64.72727272727272*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             19.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(x,3)*(64.72727272727272*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             64.72727272727272*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             19.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(x,2)*(64.72727272727272*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             64.72727272727272*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             27.272727272727273*y*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5))) + 
       Power(M,2)*(54.54545454545455*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          54.54545454545455*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          54.54545454545455*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          54.54545454545455*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
          14.90909090909091*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          0.7272727272727273*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          x*(-14.90909090909091*Power(y,2) - 14.90909090909091*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          Power(y,3)*(54.54545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             14.181818181818182*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(x,3)*(54.54545454545455*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             54.54545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             14.181818181818182*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(x,2)*(54.54545454545455*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             54.54545454545455*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             14.90909090909091*y*Power(Power(x,2) + Power(y,2) + Power(z,2),13.))) + 
       M*(16.727272727272727*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          16.727272727272727*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          16.727272727272727*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          16.727272727272727*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
          4.*y*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          x*(-4.*Power(y,2) - 4.*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          Power(y,3)*(16.727272727272727*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             4.*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(x,3)*(16.727272727272727*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             16.727272727272727*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             4.*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(x,2)*(16.727272727272727*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             16.727272727272727*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             4.*y*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),13.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),8)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
          (0.09090909090909091*Power(x,4) + 0.09090909090909091*Power(y,4) + 
            0.18181818181818182*Power(y,2)*Power(z,2) + 0.09090909090909091*Power(z,4) + 
            Power(x,2)*(0.18181818181818182*Power(y,2) + 0.18181818181818182*Power(z,2)))
           + M*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)*
          (0.8181818181818181*Power(x,4) + 0.8181818181818181*Power(y,4) + 
            1.6363636363636362*Power(y,2)*Power(z,2) + 0.8181818181818181*Power(z,4) + 
            Power(x,2)*(1.6363636363636362*Power(y,2) + 1.6363636363636362*Power(z,2)))\
          + Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.5)*
          (3.2727272727272725*Power(x,4) + 3.2727272727272725*Power(y,4) + 
            6.545454545454545*Power(y,2)*Power(z,2) + 3.2727272727272725*Power(z,4) + 
            Power(x,2)*(6.545454545454545*Power(y,2) + 6.545454545454545*Power(z,2))) + 
         Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.)*
          (7.636363636363637*Power(x,4) + 7.636363636363637*Power(y,4) + 
            15.272727272727273*Power(y,2)*Power(z,2) + 7.636363636363637*Power(z,4) + 
            Power(x,2)*(15.272727272727273*Power(y,2) + 15.272727272727273*Power(z,2)))\
          + Power(M,9)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)*
          (-1.3636363636363635*Power(x,4) + 1.4545454545454546*Power(x,6) - 
            1.3636363636363635*Power(y,4) + 1.4545454545454546*Power(y,6) - 
            2.727272727272727*Power(y,2)*Power(z,2) + 
            2.909090909090909*Power(y,3)*Power(z,3) - 1.3636363636363635*Power(z,4) + 
            1.4545454545454546*Power(z,6) + 
            Power(x,2)*(-2.727272727272727*Power(y,2) - 2.727272727272727*Power(z,2)) + 
            Power(x,3)*(2.909090909090909*Power(y,3) + 2.909090909090909*Power(z,3))) + 
         Power(M,8)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)*
          (-6.454545454545454*Power(x,4) + 7.2727272727272725*Power(x,6) - 
            6.454545454545454*Power(y,4) + 7.2727272727272725*Power(y,6) - 
            12.909090909090908*Power(y,2)*Power(z,2) + 
            14.545454545454545*Power(y,3)*Power(z,3) - 6.454545454545454*Power(z,4) + 
            7.2727272727272725*Power(z,6) + 
            Power(x,2)*(-12.909090909090908*Power(y,2) - 
               12.909090909090908*Power(z,2)) + 
            Power(x,3)*(14.545454545454545*Power(y,3) + 14.545454545454545*Power(z,3)))\
          + Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.5)*
          (1.4545454545454546*Power(x,6) + 1.4545454545454546*Power(y,6) + 
            10.*Power(z,4) + 1.4545454545454546*Power(z,6) + 
            Power(y,4)*(10. + 1.4545454545454546*Power(z,2)) + 
            Power(x,4)*(10. + 1.4545454545454546*Power(y,2) + 
               1.4545454545454546*Power(z,2)) + 
            Power(y,2)*(20.*Power(z,2) + 1.4545454545454546*Power(z,4)) + 
            Power(x,2)*(20.*Power(y,2) + 1.4545454545454546*Power(y,4) + 
               20.*Power(z,2) + 1.4545454545454546*Power(z,4))) + 
         Power(M,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
          (7.2727272727272725*Power(x,6) + 7.2727272727272725*Power(y,6) + 
            4.181818181818182*Power(z,4) + 7.2727272727272725*Power(z,6) + 
            Power(y,4)*(4.181818181818182 + 7.2727272727272725*Power(z,2)) + 
            Power(x,4)*(4.181818181818182 + 7.2727272727272725*Power(y,2) + 
               7.2727272727272725*Power(z,2)) + 
            Power(y,2)*(8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4)) + 
            Power(x,2)*(8.363636363636363*Power(y,2) + 7.2727272727272725*Power(y,4) + 
               8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4))) + 
         Power(M,6)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            2.90909090909091*Power(y,3)*Power(z,3) - 6.909090909090909*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-6.909090909090909 + 16.*Power(z,2)) + 
            Power(x,4)*(-6.909090909090909 + 16.*Power(y,2) + 16.*Power(z,2)) + 
            Power(x,3)*(-2.90909090909091*Power(y,3) - 2.90909090909091*Power(z,3)) + 
            Power(y,2)*(-13.818181818181818*Power(z,2) + 16.*Power(z,4)) + 
            Power(x,2)*(-13.818181818181818*Power(y,2) + 16.*Power(y,4) - 
               13.818181818181818*Power(z,2) + 16.*Power(z,4))) + 
         Power(M,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            14.545454545454543*Power(y,3)*Power(z,3) - 11.272727272727273*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-11.272727272727273 + 21.818181818181817*Power(z,2)) + 
            Power(x,4)*(-11.272727272727273 + 21.818181818181817*Power(y,2) + 
               21.818181818181817*Power(z,2)) + 
            Power(x,3)*(-14.545454545454543*Power(y,3) - 
               14.545454545454543*Power(z,3)) + 
            Power(y,2)*(-22.545454545454547*Power(z,2) + 
               21.818181818181817*Power(z,4)) + 
            Power(x,2)*(-22.545454545454547*Power(y,2) + 21.818181818181817*Power(y,4) - 
               22.545454545454547*Power(z,2) + 21.818181818181817*Power(z,4))))/
       (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
         Power(1.*M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6)*
         (1.3636363636363635*Power(M,3) + 
           1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
           0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
           0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))),0.5));
}

static void
func_Pi13(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*x*z*(1.4545454545454546*Power(x,5)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
       1.4545454545454546*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
       1.4545454545454546*Power(y,3)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
       1.4545454545454546*Power(y,2)*Power(z,3)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
       1.4545454545454546*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
       0.36363636363636365*Power(y,2)*z*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
       0.36363636363636365*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       x*(-0.36363636363636365*Power(y,2) - 0.36363636363636365*Power(z,2))*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       Power(x,3)*(1.4545454545454546*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          1.4545454545454546*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
          0.36363636363636365*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)) + 
       Power(x,2)*(1.4545454545454546*Power(y,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          1.4545454545454546*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
          0.36363636363636365*z*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)) + 
       Power(M,9)*(-1.4545454545454546*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          1.4545454545454546*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          1.4545454545454546*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          1.4545454545454546*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          0.7272727272727273*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*
           (-1.4545454545454546*Power(y,3) - 1.4545454545454546*Power(z,3)) + 
          Power(y,3)*(-1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
          Power(x,3)*(-1.4545454545454546*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5))) + 
       Power(M,8)*(-15.272727272727272*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          15.272727272727272*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          15.272727272727272*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          15.272727272727272*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          8.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (-15.272727272727272*Power(y,3) - 15.272727272727272*Power(z,3)) + 
          Power(y,3)*(-15.272727272727272*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
          Power(x,3)*(-15.272727272727272*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             15.272727272727272*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))) + 
       Power(M,7)*(-61.81818181818181*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          61.81818181818181*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          61.81818181818181*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          23.636363636363637*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          x*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
           (5.454545454545455*Power(y,2) + 5.454545454545455*Power(z,2)) + 
          Power(y,3)*(-61.81818181818181*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             18.181818181818183*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(x,3)*(-61.81818181818181*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             61.81818181818181*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             23.636363636363637*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(y,2)*(-61.81818181818181*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             5.454545454545455*z*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(x,2)*(-61.81818181818181*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             61.81818181818181*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             5.454545454545455*z*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5))) + 
       Power(M,6)*(-121.45454545454545*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          121.45454545454545*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          121.45454545454545*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          9.818181818181818*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          x*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
           (36.72727272727273*Power(y,2) + 36.72727272727273*Power(z,2)) + 
          Power(y,3)*(-121.45454545454545*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             26.909090909090907*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(x,3)*(-121.45454545454545*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             121.45454545454545*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             9.818181818181818*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(y,2)*(-121.45454545454545*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             36.72727272727273*z*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(x,2)*(-121.45454545454545*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             121.45454545454545*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             36.72727272727273*z*Power(Power(x,2) + Power(y,2) + Power(z,2),11.))) + 
       Power(M,5)*(-111.27272727272727*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          111.27272727272727*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          111.27272727272727*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          1.090909090909091*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          x*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)*
           (25.818181818181817*Power(y,2) + 25.818181818181817*Power(z,2)) + 
          Power(y,3)*(-111.27272727272727*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             26.909090909090907*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(x,3)*(-111.27272727272727*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             111.27272727272727*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             1.090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(y,2)*(-111.27272727272727*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             25.818181818181817*z*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(x,2)*(-111.27272727272727*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             111.27272727272727*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             25.818181818181817*z*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5))) + 
       Power(M,4)*(-12.363636363636363*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          12.363636363636363*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          12.363636363636363*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          3.272727272727273*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          x*(-21.454545454545453*Power(y,2) - 21.454545454545453*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          Power(x,3)*(-12.363636363636363*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             12.363636363636363*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             3.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(y,3)*(-12.363636363636363*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             18.181818181818183*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(y,2)*(-12.363636363636363*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             21.454545454545453*z*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(x,2)*(-12.363636363636363*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             12.363636363636363*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             21.454545454545453*z*Power(Power(x,2) + Power(y,2) + Power(z,2),12.))) + 
       Power(M,3)*(64.72727272727272*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          64.72727272727272*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          64.72727272727272*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
          19.272727272727273*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          x*(-27.272727272727273*Power(y,2) - 27.272727272727273*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          Power(x,3)*(64.72727272727272*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             64.72727272727272*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             19.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(y,3)*(64.72727272727272*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(y,2)*(64.72727272727272*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             27.272727272727273*z*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(x,2)*(64.72727272727272*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             64.72727272727272*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             27.272727272727273*z*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5))) + 
       Power(M,2)*(54.54545454545455*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          54.54545454545455*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          54.54545454545455*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
          14.181818181818182*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          x*(-14.90909090909091*Power(y,2) - 14.90909090909091*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          Power(x,3)*(54.54545454545455*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             54.54545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             14.181818181818182*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(y,3)*(54.54545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(y,2)*(54.54545454545455*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             14.90909090909091*z*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(x,2)*(54.54545454545455*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             54.54545454545455*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             14.90909090909091*z*Power(Power(x,2) + Power(y,2) + Power(z,2),13.))) + 
       M*(16.727272727272727*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          16.727272727272727*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          16.727272727272727*Power(y,3)*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          16.727272727272727*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
          4.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          x*(-4.*Power(y,2) - 4.*Power(z,2))*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          Power(x,3)*(16.727272727272727*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             16.727272727272727*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             4.*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(y,2)*(16.727272727272727*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             4.*z*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(x,2)*(16.727272727272727*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             16.727272727272727*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             4.*z*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),13.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),8)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
          (0.09090909090909091*Power(x,4) + 0.09090909090909091*Power(y,4) + 
            0.18181818181818182*Power(y,2)*Power(z,2) + 0.09090909090909091*Power(z,4) + 
            Power(x,2)*(0.18181818181818182*Power(y,2) + 0.18181818181818182*Power(z,2)))
           + M*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)*
          (0.8181818181818181*Power(x,4) + 0.8181818181818181*Power(y,4) + 
            1.6363636363636362*Power(y,2)*Power(z,2) + 0.8181818181818181*Power(z,4) + 
            Power(x,2)*(1.6363636363636362*Power(y,2) + 1.6363636363636362*Power(z,2)))\
          + Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.5)*
          (3.2727272727272725*Power(x,4) + 3.2727272727272725*Power(y,4) + 
            6.545454545454545*Power(y,2)*Power(z,2) + 3.2727272727272725*Power(z,4) + 
            Power(x,2)*(6.545454545454545*Power(y,2) + 6.545454545454545*Power(z,2))) + 
         Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.)*
          (7.636363636363637*Power(x,4) + 7.636363636363637*Power(y,4) + 
            15.272727272727273*Power(y,2)*Power(z,2) + 7.636363636363637*Power(z,4) + 
            Power(x,2)*(15.272727272727273*Power(y,2) + 15.272727272727273*Power(z,2)))\
          + Power(M,9)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)*
          (-1.3636363636363635*Power(x,4) + 1.4545454545454546*Power(x,6) - 
            1.3636363636363635*Power(y,4) + 1.4545454545454546*Power(y,6) - 
            2.727272727272727*Power(y,2)*Power(z,2) + 
            2.909090909090909*Power(y,3)*Power(z,3) - 1.3636363636363635*Power(z,4) + 
            1.4545454545454546*Power(z,6) + 
            Power(x,2)*(-2.727272727272727*Power(y,2) - 2.727272727272727*Power(z,2)) + 
            Power(x,3)*(2.909090909090909*Power(y,3) + 2.909090909090909*Power(z,3))) + 
         Power(M,8)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)*
          (-6.454545454545454*Power(x,4) + 7.2727272727272725*Power(x,6) - 
            6.454545454545454*Power(y,4) + 7.2727272727272725*Power(y,6) - 
            12.909090909090908*Power(y,2)*Power(z,2) + 
            14.545454545454545*Power(y,3)*Power(z,3) - 6.454545454545454*Power(z,4) + 
            7.2727272727272725*Power(z,6) + 
            Power(x,2)*(-12.909090909090908*Power(y,2) - 
               12.909090909090908*Power(z,2)) + 
            Power(x,3)*(14.545454545454545*Power(y,3) + 14.545454545454545*Power(z,3)))\
          + Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.5)*
          (1.4545454545454546*Power(x,6) + 1.4545454545454546*Power(y,6) + 
            10.*Power(z,4) + 1.4545454545454546*Power(z,6) + 
            Power(y,4)*(10. + 1.4545454545454546*Power(z,2)) + 
            Power(x,4)*(10. + 1.4545454545454546*Power(y,2) + 
               1.4545454545454546*Power(z,2)) + 
            Power(y,2)*(20.*Power(z,2) + 1.4545454545454546*Power(z,4)) + 
            Power(x,2)*(20.*Power(y,2) + 1.4545454545454546*Power(y,4) + 
               20.*Power(z,2) + 1.4545454545454546*Power(z,4))) + 
         Power(M,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
          (7.2727272727272725*Power(x,6) + 7.2727272727272725*Power(y,6) + 
            4.181818181818182*Power(z,4) + 7.2727272727272725*Power(z,6) + 
            Power(y,4)*(4.181818181818182 + 7.2727272727272725*Power(z,2)) + 
            Power(x,4)*(4.181818181818182 + 7.2727272727272725*Power(y,2) + 
               7.2727272727272725*Power(z,2)) + 
            Power(y,2)*(8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4)) + 
            Power(x,2)*(8.363636363636363*Power(y,2) + 7.2727272727272725*Power(y,4) + 
               8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4))) + 
         Power(M,6)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            2.90909090909091*Power(y,3)*Power(z,3) - 6.909090909090909*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-6.909090909090909 + 16.*Power(z,2)) + 
            Power(x,4)*(-6.909090909090909 + 16.*Power(y,2) + 16.*Power(z,2)) + 
            Power(x,3)*(-2.90909090909091*Power(y,3) - 2.90909090909091*Power(z,3)) + 
            Power(y,2)*(-13.818181818181818*Power(z,2) + 16.*Power(z,4)) + 
            Power(x,2)*(-13.818181818181818*Power(y,2) + 16.*Power(y,4) - 
               13.818181818181818*Power(z,2) + 16.*Power(z,4))) + 
         Power(M,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            14.545454545454543*Power(y,3)*Power(z,3) - 11.272727272727273*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-11.272727272727273 + 21.818181818181817*Power(z,2)) + 
            Power(x,4)*(-11.272727272727273 + 21.818181818181817*Power(y,2) + 
               21.818181818181817*Power(z,2)) + 
            Power(x,3)*(-14.545454545454543*Power(y,3) - 
               14.545454545454543*Power(z,3)) + 
            Power(y,2)*(-22.545454545454547*Power(z,2) + 
               21.818181818181817*Power(z,4)) + 
            Power(x,2)*(-22.545454545454547*Power(y,2) + 21.818181818181817*Power(y,4) - 
               22.545454545454547*Power(z,2) + 21.818181818181817*Power(z,4))))/
       (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
         Power(1.*M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6)*
         (1.3636363636363635*Power(M,3) + 
           1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
           0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
           0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))),0.5));
}

static void
func_Pi22(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*(0.7272727272727273*Power(x,7)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       0.7272727272727273*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       1.4545454545454546*Power(y,5)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       0.7272727272727273*Power(y,4)*Power(z,3)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       0.7272727272727273*Power(y,3)*Power(z,4)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       1.4545454545454546*Power(y,2)*Power(z,5)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       0.7272727272727273*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
       0.7272727272727273*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),15.) - 
       0.7272727272727273*Power(y,3)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),15.) - 
       0.7272727272727273*Power(y,2)*Power(z,3)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),15.) - 
       0.7272727272727273*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),15.) + 
       Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)*
        (0.7272727272727273*Power(y,3) + 0.7272727272727273*Power(z,3)) + 
       Power(x,5)*(1.4545454545454546*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          1.4545454545454546*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
          0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),15.)) + 
       Power(x,3)*(0.7272727272727273*Power(y,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          0.7272727272727273*Power(z,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
          0.7272727272727273*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),15.) + 
          Power(y,2)*(1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),15.))) + 
       Power(x,2)*(1.4545454545454546*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          1.4545454545454546*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          1.4545454545454546*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
          0.7272727272727273*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),15.) + 
          Power(y,3)*(1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),15.))) + 
       Power(M,11)*(1.4545454545454546*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*
           (1.4545454545454546*Power(y,3) + 1.4545454545454546*Power(z,3)) + 
          Power(x,3)*Power(z,2)*(1.4545454545454546*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
          Power(x,5)*(1.4545454545454546*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             2.909090909090909*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
          Power(z,2)*(1.4545454545454546*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             1.4545454545454546*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             1.4545454545454546*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             0.7272727272727273*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             Power(y,3)*(1.4545454545454546*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
                0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5))) + 
          Power(x,2)*(1.4545454545454546*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             1.4545454545454546*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             2.909090909090909*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             0.7272727272727273*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             Power(y,3)*(2.909090909090909*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
                0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)))) + 
       Power(M,10)*(15.272727272727272*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (15.272727272727272*Power(y,3) + 15.272727272727272*Power(z,3)) + 
          Power(x,3)*Power(z,2)*(15.272727272727272*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             15.272727272727272*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
          Power(x,5)*(15.272727272727272*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             30.545454545454543*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
          Power(z,2)*(15.272727272727272*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             15.272727272727272*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             15.272727272727272*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             8.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             Power(y,3)*(15.272727272727272*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
                8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))) + 
          Power(x,2)*(15.272727272727272*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             15.272727272727272*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             30.545454545454543*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             8.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             Power(y,3)*(30.545454545454543*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
                8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)))) + 
       Power(M,9)*(72.72727272727273*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          10.90909090909091*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          10.90909090909091*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          72.72727272727273*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          40.*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*
           (72.72727272727273*Power(y,3) + 72.72727272727273*Power(z,3)) + 
          Power(x,5)*(83.63636363636364*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             145.45454545454547*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             40.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(y,5)*(83.63636363636364*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             10.90909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(y,3)*(72.72727272727273*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             29.090909090909093*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(y,2)*(83.63636363636364*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             21.81818181818182*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(x,2)*(83.63636363636364*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             83.63636363636364*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             145.45454545454547*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             40.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             Power(y,3)*(145.45454545454547*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
                29.090909090909093*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5))) + 
          Power(x,3)*(10.90909090909091*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             72.72727272727273*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             40.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             Power(y,2)*(83.63636363636364*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
                21.81818181818182*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)))) + 
       Power(M,8)*(207.27272727272728*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          85.81818181818183*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          85.81818181818183*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          207.27272727272728*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          120.*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)*
           (207.27272727272728*Power(y,3) + 207.27272727272728*Power(z,3)) + 
          Power(x,5)*(293.09090909090907*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             414.54545454545456*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             120.*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(y,5)*(293.09090909090907*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             73.45454545454545*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(y,3)*(207.27272727272728*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             46.54545454545455*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)
             ) + Power(y,2)*(293.09090909090907*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             146.9090909090909*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)
             ) + Power(x,3)*(85.81818181818183*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             207.27272727272728*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             120.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             Power(y,2)*(293.09090909090907*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
                146.9090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),11.))) + 
          Power(x,2)*(293.09090909090907*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             293.09090909090907*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             414.54545454545456*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             120.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             Power(y,3)*(414.54545454545456*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
                46.54545454545455*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)))) + 
       Power(M,7)*(392.72727272727275*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          281.45454545454544*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          281.45454545454544*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          392.72727272727275*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          240.*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
           (392.72727272727275*Power(y,3) + 392.72727272727275*Power(z,3)) + 
          Power(x,5)*(674.1818181818181*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             785.4545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             240.*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(y,5)*(674.1818181818181*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             215.27272727272728*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(y,3)*(392.72727272727275*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             188.36363636363635*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(y,2)*(674.1818181818181*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             266.9090909090909*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(x,3)*(281.45454545454544*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             392.72727272727275*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             240.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             Power(y,2)*(674.1818181818181*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
                266.9090909090909*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5))) + 
          Power(x,2)*(674.1818181818181*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             674.1818181818181*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             785.4545454545455*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             240.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             Power(y,3)*(785.4545454545455*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
                188.36363636363635*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)))) + 
       Power(M,6)*(519.2727272727273*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          506.90909090909093*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          506.90909090909093*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          519.2727272727273*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          336.*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
           (519.2727272727273*Power(y,3) + 519.2727272727273*Power(z,3)) + 
          Power(y,5)*(1026.181818181818*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             360.7272727272727*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(x,5)*(1026.181818181818*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             1038.5454545454545*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             336.*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(y,3)*(519.2727272727273*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             378.90909090909093*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(y,2)*(1026.181818181818*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             317.8181818181818*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)
             ) + Power(x,2)*(1026.181818181818*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             1026.181818181818*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             1038.5454545454545*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             336.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             Power(y,3)*(1038.5454545454545*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
                378.90909090909093*Power(Power(x,2) + Power(y,2) + Power(z,2),12.))) + 
          Power(x,3)*(506.90909090909093*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             519.2727272727273*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             336.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             Power(y,2)*(1026.181818181818*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
                317.8181818181818*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)))) + 
       Power(M,5)*(488.7272727272727*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          553.4545454545455*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          553.4545454545455*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          488.7272727272727*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
          336.*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)*
           (488.7272727272727*Power(y,3) + 488.7272727272727*Power(z,3)) + 
          Power(y,5)*(1042.1818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             382.54545454545456*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(x,5)*(1042.1818181818182*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             977.4545454545454*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             336.*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(y,3)*(488.7272727272727*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             390.5454545454545*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(y,2)*(1042.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             328.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(x,2)*(1042.1818181818182*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             1042.1818181818182*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             977.4545454545454*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             336.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             Power(y,3)*(977.4545454545454*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
                390.5454545454545*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5))) + 
          Power(x,3)*(553.4545454545455*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             488.7272727272727*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             336.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             Power(y,2)*(1042.1818181818182*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
                328.*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)))) + 
       Power(M,4)*(327.27272727272725*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          381.8181818181818*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          381.8181818181818*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          327.27272727272725*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
          240.*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)*
           (327.27272727272725*Power(y,3) + 327.27272727272725*Power(z,3)) + 
          Power(y,5)*(709.0909090909091*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             269.09090909090907*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(x,5)*(709.0909090909091*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             654.5454545454545*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             240.*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(y,3)*(327.27272727272725*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             269.8181818181818*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)
             ) + Power(y,2)*(709.0909090909091*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             239.27272727272728*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(x,2)*(709.0909090909091*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             709.0909090909091*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             654.5454545454545*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             240.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             Power(y,3)*(654.5454545454545*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
                269.8181818181818*Power(Power(x,2) + Power(y,2) + Power(z,2),13.))) + 
          Power(x,3)*(381.8181818181818*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             327.27272727272725*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             240.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             Power(y,2)*(709.0909090909091*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
                239.27272727272728*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)))) + 
       Power(M,3)*(152.72727272727272*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          169.45454545454544*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          169.45454545454544*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          152.72727272727272*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
          120.*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)*
           (152.72727272727272*Power(y,3) + 152.72727272727272*Power(z,3)) + 
          Power(y,5)*(322.1818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             128.*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(x,5)*(322.1818181818182*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             305.45454545454544*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             120.*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(y,3)*(152.72727272727272*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             128.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(y,2)*(322.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             120.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(x,2)*(322.1818181818182*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             322.1818181818182*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             305.45454545454544*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             120.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             Power(y,3)*(305.45454545454544*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
                128.*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5))) + 
          Power(x,3)*(169.45454545454544*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             152.72727272727272*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             120.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             Power(y,2)*(322.1818181818182*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
                120.*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)))) + 
       Power(M,2)*(47.27272727272727*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          48.72727272727273*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          48.72727272727273*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          47.27272727272727*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
          40.*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)*
           (47.27272727272727*Power(y,3) + 47.27272727272727*Power(z,3)) + 
          Power(y,5)*(96.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.72727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)) + 
          Power(x,5)*(96.*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             94.54545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)) + 
          Power(y,3)*(47.27272727272727*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.72727272727273*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)
             ) + Power(y,2)*(96.*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)) + 
          Power(x,2)*(96.*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             96.*Power(y,2)*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             94.54545454545455*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
             Power(y,3)*(94.54545454545455*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
                40.72727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),14.))) + 
          Power(x,3)*(48.72727272727273*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             47.27272727272727*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
             Power(y,2)*(96.*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
                40.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)))) + 
       M*(8.727272727272727*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          8.727272727272727*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          8.727272727272727*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          8.727272727272727*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
          8.*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)*
           (8.727272727272727*Power(y,3) + 8.727272727272727*Power(z,3)) + 
          Power(y,5)*(17.454545454545453*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5)) + 
          Power(x,5)*(17.454545454545453*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             17.454545454545453*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5)) + 
          Power(y,3)*(8.727272727272727*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5)) + 
          Power(y,2)*(17.454545454545453*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5)) + 
          Power(x,3)*(8.727272727272727*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             8.727272727272727*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5) + 
             Power(y,2)*(17.454545454545453*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
                8.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5))) + 
          Power(x,2)*(17.454545454545453*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             17.454545454545453*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             17.454545454545453*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5) + 
             Power(y,3)*(17.454545454545453*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
                8.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5))))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),13.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),8)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
          (0.09090909090909091*Power(x,4) + 0.09090909090909091*Power(y,4) + 
            0.18181818181818182*Power(y,2)*Power(z,2) + 0.09090909090909091*Power(z,4) + 
            Power(x,2)*(0.18181818181818182*Power(y,2) + 0.18181818181818182*Power(z,2)))
           + M*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)*
          (0.8181818181818181*Power(x,4) + 0.8181818181818181*Power(y,4) + 
            1.6363636363636362*Power(y,2)*Power(z,2) + 0.8181818181818181*Power(z,4) + 
            Power(x,2)*(1.6363636363636362*Power(y,2) + 1.6363636363636362*Power(z,2)))\
          + Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.5)*
          (3.2727272727272725*Power(x,4) + 3.2727272727272725*Power(y,4) + 
            6.545454545454545*Power(y,2)*Power(z,2) + 3.2727272727272725*Power(z,4) + 
            Power(x,2)*(6.545454545454545*Power(y,2) + 6.545454545454545*Power(z,2))) + 
         Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.)*
          (7.636363636363637*Power(x,4) + 7.636363636363637*Power(y,4) + 
            15.272727272727273*Power(y,2)*Power(z,2) + 7.636363636363637*Power(z,4) + 
            Power(x,2)*(15.272727272727273*Power(y,2) + 15.272727272727273*Power(z,2)))\
          + Power(M,9)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)*
          (-1.3636363636363635*Power(x,4) + 1.4545454545454546*Power(x,6) - 
            1.3636363636363635*Power(y,4) + 1.4545454545454546*Power(y,6) - 
            2.727272727272727*Power(y,2)*Power(z,2) + 
            2.909090909090909*Power(y,3)*Power(z,3) - 1.3636363636363635*Power(z,4) + 
            1.4545454545454546*Power(z,6) + 
            Power(x,2)*(-2.727272727272727*Power(y,2) - 2.727272727272727*Power(z,2)) + 
            Power(x,3)*(2.909090909090909*Power(y,3) + 2.909090909090909*Power(z,3))) + 
         Power(M,8)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)*
          (-6.454545454545454*Power(x,4) + 7.2727272727272725*Power(x,6) - 
            6.454545454545454*Power(y,4) + 7.2727272727272725*Power(y,6) - 
            12.909090909090908*Power(y,2)*Power(z,2) + 
            14.545454545454545*Power(y,3)*Power(z,3) - 6.454545454545454*Power(z,4) + 
            7.2727272727272725*Power(z,6) + 
            Power(x,2)*(-12.909090909090908*Power(y,2) - 
               12.909090909090908*Power(z,2)) + 
            Power(x,3)*(14.545454545454545*Power(y,3) + 14.545454545454545*Power(z,3)))\
          + Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.5)*
          (1.4545454545454546*Power(x,6) + 1.4545454545454546*Power(y,6) + 
            10.*Power(z,4) + 1.4545454545454546*Power(z,6) + 
            Power(y,4)*(10. + 1.4545454545454546*Power(z,2)) + 
            Power(x,4)*(10. + 1.4545454545454546*Power(y,2) + 
               1.4545454545454546*Power(z,2)) + 
            Power(y,2)*(20.*Power(z,2) + 1.4545454545454546*Power(z,4)) + 
            Power(x,2)*(20.*Power(y,2) + 1.4545454545454546*Power(y,4) + 
               20.*Power(z,2) + 1.4545454545454546*Power(z,4))) + 
         Power(M,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
          (7.2727272727272725*Power(x,6) + 7.2727272727272725*Power(y,6) + 
            4.181818181818182*Power(z,4) + 7.2727272727272725*Power(z,6) + 
            Power(y,4)*(4.181818181818182 + 7.2727272727272725*Power(z,2)) + 
            Power(x,4)*(4.181818181818182 + 7.2727272727272725*Power(y,2) + 
               7.2727272727272725*Power(z,2)) + 
            Power(y,2)*(8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4)) + 
            Power(x,2)*(8.363636363636363*Power(y,2) + 7.2727272727272725*Power(y,4) + 
               8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4))) + 
         Power(M,6)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            2.90909090909091*Power(y,3)*Power(z,3) - 6.909090909090909*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-6.909090909090909 + 16.*Power(z,2)) + 
            Power(x,4)*(-6.909090909090909 + 16.*Power(y,2) + 16.*Power(z,2)) + 
            Power(x,3)*(-2.90909090909091*Power(y,3) - 2.90909090909091*Power(z,3)) + 
            Power(y,2)*(-13.818181818181818*Power(z,2) + 16.*Power(z,4)) + 
            Power(x,2)*(-13.818181818181818*Power(y,2) + 16.*Power(y,4) - 
               13.818181818181818*Power(z,2) + 16.*Power(z,4))) + 
         Power(M,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            14.545454545454543*Power(y,3)*Power(z,3) - 11.272727272727273*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-11.272727272727273 + 21.818181818181817*Power(z,2)) + 
            Power(x,4)*(-11.272727272727273 + 21.818181818181817*Power(y,2) + 
               21.818181818181817*Power(z,2)) + 
            Power(x,3)*(-14.545454545454543*Power(y,3) - 
               14.545454545454543*Power(z,3)) + 
            Power(y,2)*(-22.545454545454547*Power(z,2) + 
               21.818181818181817*Power(z,4)) + 
            Power(x,2)*(-22.545454545454547*Power(y,2) + 21.818181818181817*Power(y,4) - 
               22.545454545454547*Power(z,2) + 21.818181818181817*Power(z,4))))/
       (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
         Power(1.*M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6)*
         (1.3636363636363635*Power(M,3) + 
           1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
           0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
           0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))),0.5));
}

static void
func_Pi23(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,4)*y*z*(1.4545454545454546*Power(x,5)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
       1.4545454545454546*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
       1.4545454545454546*Power(y,3)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
       1.4545454545454546*Power(y,2)*Power(z,3)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
       1.4545454545454546*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
       0.36363636363636365*Power(y,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
       0.36363636363636365*Power(y,2)*z*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
       0.36363636363636365*y*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
       0.36363636363636365*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       Power(x,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)*
        (1.4545454545454546*Power(y,2) + 1.4545454545454546*Power(z,2)) + 
       Power(x,2)*(1.4545454545454546*Power(y,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          1.4545454545454546*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
          0.36363636363636365*y*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
          0.36363636363636365*z*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)) + 
       Power(M,9)*(-1.4545454545454546*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          1.4545454545454546*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          1.4545454545454546*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
          1.4545454545454546*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          0.7272727272727273*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*
           (-1.4545454545454546*Power(y,3) - 1.4545454545454546*Power(z,3)) + 
          Power(y,3)*(-1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
          Power(x,3)*(-1.4545454545454546*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5))) + 
       Power(M,8)*(-15.272727272727273*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          15.272727272727273*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          15.272727272727273*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
          15.272727272727273*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          8.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (-15.272727272727273*Power(y,3) - 15.272727272727273*Power(z,3)) + 
          Power(y,3)*(-15.272727272727273*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
          Power(x,3)*(-15.272727272727273*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             15.272727272727273*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))) + 
       Power(M,7)*(-61.81818181818182*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          61.81818181818182*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          61.81818181818182*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          5.454545454545455*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          23.636363636363637*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          Power(x,3)*(-61.81818181818182*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             61.81818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             18.181818181818183*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(y,3)*(-61.81818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             23.636363636363637*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(y,2)*(-61.81818181818182*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             5.454545454545455*z*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(x,2)*(-61.81818181818182*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             61.81818181818182*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             5.454545454545455*y*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             5.454545454545455*z*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5))) + 
       Power(M,6)*(-121.45454545454545*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          121.45454545454545*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          121.45454545454545*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          36.72727272727273*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          9.818181818181818*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          Power(x,3)*(-121.45454545454545*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             121.45454545454545*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             26.90909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(y,3)*(-121.45454545454545*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             9.818181818181818*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(y,2)*(-121.45454545454545*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             36.72727272727273*z*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(x,2)*(-121.45454545454545*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             121.45454545454545*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             36.72727272727273*y*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             36.72727272727273*z*Power(Power(x,2) + Power(y,2) + Power(z,2),11.))) + 
       Power(M,5)*(-111.27272727272728*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          111.27272727272728*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          111.27272727272728*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          25.81818181818182*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
          1.0909090909090908*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          Power(x,3)*(-111.27272727272728*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             111.27272727272728*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             26.90909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(y,3)*(-111.27272727272728*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             1.0909090909090908*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(y,2)*(-111.27272727272728*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             25.81818181818182*z*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(x,2)*(-111.27272727272728*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             111.27272727272728*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             25.81818181818182*y*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             25.81818181818182*z*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5))) + 
       Power(M,4)*(-12.363636363636363*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          12.363636363636363*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          12.363636363636363*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          21.454545454545457*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
          3.272727272727273*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          Power(y,3)*(-12.363636363636363*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             3.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(x,3)*(-12.363636363636363*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             12.363636363636363*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             18.181818181818183*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(y,2)*(-12.363636363636363*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             21.454545454545457*z*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(x,2)*(-12.363636363636363*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             12.363636363636363*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             21.454545454545457*y*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             21.454545454545457*z*Power(Power(x,2) + Power(y,2) + Power(z,2),12.))) + 
       Power(M,3)*(64.72727272727273*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          64.72727272727273*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          64.72727272727273*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
          27.272727272727273*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
          19.272727272727273*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          Power(y,3)*(64.72727272727273*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             19.272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(x,3)*(64.72727272727273*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             64.72727272727273*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(y,2)*(64.72727272727273*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             27.272727272727273*z*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(x,2)*(64.72727272727273*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             64.72727272727273*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             27.272727272727273*y*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             27.272727272727273*z*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5))) + 
       Power(M,2)*(54.54545454545455*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          54.54545454545455*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          54.54545454545455*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
          14.90909090909091*y*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
          14.181818181818182*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          Power(y,3)*(54.54545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             14.181818181818182*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(x,3)*(54.54545454545455*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             54.54545454545455*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(y,2)*(54.54545454545455*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             14.90909090909091*z*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(x,2)*(54.54545454545455*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             54.54545454545455*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             14.90909090909091*y*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             14.90909090909091*z*Power(Power(x,2) + Power(y,2) + Power(z,2),13.))) + 
       M*(16.727272727272727*Power(x,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          16.727272727272727*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          16.727272727272727*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
          4.*y*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
          4.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          Power(x,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)*
           (16.727272727272727*Power(y,2) + 16.727272727272727*Power(z,2)) + 
          Power(y,3)*(16.727272727272727*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             4.*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(y,2)*(16.727272727272727*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             4.*z*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(x,2)*(16.727272727272727*Power(y,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             16.727272727272727*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             4.*y*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             4.*z*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),13.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),8)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
          (0.09090909090909091*Power(x,4) + 0.09090909090909091*Power(y,4) + 
            0.18181818181818182*Power(y,2)*Power(z,2) + 0.09090909090909091*Power(z,4) + 
            Power(x,2)*(0.18181818181818182*Power(y,2) + 0.18181818181818182*Power(z,2)))
           + M*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)*
          (0.8181818181818181*Power(x,4) + 0.8181818181818181*Power(y,4) + 
            1.6363636363636362*Power(y,2)*Power(z,2) + 0.8181818181818181*Power(z,4) + 
            Power(x,2)*(1.6363636363636362*Power(y,2) + 1.6363636363636362*Power(z,2)))\
          + Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.5)*
          (3.2727272727272725*Power(x,4) + 3.2727272727272725*Power(y,4) + 
            6.545454545454545*Power(y,2)*Power(z,2) + 3.2727272727272725*Power(z,4) + 
            Power(x,2)*(6.545454545454545*Power(y,2) + 6.545454545454545*Power(z,2))) + 
         Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.)*
          (7.636363636363637*Power(x,4) + 7.636363636363637*Power(y,4) + 
            15.272727272727273*Power(y,2)*Power(z,2) + 7.636363636363637*Power(z,4) + 
            Power(x,2)*(15.272727272727273*Power(y,2) + 15.272727272727273*Power(z,2)))\
          + Power(M,9)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)*
          (-1.3636363636363635*Power(x,4) + 1.4545454545454546*Power(x,6) - 
            1.3636363636363635*Power(y,4) + 1.4545454545454546*Power(y,6) - 
            2.727272727272727*Power(y,2)*Power(z,2) + 
            2.909090909090909*Power(y,3)*Power(z,3) - 1.3636363636363635*Power(z,4) + 
            1.4545454545454546*Power(z,6) + 
            Power(x,2)*(-2.727272727272727*Power(y,2) - 2.727272727272727*Power(z,2)) + 
            Power(x,3)*(2.909090909090909*Power(y,3) + 2.909090909090909*Power(z,3))) + 
         Power(M,8)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)*
          (-6.454545454545454*Power(x,4) + 7.2727272727272725*Power(x,6) - 
            6.454545454545454*Power(y,4) + 7.2727272727272725*Power(y,6) - 
            12.909090909090908*Power(y,2)*Power(z,2) + 
            14.545454545454545*Power(y,3)*Power(z,3) - 6.454545454545454*Power(z,4) + 
            7.2727272727272725*Power(z,6) + 
            Power(x,2)*(-12.909090909090908*Power(y,2) - 
               12.909090909090908*Power(z,2)) + 
            Power(x,3)*(14.545454545454545*Power(y,3) + 14.545454545454545*Power(z,3)))\
          + Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.5)*
          (1.4545454545454546*Power(x,6) + 1.4545454545454546*Power(y,6) + 
            10.*Power(z,4) + 1.4545454545454546*Power(z,6) + 
            Power(y,4)*(10. + 1.4545454545454546*Power(z,2)) + 
            Power(x,4)*(10. + 1.4545454545454546*Power(y,2) + 
               1.4545454545454546*Power(z,2)) + 
            Power(y,2)*(20.*Power(z,2) + 1.4545454545454546*Power(z,4)) + 
            Power(x,2)*(20.*Power(y,2) + 1.4545454545454546*Power(y,4) + 
               20.*Power(z,2) + 1.4545454545454546*Power(z,4))) + 
         Power(M,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
          (7.2727272727272725*Power(x,6) + 7.2727272727272725*Power(y,6) + 
            4.181818181818182*Power(z,4) + 7.2727272727272725*Power(z,6) + 
            Power(y,4)*(4.181818181818182 + 7.2727272727272725*Power(z,2)) + 
            Power(x,4)*(4.181818181818182 + 7.2727272727272725*Power(y,2) + 
               7.2727272727272725*Power(z,2)) + 
            Power(y,2)*(8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4)) + 
            Power(x,2)*(8.363636363636363*Power(y,2) + 7.2727272727272725*Power(y,4) + 
               8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4))) + 
         Power(M,6)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            2.90909090909091*Power(y,3)*Power(z,3) - 6.909090909090909*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-6.909090909090909 + 16.*Power(z,2)) + 
            Power(x,4)*(-6.909090909090909 + 16.*Power(y,2) + 16.*Power(z,2)) + 
            Power(x,3)*(-2.90909090909091*Power(y,3) - 2.90909090909091*Power(z,3)) + 
            Power(y,2)*(-13.818181818181818*Power(z,2) + 16.*Power(z,4)) + 
            Power(x,2)*(-13.818181818181818*Power(y,2) + 16.*Power(y,4) - 
               13.818181818181818*Power(z,2) + 16.*Power(z,4))) + 
         Power(M,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            14.545454545454543*Power(y,3)*Power(z,3) - 11.272727272727273*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-11.272727272727273 + 21.818181818181817*Power(z,2)) + 
            Power(x,4)*(-11.272727272727273 + 21.818181818181817*Power(y,2) + 
               21.818181818181817*Power(z,2)) + 
            Power(x,3)*(-14.545454545454543*Power(y,3) - 
               14.545454545454543*Power(z,3)) + 
            Power(y,2)*(-22.545454545454547*Power(z,2) + 
               21.818181818181817*Power(z,4)) + 
            Power(x,2)*(-22.545454545454547*Power(y,2) + 21.818181818181817*Power(y,4) - 
               22.545454545454547*Power(z,2) + 21.818181818181817*Power(z,4))))/
       (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
         Power(1.*M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6)*
         (1.3636363636363635*Power(M,3) + 
           1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
           0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
           0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))),0.5));
}

static void
func_Pi33(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*(0.7272727272727273*Power(x,7)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       0.7272727272727273*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       1.4545454545454546*Power(y,5)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       0.727272727272727*Power(y,4)*Power(z,3)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       0.7272727272727273*Power(y,3)*Power(z,4)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       1.4545454545454546*Power(y,2)*Power(z,5)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
       0.7272727272727273*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
       0.7272727272727273*Power(y,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),15.) - 
       0.7272727272727273*Power(y,3)*Power(z,2)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),15.) - 
       0.7272727272727273*Power(y,2)*Power(z,3)*
        Power(Power(x,2) + Power(y,2) + Power(z,2),15.) - 
       0.7272727272727273*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),15.) + 
       Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)*
        (0.7272727272727273*Power(y,3) + 0.727272727272727*Power(z,3)) + 
       Power(x,5)*(1.4545454545454546*Power(y,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          1.4545454545454546*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
          0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),15.)) + 
       Power(x,3)*(0.7272727272727273*Power(y,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          0.7272727272727273*Power(z,4)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
          0.7272727272727273*Power(z,2)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),15.) + 
          Power(y,2)*(1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),15.))) + 
       Power(x,2)*(1.4545454545454546*Power(y,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          1.454545454545454*Power(y,2)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          1.4545454545454546*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
          0.7272727272727273*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),15.) + 
          Power(y,3)*(1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.) - 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),15.))) + 
       Power(M,11)*(1.4545454545454546*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*
           (1.4545454545454546*Power(y,3) + 1.4545454545454546*Power(z,3)) + 
          Power(x,3)*Power(y,2)*(1.4545454545454546*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
          Power(x,5)*(2.909090909090909*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             1.4545454545454546*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)) + 
          Power(y,2)*(1.4545454545454546*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             1.4545454545454546*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             1.4545454545454546*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             0.7272727272727273*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             Power(y,3)*(1.4545454545454546*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
                0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5))) + 
          Power(x,2)*(2.909090909090909*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             2.909090909090909*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) + 
             1.4545454545454546*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
             0.7272727272727273*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             Power(y,3)*(1.4545454545454546*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),8.5) - 
                0.7272727272727273*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)))) + 
       Power(M,10)*(15.272727272727273*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*
           (15.272727272727273*Power(y,3) + 15.272727272727273*Power(z,3)) + 
          Power(x,3)*Power(y,2)*(15.272727272727273*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             15.272727272727273*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
          Power(x,5)*(30.545454545454547*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             15.272727272727273*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)) + 
          Power(y,2)*(15.272727272727273*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             15.272727272727273*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             15.272727272727273*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             8.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             Power(y,3)*(15.272727272727273*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
                8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.))) + 
          Power(x,2)*(30.545454545454547*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             30.545454545454547*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) + 
             15.272727272727273*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
             8.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             Power(y,3)*(15.272727272727273*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),9.) - 
                8.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)))) + 
       Power(M,9)*(72.72727272727272*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          72.72727272727272*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          72.72727272727272*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
          10.90909090909091*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
          10.90909090909091*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*
           (72.72727272727272*Power(y,3) + 72.72727272727272*Power(z,3)) + 
          Power(y,5)*(83.63636363636363*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             40.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(x,5)*(145.45454545454544*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             83.63636363636363*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             40.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(y,3)*(10.90909090909091*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             21.81818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(y,2)*(83.63636363636363*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             29.090909090909086*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)) + 
          Power(x,3)*(72.72727272727273*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             10.90909090909091*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             21.81818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             Power(y,2)*(83.63636363636363*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
                40.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5))) + 
          Power(x,2)*(145.45454545454544*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             145.45454545454544*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) + 
             83.63636363636363*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
             29.090909090909086*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             Power(y,3)*(83.63636363636363*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),9.5) - 
                40.*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)))) + 
       Power(M,8)*(207.27272727272728*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          207.27272727272728*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          207.27272727272728*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
          85.81818181818183*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
          73.45454545454545*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.)*
           (207.27272727272728*Power(y,3) + 207.27272727272728*Power(z,3)) + 
          Power(y,5)*(293.0909090909091*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             120.00000000000001*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(x,5)*(414.54545454545456*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             293.0909090909091*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             120.00000000000001*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)) + 
          Power(y,3)*(85.81818181818183*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             146.9090909090909*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)
             ) + Power(y,2)*(293.0909090909091*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             46.54545454545455*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)
             ) + Power(x,3)*(207.27272727272728*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             85.81818181818183*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             146.9090909090909*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             Power(y,2)*(293.0909090909091*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
                120.00000000000001*Power(Power(x,2) + Power(y,2) + Power(z,2),11.))) + 
          Power(x,2)*(414.54545454545456*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             414.54545454545456*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) + 
             293.0909090909091*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
             46.54545454545455*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             Power(y,3)*(293.0909090909091*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),10.) - 
                120.00000000000001*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)))) + 
       Power(M,7)*(392.72727272727263*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          392.72727272727263*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          392.72727272727263*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
          281.45454545454544*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
          215.2727272727273*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),10.5)*
           (392.72727272727263*Power(y,3) + 392.72727272727263*Power(z,3)) + 
          Power(y,5)*(674.1818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             239.99999999999997*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(x,5)*(785.4545454545453*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             674.1818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             239.99999999999997*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(y,3)*(281.45454545454544*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             266.90909090909093*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(y,2)*(674.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             188.36363636363635*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)) + 
          Power(x,3)*(392.72727272727263*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             281.45454545454544*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             266.90909090909093*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             Power(y,2)*(674.1818181818182*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
                239.99999999999997*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5))) + 
          Power(x,2)*(785.4545454545453*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             785.4545454545453*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) + 
             674.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
             188.36363636363635*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             Power(y,3)*(674.1818181818182*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),10.5) - 
                239.99999999999997*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)))) + 
       Power(M,6)*(519.2727272727274*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          519.2727272727274*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          519.2727272727273*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
          506.909090909091*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
          360.72727272727275*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.)*
           (519.2727272727274*Power(y,3) + 519.2727272727273*Power(z,3)) + 
          Power(y,5)*(1026.1818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             336.00000000000006*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(x,5)*(1038.5454545454547*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             1026.1818181818182*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             336.00000000000006*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(y,3)*(506.909090909091*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             317.81818181818187*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(y,2)*(1026.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             378.90909090909093*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.)) + 
          Power(x,3)*(519.2727272727274*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             506.909090909091*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             317.81818181818187*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             Power(y,2)*(1026.1818181818182*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
                336.00000000000006*Power(Power(x,2) + Power(y,2) + Power(z,2),12.))) + 
          Power(x,2)*(1038.5454545454547*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             1038.5454545454545*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) + 
             1026.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
             378.90909090909093*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             Power(y,3)*(1026.1818181818182*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),11.) - 
                336.00000000000006*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)))) + 
       Power(M,5)*(488.7272727272728*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          488.7272727272728*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          488.7272727272728*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
          553.4545454545455*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
          382.5454545454546*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),11.5)*
           (488.7272727272728*Power(y,3) + 488.7272727272728*Power(z,3)) + 
          Power(y,5)*(1042.1818181818185*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             336.00000000000006*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(x,5)*(977.4545454545456*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             1042.1818181818185*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             336.00000000000006*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(y,3)*(553.4545454545455*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             327.99999999999994*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(y,2)*(1042.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             390.54545454545456*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)) + 
          Power(x,3)*(488.7272727272728*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             553.4545454545455*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             327.99999999999994*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             Power(y,2)*(1042.1818181818185*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
                336.00000000000006*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5))) + 
          Power(x,2)*(977.4545454545456*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             977.4545454545456*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) + 
             1042.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
             390.54545454545456*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             Power(y,3)*(1042.1818181818185*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),11.5) - 
                336.00000000000006*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)))) + 
       Power(M,4)*(327.27272727272725*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          327.27272727272725*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          327.2727272727274*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
          381.81818181818187*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
          269.0909090909091*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.)*
           (327.27272727272725*Power(y,3) + 327.2727272727274*Power(z,3)) + 
          Power(y,5)*(709.0909090909088*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             239.99999999999997*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(x,5)*(654.5454545454544*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             709.0909090909088*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             239.99999999999997*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(y,3)*(381.81818181818187*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             239.27272727272728*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.)) + 
          Power(y,2)*(709.0909090909091*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             269.8181818181819*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)
             ) + Power(x,2)*(654.5454545454545*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             654.5454545454548*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             709.0909090909091*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             269.8181818181819*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             Power(y,3)*(709.090909090909*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
                240.00000000000003*Power(Power(x,2) + Power(y,2) + Power(z,2),13.))) + 
          Power(x,3)*(327.2727272727272*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) + 
             381.81818181818187*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
             239.27272727272728*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             Power(y,2)*(709.090909090909*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),12.) - 
                239.99999999999994*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)))) + 
       Power(M,3)*(152.72727272727272*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          152.72727272727272*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          152.72727272727272*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
          169.4545454545455*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
          128.*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),12.5)*
           (152.72727272727272*Power(y,3) + 152.72727272727272*Power(z,3)) + 
          Power(y,5)*(322.18181818181813*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             119.99999999999999*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(x,5)*(305.45454545454544*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             322.18181818181813*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             119.99999999999999*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(y,3)*(169.4545454545455*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             120.00000000000001*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(y,2)*(322.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             128.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)) + 
          Power(x,3)*(152.72727272727272*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             169.4545454545455*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             120.00000000000001*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             Power(y,2)*(322.18181818181813*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
                119.99999999999999*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5))) + 
          Power(x,2)*(305.45454545454544*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             305.45454545454544*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) + 
             322.1818181818182*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
             128.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             Power(y,3)*(322.18181818181813*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),12.5) - 
                119.99999999999999*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)))) + 
       Power(M,2)*(47.27272727272727*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          47.27272727272727*Power(y,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          47.272727272727245*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
          48.72727272727273*Power(z,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
          40.727272727272734*Power(z,5)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.)*
           (47.27272727272727*Power(y,3) + 47.272727272727245*Power(z,3)) + 
          Power(y,5)*(96.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)) + 
          Power(x,5)*(94.54545454545455*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             96.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)) + 
          Power(y,3)*(48.72727272727273*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.00000000000001*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)
             ) + Power(y,2)*(96.00000000000001*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.72727272727272*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)
             ) + Power(x,3)*(47.27272727272727*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             48.72727272727273*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.00000000000001*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
             Power(y,2)*(96.*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
                40.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.))) + 
          Power(x,2)*(94.54545454545455*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             94.54545454545449*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) + 
             96.00000000000001*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
             40.72727272727272*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),14.) + 
             Power(y,3)*(96.*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),13.) - 
                40.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.)))) + 
       M*(8.727272727272727*Power(x,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          8.727272727272727*Power(y,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          8.727272727272725*Power(y,4)*Power(z,3)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
          8.727272727272727*Power(z,7)*
           Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
          8.*Power(z,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5) + 
          Power(x,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),13.5)*
           (8.727272727272727*Power(y,3) + 8.727272727272725*Power(z,3)) + 
          Power(y,5)*(17.454545454545453*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5)) + 
          Power(x,5)*(17.454545454545453*Power(y,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             17.454545454545453*Power(z,2)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5)) + 
          Power(y,3)*(8.727272727272727*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5)) + 
          Power(y,2)*(17.454545454545453*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5)) + 
          Power(x,3)*(8.727272727272727*Power(y,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             8.727272727272727*Power(z,4)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5) + 
             Power(y,2)*(17.454545454545453*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
                8.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5))) + 
          Power(x,2)*(17.454545454545453*Power(y,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             17.45454545454545*Power(y,2)*Power(z,3)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) + 
             17.454545454545453*Power(z,5)*
              Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
             8.*Power(z,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5) + 
             Power(y,3)*(17.454545454545453*Power(z,2)*
                 Power(Power(x,2) + Power(y,2) + Power(z,2),13.5) - 
                8.*Power(Power(x,2) + Power(y,2) + Power(z,2),14.5))))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),13.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),8)*
     (1.3636363636363635*Power(M,3) + 
       1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
       0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))*
     Power((Power(Power(x,2) + Power(y,2) + Power(z,2),5.5)*
          (0.09090909090909091*Power(x,4) + 0.09090909090909091*Power(y,4) + 
            0.18181818181818182*Power(y,2)*Power(z,2) + 0.09090909090909091*Power(z,4) + 
            Power(x,2)*(0.18181818181818182*Power(y,2) + 0.18181818181818182*Power(z,2)))
           + M*Power(Power(x,2) + Power(y,2) + Power(z,2),5.)*
          (0.8181818181818181*Power(x,4) + 0.8181818181818181*Power(y,4) + 
            1.6363636363636362*Power(y,2)*Power(z,2) + 0.8181818181818181*Power(z,4) + 
            Power(x,2)*(1.6363636363636362*Power(y,2) + 1.6363636363636362*Power(z,2)))\
          + Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.5)*
          (3.2727272727272725*Power(x,4) + 3.2727272727272725*Power(y,4) + 
            6.545454545454545*Power(y,2)*Power(z,2) + 3.2727272727272725*Power(z,4) + 
            Power(x,2)*(6.545454545454545*Power(y,2) + 6.545454545454545*Power(z,2))) + 
         Power(M,3)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.)*
          (7.636363636363637*Power(x,4) + 7.636363636363637*Power(y,4) + 
            15.272727272727273*Power(y,2)*Power(z,2) + 7.636363636363637*Power(z,4) + 
            Power(x,2)*(15.272727272727273*Power(y,2) + 15.272727272727273*Power(z,2)))\
          + Power(M,9)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.)*
          (-1.3636363636363635*Power(x,4) + 1.4545454545454546*Power(x,6) - 
            1.3636363636363635*Power(y,4) + 1.4545454545454546*Power(y,6) - 
            2.727272727272727*Power(y,2)*Power(z,2) + 
            2.909090909090909*Power(y,3)*Power(z,3) - 1.3636363636363635*Power(z,4) + 
            1.4545454545454546*Power(z,6) + 
            Power(x,2)*(-2.727272727272727*Power(y,2) - 2.727272727272727*Power(z,2)) + 
            Power(x,3)*(2.909090909090909*Power(y,3) + 2.909090909090909*Power(z,3))) + 
         Power(M,8)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)*
          (-6.454545454545454*Power(x,4) + 7.2727272727272725*Power(x,6) - 
            6.454545454545454*Power(y,4) + 7.2727272727272725*Power(y,6) - 
            12.909090909090908*Power(y,2)*Power(z,2) + 
            14.545454545454545*Power(y,3)*Power(z,3) - 6.454545454545454*Power(z,4) + 
            7.2727272727272725*Power(z,6) + 
            Power(x,2)*(-12.909090909090908*Power(y,2) - 
               12.909090909090908*Power(z,2)) + 
            Power(x,3)*(14.545454545454545*Power(y,3) + 14.545454545454545*Power(z,3)))\
          + Power(M,4)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.5)*
          (1.4545454545454546*Power(x,6) + 1.4545454545454546*Power(y,6) + 
            10.*Power(z,4) + 1.4545454545454546*Power(z,6) + 
            Power(y,4)*(10. + 1.4545454545454546*Power(z,2)) + 
            Power(x,4)*(10. + 1.4545454545454546*Power(y,2) + 
               1.4545454545454546*Power(z,2)) + 
            Power(y,2)*(20.*Power(z,2) + 1.4545454545454546*Power(z,4)) + 
            Power(x,2)*(20.*Power(y,2) + 1.4545454545454546*Power(y,4) + 
               20.*Power(z,2) + 1.4545454545454546*Power(z,4))) + 
         Power(M,5)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
          (7.2727272727272725*Power(x,6) + 7.2727272727272725*Power(y,6) + 
            4.181818181818182*Power(z,4) + 7.2727272727272725*Power(z,6) + 
            Power(y,4)*(4.181818181818182 + 7.2727272727272725*Power(z,2)) + 
            Power(x,4)*(4.181818181818182 + 7.2727272727272725*Power(y,2) + 
               7.2727272727272725*Power(z,2)) + 
            Power(y,2)*(8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4)) + 
            Power(x,2)*(8.363636363636363*Power(y,2) + 7.2727272727272725*Power(y,4) + 
               8.363636363636363*Power(z,2) + 7.2727272727272725*Power(z,4))) + 
         Power(M,6)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            2.90909090909091*Power(y,3)*Power(z,3) - 6.909090909090909*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-6.909090909090909 + 16.*Power(z,2)) + 
            Power(x,4)*(-6.909090909090909 + 16.*Power(y,2) + 16.*Power(z,2)) + 
            Power(x,3)*(-2.90909090909091*Power(y,3) - 2.90909090909091*Power(z,3)) + 
            Power(y,2)*(-13.818181818181818*Power(z,2) + 16.*Power(z,4)) + 
            Power(x,2)*(-13.818181818181818*Power(y,2) + 16.*Power(y,4) - 
               13.818181818181818*Power(z,2) + 16.*Power(z,4))) + 
         Power(M,7)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.)*
          (14.545454545454545*Power(x,6) + 14.545454545454545*Power(y,6) - 
            14.545454545454543*Power(y,3)*Power(z,3) - 11.272727272727273*Power(z,4) + 
            14.545454545454545*Power(z,6) + 
            Power(y,4)*(-11.272727272727273 + 21.818181818181817*Power(z,2)) + 
            Power(x,4)*(-11.272727272727273 + 21.818181818181817*Power(y,2) + 
               21.818181818181817*Power(z,2)) + 
            Power(x,3)*(-14.545454545454543*Power(y,3) - 
               14.545454545454543*Power(z,3)) + 
            Power(y,2)*(-22.545454545454547*Power(z,2) + 
               21.818181818181817*Power(z,4)) + 
            Power(x,2)*(-22.545454545454547*Power(y,2) + 21.818181818181817*Power(y,4) - 
               22.545454545454547*Power(z,2) + 21.818181818181817*Power(z,4))))/
       (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
         Power(1.*M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),6)*
         (1.3636363636363635*Power(M,3) + 
           1.*Power(M,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5) + 
           0.45454545454545453*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
           0.09090909090909091*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5))),0.5));
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
    *value = (Power(M,2)*x*(-4.*M*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
       12.*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
       8.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
       8.*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3)); 
}

static void
func_Phi102(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4.*Power(M,2)*x*Power(y,2)*(1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       3.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3)); 
}

static void
func_Phi103(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4.*Power(M,2)*x*Power(z,2)*(1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
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
    *value = (-4.*Power(M,2)*Power(x,2)*y*(1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       3.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3));
}

static void
func_Phi202(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*y*(-4.*M*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
       12.*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
       8.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
       8.*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3));
}

static void
func_Phi203(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4.*Power(M,2)*y*Power(z,2)*(1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
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
    *value = (-4.*Power(M,2)*Power(x,2)*z*(1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       3.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3));
}

static void
func_Phi302(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-4.*Power(M,2)*Power(y,2)*z*(1.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),1.) + 
       3.*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*
     Power(M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),3));
}

static void
func_Phi303(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (Power(M,2)*z*(-4.*M*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) - 
       12.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),2.) + 
       8.*M*Power(Power(x,2) + Power(y,2) + Power(z,2),2.5) + 
       8.*Power(Power(x,2) + Power(y,2) + Power(z,2),3.)))/
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

