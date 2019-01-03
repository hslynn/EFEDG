#define M 0.2
#define R (Pow(x*x + y*y + z*z, 0.5))

static void
func_zero(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0; 
}

static void
func_Psi00(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -(1.0-M/R)/(1.0+M/R);
}

static void
func_Psi11(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = Pow(1.0+M/R, 2.0) + (1.0+M/R)/(1-M/R)*Pow(M, 2.0)/Pow(R, 4.0)*Pow(x, 2.0);
}

static void
func_Psi12(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (1.0+M/R)/(1-M/R)*Pow(M, 2.0)/Pow(R, 4.0)*x*y;
}

static void
func_Psi13(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (1.0+M/R)/(1-M/R)*Pow(M, 2.0)/Pow(R, 4.0)*x*z;
}


static void
func_Psi22(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = Pow(1.0+M/R, 2.0) + (1.0+M/R)/(1-M/R)*Pow(M, 2.0)/Pow(R, 4.0)*Pow(y, 2.0);
}

static void
func_Psi23(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (1.0+M/R)/(1-M/R)*Pow(M, 2.0)/Pow(R, 4.0)*y*z;
}

static void
func_Psi33(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = Pow(1.0+M/R, 2.0) + (1.0+M/R)/(1-M/R)*Pow(M, 2.0)/Pow(R, 4.0)*Pow(z, 2.0);
}

//nonzero funcs of Phi
static void
func_Phi100(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2.*M*x)/(Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)*Power(M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi111(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2.*Power(M,3.)*Power(x,3.)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 2.*Power(M,3)*x*Power(Power(x,2) + Power(y,2)
     + Power(z,2),8.5)*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.) + 
     4.*Power(M,2)*x*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.) - 
     2.*M*x*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.) + 
     Power(M,4.)*(4.*Power(x,3.)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) - 2.*Power(x,1.)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.)) + 
     Power(M,2.)*(-4.*Power(x,3.)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 2.*Power(x,1.)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
} 

static void
func_Phi112(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
   *value = (y*(-2.*Power(M,3.)*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + Power(M,4.)*(4.*Power(x,2)*Power(Power(x,2) + Power(y,2) 
    + Power(z,2),5.5) - 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)) + 
       Power(M,2.)*(-4.*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)); 
}

static void
func_Phi113(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (z*(-2.*Power(M,3.)*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + Power(M,4.)*(4.*Power(x,2)*Power(Power(x,2) + Power(y,2) 
    + Power(z,2),5.5) - 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)) + 
       Power(M,2.)*(-4.*Power(x,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi122(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = x*((-2.*M*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
     (4.*Power(M,4.)*Power(y,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)) - 
     (6.*Power(M,3.)*Power(y,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)) - 
     (4.*Power(M,2.)*Power(y,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5))));
}

static void
func_Phi123(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (x*y*z*(4.*Power(M,4.)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.5) - 2.*Power(M,3.)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.)
     - 4.*Power(M,2.)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi133(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = x*((-2.*M*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
     (4.*Power(M,4.)*Power(z,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)) - 
     (6.*Power(M,3.)*Power(z,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)) - 
     (4.*Power(M,2.)*Power(z,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5))));
}

static void
func_Phi200(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2.*M*y)/(Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)*Power(M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi211(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = y*((-2.*M*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
     (4.*Power(M,4.)*Power(x,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)) - 
     (6.*Power(M,3.)*Power(x,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)) - 
     (4.*Power(M,2.)*Power(x,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)))); 
}

static void
func_Phi212(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (x*(-2.*Power(M,3.)*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + Power(M,4.)*(4.*Power(y,2)*Power(Power(x,2) + 
      Power(y,2) + Power(z,2),5.5) - 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)) + 
       Power(M,2.)*(-4.*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi213(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (x*y*z*(4.*Power(M,4.)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.5) - 2.*Power(M,3.)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.) -
     4.*Power(M,2.)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi222(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2.*Power(M,3.)*Power(y,3.)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 2.*Power(M,3)*y*Power(Power(x,2) + Power(y,2) + 
    Power(z,2),8.5)*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.) + 
     4.*Power(M,2)*y*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.) - 
     2.*M*y*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.) + 
     Power(M,4.)*(4.*Power(y,3.)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) - 2.*Power(y,1.)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.)) + 
     Power(M,2.)*(-4.*Power(y,3.)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 2.*Power(y,1.)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi223(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (z*(-2.*Power(M,3.)*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + Power(M,4.)*(4.*Power(y,2)*Power(Power(x,2) + 
      Power(y,2) + Power(z,2),5.5) - 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)) + 
       Power(M,2.)*(-4.*Power(y,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi233(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = y*((-2.*M*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
     (4.*Power(M,4.)*Power(z,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)) - 
     (6.*Power(M,3.)*Power(z,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)) - 
     (4.*Power(M,2.)*Power(z,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5))));
}

static void
func_Phi300(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2.*M*z)/(Power(Power(x,2) + Power(y,2) + Power(z,2),0.5)*Power(M + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi311(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = z*((-2.*M*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
     (4.*Power(M,4.)*Power(x,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)) - 
     (6.*Power(M,3.)*Power(x,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)) - 
     (4.*Power(M,2.)*Power(x,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5))));
}

static void
func_Phi312(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (x*y*z*(4.*Power(M,4.)*Power(Power(x,2) + Power(y,2) + Power(z,2),3.5) - 2.*Power(M,3.)*Power(Power(x,2) + Power(y,2) + 
    Power(z,2),4.) - 4.*Power(M,2.)*Power(Power(x,2) + Power(y,2) + Power(z,2),4.5)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi313(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (x*(-2.*Power(M,3.)*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + Power(M,4.)*(4.*Power(z,2)*Power(Power(x,2) +
     Power(y,2) + Power(z,2),5.5) - 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)) + 
       Power(M,2.)*(-4.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)); 
}

static void
func_Phi322(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = z*((-2.*M*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.))/Power(Power(x,2) + Power(y,2) + Power(z,2),1.5) + 
     (4.*Power(M,4.)*Power(y,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),3.)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)) - 
     (6.*Power(M,3.)*Power(y,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2)) - 
     (4.*Power(M,2.)*Power(y,2.))/(Power(Power(x,2) + Power(y,2) + Power(z,2),2.5)*(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5))));
}

static void
func_Phi323(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (y*(-2.*Power(M,3.)*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.) + Power(M,4.)*(4.*Power(z,2)*Power(Power(x,2) + 
    Power(y,2) + Power(z,2),5.5) - 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5)) + 
       Power(M,2.)*(-4.*Power(z,2)*Power(Power(x,2) + Power(y,2) + Power(z,2),6.5) + 1.*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5))))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),8.5)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

static void
func_Phi333(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (-2.*Power(M,3.)*Power(z,3.)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.5) - 2.*Power(M,3)*z*Power(Power(x,2) + Power(y,2) + 
    Power(z,2),8.5)*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.) + 
     4.*Power(M,2)*z*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.) - 
     2.*M*z*Power(Power(x,2) + Power(y,2) + Power(z,2),9.5)*Power(1. + M/Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),1.) + 
     Power(M,4.)*(4.*Power(z,3.)*Power(Power(x,2) + Power(y,2) + Power(z,2),7.) - 2.*Power(z,1.)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.)) + 
     Power(M,2.)*(-4.*Power(z,3.)*Power(Power(x,2) + Power(y,2) + Power(z,2),8.) + 2.*Power(z,1.)*Power(Power(x,2) + Power(y,2) + Power(z,2),9.)))/
   (Power(Power(x,2) + Power(y,2) + Power(z,2),10.)*Power(-1.*M + Power(Power(x,2) + Power(y,2) + Power(z,2),0.5),2));
}

/*set dof data using above functions*/
static void
set_data_dofs(DOF **dofs_var)
{
    int i;
    DOF_USER_FUNC funcs_Psi[10] = {func_Psi00, func_zero, func_zero, func_zero,
                                                    func_Psi11, func_Psi12, func_Psi13,
                                                                func_Psi22, func_Psi23,
                                                                            func_Psi33};
    DOF_USER_FUNC funcs_Pi[10];
    for(i=0; i<10; i++){
        funcs_Pi[i] = func_zero;
    }

    DOF_USER_FUNC funcs_Phi[30] = {func_Phi100, func_zero, func_zero, func_zero,        
                                               func_Phi111, func_Phi112, func_Phi113,
                                                          func_Phi122, func_Phi123,
                                                                      func_Phi133,
                                   func_Phi200, func_zero, func_zero, func_zero,        
                                               func_Phi211, func_Phi212, func_Phi213,
                                                          func_Phi222, func_Phi223,
                                                                      func_Phi233,
                                   func_Phi300, func_zero, func_zero, func_zero,        
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

