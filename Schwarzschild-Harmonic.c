#define M 0.5 ,
#define R pow(x*x + y*y + z*z, 0.5)

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
    *value = pow(1.0+M/R, 2.0) + (1.0+M/R)/(1-M/R)*pow(M, 2.0)/pow(R, 4.0)*pow(x, 2.0);
}

static void
func_Psi12(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (1.0+M/R)/(1-M/R)*pow(M, 2.0)/pow(R, 4.0)*x*y;
}

static void
func_Psi13(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (1.0+M/R)/(1-M/R)*pow(M, 2.0)/pow(R, 4.0)*x*z;
}


static void
func_Psi22(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = pow(1.0+M/R, 2.0) + (1.0+M/R)/(1-M/R)*pow(M, 2.0)/pow(R, 4.0)*pow(y, 2.0);
}

static void
func_Psi23(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (1.0+M/R)/(1-M/R)*pow(M, 2.0)/pow(R, 4.0)*y*z;
}

static void
func_Psi33(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = pow(1.0+M/R, 2.0) + (1.0+M/R)/(1-M/R)*pow(M, 2.0)/pow(R, 4.0)*pow(z, 2.0);
}

static void
init_bdry_funcs(DOF_USER_FUNC *func_bdrys_Psi, DOF_USER_FUNC *func_bdrys_Pi, DOF_USER_FUNC *func_bdrys_Phi)
{
    int i;
    func_bdrys_Psi[10] = {func_Psi00, func_zero, func_zero, func_zero,
                                                func_Psi11, func_Psi12, func_Psi13,
                                                           func_Psi22, func_Psi23,
                                                                        func_Psi33};
    for(i=0; i<10; i++){
        *(func_bdrys_Pi++) = func_zero;
    }

    //phi
}

