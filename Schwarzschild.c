#define M 0.5 
#define R pow(x*x + y*y + z*z, 0.5)

#define func_bdry_Psi00(x, y, z, value)\
        func_Psi00(x, y, z, value)
#define func_bdry_Psi11(x, y, z, value)\
        func_Psi00(x, y, z, value)
#define func_bdry_Psi12(x, y, z, value)\
        func_Psi00(x, y, z, value)
#define func_bdry_Psi13(x, y, z, value)\
        func_Psi00(x, y, z, value)
#define func_bdry_Psi22(x, y, z, value)\
        func_Psi00(x, y, z, value)
#define func_bdry_Psi23(x, y, z, value)\
        func_Psi00(x, y, z, value)
#define func_bdry_Psi33(x, y, z, value)\
        func_Psi00(x, y, z, value)

static void
func_bdry_zero(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
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
