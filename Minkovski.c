#define M 0.5
#define R (Pow(x*x + y*y + z*z, 0.5))

static void
func_zero(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0; 
}

static void
func_one(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 1; 
}

static void
func_minus_one(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -1; 
}

/*set dof data using above functions*/
static void
set_data_var(DOF **dofs_var)
{
    int i;
    DOF_USER_FUNC funcs_Psi[10] = {func_minus_one, func_zero, func_zero, func_zero,
                                                    func_one, func_zero, func_zero,
                                                                func_one, func_zero,
                                                                            func_one};
    DOF_USER_FUNC funcs_Pi[10];
    for(i=0; i<10; i++){
        funcs_Pi[i] = func_zero;
    }

    DOF_USER_FUNC funcs_Phi[30];
    for(i=0;i<30;i++){
        funcs_Phi[i] = func_zero;
    }
    for(i=0; i<10; i++){
        phgDofSetDataByFunction(dofs_var[i], funcs_Psi[i]);
        phgDofSetDataByFunction(dofs_var[10 + i], funcs_Pi[i]);
        
        phgDofSetDataByFunction(dofs_var[20 + i], funcs_Phi[i]);
        phgDofSetDataByFunction(dofs_var[30 + i], funcs_Phi[10 + i]);
        phgDofSetDataByFunction(dofs_var[40 + i], funcs_Phi[20 + i]);
    }
}

