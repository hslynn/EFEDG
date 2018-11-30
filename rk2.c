static void
get_dofs_rhs(DOF **dofs_var, DOF **dofs_g, DOF **dofs_N, DOF **dofs_src,
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF **dofs_Hhat, DOF **dofs_rhs)
{
    get_dofs_auxi(dofs_var, dofs_g, dofs_N, dofs_src); 
    get_dofs_Hhat(dofs_var, dofs_g, dofs_N, 
            dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat);
    for(int i=0;i<NVAR;i++){
        dofs_rhs[i] = phgDofAXPBY(-1.0, dofs_Hhat[i], 1.0, dofs_src + i);
    }
}
static void
ssp_rk2(FLOAT dt, DOF **dofs_var, DOF **dofs_g, DOF **dofs_N, DOF **dofs_src,
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF ** dofs_Hhat, DOF **dofs_rhs)
{
    INT i; 
    DOF *dofs_var_temp[NVAR];
    create_dofs(dofs_var[0]->g, dofs_var[0]->type, 1, dofs_var_temp, "temp_var", NVAR); 
    for(i=0;i<NVAR;i++){
        phgDofCopy(dofs_var[i], dofs_var_temp + i, NULL, "temp_var");
        phgDofAXPBY(dt, dofs_rhs[i], 1.0, dofs_var_temp + i);
    }
    
    get_dofs_rhs(dofs_var_temp, dofs_g, dofs_N, dofs_src,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs); 
    for(i=0;i<NVAR;i++){
       phgDofAXPBY(dt, dofs_rhs[i], 1.0, dofs_var_temp + i); 
       phgDofAXPBY(0.5, dofs_var_temp[i], 0.5, dofs_var + i);
       phgDofFree(dofs_var_temp + i);  
    }
}
