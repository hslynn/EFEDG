static void
rk2(FLOAT dt, DOF **dofs_var, DOF **dofs_bdry, DOF **dofs_g, DOF **dofs_N, DOF **dofs_H, DOF **dofs_deriH, DOF **dofs_src,
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF **dofs_Hhat, DOF **dofs_rhs)
{
    INT i; 
    DOF *dofs_var_temp[NVAR];
    get_dofs_rhs(dofs_var, dofs_bdry, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
    create_dofs(dofs_var[0]->g, dofs_var[0]->type, 1, dofs_var_temp, "temp_var", NVAR); 
    copy_dofs(dofs_var, dofs_var_temp, "temp_var", NVAR);
    for(i=0;i<NVAR;i++){
        phgDofAXPBY(dt, dofs_rhs[i], 1.0, dofs_var_temp + i);
    }
    get_dofs_rhs(dofs_var_temp, dofs_bdry, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs); 
    for(i=0;i<NVAR;i++){
        phgDofAXPBY(dt, dofs_rhs[i], 1.0, dofs_var_temp + i); 
        phgDofAXPBY(0.5, dofs_var_temp[i], 0.5, dofs_var + i);
    }
    
    free_dofs(dofs_var_temp, NVAR);
}
