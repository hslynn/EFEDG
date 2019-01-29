static void
rk4(FLOAT dt, DOF **dofs_var, DOF **dofs_bdry, DOF **dofs_g, DOF **dofs_N, DOF **dofs_H, DOF **dofs_deriH, DOF **dofs_src,
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF **dofs_Hhat, DOF **dofs_rhs)
{
    INT i; 
    
    DOF *dofs_k[4][NVAR];
    DOF *dofs_var_temp[NVAR];
    for(i=0;i<4;i++){
        create_dofs(dofs_var[0]->g, dofs_var[0]->type, 1, dofs_k[i], "dofs_k", NVAR);
    }
    create_dofs(dofs_var[0]->g, dofs_var[0]->type, 1, dofs_var_temp, "temp_var", NVAR);
    copy_dofs(dofs_var, dofs_var_temp, "temp_var", NVAR);

    get_dofs_rhs(dofs_var_temp, dofs_bdry, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_k[0]);
    for(i=0;i<NVAR;i++){
        phgDofAXPBY(0.5*dt, dofs_k[0][i], 1.0, dofs_var_temp + i);
    }
    get_dofs_rhs(dofs_var_temp, dofs_bdry, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_k[1]);

    copy_dofs(dofs_var, dofs_var_temp, "temp_var", NVAR);
    for(i=0;i<NVAR;i++){
        phgDofAXPBY(0.5*dt, dofs_k[1][i], 1.0, dofs_var_temp + i);
    }
    get_dofs_rhs(dofs_var_temp, dofs_bdry, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_k[2]);

    copy_dofs(dofs_var, dofs_var_temp, "temp_var", NVAR);
    for(i=0;i<NVAR;i++){
        phgDofAXPBY(dt, dofs_k[2][i], 1.0, dofs_var_temp + i);
    }
    get_dofs_rhs(dofs_var_temp, dofs_bdry, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_k[3]);

    for(i=0;i<NVAR;i++){
        phgDofAXPBY(1.0/6.0 * dt, dofs_k[0][i], 1.0, dofs_var + i);
        phgDofAXPBY(1.0/3.0 * dt, dofs_k[1][i], 1.0, dofs_var + i);
        phgDofAXPBY(1.0/3.0 * dt, dofs_k[2][i], 1.0, dofs_var + i);
        phgDofAXPBY(1.0/6.0 * dt, dofs_k[3][i], 1.0, dofs_var + i);
    }

    for(i=0;i<4;i++){
        free_dofs(dofs_k[i], NVAR);
    }
    free_dofs(dofs_var_temp, NVAR);
}
