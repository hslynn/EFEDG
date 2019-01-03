static void
get_dofs_rhs(DOF **dofs_var, DOF **dofs_bdry, DOF **dofs_g, DOF **dofs_N, DOF **dofs_src,
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF **dofs_Hhat, DOF **dofs_rhs)
{
    short i;
    get_dofs_auxi(dofs_var, dofs_g, dofs_N, dofs_src); 
    //phgExportVTKn(dofs_src[0]->g, "src.vtk", NVAR, dofs_src);
    //phgExportVTKn(dofs_g[0]->g, "g.vtk", 6, dofs_g);
    //phgExportVTKn(dofs_N[0]->g, "N.vtk", 4, dofs_N);

    get_dofs_Hhat(dofs_var, dofs_bdry, dofs_g, dofs_N, 
            dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat);

    //phgExportVTKn(dofs_gradPsi[0]->g, "gradPsi.vtk", 30, dofs_gradPsi);
    //phgExportVTKn(dofs_gradPsi[0]->g, "gradPi.vtk", 30, dofs_gradPi);
    //phgExportVTKn(dofs_gradPsi[0]->g, "gradxPhi.vtk", 30, dofs_gradPhi);
    //phgExportVTKn(dofs_gradPsi[0]->g, "gradyPhi.vtk", 30, dofs_gradPhi + 30);
    //phgExportVTKn(dofs_gradPsi[0]->g, "gradzPhi.vtk", 30, dofs_gradPhi + 60);

    for(i=0;i<NVAR;i++){
        dofs_rhs[i] = phgDofAXPBY(-1.0, dofs_Hhat[i], 1.0, dofs_src + i);
    }
}
static void
ssp_rk2(FLOAT dt, DOF **dofs_var, DOF **dofs_bdry, DOF **dofs_g, DOF **dofs_N, DOF **dofs_src,
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF **dofs_Hhat, DOF **dofs_rhs)
{
    INT i; 
    DOF *dofs_var_temp[NVAR];
    get_dofs_rhs(dofs_var, dofs_bdry, dofs_g, dofs_N, dofs_src,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
    create_dofs(dofs_var[0]->g, dofs_var[0]->type, 1, dofs_var_temp, "temp_var", NVAR); 
    copy_dofs(dofs_var, dofs_var_temp, "temp_var", NVAR);
    for(i=0;i<NVAR;i++){
        phgDofAXPBY(dt, dofs_rhs[i], 1.0, dofs_var_temp + i);
    }
    get_dofs_rhs(dofs_var_temp, dofs_bdry, dofs_g, dofs_N, dofs_src,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs); 
    for(i=0;i<NVAR;i++){
       phgDofAXPBY(dt, dofs_rhs[i], 1.0, dofs_var_temp + i); 
       phgDofAXPBY(0.5, dofs_var_temp[i], 0.5, dofs_var + i);
    }
    
    free_dofs(dofs_var_temp, NVAR);
}
