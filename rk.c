#include "global_def.h"
#include "rhs.h"
#include "hdw.h"

void
rk2(FLOAT dt, DOF **dofs_var, DOF **dofs_exact, DOF **dofs_g, DOF **dofs_N, DOF **dofs_H, DOF **dofs_deriH, DOF **dofs_src, DOF **dofs_C,
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF **dofs_Hhat, DOF **dofs_rhs)
{
    INT i; 
    DOF *dofs_var_temp[NVAR];
    get_dofs_rhs(dofs_var, dofs_var, dofs_exact, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src, dofs_C,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
    //phgExportVTKn(dofs_src[0]->g, "src.vtk", NVAR, dofs_src);
    //phgExportVTKn(dofs_src[0]->g, "Hhat.vtk", NVAR, dofs_Hhat);
    //phgExportVTKn(dofs_src[0]->g, "rhs.vtk", NVAR, dofs_rhs);
    create_dofs(dofs_var[0]->g, dofs_var[0]->type, 1, dofs_var_temp, "temp_var", NVAR); 
    copy_dofs(dofs_var, dofs_var_temp, "temp_var", NVAR);
    for(i=0;i<NVAR;i++){
        phgDofAXPBY(dt, dofs_rhs[i], 1.0, dofs_var_temp + i);
    }
    get_dofs_rhs(dofs_var_temp, dofs_var_temp, dofs_exact, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src, dofs_C,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs); 
    for(i=0;i<NVAR;i++){
        phgDofAXPBY(dt, dofs_rhs[i], 1.0, dofs_var_temp + i); 
        phgDofAXPBY(0.5, dofs_var_temp[i], 0.5, dofs_var + i);
    }
    
    free_dofs(dofs_var_temp, NVAR);
}

void
rk3(FLOAT dt, DOF **dofs_var, DOF **dofs_exact, DOF **dofs_g, DOF **dofs_N, DOF **dofs_H, DOF **dofs_deriH, DOF **dofs_src, DOF **dofs_C,
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF **dofs_Hhat, DOF **dofs_rhs)
{
    INT i; 
    DOF *dofs_var_temp[NVAR];
    create_dofs(dofs_var[0]->g, dofs_var[0]->type, 1, dofs_var_temp, "temp_var", NVAR); 
    copy_dofs(dofs_var, dofs_var_temp, "temp_var", NVAR);
    
    /*u_1 = u_n + dt*rhs(u_n)*/
    get_dofs_rhs(dofs_var, dofs_exact, dofs_exact, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src, dofs_C,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
    //phgExportVTKn(dofs_src[0]->g, "src.vtk", NVAR, dofs_src);
    //phgExportVTKn(dofs_src[0]->g, "Hhat.vtk", NVAR, dofs_Hhat);
    //phgExportVTKn(dofs_src[0]->g, "rhs.vtk", NVAR, dofs_rhs);
    for(i=0;i<NVAR;i++){
        phgDofAXPBY(dt, dofs_rhs[i], 1.0, dofs_var_temp + i);
    }

    /*u_2 = 3/4*u_n + 1/4*u_1 + 1/4*dt*rhs(u_1)*/
    get_dofs_rhs(dofs_var_temp, dofs_exact, dofs_exact, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src, dofs_C,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs); 
    for(i=0;i<NVAR;i++){
        phgDofAXPBY(dt, dofs_rhs[i], 1.0, dofs_var_temp + i); 
        phgDofAXPBY(3./4., dofs_var[i], 1./4., dofs_var_temp + i);
    }

    /*u_{n+1} = 1/3*u_n + 2/3*u_2 + 2/3*dt*rhs(u_2)*/
    get_dofs_rhs(dofs_var_temp, dofs_exact, dofs_exact, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src, dofs_C,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
    for(i=0;i<NVAR;i++){
        phgDofAXPBY(dt, dofs_rhs[i], 1.0, dofs_var_temp + i); 
        phgDofAXPBY(2./3., dofs_var_temp[i], 1./3., dofs_var + i);
    }
    
    free_dofs(dofs_var_temp, NVAR);
}

void
rk1(FLOAT dt, DOF **dofs_var, DOF **dofs_exact, DOF **dofs_g, DOF **dofs_N, DOF **dofs_H, DOF **dofs_deriH, DOF **dofs_src, DOF **dofs_C,
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF **dofs_Hhat, DOF **dofs_rhs)
{
    INT i; 
    /*u_n = u_n + dt*rhs(u_n)*/
    get_dofs_rhs(dofs_var, dofs_exact, dofs_exact, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src, dofs_C,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
    for(i=0;i<NVAR;i++){
        phgDofAXPBY(dt, dofs_rhs[i], 1.0, dofs_var + i);
    }
}
