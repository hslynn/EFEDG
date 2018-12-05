#include "phg.h"
#include <string.h>
#include <math.h>

#define NVAR 50

#include "hdw.c"
#include "auxi_dofs.c"
#include "Hhat.c"
#include "rk2.c"
#include "Schwarzschild_Harmonic.c"

int 
main(int argc, char * argv[])
{
    char *meshfile ="./mesh/hollowed_icosahedron.mesh";
    GRID *g; 
    DOF_TYPE *dof_tp = DOF_DG2;
    FLOAT t0 = 0.0, t1 = 0.0;
    FLOAT dt = 0.01, max_time = 1.0;
    INT max_step = ceil(max_time/dt);
    INT i;

    phgInit(&argc, &argv);	

    g = phgNewGrid(-1); 
    phgImport(g, meshfile, FALSE);
    phgBalanceGrid(g, 1.2, 1, NULL, 0.);

    /*creat dofs for all the functions to be solved*/ 
    DOF *dofs_Psi[10], *dofs_Pi[10], *dofs_Phi[30];
    DOF *dofs_ini_var[50], *dofs_bdry[50];
    create_dofs(g, dof_tp, 1, dofs_Psi, "Psi", 10);
    create_dofs(g, dof_tp, 1, dofs_Pi, "Pi", 10);
    create_dofs(g, dof_tp, 1, dofs_Phi, "Phi", 30);
    create_dofs(g, dof_tp, 1, dofs_ini_var, "ini_var", 50);

    /*create dofs for all source terms*/
    DOF *dofs_srcPsi[10], *dofs_srcPi[10], *dofs_srcPhi[30];
    create_dofs(g, dof_tp, 1, dofs_srcPsi, "srcPsi", 10);
    create_dofs(g, dof_tp, 1, dofs_srcPi, "srcPi", 10);
    create_dofs(g, dof_tp, 1, dofs_srcPhi, "srcPhi", 30);
    
    /*create dofs for derivatives of vars*/ 
    DOF *dofs_gradPsi[30], *dofs_gradPi[30], *dofs_gradPhi[90];
    create_dofs(g, dof_tp, 2, dofs_gradPsi, "gradPsi", 30);
    create_dofs(g, dof_tp, 2, dofs_gradPi, "gradPi", 30);
    create_dofs(g, dof_tp, 2, dofs_gradPhi, "gradPhi", 90);

    //create lists to store dofs of var, src
    DOF *dofs_var[NVAR], *dofs_src[NVAR];
    for(i =0;i<10;i++){
        dofs_var[i] = dofs_Psi[i];
        dofs_var[10+i] = dofs_Pi[i];
        dofs_var[20+i] = dofs_Phi[i]; 
        dofs_var[30+i] = dofs_Phi[10 + i];
        dofs_var[40+i] = dofs_Phi[20 + i];

        dofs_src[i] = dofs_srcPsi[i];
        dofs_src[10+i] = dofs_srcPi[i];
        dofs_src[20+i] = dofs_srcPhi[i]; 
        dofs_src[30+i] = dofs_srcPhi[10 + i];
        dofs_src[40+i] = dofs_srcPhi[20 + i];
    }

    /*create auxi dofs*/
    DOF *dofs_g[6], *dofs_N[4], *dofs_Hhat[NVAR], *dofs_rhs[NVAR];
    create_dofs(g, dof_tp, 1, dofs_g, "g", 6);
    create_dofs(g, dof_tp, 1, dofs_N, "N", 4);
    create_dofs(g, dof_tp, 1, dofs_Hhat, "Hhat", NVAR);  

    //for testing
    for(i=0;i<10;i++){
        phgDofSetDataByValue(dofs_Pi[i], 0);
        phgDofSetDataByValue(dofs_Phi[i], 0);
        phgDofSetDataByValue(dofs_Phi[10 + i], 0);
        phgDofSetDataByValue(dofs_Phi[20 + i], 0);
    } 

    //phgDofSetDataByFunction(dofs_Pi[1], func_u);
    set_data_dofs(dofs_var);
    copy_dofs(dofs_var, dofs_ini_var, "ini_var", NVAR);
    copy_dofs(dofs_var, dofs_bdry, "bdry", NVAR);
    
    
    if(phgRank == 0){    
        t0 = phgGetTime(NULL);
    }   

    get_dofs_rhs(dofs_var, dofs_bdry, dofs_g, dofs_N, dofs_src,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);

    char vtk_name[20]; 
    for(i=0;i<max_step;i++){
        sprintf(vtk_name, "%f", i*dt);
        strcat(vtk_name, ".vtk");
        printf("step %d completed\n", i);
        printf("time length: %f\n\n", i*dt);
        phgExportVTKn(g, vtk_name, 10, dofs_Psi);
        ssp_rk2(dt, dofs_var, dofs_bdry, dofs_g, dofs_N, dofs_src,
                dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
    }
      
    if(phgRank == 0){ 
        t1 = phgGetTime(NULL);
        phgPrintf("Total time cost = %f\n\n", t1 - t0);
    }

     
    /*release the mem*/
    for(i=0;i<NVAR;i++){
        phgDofFree(dofs_var + i);
        phgDofFree(dofs_src + i);
        phgDofFree(dofs_Hhat + i);
        }
    for(i=0;i<30;i++){
        phgDofFree(dofs_gradPsi + i);
        phgDofFree(dofs_gradPi + i);
        phgDofFree(dofs_gradPhi + i);
        phgDofFree(dofs_gradPhi + 30 + i);
        phgDofFree(dofs_gradPhi + 60 + i);
        }
    for(i=0; i<6; i++){
        phgDofFree(dofs_g + i);
    }
    for(i=0; i<4; i++){
        phgDofFree(dofs_N + i);
    }
    phgFreeGrid(&g); 
 
    phgFinalize(); 

    return 0;
}

