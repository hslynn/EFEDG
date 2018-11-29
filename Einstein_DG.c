#include "phg.h"

#include <string.h>
#include <math.h>


#include "hdw.c"
#include "auxi_dofs.c"
#include "Hhat.c"

static void 
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) 
{   
    *value = x*x + y;
}

//static void
//func_zero(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
//{
//    *value = 0.;
//}
//
//static void
//func_one(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
//{
//    *value = 1.;
//}

int 
main(int argc, char * argv[])
{
    char *meshfile ="icosahedron.mesh";
    //const char *matrix_fn = "mat.m", *var_name = "mat_M";
    GRID *g; 
    DOF_TYPE *dof_tp = DOF_DG2;
    FLOAT t0 = 0.0, t1 = 0.0;
    //MAT *mat_M;

    phgInit(&argc, &argv);	


    g = phgNewGrid(-1); 
    phgImport(g, meshfile, FALSE);
    phgBalanceGrid(g, 1.2, 1, NULL, 0.);

    /*creat dofs for all the functions to be solved*/ 
    DOF *dofs_Psi[10], *dofs_Pi[10], *dofs_Phi[30];
    create_dofs(g, dof_tp, 1, dofs_Psi, "Psi", 10);
    create_dofs(g, dof_tp, 1, dofs_Pi, "Pi", 10);
    create_dofs(g, dof_tp, 1, dofs_Phi, "Phi", 30);

    /*create dofs for all RHS terms*/
    DOF *dofs_rhsPsi[10], *dofs_rhsPi[10], *dofs_rhsPhi[30];
    create_dofs(g, dof_tp, 1, dofs_rhsPsi, "rhsPsi", 10);
    create_dofs(g, dof_tp, 1, dofs_rhsPi, "rhsPi", 10);
    create_dofs(g, dof_tp, 1, dofs_rhsPhi, "rhsPhi", 30);
    
    /*create dofs for derivatives of vars*/ 
    DOF *dofs_gradPsi[30], *dofs_gradPi[30], *dofs_gradPhi[90];
    create_dofs(g, dof_tp, 2, dofs_gradPsi, "gradPsi", 30);
    create_dofs(g, dof_tp, 2, dofs_gradPi, "gradPi", 30);
    create_dofs(g, dof_tp, 2, dofs_gradPhi, "gradPhi", 90);

    //create lists to store dofs of var, rhs and grad
    DOF *dofs_var[50], *dofs_rhs[50];
    for(int i =0;i<10;i++){
        dofs_var[i] = dofs_Psi[i];
        dofs_var[10+i] = dofs_Pi[i];
        dofs_var[20+i] = dofs_Phi[i]; 
        dofs_var[30+i] = dofs_Phi[10 + i];
        dofs_var[40+i] = dofs_Phi[20 + i];

        dofs_rhs[i] = dofs_rhsPsi[i];
        dofs_rhs[10+i] = dofs_rhsPi[i];
        dofs_rhs[20+i] = dofs_rhsPhi[i]; 
        dofs_rhs[30+i] = dofs_rhsPhi[10 + i];
        dofs_rhs[40+i] = dofs_rhsPhi[20 + i];
    }

    /*create auxi dofs*/
    DOF *dofs_g[6], *dofs_N[4], *dofs_Hhat[50];
    create_dofs(g, dof_tp, 1, dofs_g, "g", 6);
    create_dofs(g, dof_tp, 1, dofs_N, "N", 4);
    create_dofs(g, dof_tp, 1, dofs_Hhat, "Hhat", 50);  

    //for testing
    for(int i=0;i<10;i++){
        phgDofSetDataByValue(dofs_Pi[i], 0);
        phgDofSetDataByValue(dofs_Phi[i], 0);
        phgDofSetDataByValue(dofs_Phi[10 + i], 0);
        phgDofSetDataByValue(dofs_Phi[20 + i], 0);
    } 

    //phgDofSetDataByFunction(dofs_Pi[1], func_u);
    phgDofSetDataByValue(dofs_Psi[0], -1);
    phgDofSetDataByValue(dofs_Psi[4], 1);
    phgDofSetDataByValue(dofs_Psi[7], 1);
    phgDofSetDataByValue(dofs_Psi[9], 1);

    
    
    if(phgRank == 0){    
        t0 = phgGetTime(NULL);
    }   

    get_dofs_auxi(dofs_var, dofs_g, dofs_N, dofs_rhs); 
    get_dofs_Hhat(dofs_var, dofs_g, dofs_N, dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat);
    
    if(phgRank == 0){ 
        t1 = phgGetTime(NULL);
        phgPrintf("Total time cost = %f\n\n", t1 - t0);
    }

    /*Export VTK*/
    //phgExportVTKn(g, "gradPsi.vtk", 10, dofs_gradPsi);
    phgDofDump(dofs_gradPsi[20]);
    phgDofDump(dofs_gradPi[1]);
    phgDofDump(dofs_gradPi[11]);
    phgDofDump(dofs_gradPi[21]);

    //phgExportVTKn(g, "gradPhi.vtk", 30, dofs_gradPhi);
    
    /*release the mem*/
    for(int i=0;i<50;i++){
        phgDofFree(dofs_var + i);
        phgDofFree(dofs_rhs + i);
        phgDofFree(dofs_Hhat + i);
        }
    for(int i=0;i<30;i++){
        phgDofFree(dofs_gradPsi + i);
        phgDofFree(dofs_gradPi + i);
        phgDofFree(dofs_gradPhi + i);
        phgDofFree(dofs_gradPhi + 30 + i);
        phgDofFree(dofs_gradPhi + 60 + i);
        }
    for(int i=0; i<6; i++){
        phgDofFree(dofs_g + i);
    }
    for(int i=0; i<4; i++){
        phgDofFree(dofs_N + i);
    }
    phgFreeGrid(&g); 
 
    phgFinalize(); 

    return 0;
}

