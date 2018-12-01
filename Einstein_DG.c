#include "phg.h"
#include <string.h>
#include <math.h>

#define NVAR 50

#include "hdw.c"
#include "auxi_dofs.c"
#include "Hhat.c"
#include "rk2.c"
static void 
func_Psi00(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) 
{   
    *value = x*x + y;
}

int 
main(int argc, char * argv[])
{
    char *meshfile ="./mesh/hollowed_icosahedron.mesh";
    //const char *matrix_fn = "mat.m", *var_name = "mat_M";
    GRID *g; 
    DOF_TYPE *dof_tp = DOF_DG2;
    FLOAT t0 = 0.0, t1 = 0.0;
    FLOAT dt = 0.01;
    INT i;
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
    phgDofSetDataByValue(dofs_Psi[0], -1);
    phgDofSetDataByValue(dofs_Psi[4], 1);
    phgDofSetDataByValue(dofs_Psi[7], 1);
    phgDofSetDataByValue(dofs_Psi[9], 1);

    
    printf("test0");
    
    if(phgRank == 0){    
        t0 = phgGetTime(NULL);
    }   
    get_dofs_rhs(dofs_var, dofs_g, dofs_N, dofs_src,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
    printf("test1");
    ssp_rk2(dt, dofs_var, dofs_g, dofs_N, dofs_src,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
      
    printf("test2");
    if(phgRank == 0){ 
        t1 = phgGetTime(NULL);
        phgPrintf("Total time cost = %f\n\n", t1 - t0);
    }

    /*Export VTK*/
    //phgExportVTKn(g, "gradPsi.vtk", 10, dofs_gradPsi);
    for(i=0;i<10;i++){
        phgDofDump(dofs_Psi[i]);
    }

    phgExportVTK(g, "Psi.vtk",dofs_Psi[0], NULL);
    
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

