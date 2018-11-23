#include "phg.h"

#include <string.h>
#include <math.h>

#include "auxi_dofs.c"
#include "dofs.c"
//static void 
//func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) 
//{   
//    *value = (x+1.)*(x+1.);
//}
//
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
    char *names_Psi[10] = {"Psi_00", "Psi_01", "Psi_02", "Psi_03", 
                           "Psi_11", "Psi_12", "Psi_13",
                           "Psi_22", "Psi_23",
                           "Psi_33"};

    char *names_Pi[10] = {"Pi_00", "Pi_01", "Pi_02", "Pi_03", 
                          "Pi_11", "Pi_12", "Pi_13",
                          "Pi_22", "Pi_23",
                          "Pi_33"};
    
    char *names_Phi[30] = {"Phi_100","Phi_101","Phi_102","Phi_103",
                           "Phi_111","Phi_112","Phi_113",
                           "Phi_122","Phi_123",
                           "Phi_133",
                           "Phi_200","Phi_201","Phi_202","Phi_203",
                           "Phi_211","Phi_212","Phi_213",
                           "Phi_222","Phi_223",
                           "Phi_233",
                           "Phi_300","Phi_301","Phi_302","Phi_303",
                           "Phi_311","Phi_312","Phi_313",
                           "Phi_322","Phi_323",
                           "Phi_333"};

    create_dofs(g, dof_tp, dofs_Psi, names_Psi, 10);
    create_dofs(g, dof_tp, dofs_Pi, names_Pi, 10);
    create_dofs(g, dof_tp, dofs_Phi, names_Phi, 30);

    /*create dofs for all RHS terms*/
    DOF *dofs_rhsPsi[10], *dofs_rhsPi[10], *dofs_rhsPhi[30];
    char *names_rhsPsi[10] = {"rhsPsi_00", "rhsPsi_01", "rhsPsi_02", "rhsPsi_03", 
                              "rhsPsi_11", "rhsPsi_12", "rhsPsi_13",
                              "rhsPsi_22", "rhsPsi_23",
                              "rhsPsi_33"};

    char *names_rhsPi[10] = {"rhsPi_00", "rhsPi_01", "rhsPi_02", "rhsPi_03", 
                             "rhsPi_11", "rhsPi_12", "rhsPi_13",
                             "rhsPi_22", "rhsPi_23",
                             "rhsPi_33"};
    
    char *names_rhsPhi[30] = {"rhsPhi_100","rhsPhi_101","rhsPhi_102","rhsPhi_103",
                              "rhsPhi_111","rhsPhi_112","rhsPhi_113",
                              "rhsPhi_122","rhsPhi_123",
                              "rhsPhi_133",
                              "rhsPhi_200","rhsPhi_201","rhsPhi_202","rhsPhi_203",
                              "rhsPhi_211","rhsPhi_212","rhsPhi_213",
                              "rhsPhi_222","rhsPhi_223",
                              "rhsPhi_233",
                              "rhsPhi_300","rhsPhi_301","rhsPhi_302","rhsPhi_303",
                              "rhsPhi_311","rhsPhi_312","rhsPhi_313",
                              "rhsPhi_322","rhsPhi_323",
                              "rhsPhi_333"};

    create_dofs(g, dof_tp, dofs_rhsPsi, names_rhsPsi, 10);
    create_dofs(g, dof_tp, dofs_rhsPi, names_rhsPi, 10);
    create_dofs(g, dof_tp, dofs_rhsPhi, names_rhsPhi, 30);

    /*create auxi dofs*/
    DOF *dofs_g[6], *dofs_N[4];
    char *names_g[6] = {"g11","g12","g13","g22","g23","g33"};
    char *names_N[4] = {"N","N1","N2","N3"};
    create_dofs(g, dof_tp, dofs_g, names_g, 6);
    create_dofs(g, dof_tp, dofs_N, names_N, 4);
    
    //for testing
    for(int i=0;i<10;i++){
        phgDofSetDataByValue(dofs_Pi[i], 0);
        phgDofSetDataByValue(dofs_Phi[i], 0);
        phgDofSetDataByValue(dofs_Phi[10 + i], 0);
        phgDofSetDataByValue(dofs_Phi[20 + i], 0);
    } 
    phgDofSetDataByValue(dofs_Psi[0], -1);
    phgDofSetDataByValue(dofs_Psi[4], 1);
    phgDofSetDataByValue(dofs_Psi[7], 1);
    phgDofSetDataByValue(dofs_Psi[9], 1);
    
    if(phgRank == 0){    
        t0 = phgGetTime(NULL);
    }   

    get_dofs_g(dofs_Psi, dofs_g);
    get_dofs_N(dofs_Psi, dofs_g, dofs_N);
   // get_dofs_rhs(dofs_Psi, dofs_Pi, dofs_Phi, dofs_g, dofs_N, 
   //                 dofs_rhsPsi, dofs_rhsPi, dofs_rhsPhi);

    if(phgRank == 0){ 
        t1 = phgGetTime(NULL);
        phgPrintf("Total time cost = %f\n\n", t1 - t0);
    }

    /*Export VTK*/
    phgExportVTK(g, "vtk.vtk", dofs_g[1], NULL);
    /*release the mem*/
    for(int i=0;i<10;i++){
        phgDofFree(dofs_Psi + i);
        phgDofFree(dofs_Pi + i);
        phgDofFree(dofs_Phi + i);
        phgDofFree(dofs_Phi + 10 + i);
        phgDofFree(dofs_Phi + 20 + i);
        phgDofFree(dofs_rhsPsi + i);
        phgDofFree(dofs_rhsPi + i);
        phgDofFree(dofs_rhsPhi + i);
        phgDofFree(dofs_rhsPhi + 10 + i);
        phgDofFree(dofs_rhsPhi + 20 + i);
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

