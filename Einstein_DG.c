#include "phg.h"

#include <string.h>
#include <math.h>

#include "auxi_dofs.c"
#include "dofs.c"
#include "Hhat_dofs.c"

//static void 
//func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) 
//{   
//    *value = x*x + y;
//}

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
    //char *names_Psi[10] = {"Psi_00", "Psi_01", "Psi_02", "Psi_03", 
    //                       "Psi_11", "Psi_12", "Psi_13",
    //                       "Psi_22", "Psi_23",
    //                       "Psi_33"};

    //char *names_Pi[10] = {"Pi_00", "Pi_01", "Pi_02", "Pi_03", 
    //                      "Pi_11", "Pi_12", "Pi_13",
    //                      "Pi_22", "Pi_23",
    //                      "Pi_33"};
    //
    //char *names_Phi[30] = {"Phi_100","Phi_101","Phi_102","Phi_103",
    //                       "Phi_111","Phi_112","Phi_113",
    //                       "Phi_122","Phi_123",
    //                       "Phi_133",
    //                       "Phi_200","Phi_201","Phi_202","Phi_203",
    //                       "Phi_211","Phi_212","Phi_213",
    //                       "Phi_222","Phi_223",
    //                       "Phi_233",
    //                       "Phi_300","Phi_301","Phi_302","Phi_303",
    //                       "Phi_311","Phi_312","Phi_313",
    //                       "Phi_322","Phi_323",
    //                       "Phi_333"};

    create_dofs(g, dof_tp, 1, dofs_Psi, "Psi", 10);
    create_dofs(g, dof_tp, 1, dofs_Pi, "Pi", 10);
    create_dofs(g, dof_tp, 1, dofs_Phi, "Phi", 30);

    /*create dofs for all RHS terms*/
    DOF *dofs_rhsPsi[10], *dofs_rhsPi[10], *dofs_rhsPhi[30];
    //char *names_rhsPsi[10] = {"rhsPsi_00", "rhsPsi_01", "rhsPsi_02", "rhsPsi_03", 
    //                          "rhsPsi_11", "rhsPsi_12", "rhsPsi_13",
    //                          "rhsPsi_22", "rhsPsi_23",
    //                          "rhsPsi_33"};

    //char *names_rhsPi[10] = {"rhsPi_00", "rhsPi_01", "rhsPi_02", "rhsPi_03", 
    //                         "rhsPi_11", "rhsPi_12", "rhsPi_13",
    //                         "rhsPi_22", "rhsPi_23",
    //                         "rhsPi_33"};
    //
    //char *names_rhsPhi[30] = {"rhsPhi_100","rhsPhi_101","rhsPhi_102","rhsPhi_103",
    //                          "rhsPhi_111","rhsPhi_112","rhsPhi_113",
    //                          "rhsPhi_122","rhsPhi_123",
    //                          "rhsPhi_133",
    //                          "rhsPhi_200","rhsPhi_201","rhsPhi_202","rhsPhi_203",
    //                          "rhsPhi_211","rhsPhi_212","rhsPhi_213",
    //                          "rhsPhi_222","rhsPhi_223",
    //                          "rhsPhi_233",
    //                          "rhsPhi_300","rhsPhi_301","rhsPhi_302","rhsPhi_303",
    //                          "rhsPhi_311","rhsPhi_312","rhsPhi_313",
    //                          "rhsPhi_322","rhsPhi_323",
    //                          "rhsPhi_333"};

    create_dofs(g, dof_tp, 1, dofs_rhsPsi, "rhsPsi", 10);
    create_dofs(g, dof_tp, 1, dofs_rhsPi, "rhsPi", 10);
    create_dofs(g, dof_tp, 1, dofs_rhsPhi, "rhsPhi", 30);
    
    /*create dofs for derivatives of vars*/ 
    DOF *dofs_gradPsi[10], *dofs_gradPi[10], *dofs_gradPhi[30];
    //char *names_gradPsi[10] = {"gradPsi_00", "gradPsi_01", "gradPsi_02", "gradPsi_03", 
    //                          "gradPsi_11", "gradPsi_12", "gradPsi_13",
    //                          "gradPsi_22", "gradPsi_23",
    //                          "gradPsi_33"};

    //char *names_gradPi[10] = {"gradPi_00", "gradPi_01", "gradPi_02", "gradPi_03", 
    //                         "gradPi_11", "gradPi_12", "gradPi_13",
    //                         "gradPi_22", "gradPi_23",
    //                         "gradPi_33"};
    //
    //char *names_gradPhi[30] = {"gradPhi_100","gradPhi_101","gradPhi_102","gradPhi_103",
    //                          "gradPhi_111","gradPhi_112","gradPhi_113",
    //                          "gradPhi_122","gradPhi_123",
    //                          "gradPhi_133",
    //                          "gradPhi_200","gradPhi_201","gradPhi_202","gradPhi_203",
    //                          "gradPhi_211","gradPhi_212","gradPhi_213",
    //                          "gradPhi_222","gradPhi_223",
    //                          "gradPhi_233",
    //                          "gradPhi_300","gradPhi_301","gradPhi_302","gradPhi_303",
    //                          "gradPhi_311","gradPhi_312","gradPhi_313",
    //                          "gradPhi_322","gradPhi_323",
    //                          "gradPhi_333"};

   
    create_dofs(g, dof_tp, 6, dofs_gradPsi, "gradPsi", 10);
    create_dofs(g, dof_tp, 6, dofs_gradPi, "gradPi", 10);
    create_dofs(g, dof_tp, 6, dofs_gradPhi, "gradPhi", 30);


    //create lists to store dofs of var, rhs and grad
    DOF *dofs_var[50], *dofs_rhs[50], *dofs_grad_var[50];
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

        dofs_grad_var[i] = dofs_gradPsi[i];
        dofs_grad_var[10+i] = dofs_gradPi[i];
        dofs_grad_var[20+i] = dofs_gradPhi[i]; 
        dofs_grad_var[30+i] = dofs_gradPhi[10 + i];
        dofs_grad_var[40+i] = dofs_gradPhi[20 + i];
    }

    /*create auxi dofs*/
    DOF *dofs_g[6], *dofs_N[4];
    //char *names_g[6] = {"g11","g12","g13","g22","g23","g33"};
    //char *names_N[4] = {"N","N1","N2","N3"};
    create_dofs(g, dof_tp, 1, dofs_g, "g_", 6);
    create_dofs(g, dof_tp, 1, dofs_N, "N_", 4);
    
    //for testing
    for(int i=0;i<10;i++){
        phgDofSetDataByValue(dofs_Phi[i], 0);
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

    get_dofs_auxi(dofs_var, dofs_g, dofs_N, dofs_rhs); 
    get_dofs_grad(dofs_var, dofs_grad_var, 50);
    DOF *dofs_Hhat[50];
    create_dofs(g, dof_tp, 1, dofs_Hhat, "Hhat", 50); 
    get_dofs_Hhat(dofs_grad_var, dofs_g, dofs_N, dofs_Hhat);

    if(phgRank == 0){ 
        t1 = phgGetTime(NULL);
        phgPrintf("Total time cost = %f\n\n", t1 - t0);
    }

    /*Export VTK*/
    //phgExportVTKn(g, "gradPsi.vtk", 10, dofs_gradPsi);
    phgDofDump(dofs_Hhat[0]);
    //phgExportVTKn(g, "gradPhi.vtk", 30, dofs_gradPhi);
    
    /*release the mem*/
    for(int i=0;i<50;i++){
        phgDofFree(dofs_var + i);
        phgDofFree(dofs_rhs + i);
        phgDofFree(dofs_grad_var + i);
        phgDofFree(dofs_Hhat + i);
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

