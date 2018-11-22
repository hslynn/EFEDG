#include "phg.h"

#include <string.h>
#include <math.h>

#include "rhs_dofs.c"
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
    FLOAT t0, t1;
    //MAT *mat_M;

    phgInit(&argc, &argv);	


    g = phgNewGrid(-1); 
    phgImport(g, meshfile, FALSE);
    phgBalanceGrid(g, 1.2, 1, NULL, 0.);

    //phgRefineAllElements(g, 1);
   
    /*creat dofs for all the functions to be solved*/ 
    DOF *Psi_00 = phgDofNew(g, dof_tp, 1, "Psi_00", DofInterpolation);
    DOF *Psi_01 = phgDofNew(g, dof_tp, 1, "Psi_01", DofInterpolation);
    DOF *Psi_02 = phgDofNew(g, dof_tp, 1, "Psi_02", DofInterpolation);
    DOF *Psi_03 = phgDofNew(g, dof_tp, 1, "Psi_03", DofInterpolation);
    DOF *Psi_11 = phgDofNew(g, dof_tp, 1, "Psi_11", DofInterpolation);
    DOF *Psi_12 = phgDofNew(g, dof_tp, 1, "Psi_12", DofInterpolation);
    DOF *Psi_13 = phgDofNew(g, dof_tp, 1, "Psi_13", DofInterpolation);
    DOF *Psi_22 = phgDofNew(g, dof_tp, 1, "Psi_22", DofInterpolation);
    DOF *Psi_23 = phgDofNew(g, dof_tp, 1, "Psi_23", DofInterpolation);
    DOF *Psi_33 = phgDofNew(g, dof_tp, 1, "Psi_33", DofInterpolation);

    //DOF *dofs_Psi[10] = {Psi_00, Psi_01, Psi_02, Psi_03, 
    //                             Psi_11, Psi_12, Psi_13,
    //                                     Psi_22, Psi_23,
    //                                             Psi_33};

    DOF *Pi_00 = phgDofNew(g, dof_tp, 1, "Pi_00", DofInterpolation);
    DOF *Pi_01 = phgDofNew(g, dof_tp, 1, "Pi_01", DofInterpolation);
    DOF *Pi_02 = phgDofNew(g, dof_tp, 1, "Pi_02", DofInterpolation);
    DOF *Pi_03 = phgDofNew(g, dof_tp, 1, "Pi_03", DofInterpolation);
    DOF *Pi_11 = phgDofNew(g, dof_tp, 1, "Pi_11", DofInterpolation);
    DOF *Pi_12 = phgDofNew(g, dof_tp, 1, "Pi_12", DofInterpolation);
    DOF *Pi_13 = phgDofNew(g, dof_tp, 1, "Pi_13", DofInterpolation);
    DOF *Pi_22 = phgDofNew(g, dof_tp, 1, "Pi_22", DofInterpolation);
    DOF *Pi_23 = phgDofNew(g, dof_tp, 1, "Pi_23", DofInterpolation);
    DOF *Pi_33 = phgDofNew(g, dof_tp, 1, "Pi_33", DofInterpolation);
    
    DOF *dofs_Pi[10] = {Pi_00, Pi_01, Pi_02, Pi_03, 
                               Pi_11, Pi_12, Pi_13,
                                      Pi_22, Pi_23,
                                             Pi_33};
   
    DOF *Phi_100 = phgDofNew(g, dof_tp, 1, "Phi_100", DofInterpolation);
    DOF *Phi_101 = phgDofNew(g, dof_tp, 1, "Phi_101", DofInterpolation);
    DOF *Phi_102 = phgDofNew(g, dof_tp, 1, "Phi_102", DofInterpolation);
    DOF *Phi_103 = phgDofNew(g, dof_tp, 1, "Phi_103", DofInterpolation);
    DOF *Phi_111 = phgDofNew(g, dof_tp, 1, "Phi_111", DofInterpolation);
    DOF *Phi_112 = phgDofNew(g, dof_tp, 1, "Phi_112", DofInterpolation);
    DOF *Phi_113 = phgDofNew(g, dof_tp, 1, "Phi_113", DofInterpolation);
    DOF *Phi_122 = phgDofNew(g, dof_tp, 1, "Phi_122", DofInterpolation);
    DOF *Phi_123 = phgDofNew(g, dof_tp, 1, "Phi_123", DofInterpolation);
    DOF *Phi_133 = phgDofNew(g, dof_tp, 1, "Phi_133", DofInterpolation);

    DOF *Phi_200 = phgDofNew(g, dof_tp, 1, "Phi_200", DofInterpolation);
    DOF *Phi_201 = phgDofNew(g, dof_tp, 1, "Phi_201", DofInterpolation);
    DOF *Phi_202 = phgDofNew(g, dof_tp, 1, "Phi_202", DofInterpolation);
    DOF *Phi_203 = phgDofNew(g, dof_tp, 1, "Phi_203", DofInterpolation);
    DOF *Phi_211 = phgDofNew(g, dof_tp, 1, "Phi_211", DofInterpolation);
    DOF *Phi_212 = phgDofNew(g, dof_tp, 1, "Phi_212", DofInterpolation);
    DOF *Phi_213 = phgDofNew(g, dof_tp, 1, "Phi_213", DofInterpolation);
    DOF *Phi_222 = phgDofNew(g, dof_tp, 1, "Phi_222", DofInterpolation);
    DOF *Phi_223 = phgDofNew(g, dof_tp, 1, "Phi_223", DofInterpolation);
    DOF *Phi_233 = phgDofNew(g, dof_tp, 1, "Phi_233", DofInterpolation);

    DOF *Phi_300 = phgDofNew(g, dof_tp, 1, "Phi_300", DofInterpolation);
    DOF *Phi_301 = phgDofNew(g, dof_tp, 1, "Phi_301", DofInterpolation);
    DOF *Phi_302 = phgDofNew(g, dof_tp, 1, "Phi_302", DofInterpolation);
    DOF *Phi_303 = phgDofNew(g, dof_tp, 1, "Phi_303", DofInterpolation);
    DOF *Phi_311 = phgDofNew(g, dof_tp, 1, "Phi_311", DofInterpolation);
    DOF *Phi_312 = phgDofNew(g, dof_tp, 1, "Phi_312", DofInterpolation);
    DOF *Phi_313 = phgDofNew(g, dof_tp, 1, "Phi_313", DofInterpolation);
    DOF *Phi_322 = phgDofNew(g, dof_tp, 1, "Phi_322", DofInterpolation);
    DOF *Phi_323 = phgDofNew(g, dof_tp, 1, "Phi_323", DofInterpolation);
    DOF *Phi_333 = phgDofNew(g, dof_tp, 1, "Phi_333", DofInterpolation);

    DOF *dofs_Phi[30] ={Phi_100, Phi_101, Phi_102, Phi_103, 
                   Phi_111, Phi_112, Phi_113,
                   Phi_122, Phi_123,
                   Phi_133, 
                   Phi_200, Phi_201, Phi_202, Phi_203,
                   Phi_211, Phi_212, Phi_213,
                   Phi_222, Phi_223,
                   Phi_233,
                   Phi_300, Phi_301, Phi_302, Phi_303,
                   Phi_311, Phi_312, Phi_313,
                   Phi_322, Phi_323,
                   Phi_333};
    DOF *dofs_var[50] ={ Psi_00, Psi_01, Psi_02, Psi_03,
                                 Psi_11, Psi_12, Psi_13,
                                         Psi_22, Psi_23,
                                                 Psi_33,
                         Pi_00, Pi_01, Pi_02, Pi_03,
                                Pi_11, Pi_12, Pi_13,
                                       Pi_22, Pi_23,
                                              Pi_33,
                   Phi_100, Phi_101, Phi_102, Phi_103, 
                   Phi_111, Phi_112, Phi_113,
                   Phi_122, Phi_123,
                   Phi_133, 
                   Phi_200, Phi_201, Phi_202, Phi_203,
                   Phi_211, Phi_212, Phi_213,
                   Phi_222, Phi_223,
                   Phi_233,
                   Phi_300, Phi_301, Phi_302, Phi_303,
                   Phi_311, Phi_312, Phi_313,
                   Phi_322, Phi_323,
                   Phi_333};

    //define rhs dofs
    DOF *rhsPsi_00 = phgDofNew(g, dof_tp, 1, "rhsPsi_00", DofInterpolation);
    DOF *rhsPsi_01 = phgDofNew(g, dof_tp, 1, "rhsPsi_01", DofInterpolation);
    DOF *rhsPsi_02 = phgDofNew(g, dof_tp, 1, "rhsPsi_02", DofInterpolation);
    DOF *rhsPsi_03 = phgDofNew(g, dof_tp, 1, "rhsPsi_03", DofInterpolation);
    DOF *rhsPsi_11 = phgDofNew(g, dof_tp, 1, "rhsPsi_11", DofInterpolation);
    DOF *rhsPsi_12 = phgDofNew(g, dof_tp, 1, "rhsPsi_12", DofInterpolation);
    DOF *rhsPsi_13 = phgDofNew(g, dof_tp, 1, "rhsPsi_13", DofInterpolation);
    DOF *rhsPsi_22 = phgDofNew(g, dof_tp, 1, "rhsPsi_22", DofInterpolation);
    DOF *rhsPsi_23 = phgDofNew(g, dof_tp, 1, "rhsPsi_23", DofInterpolation);
    DOF *rhsPsi_33 = phgDofNew(g, dof_tp, 1, "rhsPsi_33", DofInterpolation);

    DOF *rhsPi_00 = phgDofNew(g, dof_tp, 1, "rhsPi_00", DofInterpolation);
    DOF *rhsPi_01 = phgDofNew(g, dof_tp, 1, "rhsPi_01", DofInterpolation);
    DOF *rhsPi_02 = phgDofNew(g, dof_tp, 1, "rhsPi_02", DofInterpolation);
    DOF *rhsPi_03 = phgDofNew(g, dof_tp, 1, "rhsPi_03", DofInterpolation);
    DOF *rhsPi_11 = phgDofNew(g, dof_tp, 1, "rhsPi_11", DofInterpolation);
    DOF *rhsPi_12 = phgDofNew(g, dof_tp, 1, "rhsPi_12", DofInterpolation);
    DOF *rhsPi_13 = phgDofNew(g, dof_tp, 1, "rhsPi_13", DofInterpolation);
    DOF *rhsPi_22 = phgDofNew(g, dof_tp, 1, "rhsPi_22", DofInterpolation);
    DOF *rhsPi_23 = phgDofNew(g, dof_tp, 1, "rhsPi_23", DofInterpolation);
    DOF *rhsPi_33 = phgDofNew(g, dof_tp, 1, "rhsPi_33", DofInterpolation);
    
    DOF *rhsPhi_100 = phgDofNew(g, dof_tp, 1, "rhsPhi_100", DofInterpolation);
    DOF *rhsPhi_101 = phgDofNew(g, dof_tp, 1, "rhsPhi_101", DofInterpolation);
    DOF *rhsPhi_102 = phgDofNew(g, dof_tp, 1, "rhsPhi_102", DofInterpolation);
    DOF *rhsPhi_103 = phgDofNew(g, dof_tp, 1, "rhsPhi_103", DofInterpolation);
    DOF *rhsPhi_111 = phgDofNew(g, dof_tp, 1, "rhsPhi_111", DofInterpolation);
    DOF *rhsPhi_112 = phgDofNew(g, dof_tp, 1, "rhsPhi_112", DofInterpolation);
    DOF *rhsPhi_113 = phgDofNew(g, dof_tp, 1, "rhsPhi_113", DofInterpolation);
    DOF *rhsPhi_122 = phgDofNew(g, dof_tp, 1, "rhsPhi_122", DofInterpolation);
    DOF *rhsPhi_123 = phgDofNew(g, dof_tp, 1, "rhsPhi_123", DofInterpolation);
    DOF *rhsPhi_133 = phgDofNew(g, dof_tp, 1, "rhsPhi_133", DofInterpolation);

    DOF *rhsPhi_200 = phgDofNew(g, dof_tp, 1, "rhsPhi_200", DofInterpolation);
    DOF *rhsPhi_201 = phgDofNew(g, dof_tp, 1, "rhsPhi_201", DofInterpolation);
    DOF *rhsPhi_202 = phgDofNew(g, dof_tp, 1, "rhsPhi_202", DofInterpolation);
    DOF *rhsPhi_203 = phgDofNew(g, dof_tp, 1, "rhsPhi_203", DofInterpolation);
    DOF *rhsPhi_211 = phgDofNew(g, dof_tp, 1, "rhsPhi_211", DofInterpolation);
    DOF *rhsPhi_212 = phgDofNew(g, dof_tp, 1, "rhsPhi_212", DofInterpolation);
    DOF *rhsPhi_213 = phgDofNew(g, dof_tp, 1, "rhsPhi_213", DofInterpolation);
    DOF *rhsPhi_222 = phgDofNew(g, dof_tp, 1, "rhsPhi_222", DofInterpolation);
    DOF *rhsPhi_223 = phgDofNew(g, dof_tp, 1, "rhsPhi_223", DofInterpolation);
    DOF *rhsPhi_233 = phgDofNew(g, dof_tp, 1, "rhsPhi_233", DofInterpolation);

    DOF *rhsPhi_300 = phgDofNew(g, dof_tp, 1, "rhsPhi_300", DofInterpolation);
    DOF *rhsPhi_301 = phgDofNew(g, dof_tp, 1, "rhsPhi_301", DofInterpolation);
    DOF *rhsPhi_302 = phgDofNew(g, dof_tp, 1, "rhsPhi_302", DofInterpolation);
    DOF *rhsPhi_303 = phgDofNew(g, dof_tp, 1, "rhsPhi_303", DofInterpolation);
    DOF *rhsPhi_311 = phgDofNew(g, dof_tp, 1, "rhsPhi_311", DofInterpolation);
    DOF *rhsPhi_312 = phgDofNew(g, dof_tp, 1, "rhsPhi_312", DofInterpolation);
    DOF *rhsPhi_313 = phgDofNew(g, dof_tp, 1, "rhsPhi_313", DofInterpolation);
    DOF *rhsPhi_322 = phgDofNew(g, dof_tp, 1, "rhsPhi_322", DofInterpolation);
    DOF *rhsPhi_323 = phgDofNew(g, dof_tp, 1, "rhsPhi_323", DofInterpolation);
    DOF *rhsPhi_333 = phgDofNew(g, dof_tp, 1, "rhsPhi_333", DofInterpolation);

    DOF *dofs_rhs[50] ={rhsPsi_00, rhsPsi_01, rhsPsi_02, rhsPsi_03,
                                 rhsPsi_11, rhsPsi_12, rhsPsi_13,
                                         rhsPsi_22, rhsPsi_23,
                                                 rhsPsi_33,
                         rhsPi_00, rhsPi_01, rhsPi_02, rhsPi_03,
                                rhsPi_11, rhsPi_12, rhsPi_13,
                                       rhsPi_22, rhsPi_23,
                                              rhsPi_33,
                   rhsPhi_100, rhsPhi_101, rhsPhi_102, rhsPhi_103, 
                   rhsPhi_111, rhsPhi_112, rhsPhi_113,
                   rhsPhi_122, rhsPhi_123,
                   rhsPhi_133, 
                   rhsPhi_200, rhsPhi_201, rhsPhi_202, rhsPhi_203,
                   rhsPhi_211, rhsPhi_212, rhsPhi_213,
                   rhsPhi_222, rhsPhi_223,
                   rhsPhi_233,
                   rhsPhi_300, rhsPhi_301, rhsPhi_302, rhsPhi_303,
                   rhsPhi_311, rhsPhi_312, rhsPhi_313,
                   rhsPhi_322, rhsPhi_323,
                   rhsPhi_333};
 

    //for testing
    for(int i=0;i<10;i++){
        phgDofSetDataByValue(dofs_Pi[i], 0);
        phgDofSetDataByValue(dofs_Phi[i], 0);
        phgDofSetDataByValue(dofs_Phi[10 + i], 0);
        phgDofSetDataByValue(dofs_Phi[20 + i], 0);
    } 
    phgDofSetDataByValue(Psi_00, -1);
    phgDofSetDataByValue(Psi_11, 1);
    phgDofSetDataByValue(Psi_22, 1);
    phgDofSetDataByValue(Psi_33, 1);
    
    if(phgRank == 0){    
        t0 = phgGetTime(NULL);
    }   
    get_rhs_dof(dofs_var, dofs_rhs);

    if(phgRank == 0){ 
        t1 = phgGetTime(NULL);
        phgPrintf("Total time cost = %f\n\n", t1 - t0);
    }

    /*Export VTK*/

    /*release the mem*/
    for(int i=0;i<50;i++){
        phgDofFree(dofs_var + i);
        phgDofFree(dofs_rhs + i);
        }
    phgFreeGrid(&g); 
 
    phgFinalize(); 

    return 0;
}

