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
    char *meshfile ="cube.dat";
    //char *meshfile = "icosahedron.mesh"; 
    //const char *matrix_fn = "mat.m", *var_name = "mat_M";
    GRID *g; 
    DOF_TYPE *dof_tp = DOF_DG2;
    DOF *Psi_00, *Psi_01, *Psi_02, *Psi_03, *Psi_11, *Psi_12, *Psi_13, *Psi_22, *Psi_23, *Psi_33;
    DOF *Pi_00, *Pi_01, *Pi_02, *Pi_03, *Pi_11, *Pi_12, *Pi_13, *Pi_22, *Pi_23, *Pi_33;
    DOF *Phi_100, *Phi_101, *Phi_102, *Phi_103, *Phi_111, *Phi_112, *Phi_113, *Phi_122, *Phi_123, *Phi_133;
    DOF *Phi_200, *Phi_201, *Phi_202, *Phi_203, *Phi_211, *Phi_212, *Phi_213, *Phi_222, *Phi_223, *Phi_233;
    DOF *Phi_300, *Phi_301, *Phi_302, *Phi_303, *Phi_311, *Phi_312, *Phi_313, *Phi_322, *Phi_323, *Phi_333;

    //MAT *mat_M;

    phgInit(&argc, &argv);	


    g = phgNewGrid(-1); 
    phgImport(g, meshfile, FALSE);

    //phgRefineAllElements(g, 1);
   
    /*creat dofs for all the functions to be solved*/ 
    Psi_00 = phgDofNew(g, dof_tp, 1, "Psi_00", DofInterpolation);
    Psi_01 = phgDofNew(g, dof_tp, 1, "Psi_01", DofInterpolation);
    Psi_02 = phgDofNew(g, dof_tp, 1, "Psi_02", DofInterpolation);
    Psi_03 = phgDofNew(g, dof_tp, 1, "Psi_03", DofInterpolation);
    Psi_11 = phgDofNew(g, dof_tp, 1, "Psi_11", DofInterpolation);
    Psi_12 = phgDofNew(g, dof_tp, 1, "Psi_12", DofInterpolation);
    Psi_13 = phgDofNew(g, dof_tp, 1, "Psi_13", DofInterpolation);
    Psi_22 = phgDofNew(g, dof_tp, 1, "Psi_22", DofInterpolation);
    Psi_23 = phgDofNew(g, dof_tp, 1, "Psi_23", DofInterpolation);
    Psi_33 = phgDofNew(g, dof_tp, 1, "Psi_33", DofInterpolation);

    //DOF *dofs_Psi[10] = {Psi_00, Psi_01, Psi_02, Psi_03, 
    //                             Psi_11, Psi_12, Psi_13,
    //                                     Psi_22, Psi_23,
    //                                             Psi_33};

    Pi_00 = phgDofNew(g, dof_tp, 1, "Pi_00", DofInterpolation);
    Pi_01 = phgDofNew(g, dof_tp, 1, "Pi_01", DofInterpolation);
    Pi_02 = phgDofNew(g, dof_tp, 1, "Pi_02", DofInterpolation);
    Pi_03 = phgDofNew(g, dof_tp, 1, "Pi_03", DofInterpolation);
    Pi_11 = phgDofNew(g, dof_tp, 1, "Pi_11", DofInterpolation);
    Pi_12 = phgDofNew(g, dof_tp, 1, "Pi_12", DofInterpolation);
    Pi_13 = phgDofNew(g, dof_tp, 1, "Pi_13", DofInterpolation);
    Pi_22 = phgDofNew(g, dof_tp, 1, "Pi_22", DofInterpolation);
    Pi_23 = phgDofNew(g, dof_tp, 1, "Pi_23", DofInterpolation);
    Pi_33 = phgDofNew(g, dof_tp, 1, "Pi_33", DofInterpolation);
    
    DOF *dofs_Pi[10] = {Pi_00, Pi_01, Pi_02, Pi_03, 
                               Pi_11, Pi_12, Pi_13,
                                      Pi_22, Pi_23,
                                             Pi_33};
   
    Phi_100 = phgDofNew(g, dof_tp, 1, "Phi_100", DofInterpolation);
    Phi_101 = phgDofNew(g, dof_tp, 1, "Phi_101", DofInterpolation);
    Phi_102 = phgDofNew(g, dof_tp, 1, "Phi_102", DofInterpolation);
    Phi_103 = phgDofNew(g, dof_tp, 1, "Phi_103", DofInterpolation);
    Phi_111 = phgDofNew(g, dof_tp, 1, "Phi_111", DofInterpolation);
    Phi_112 = phgDofNew(g, dof_tp, 1, "Phi_112", DofInterpolation);
    Phi_113 = phgDofNew(g, dof_tp, 1, "Phi_113", DofInterpolation);
    Phi_122 = phgDofNew(g, dof_tp, 1, "Phi_122", DofInterpolation);
    Phi_123 = phgDofNew(g, dof_tp, 1, "Phi_123", DofInterpolation);
    Phi_133 = phgDofNew(g, dof_tp, 1, "Phi_133", DofInterpolation);

    Phi_200 = phgDofNew(g, dof_tp, 1, "Phi_200", DofInterpolation);
    Phi_201 = phgDofNew(g, dof_tp, 1, "Phi_201", DofInterpolation);
    Phi_202 = phgDofNew(g, dof_tp, 1, "Phi_202", DofInterpolation);
    Phi_203 = phgDofNew(g, dof_tp, 1, "Phi_203", DofInterpolation);
    Phi_211 = phgDofNew(g, dof_tp, 1, "Phi_211", DofInterpolation);
    Phi_212 = phgDofNew(g, dof_tp, 1, "Phi_212", DofInterpolation);
    Phi_213 = phgDofNew(g, dof_tp, 1, "Phi_213", DofInterpolation);
    Phi_222 = phgDofNew(g, dof_tp, 1, "Phi_222", DofInterpolation);
    Phi_223 = phgDofNew(g, dof_tp, 1, "Phi_223", DofInterpolation);
    Phi_233 = phgDofNew(g, dof_tp, 1, "Phi_233", DofInterpolation);

    Phi_300 = phgDofNew(g, dof_tp, 1, "Phi_300", DofInterpolation);
    Phi_301 = phgDofNew(g, dof_tp, 1, "Phi_301", DofInterpolation);
    Phi_302 = phgDofNew(g, dof_tp, 1, "Phi_302", DofInterpolation);
    Phi_303 = phgDofNew(g, dof_tp, 1, "Phi_303", DofInterpolation);
    Phi_311 = phgDofNew(g, dof_tp, 1, "Phi_311", DofInterpolation);
    Phi_312 = phgDofNew(g, dof_tp, 1, "Phi_312", DofInterpolation);
    Phi_313 = phgDofNew(g, dof_tp, 1, "Phi_313", DofInterpolation);
    Phi_322 = phgDofNew(g, dof_tp, 1, "Phi_322", DofInterpolation);
    Phi_323 = phgDofNew(g, dof_tp, 1, "Phi_323", DofInterpolation);
    Phi_333 = phgDofNew(g, dof_tp, 1, "Phi_333", DofInterpolation);

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
        
    DOF *dof_rhs = phgDofNew(g, dof_tp, 50, "dof_rhs", DofInterpolation);

    get_rhs_dof(dofs_var, dof_rhs);
    //printf("test\n");
    /*Export VTK*/
    phgExportVTK(g, "tmp.vtk", dof_rhs, NULL);

    /*release the mem*/
    phgDofFree(&dof_rhs);
    for(int i=0;i<50;i++)
        phgDofFree(dofs_var + i);
    phgFreeGrid(&g); 
 
    phgFinalize(); 

    return 0;
}

