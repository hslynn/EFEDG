#include "phg.h"
#include <string.h>
#include <math.h>

#define NVAR 50

#include "hdw.c"
#include "auxi_dofs.c"
#include "Hhat.c"
#include "rk2.c"
#include "Schwarzschild_Harmonic.c"
#include "error.c"

int 
main(int argc, char * argv[])
{
    char *meshfile ="./mesh/hollowed_icosahedron.mesh";
    GRID *g; 
    FLOAT t0 = 0.0, t1 = 0.0;
    FLOAT dt = 0.01, max_time = 0.1;
    INT max_step = ceil(max_time/dt);
    DOF_TYPE *dg_type;
    INT i, p_order = 2, refine_time = 0;
    
    //command line options
    phgOptionsRegisterInt("p", "polynomial order of DG basis, default is 2", &p_order);
    phgOptionsRegisterString("m", "name of the mesh file, default is \"./mesh/hollowed_icsahedron.mesh\"", 
            &meshfile);
    phgOptionsRegisterInt("r", "mesh refine times, default is 0", &refine_time);

    phgInit(&argc, &argv);	
    switch(p_order){
        case 0: dg_type = DOF_DG0;
            break;
        case 1: dg_type = DOF_DG1;
            break;
        case 2: dg_type = DOF_DG2; 
            break;                
        case 3: dg_type = DOF_DG3;
            break;
        case 4: dg_type = DOF_DG4;
            break;                
        case 5: dg_type = DOF_DG5;
            break;
        case 6: dg_type = DOF_DG6;
            break;                
        case 7: dg_type = DOF_DG7;
            break;
        case 8: dg_type = DOF_DG8;
            break;                
        case 9: dg_type = DOF_DG9;
            break;
        case 10: dg_type = DOF_DG10;
            break;
        case 11: dg_type = DOF_DG11;
            break;
        case 12: dg_type = DOF_DG12;
            break;
        case 13: dg_type = DOF_DG13;
            break;
        case 14: dg_type = DOF_DG14;
            break;
        case 15: dg_type = DOF_DG15;
            break;
        default: dg_type = DOF_DG2;
            printf("Unavailable polynomial order, using default set DOF_DG2\n\n");
    }

    g = phgNewGrid(-1); 
    phgImport(g, meshfile, FALSE);
    phgRefineAllElements(g, refine_time);
    phgBalanceGrid(g, 1.2, 1, NULL, 0.);

    /*creat dofs for all the functions to be solved*/ 
    DOF *dofs_Psi[10], *dofs_Pi[10], *dofs_Phi[30];
    DOF *dofs_sol[50], *dofs_bdry[50], *dofs_diff[NVAR];
    create_dofs(g, dg_type, 1, dofs_Psi, "Psi", 10);
    create_dofs(g, dg_type, 1, dofs_Pi, "Pi", 10);
    create_dofs(g, dg_type, 1, dofs_Phi, "Phi", 30);
    create_dofs(g, dg_type, 1, dofs_sol, "sol", 50);
    create_dofs(g, dg_type, 1, dofs_bdry, "bdry", NVAR);
    create_dofs(g, dg_type, 1, dofs_diff, "diff", NVAR);

    /*create dofs for all source terms*/
    DOF *dofs_srcPsi[10], *dofs_srcPi[10], *dofs_srcPhi[30];
    create_dofs(g, dg_type, 1, dofs_srcPsi, "srcPsi", 10);
    create_dofs(g, dg_type, 1, dofs_srcPi, "srcPi", 10);
    create_dofs(g, dg_type, 1, dofs_srcPhi, "srcPhi", 30);
    
    /*create dofs for derivatives of vars*/ 
    DOF *dofs_gradPsi[30], *dofs_gradPi[30], *dofs_gradPhi[90];
    create_dofs(g, dg_type, 2, dofs_gradPsi, "gradPsi", 30);
    create_dofs(g, dg_type, 2, dofs_gradPi, "gradPi", 30);
    create_dofs(g, dg_type, 2, dofs_gradPhi, "gradPhi", 90);

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
    create_dofs(g, dg_type, 1, dofs_g, "g", 6);
    create_dofs(g, dg_type, 1, dofs_N, "N", 4);
    create_dofs(g, dg_type, 1, dofs_Hhat, "Hhat", NVAR);  

    //set dof data
    set_data_dofs(dofs_var);
    phgExportVTKn(g, "var.vtk", 50, dofs_var);
    copy_dofs(dofs_var, dofs_sol, "sol", NVAR);
    copy_dofs(dofs_var, dofs_bdry, "bdry", NVAR);
   
    t0 = phgGetTime(NULL);
    
    get_dofs_rhs(dofs_var, dofs_bdry, dofs_g, dofs_N, dofs_src,
        dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
    
    char Hhat_name[30], rhs_name[30]; 
    for(i=0;i<max_step;i++){
        
        //sprintf(Hhat_name, "Hhat_%f", i*dt);
        //sprintf(rhs_name, "rhs_%f", i*dt);
        //strcat(Hhat_name, ".vtk");
        //strcat(rhs_name, ".vtk");
        if(phgRank==0){
            printf("step %d completed\n", i);
            printf("time length: %f\n\n", i*dt);
        }
        //phgExportVTKn(g, Hhat_name, 50, dofs_Hhat + 0);
        //phgExportVTKn(g, rhs_name, 50, dofs_rhs + 0);
        
        ssp_rk2(dt, dofs_var, dofs_bdry, dofs_g, dofs_N, dofs_src,
               dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
    }
     
    phgPrintf("Using mesh file: %s\n", meshfile);
    phgPrintf("Refine times: %d\n", refine_time);
    phgPrintf("Highest polymonial order: %d\n", dg_type->order);
    phgPrintf("Total elements = %d\n", dofs_var[0]->g->nelem_global);
    phgPrintf("Total processes = %d\n", phgNProcs);
    t1 = phgGetTime(NULL);
    phgPrintf("Total time cost = %f\n\n", t1 - t0);

    get_dofs_diff(dofs_var, dofs_sol, dofs_diff); 
    phgExportVTKn(g, "diff.vtk", 10, dofs_diff);
    //compute the L2 error of dofs_diff terms
    FLOAT err[10];
    for(i=0;i<10;i++){
        err[i] = phgDofNormL2(dofs_diff[i]);
        if(phgRank == 0){
            printf("L2 error of psi[%d]: %f\n", i, err[i]);
        }
    }

    
   
    /*release the mem*/
    free_dofs(dofs_var, NVAR);
    free_dofs(dofs_src, NVAR);
    free_dofs(dofs_Hhat, NVAR);
    free_dofs(dofs_sol, NVAR);
    free_dofs(dofs_bdry, NVAR);
    free_dofs(dofs_diff, NVAR);

    free_dofs(dofs_gradPsi, 30);
    free_dofs(dofs_gradPi, 30);
    free_dofs(dofs_gradPhi, 90);

    free_dofs(dofs_g, 6);
    free_dofs(dofs_N, 4);
    phgFreeGrid(&g); 

    phgFinalize(); 

    return 0;
}

