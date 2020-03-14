#include "phg.h"
#include <string.h>
#include <math.h>

#define NVAR 50

#include "hdw.c"
#include "auxi_dofs.c"
#include "Hhat.c"
#include "rhs.c"
#include "rk2.c"
//#include "rk4.c"
#include "Sch_Kerr_Schild.c"
//#include "Schwarzschild_Harmonic.c"
//#include "Schwarzschild_Horizon_Penitrating.c"
//#include "Minkovski.c"
#include "error.c"

int 
main(int argc, char *argv[])
{
    //char *meshfile ="./mesh/hollowed_icosahedron.mesh";
    char *meshfile ="./mesh/SinS.albert";
    GRID *g; 
    ELEMENT *e;
    FLOAT ele_diam, min_diam = 1000.0, max_diam = 0.0;
    FLOAT t0 = 0.0, t1 = 0.0; 
    FLOAT dt, max_time = 100*M;
    DOF_TYPE *dg_type;
    INT i, j, p_order = 2, refine_time = 0;
    MPI_Status status; 

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
            phgPrintf("\nUnavailable polynomial order, using default set DOF_DG2\n");
    }

    g = phgNewGrid(-1); 
    phgImport(g, meshfile, FALSE);
    phgBalanceGrid(g, 1.2, 1, NULL, 0.);
    phgRefineAllElements(g, refine_time);
    phgBalanceGrid(g, 1.2, 1, NULL, 0.);

    phgPrintf("\nUsing mesh file: %s\n", meshfile);
    phgPrintf("Refine times: %d\n", refine_time);
    phgPrintf("Highest polymonial order: %d\n", dg_type->order);
    phgPrintf("Total elements = %d\n", g->nelem_global);
    phgPrintf("Total processes = %d\n", phgNProcs);


    /*find the max and min diameter of all elements*/
    ForAllElements(g, e){
        ele_diam = phgGeomGetDiameter(g ,e);
        max_diam = (max_diam > ele_diam)?max_diam:ele_diam;
        min_diam = (min_diam < ele_diam)?min_diam:ele_diam;
    }
    if(phgRank!=0){
        MPI_Send(&max_diam, 1, PHG_MPI_FLOAT, 0, phgRank, phgComm);
        MPI_Send(&min_diam, 1, PHG_MPI_FLOAT, 0, phgNProcs + phgRank, phgComm);
    }
    if(phgRank==0){
        FLOAT recv_diam;
        for(i=1;i<phgNProcs;i++){
            MPI_Recv(&recv_diam, 1, PHG_MPI_FLOAT, i, i, phgComm, &status);
            max_diam = (max_diam > recv_diam)?max_diam:recv_diam;
            
            MPI_Recv(&recv_diam, 1, PHG_MPI_FLOAT, i, phgNProcs + i, phgComm, &status);
            min_diam = (min_diam < recv_diam)?min_diam:recv_diam;
        } 
        for(i=1;i<phgNProcs;i++){
            MPI_Send(&max_diam, 1, PHG_MPI_FLOAT, i, i, phgComm);
            MPI_Send(&min_diam, 1, PHG_MPI_FLOAT, i, phgNProcs + i, phgComm);
        }
    }
    if(phgRank!=0){
        MPI_Recv(&max_diam, 1, PHG_MPI_FLOAT, 0, phgRank, phgComm, &status);
        MPI_Recv(&min_diam, 1, PHG_MPI_FLOAT, 0, phgNProcs + phgRank, phgComm, &status);
    }
    
    phgPrintf("\nmax_diam = %.16lf\n", max_diam);
    phgPrintf("min_diam = %.16lf\n\n", min_diam);

    dt = 0.3*min_diam; 
    //dt = 0.3*Pow(10000./g->nelem_global, 1./3.); 

    /*creat dofs for all the functions to be solved*/ 
    DOF *dofs_Psi[10], *dofs_Pi[10], *dofs_Phi[30];
    DOF *dofs_sol[50], *dofs_bdry[50], *dofs_err[NVAR];
    create_dofs(g, dg_type, 1, dofs_Psi, "Psi", 10);
    create_dofs(g, dg_type, 1, dofs_Pi, "Pi", 10);
    create_dofs(g, dg_type, 1, dofs_Phi, "Phi", 30);
    create_dofs(g, dg_type, 1, dofs_sol, "sol", 50);
    create_dofs(g, dg_type, 1, dofs_bdry, "bdry", NVAR);
    create_dofs(g, dg_type, 1, dofs_err, "err", NVAR);

    /*create dofs for all source terms*/
    DOF *dofs_srcPsi[10], *dofs_srcPi[10], *dofs_srcPhi[30];
    create_dofs(g, dg_type, 1, dofs_srcPsi, "srcPsi", 10);
    create_dofs(g, dg_type, 1, dofs_srcPi, "srcPi", 10);
    create_dofs(g, dg_type, 1, dofs_srcPhi, "srcPhi", 30);
    
    /*create dofs for derivatives of vars*/ 
    DOF *dofs_gradPsi[30], *dofs_gradPi[30], *dofs_gradPhi[90];
    DOF *dofs_gradPsi_ave[30], *dofs_gradPsi_diff[30], *dofs_gradPsi_err[30];
    create_dofs(g, dg_type, 2, dofs_gradPsi, "gradPsi", 30);
    create_dofs(g, dg_type, 2, dofs_gradPi, "gradPi", 30);
    create_dofs(g, dg_type, 2, dofs_gradPhi, "gradPhi", 90);
    create_dofs(g, dg_type, 1, dofs_gradPsi_ave, "gradPsi_ave", 30);
    create_dofs(g, dg_type, 1, dofs_gradPsi_diff, "gradPsi_diff", 30);
    create_dofs(g, dg_type, 1, dofs_gradPsi_err, "gradPsi_err", 30);

    //create lists to store dofs of var, src
    DOF *dofs_var[NVAR], *dofs_src[NVAR];
    for(i=0;i<10;i++){
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
    DOF *dofs_C[4];
    create_dofs(g, dg_type, 1, dofs_g, "g", 6);
    create_dofs(g, dg_type, 1, dofs_N, "N", 4);
    create_dofs(g, dg_type, 1, dofs_Hhat, "Hhat", NVAR);  
    create_dofs(g, dg_type, 1, dofs_rhs, "rhs", NVAR);  
    create_dofs(g, dg_type, 1, dofs_C, "C", 4);

    /*create dofs for gauge functions H_a and their derivatives*/
    DOF *dofs_H[4], *dofs_deriH[16];
    create_dofs(g, dg_type, 1, dofs_H, "H", 4);
    create_dofs(g, dg_type, 1, dofs_deriH, "deriH", 16);

    //set dof data
    set_data_var(dofs_var);
    set_data_H(dofs_H);
    set_data_deriH(dofs_deriH);
    //phgExportVTKn(g, "var.vtk", 50, dofs_var);
    //phgExportVTKn(g, "H.vtk", 50, dofs_H);
    //phgExportVTKn(g, "deriH.vtk", 50, dofs_deriH);
    copy_dofs(dofs_var, dofs_sol, "sol", NVAR);
    copy_dofs(dofs_var, dofs_bdry, "bdry", NVAR);

    //get_dofs_rhs(dofs_var, dofs_bdry, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src, dofs_C,
    //    dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);

    //for(j=0;j<NVAR;j++){
    //        phgPrintf("L2 norm of rhs_%d: %.16lf\n", j, phgDofNormL2(dofs_rhs[j]));
    //    }

    ////phgExportVTKn(g, "rhs.vtk", 50, dofs_rhs);
    //phgAbort(0);

    //for(i=0;i<30;i++){
    //    split_dof(dofs_gradPsi[i], dofs_gradPsi_ave[i], dofs_gradPsi_diff[i]);
    //}
    //get_dofs_diff(dofs_gradPsi_ave, dofs_var + 20, dofs_gradPsi_err, 30);
    //for(j=0;j<10;j++){
    //    phgPrintf("L2 error of gradPsi[%d] = %.16lf\n", j, phgDofNormL2(dofs_gradPsi_err[j]));
    //    phgPrintf("L2 error of gradPsi[%d] = %.16lf\n", j + 10, phgDofNormL2(dofs_gradPsi_err[j+10]));
    //    phgPrintf("L2 error of gradPsi[%d] = %.16lf\n\n", j + 20, phgDofNormL2(dofs_gradPsi_err[j+20]));
    //}

    ///*test blocks*/
    //DOF *u, *gradz_numerical, *gradz_ave, *gradz_diff, *gradz_analytic, *gradz_err, *compare_err;
    //u = phgDofNew(g, dg_type, 1, "u", DofInterpolation);
    //gradz_numerical = phgDofNew(g, dg_type, 2, "numerical", DofInterpolation);
    //gradz_ave = phgDofNew(g, dg_type, 1, "ave", DofInterpolation);
    //gradz_diff = phgDofNew(g, dg_type, 1, "diff", DofInterpolation);
    //gradz_analytic = phgDofNew(g, dg_type, 1, "analytic", DofNoAction);
    //gradz_err = phgDofNew(g, dg_type, 1, "err", DofInterpolation);
    //compare_err = phgDofNew(g, dg_type, 1, "compare_err", DofInterpolation);
    //phgDofSetDataByFunction(u, func_Psi33);
    //phgDofSetDataByFunction(gradz_analytic, func_Phi333);
    //
    //NEIGHBOUR_DATA *nd = phgDofInitNeighbourData(u, NULL);
    //dgHJGradDof(u, nd, u, gradz_numerical, 2);
    //get_dof_grad_hat(gradz_numerical);
    //
    //split_dof(gradz_numerical, gradz_ave, gradz_diff);
    ////phgDofDump(gradz_ave);
    //phgDofCopy(gradz_ave, &gradz_err, NULL, "err");
    //phgDofCopy(dofs_gradPsi_ave[29], &compare_err, NULL, "compare_err");
    //phgDofAXPBY(1.0, gradz_analytic, -1.0, &gradz_err);
    //phgDofAXPBY(1.0, gradz_ave, -1.0, &compare_err);
    //
    //phgPrintf("L2 err of gradz: %.16lf\n", phgDofNormL2(gradz_err));
    //phgPrintf("L2 err of compare: %.16lf\n", phgDofNormL2(compare_err));
    //phgAbort(0);



    t0 = phgGetTime(NULL);
    char Hhat_name[30], rhs_name[30], err_name[30]; 
    for(i=0;i*dt<max_time;i++){
        sprintf(rhs_name, "rhs_%lf", i*dt);
        sprintf(Hhat_name, "Hhat_%lf", i*dt);
        sprintf(err_name, "err_%lf", i*dt);
        strcat(rhs_name, ".vtk");
        strcat(Hhat_name, ".vtk");
        strcat(err_name, ".vtk");

        get_dofs_diff(dofs_var, dofs_sol, dofs_err, NVAR); 

        phgPrintf("\n\nStep %d completed\n", i);
        phgPrintf("Advanced in time: %lf\n", i*dt);
        
        t1 = phgGetTime(NULL);
        phgPrintf("Wall time past: %lf\n", t1 - t0);
        for(j=0;j<NVAR;j++){
            phgPrintf("L2 norm of rhs_%d: %.16lf\n", j, phgDofNormL2(dofs_rhs[j]));
            phgPrintf("L2 err of var_%d: %.16lf\n\n", j, phgDofNormL2(dofs_err[j]));
        }

        for(j=0;j<4;j++){
            phgPrintf("L2 norm of C_%d: %.16lf\n", j, phgDofNormL2(dofs_C[j]));
        }
        //phgExportVTKn(g, "Hhat.vtk", 50, dofs_Hhat + 0);
        //phgExportVTKn(g, "rhs.vtk", 50, dofs_rhs + 0);
        //phgExportVTKn(g, "src.vtk", 50, dofs_src + 0);
        //phgExportVTKn(g, "var_err.vtk", NVAR, dofs_err + 0);
            
        rk2(dt, dofs_var, dofs_bdry, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src, dofs_C,
               dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
    }
     
    t1 = phgGetTime(NULL);
    phgPrintf("Total time cost = %lf\n", t1 - t0);
   
    /*release the mem*/
    free_dofs(dofs_var, NVAR);
    free_dofs(dofs_src, NVAR);
    free_dofs(dofs_Hhat, NVAR);
    free_dofs(dofs_sol, NVAR);
    free_dofs(dofs_bdry, NVAR);
    free_dofs(dofs_err, NVAR);
    free_dofs(dofs_rhs, NVAR);
    
    free_dofs(dofs_gradPsi, 30);
    free_dofs(dofs_gradPsi_ave, 30);
    free_dofs(dofs_gradPsi_diff, 30);
    free_dofs(dofs_gradPsi_err, 30);
    free_dofs(dofs_gradPi, 30);
    free_dofs(dofs_gradPhi, 90);

    free_dofs(dofs_g, 6);
    free_dofs(dofs_N, 4);
    free_dofs(dofs_H, 4);
    free_dofs(dofs_deriH, 16);
    phgFreeGrid(&g); 

    phgFinalize(); 

    return 0;
}

