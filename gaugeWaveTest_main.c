#include <string.h>

#include "global_def.h"
#include "initial_condition.h"
#include "hdw.h"
#include "rk.h"
#include "rhs.h"
#include "runtime.h"

int 
main(int argc, char *argv[])
{
    char *meshfile ="./mesh/cube4.dat";
    GRID *g; 
    ELEMENT *e;
    FLOAT ele_diam, min_diam = 1000.0, max_diam = 0.0;
    FLOAT t0 = 0.0, t1 = 0.0; 
    FLOAT dt, max_time = 100;
    INT i, j, p_order = 1, refine_time = 0, rk_order = 3;
    INT isPeri = 0;
    DOF_TYPE *dg_type = DG_TYPE, *dg_filter = DG_FILTER;
    MPI_Status status; 

    //command line options
    phgOptionsRegisterInt("p", "polynomial order of DG basis, default is 1", &p_order);
    phgOptionsRegisterString("m", "name of the mesh file,  default is \"./mesh/hollowed_icsahedron.mesh\"", 
            &meshfile);
    phgOptionsRegisterInt("r", "mesh refine times, default is 0", &refine_time);
    phgOptionsRegisterInt("o", "Periodic or not, default is 0", &isPeri);

    phgInit(&argc, &argv);	
    switch(p_order){
        case 0: dg_type = DOF_DG0;
            break;
        case 1: dg_type = DOF_DG1;
            dg_filter = DOF_DG0;
            break;
        case 2: dg_type = DOF_DG2; 
            dg_filter = DOF_DG1;
            break;                
        case 3: dg_type = DOF_DG3;
            dg_filter = DOF_DG2;
            break;
        case 4: dg_type = DOF_DG4;
            dg_filter = DOF_DG3;
            break;                
        case 5: dg_type = DOF_DG5;
            dg_filter = DOF_DG4;
            break;
        case 6: dg_type = DOF_DG6;
            dg_filter = DOF_DG5;
            break;                
        case 7: dg_type = DOF_DG7;
            dg_filter = DOF_DG6;
            break;
        case 8: dg_type = DOF_DG8;
            dg_filter = DOF_DG7;
            break;                
        case 9: dg_type = DOF_DG9;
            dg_filter = DOF_DG8;
            break;
        case 10: dg_type = DOF_DG10;
            dg_filter = DOF_DG9;
            break;
        case 11: dg_type = DOF_DG11;
            dg_filter = DOF_DG10;
            break;
        case 12: dg_type = DOF_DG12;
            dg_filter = DOF_DG11;
            break;
        case 13: dg_type = DOF_DG13;
            dg_filter = DOF_DG12;
            break;
        case 14: dg_type = DOF_DG14;
            dg_filter = DOF_DG13;
            break;
        case 15: dg_type = DOF_DG15;
            dg_filter = DOF_DG14;
            break;
        default: dg_type = DOF_DG2;
            dg_filter = DOF_DG1;
            phgPrintf("\nUnavailable polynomial order, using default set DOF_DG2\n");
    }
    g = phgNewGrid(-1); 
    if(isPeri == 1){
        phgSetPeriodicity(g, X_MASK | Y_MASK | Z_MASK);
    }
    phgImport(g, meshfile, FALSE);
    phgBalanceGrid(g, 1.2, 1, NULL, 0.);
    if (refine_time > 0){
        phgRefineAllElements(g, refine_time);
        phgBalanceGrid(g, 1.2, 1, NULL, 0.);
    }

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

    phgPrintf("Total processes = %d\n\n", phgNProcs);

    phgPrintf("Using mesh file: %s\n", meshfile);
    phgPrintf("Refine times = %d\n", refine_time);
    phgPrintf("Total elements = %d\n", g->nelem_global);
    phgPrintf("max_diam = %.16lf\n", max_diam);
    phgPrintf("min_diam = %.16lf\n\n", min_diam);

    phgPrintf("Highest polymonial order = %d\n", dg_type->order);
    phgPrintf("Runge-Kutta order = %d\n", rk_order);
    dt = 0.05/(2.*p_order+1.)*min_diam; 
    phgPrintf("time step length = %.16lf\n\n", dt);

    phgPrintf("Spacetime type = %d\n", SPACETIME);
    phgPrintf("Black hole mass = %.2f\n", M);
    phgPrintf("Black hole spin = %.2f\n\n", A);
 
    /*creat dofs for all the functions to be solved*/ 
    DOF *dofs_Psi[10], *dofs_Pi[10], *dofs_Phi[30];
    DOF *dofs_sol[50], *dofs_exact[50], *dofs_err[NVAR];
    create_dofs(g, dg_type, 1, dofs_Psi, "Psi", 10);
    create_dofs(g, dg_type, 1, dofs_Pi, "Pi", 10);
    create_dofs(g, dg_type, 1, dofs_Phi, "Phi", 30);
    create_dofs(g, dg_type, 1, dofs_sol, "sol", 50);
    create_dofs(g, dg_type, 1, dofs_exact, "exact", NVAR);
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
    copy_dofs(dofs_var, dofs_sol, "sol", NVAR);
    copy_dofs(dofs_var, dofs_exact, "exact", NVAR);

    char fn_err[50], fn_err_inf[NVAR], fn_rhs[50], fn_C[50];
    FILE *fp_err, *fp_err_inf, *fp_rhs, *fp_C;
    sprintf(fn_err, "./data/gaugeWaveTest_err_A%.2f_r%dp%d_peri%d.data", AA, refine_time, p_order, isPeri);
    sprintf(fn_err_inf, "./data/gaugeWaveTest_err_inf_A%.2f_r%dp%d_peri%d.data", AA, refine_time, p_order, isPeri);
    sprintf(fn_rhs, "./data/gaugeWaveTest_rhs_A%.2f_r%dp%d_peri%d.data", AA, refine_time, p_order, isPeri);
    sprintf(fn_C, "./data/gaugeWaveTest_C_A%.2f_r%dp%d_peri%d.data", AA, refine_time, p_order, isPeri);
    INT steps_complished = 0;
    t0 = phgGetTime(NULL);
    FLOAT L2_err_array[NVAR], Linf_err_array[NVAR], L2_rhs_array[NVAR], L2_C_array[4];

    /* advance in time */
    for(i=steps_complished;TIME<max_time;i++){
        if(TIME+dt <= max_time){
            rk1(dt, dofs_var, dofs_exact, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src, dofs_C,
                dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
            TIME += dt;
        }
        else{
            rk1(max_time - TIME, dofs_var, dofs_exact, dofs_g, dofs_N, dofs_H, dofs_deriH, dofs_src, dofs_C,
                dofs_gradPsi, dofs_gradPi, dofs_gradPhi, dofs_Hhat, dofs_rhs);
            TIME = max_time;
        }

        phgPrintf("Step %d completed\n", i + 1);
        phgPrintf("Advanced in time: %lf\n", TIME);
        t1 = phgGetTime(NULL);
        phgPrintf("Wall time past: %lf\n\n", t1 - t0);

        set_data_var(dofs_sol);
        set_data_var(dofs_exact);
        get_dofs_diff(dofs_var, dofs_sol, dofs_err, NVAR); 
        for(j=0;j<NVAR;j++){
            L2_err_array[j] = phgDofNormL2(dofs_err[j]);
            Linf_err_array[j] = phgDofNormInftyVec(dofs_err[j]);
            L2_rhs_array[j] = phgDofNormL2(dofs_rhs[j]);
        }
        for(j=0;j<4;j++){
            L2_C_array[j] = phgDofNormL2(dofs_C[j]);
        }

        if(phgRank==0){
            fp_err=fopen(fn_err,"a");
            fp_err_inf=fopen(fn_err_inf,"a");
            fp_C=fopen(fn_C, "a");
            fp_rhs=fopen(fn_rhs, "a");

            fprintf(fp_err, "%.16lf ", TIME);
            fprintf(fp_err_inf, "%.16lf ", TIME);
            fprintf(fp_rhs, "%.16lf ", TIME);
            fprintf(fp_C, "%.16lf ", TIME);
            for(j=0;j<NVAR;j++){
                fprintf(fp_err, "%.16lf ", L2_err_array[j]);
                fprintf(fp_err_inf, "%.16lf ", Linf_err_array[j]);
                fprintf(fp_rhs, "%.16lf ", L2_rhs_array[j]);
            }
            fprintf(fp_rhs, "\n");
            fprintf(fp_err, "\n");
            fprintf(fp_err_inf, "\n");
            for(j=0;j<4;j++){
                fprintf(fp_C, "%.16lf ", L2_C_array[j]);
            }
            fprintf(fp_C, "\n");

            fclose(fp_C);
            fclose(fp_err);
            fclose(fp_err_inf);
            fclose(fp_rhs);
        }
    }
    t1 = phgGetTime(NULL);
    phgPrintf("Total time cost = %lf\n", t1 - t0);

    phgExportVTKn(g, "g00.vtk", 1, dofs_var);
    phgExportVTKn(g, "g11.vtk", 1, dofs_var + 4);
    phgExportVTKn(g, "g00_exact.vtk", 1, dofs_sol);
    phgExportVTKn(g, "g11_exact.vtk", 1, dofs_sol + 4);
   
    /*release the mem*/
    free_dofs(dofs_var, NVAR);
    free_dofs(dofs_src, NVAR);
    free_dofs(dofs_C, 4);
    free_dofs(dofs_Hhat, NVAR);
    free_dofs(dofs_sol, NVAR);
    free_dofs(dofs_exact, NVAR);
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

