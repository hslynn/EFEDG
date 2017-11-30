#include "phg.h"

#include <string.h>
#include <math.h>


static void 
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) 
{   
    *value = x + 2 * y + 0 * z;
}

void
phgQuadDofGradBas(ELEMENT *e, DOF *p, DOF *v, int m, int order, FLOAT *values)
{
    int i;
    const FLOAT *g1, *g2, *w;
    FLOAT d1, d2, d3, vol;
    QUAD *quad;

    assert(!SpecialDofType(v->type));
    assert(p->dim == 1);

    quad = phgQuadGetQuad3D(order);
    g1 = phgQuadGetDofValues(e, p, quad);
 
    g2 = phgQuadGetBasisGradient(e, v, m, quad);
    w = quad->weights;
    
    
    d1 = d2 = d3 = 0.;
    for (i = 0; i < quad->npoints; i++) {
        d1 += *(g1) * (*g2++) * (*w);
        d2 += *(g1) * (*g2++) * (*w);
        d3 += *(g1) * (*g2++) * (*w);
        g1++;
        w++;
    }

    vol = phgGeomGetVolume(p->g, e);
    values[0] = d1 * vol;
    values[1] = d2 * vol;
    values[2] = d3 * vol;
}


static void
build_linear_system(SOLVER *solver, DOF *u_h)
{
    int i, j, s, ss, coord, in_out, dof_no = 0;
    GRID *g = u_h->g;
    ELEMENT *e, *neigh_e;
    FLOAT val_p, val_m, mass_term, boundary_term, stiff_term[3];

    for(coord = 0; coord<3; coord++){
        for(in_out = 0; in_out < 2; in_out++){

            ForAllElements(g, e){ 
                int N = DofGetNBas(u_h, e); 
                int I[N];
        
                /*mass matrix*/
                for(i = 0; i<N; i++){
                    /*I[N]建立了一个映射，每个单元中的基函数在线性解法器中被给予一个整体编号，以便于生成一个大矩阵*/
                    I[i] = phgSolverMapE2L(solver, dof_no, e, i);
                    for(j = 0; j <= i; j++){
                        mass_term = phgQuadBasDotBas(e, u_h, i, u_h, j, QUAD_DEFAULT);
                        /*添加到解法器中矩阵的相应项*/
                        phgSolverAddMatrixEntry(solver, I[i], I[j], mass_term);
                        if(j < i){
                            phgSolverAddMatrixEntry(solver, I[j], I[i], mass_term);
                        }
                    }
                }
                
                /*stiffness part*/ 
                for(i = 0; i < N; i++){ 
                    
                    phgQuadDofGradBas(e, u_h, u_h, i, 10, stiff_term);
        
                    /*添加到解法器相应右端项*/
                    phgSolverAddRHSEntry(solver, I[i], - stiff_term[coord]);
                }
        
                for (s=0; s<4; s++){
                    const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);                        
                    
                    
                    /*boundary term, using numerical flux*/
                    if(e->bound_type[s] & INTERIOR){//如果s是内部面,找到邻居单元及s面在其上的编号ss
                        neigh_e = (ELEMENT *)e->neighbours[s];
                        for(ss = 0; ss < 4; ss++){ 
                            if(neigh_e->faces[ss] == e->faces[s]) break; 
                        }               
        
                        for(i = 0; i < N; i++){ 
                            val_p = phgQuadFaceDofDotBas(neigh_e, ss, u_h, DOF_PROJ_NONE, 
                                    u_h, i, QUAD_DEFAULT); //计算u^+ 
                            val_m = phgQuadFaceDofDotBas(e, s, u_h, DOF_PROJ_NONE, 
                                    u_h, i, QUAD_DEFAULT); //计算u^-
                            boundary_term = 0.5 * normal[coord] * (val_m + val_p) + 
                                    (in_out - 0.5) * fabs(normal[coord]) * (val_p - val_m);//数值通量 
                        
                            /*添加到解法器相应右端项*/
                            phgSolverAddRHSEntry(solver, I[i], boundary_term);
                        } 
                    } 
                    else{//s是边界面时,只有u^- 
                        
                        for(i = 0; i < N; i++){ 
                            boundary_term = normal[coord] * phgQuadFaceDofDotBas(e, s,
                                    u_h, DOF_PROJ_NONE, u_h, i, QUAD_DEFAULT); 
                            
                            phgSolverAddRHSEntry(solver, I[i], boundary_term);
                        } 
                    } 
                }
            }

            dof_no++;
        }
    }
}


int 
main(int argc, char * argv[])
{
    char *fn = "cube3.dat"; 
    //const char *matrix_fn = "mat.m", *var_name = "mat_M";
    GRID *g; 
    DOF_TYPE *dof_tp = DOF_DG2;
    DOF *u_h, *p_0, *p_1, *q_0, *q_1, *w_0, *w_1;
    SOLVER *solver;
    //MAT *mat_M;

    phgInit(&argc, &argv);	
    
    g = phgNewGrid(-1); 
    phgImport(g, fn, FALSE);

    //phgRefineAllElements(g, 1);

    u_h = phgDofNew(g, dof_tp, 1, "u_h", func_u); 
    p_0 = phgDofNew(g, dof_tp, 1, "p_0", DofInterpolation); 
    p_1 = phgDofNew(g, dof_tp, 1, "p_1", DofInterpolation); 
    q_0 = phgDofNew(g, dof_tp, 1, "q_0", DofInterpolation); 
    q_1 = phgDofNew(g, dof_tp, 1, "q_1", DofInterpolation); 
    w_0 = phgDofNew(g, dof_tp, 1, "w_0", DofInterpolation); 
    w_1 = phgDofNew(g, dof_tp, 1, "w_1", DofInterpolation); 
        
    /*建立解法器，生成矩阵及右端项*/
    solver = phgSolverCreate(SOLVER_DEFAULT, p_0, p_1, q_0, q_1, w_0, w_1, NULL);
    
    build_linear_system(solver, u_h);

   
    /*将解法器中的矩阵输出到文件*/ 
    //mat_M = phgSolverGetMat(solver);
    //phgMatDumpMATLAB(mat_M, var_name, matrix_fn);

    /*解法器求解*/
    phgSolverSolve(solver, TRUE, p_0, p_1, q_0, q_1, w_0, w_1, NULL);
    phgSolverDestroy(&solver);

    /*输出结果*/
    phgDofDump(p_0);
    phgDofDump(p_1);

    phgDofDump(q_0);
    phgDofDump(q_1);

    phgDofDump(w_0);
    phgDofDump(w_1);

    
    
    /*Export VTK*/
    phgExportVTK(g, "tmp.vtk", p_0, p_1, q_0, q_1, w_0, w_1, NULL);


    /*release the mem*/
    phgDofFree(&u_h); 
    phgDofFree(&p_0); 
    phgDofFree(&p_1);
    phgDofFree(&q_0);
    phgDofFree(&q_1);
    phgDofFree(&w_0);
    phgDofFree(&w_1);

    phgFreeGrid(&g); 
 
    phgFinalize(); 

    return 0;
}
