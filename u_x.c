#include "phg.h"

#include <string.h>
#include <math.h>

static void 
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) 
{   
    *value = x + 0*y + 0*z;
}

float
phgQuadDofGradxBas(ELEMENT *e, DOF *p, DOF *v, int m, int order)
{
    int i;
    const FLOAT *g1, *g2, *w;
    FLOAT d1, d2, d3, vol, value;
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
    value = d1 * vol;
    //values[1] = d2 * vol;
    //values[2] = d3 * vol;
    return value;
}


static void
build_linear_system(SOLVER *solver, DOF *u_h)
{
    int i, j, s, ss;
    GRID *g = u_h->g;
    ELEMENT *e, *neigh_e;
    float val_p, val_m, mass_term, stiff_term, boundary_term;

    ForAllElements(g, e){ 
        int N = DofGetNBas(u_h, e); 
        int I[N];

        /*mass matrix*/
        for(i = 0; i<N; i++){
            /*I[N]建立了一个映射，每个单元中的基函数在线性解法器中被给予一个整体编号，以便于生成一个大矩阵*/
            I[i] = phgSolverMapE2L(solver, 0, e, i);
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
            
            stiff_term = - phgQuadDofGradxBas(e, u_h, u_h, i, 10);

            /*添加到解法器相应右端项*/
            phgSolverAddRHSEntry(solver, I[i], stiff_term);
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
                    boundary_term = 0.5 * normal[0] * (val_m + val_p) + 
                            0.5 * fabs(normal[0]) * (val_p - val_m);//数值通量 
                
                    /*添加到解法器相应右端项*/
                    phgSolverAddRHSEntry(solver, I[i], boundary_term);
                } 
            } 
            else{//s是边界面时,只有u^- 
                
                for(i = 0; i < N; i++){ 
                    boundary_term = normal[0] * phgQuadFaceDofDotBas(e, s,
                            u_h, DOF_PROJ_NONE, u_h, i, QUAD_DEFAULT); 
                    
                    phgSolverAddRHSEntry(solver, I[i], boundary_term);
                } 
            } 
        }
    }
}


int 
main(int argc, char * argv[])
{
    char *fn = "cube1o.dat"; 
    const char *matrix_fn = "mat.m", *var_name = "mat_M";
    GRID *g; 
    DOF *u_h, *p_1;
    SOLVER *solver;
    MAT *mat_M;

    phgInit(&argc, &argv);	
    
    g = phgNewGrid(-1); 
    phgImport(g, fn, FALSE);

    //phgRefineAllElements(g, 1);


    u_h = phgDofNew(g, DOF_DG1, 1, "u_h", func_u); 
    p_1 = phgDofNew(g, DOF_DG1, 1, "p_1", DofInterpolation); 
    
    /*建立解法器，生成矩阵及右端项*/
    solver = phgSolverCreate(SOLVER_DEFAULT, p_1, NULL);
    build_linear_system(solver, u_h);
   
    /*将解法器中的矩阵输出到文件*/ 
    //mat_M = phgSolverGetMat(solver);
    //phgMatDumpMATLAB(mat_M, var_name, matrix_fn);

    /*解法器求解*/
    phgSolverSolve(solver, TRUE, p_1, NULL);
    phgSolverDestroy(&solver);

    /*输出结果*/
    phgDofDump(u_h);
    phgDofDump(p_1);
    
    /*Export VTK*/
    //phgExportVTK(g, "tmp.vtk", p_1, u_h,  NULL);


    /*release the mem*/
    phgDofFree(&u_h); 
    phgDofFree(&p_1); 
 
    phgFreeGrid(&g); 
 
    phgFinalize(); 

    return 0;
}
