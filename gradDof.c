#include "quad.c"

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
                    /*I[N]建立了一个映射，每个自由度对象的每个单元的基函数在线性解法器中
                        被给予一个整体编号，以便于生成一个整体矩阵*/
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
                            val_p = dgQuadFaceNeighDofDotBas(e, u_h, s, i, neigh_e, 
                                    u_h, ss, QUAD_DEFAULT);//计算u^+ 
                            val_m = phgQuadFaceDofDotBas(e, s, u_h, DOF_PROJ_NONE, 
                                    u_h, i, QUAD_DEFAULT); //计算u^-
                            boundary_term = 0.5 * normal[coord] * (val_p + val_m) + 
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
            dof_no %= 6; 
        }
    }
}

static void
dgHJGradDof(DOF *u_h, DOF *p_0, DOF *p_1, DOF *q_0, DOF *q_1, DOF *w_0, DOF *w_1)
{
    SOLVER *solver;
    solver = phgSolverCreate(SOLVER_DEFAULT, p_0, p_1, q_0, q_1, w_0, w_1, NULL);
    
    build_linear_system(solver, u_h);
    
    phgSolverSolve(solver, TRUE, p_0, p_1, q_0, q_1, w_0, w_1, NULL);
    phgSolverDestroy(&solver);

}



















