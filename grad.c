#include "quad.c"

static void
build_linear_system(SOLVER *solver, DOF *u_h)
{
    int n, k, i, j, s, ss, coord, in_out;
    GRID *g = u_h->g;
    ELEMENT *e, *neigh_e;
    FLOAT val_ext, val_int, mass_term, boundary_term, stiff_term[3];

    ForAllElements(g, e){ 
    int nbas = DofGetNBas(u_h, e); 
    int I[6*nbas];
        
        for(n = 0; n<nbas; n++){
            
            /*stiffness part*/ 
            phgQuadDofGradBas(e, u_h, u_h, n, 10, stiff_term);

            /*mass matrix*/
            for(coord = 0; coord<3; coord++){
                for(in_out = 0; in_out < 2; in_out++){
                    /*I[]建立了一个映射，每个单元的自由度在线性解法器中被给予一个编号*/
                    i = 6*n + 2*coord + in_out;//i是e单元中基函数编号为n、坐标编号为coord、flux编号为in_out对应的自由度编号 

                    I[i] = phgSolverMapE2L(solver, 0, e, i);

                    /*mass part*/
                    for(k = 0; k <= n; k++){
                        mass_term = phgQuadBasDotBas(e, u_h, n, u_h, k, QUAD_DEFAULT);
                        j = 6*k + 2*coord + in_out;
                        /*添加到解法器中矩阵的相应项*/
                        phgSolverAddMatrixEntry(solver, I[i], I[j], mass_term);
                        if(k < n){
                            phgSolverAddMatrixEntry(solver, I[j], I[i], mass_term);
                        }
                    }

                    /*添加stiffness到解法器相应右端项*/
                    phgSolverAddRHSEntry(solver, I[i], - stiff_term[coord]);

                
                    /*boundary term, using numerical flux*/ 
                    for (s=0; s<4; s++){
                        const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);                        
                        
                        if(e->bound_type[s] & INTERIOR){//如果s是内部面,找到邻居单元及s面在其上的编号ss
                            neigh_e = (ELEMENT *)e->neighbours[s];
                            for(ss = 0; ss < 4; ss++){ 
                                if(neigh_e->faces[ss] == e->faces[s]) break; 
                            }               
        
                            val_ext = dgQuadFaceNeighDofDotBas(e, u_h, s, n, neigh_e, 
                                    u_h, ss, QUAD_DEFAULT);//计算u^+ 
                            val_int = phgQuadFaceDofDotBas(e, s, u_h, DOF_PROJ_NONE, 
                                    u_h, n, QUAD_DEFAULT); //计算u^-
                            boundary_term = 0.5 * normal[coord] * (val_ext + val_int) + 
                                    (0.5 - in_out) * fabs(normal[coord]) * (val_ext - val_int);//数值通量, in_out = 0 对应于取plus 
                            
                            /*添加到解法器相应右端项*/
                            phgSolverAddRHSEntry(solver, I[i], boundary_term);
                        } 
                        else{//s是边界面时,只有u^- 
                            boundary_term = normal[coord] * phgQuadFaceDofDotBas(e, s,
                                    u_h, DOF_PROJ_NONE, u_h, n, QUAD_DEFAULT); 
                                
                            phgSolverAddRHSEntry(solver, I[i], boundary_term);
                        } 
                    } 
                }
            }
        }
    }
}

static void
dgHJGradDof(DOF *u_h, DOF *dof_grad)
{
    SOLVER *solver;
    solver = phgSolverCreate(SOLVER_DEFAULT, dof_grad, NULL);
    
    build_linear_system(solver, u_h);
    
    phgSolverSolve(solver, TRUE, dof_grad, NULL);
    phgSolverDestroy(&solver);

}

static void
get_dofs_grad(DOF **dofs_var, DOF ** dofs_grad, int ndof)
{
    for(INT i=0;i<ndof;i++)
        { 
            dgHJGradDof(dofs_var[i], dofs_grad[i]);
        }
}

















