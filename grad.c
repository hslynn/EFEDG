#include "quad.c"

static void
build_linear_system(SOLVER *solver, DOF *u_h, int coord)
{
    int n, k, i, j, s, ss, in_out;
    GRID *g = u_h->g;
    ELEMENT *e, *neigh_e;
    FLOAT val_ext, val_int, mass_term, boundary_term, stiff_term[3];

    ForAllElements(g, e){ 
    int nbas = DofGetNBas(u_h, e); 
    int I[2*nbas];
        
        for(n = 0; n<nbas; n++){
            
            /*stiffness part*/ 
            phgQuadDofGradBas(e, u_h, u_h, n, 10, stiff_term);

            /*mass matrix*/
            for(in_out = 0; in_out < 2; in_out++){
                /*I[]建立了一个映射，每个单元的自由度在线性解法器中被给予一个编号*/
                i = n*2 + in_out;//i是e单元中基函数编号为n、flux编号为in_out对应的自由度编号 

                I[i] = phgSolverMapE2L(solver, 0, e, i);

                /*mass part*/
                for(k = 0; k <= n; k++){
                    mass_term = phgQuadBasDotBas(e, u_h, n, u_h, k, QUAD_DEFAULT);
                    j = 2*k + in_out;
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
                    else{//s是边界面时,施加边界条件 
                        boundary_term = normal[coord] * phgQuadFaceDofDotBas(e, s,
                                u_h, DOF_PROJ_NONE, u_h, n, QUAD_DEFAULT); 
                            
                        phgSolverAddRHSEntry(solver, I[i], boundary_term);
                    } 
                } 
            }
        }
    }
}

static void
dgHJGradDof(DOF *u_h, DOF *dof_grad, int coord)
{
    SOLVER *solver;
    solver = phgSolverCreate(SOLVER_DEFAULT, dof_grad, NULL);
    
    build_linear_system(solver, u_h, coord);
    
    phgSolverSolve(solver, TRUE, dof_grad, NULL);
    phgSolverDestroy(&solver);
}

static void
get_dofs_grad(DOF **dofs, DOF **dofs_grad, int ndof)
{
    for(INT i=0;i<ndof;i++){ 
        dgHJGradDof(dofs[i], dofs_grad[i], 0);
        dgHJGradDof(dofs[i], dofs_grad[ndof + i], 1);
        dgHJGradDof(dofs[i], dofs_grad[2*ndof + i], 2);
    }
}

static void
get_dof_grad_hat(DOF *dof_grad)
{
    
    INT i, npair = DofGetDataCount(dof_grad)/2;
    FLOAT grad_plus, grad_minus;
    FLOAT *p = DofData(dof_grad);
    
    for(i=0;i<npair;i++){
        grad_plus = *p, grad_minus = *(p+1); 
        *p = 0.5 * (grad_plus + grad_minus);
        *(++p) = 0.5 * (grad_plus - grad_minus); 
        p++;
    } 
}

static void
get_dofs_grad_hat(DOF **dofs, DOF **dofs_grad_hat, INT ndof) 
{
    get_dofs_grad(dofs, dofs_grad_hat, ndof);
    
    for(int i=0;i<3*ndof;i++){
        get_dof_grad_hat(dofs_grad_hat[i]);
    }    
}
