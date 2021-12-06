#include "phg.h"
#include "filter.h"
#include "global_def.h"
#include <math.h>

void
phgQuadDofGradBas(ELEMENT *e, DOF *p, DOF *v, int m, int order, FLOAT *values)
{
    int i;
    const FLOAT *g1, *g2, *w;
    FLOAT d1, d2, d3, vol;
    QUAD *quad;

    assert(!SpecialDofType(v->type));
    assert(p->dim == 1);

    if (order < 0){
	    order = DofTypeOrder(p, e) + DofTypeOrder(v, e);
	}
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

void
set_neighbour_data(DOF *u, ELEMENT *e, INT s, DOF *neigh_u, NEIGHBOUR_DATA *nd)
{ 
    assert(e->bound_type[s] & INTERIOR); 
    SHORT i, nbasface;
    FLOAT *ele_data;

    ele_data = DofElementData(neigh_u, e->index);
    nbasface = phgDofNeighbourNBas(nd, e, s, NULL);
    SHORT bases[nbasface];
    phgDofGetBasesOnFace(u, e, s, bases);
    for(i=0;i<nbasface;i++){
        ele_data[bases[i]] = *phgDofNeighbourData(nd, e, s, i, NULL);
    }
}

void
reset_neighbour_data(DOF *neigh_u, ELEMENT *e, INT s)
{
    assert(e->bound_type[s] & INTERIOR); 
    SHORT i, nbasface, bases[DofGetNBas(neigh_u, e)];
    FLOAT *ele_data;

    ele_data = DofElementData(neigh_u, e->index);
    nbasface = phgDofGetBasesOnFace(neigh_u, e, s, bases);
    for(i=0;i<nbasface;i++){
        ele_data[bases[i]] = 0.;   
    }
}

void
build_linear_system(SOLVER *solver, DOF *u_h, NEIGHBOUR_DATA *nd, DOF *dof_bdry_in, DOF *dof_bdry_out, INT coord)
{
    int n, k, i, j, s, in_out;
    GRID *g = u_h->g;
    DOF *neigh_u;
    ELEMENT *e;
    FLOAT val_ext, val_int, mass_term, boundary_term, stiff_term[3];
   
    neigh_u = phgDofNew(g, u_h->type, 1, "neigh_u", DofInterpolation); 

    ForAllElements(g, e){ 
    int nbas = DofGetNBas(u_h, e); 
    int I[2*nbas];
        
        for(n = 0; n<nbas; n++){
            
            /*stiffness part*/ 
            phgQuadDofGradBas(e, u_h, u_h, n, QUAD_DEFAULT, stiff_term);

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
                    FLOAT normal[Dim];
                    phgGeomGetFaceOutNormal(g, e, s, normal);                        

                    val_int = phgQuadFaceDofDotBas(e, s, u_h, DOF_PROJ_NONE, 
                                u_h, n, QUAD_DEFAULT); //计算u^-
                   
                    /*分是否边界面来计算u^+*/ 
                    if(e->bound_type[s] & INTERIOR){//如果s是内部面,则利用neigh_u来计算u^+
                        set_neighbour_data(u_h, e, s, neigh_u, nd); 
                        val_ext = phgQuadFaceDofDotBas(e, s, neigh_u, DOF_PROJ_NONE, 
                                u_h, n, QUAD_DEFAULT); 
                        reset_neighbour_data(neigh_u, e, s);
                    } 
                    else if(e->bound_type[s] & DIRICHLET){//s是外边界面,tetgen中设置边界类型为１,PHG处理为DIRICHLET
                        val_ext = phgQuadFaceDofDotBas(e, s, dof_bdry_out, DOF_PROJ_NONE, 
                                u_h, n, QUAD_DEFAULT);
                    }
                    else{
                        val_ext = phgQuadFaceDofDotBas(e, s, dof_bdry_in, DOF_PROJ_NONE, 
                                u_h, n, QUAD_DEFAULT);
                    }

                    boundary_term = 0.5 * normal[coord] * (val_ext + val_int) + (0.5 - in_out) * 
                            Fabs(normal[coord]) * (val_ext - val_int);//数值通量, in_out = 0 对应于取plus  
                    /*添加到解法器相应右端项*/
                    phgSolverAddRHSEntry(solver, I[i], boundary_term);
                } 
            }
        }
    }
    phgDofFree(&neigh_u);
}

void
dgHJGradDof(DOF *u_h, NEIGHBOUR_DATA *nd, DOF *dof_bdry_in, DOF *dof_bdry_out, DOF *dof_grad, int coord)
{
    SOLVER *solver;
    solver = phgSolverCreate(SOLVER_DEFAULT, dof_grad, NULL);
    
    build_linear_system(solver, u_h, nd, dof_bdry_in, dof_bdry_out, coord);
    
    phgSolverSolve(solver, TRUE, dof_grad, NULL);
    phgSolverDestroy(&solver);
}

void
get_dofs_grad(DOF **dofs, DOF **dofs_bdry_in, DOF **dofs_bdry_out, DOF **dofs_grad, int ndof)
{
    short i;
    NEIGHBOUR_DATA *nd;
    for(i=0;i<ndof;i++){ 
        nd = phgDofInitNeighbourData(dofs[i], NULL);
        dgHJGradDof(dofs[i], nd, dofs_bdry_in[i], dofs_bdry_out[i], dofs_grad[i], 0);
        dgHJGradDof(dofs[i], nd, dofs_bdry_in[i], dofs_bdry_out[i], dofs_grad[ndof + i], 1);
        dgHJGradDof(dofs[i], nd, dofs_bdry_in[i], dofs_bdry_out[i], dofs_grad[2*ndof + i], 2);
        phgDofReleaseNeighbourData(&nd);
    }
}

void
get_dof_grad_hat(DOF *dof_grad)
{
    
    FLOAT grad_plus, grad_minus;
    FLOAT *p;
    GRID *g = dof_grad->g;
    ELEMENT *e;
    int i, idx;
    
    int np = dof_grad->type->np_elem;
    ForAllElements(g, e){
        idx = e->index;
        p = DofElementData(dof_grad, idx);
        for(i=0;i<np;i++){
            grad_plus = *p, grad_minus = *(p+1); 
            *p = 0.5 * (grad_plus + grad_minus);
            *(++p) = 0.5 * (grad_plus - grad_minus); 
            p++;
        } 
    }
}

void
get_dofs_grad_hat(DOF **dofs, DOF **dofs_bdry_in, DOF **dofs_bdry_out, DOF **dofs_grad_hat, INT ndof) 
{
    get_dofs_grad(dofs, dofs_bdry_in, dofs_bdry_out, dofs_grad_hat, ndof);
    //phgExportVTKn(dofs[0]->g, "grad.vtk", 3*ndof, dofs_grad_hat);
    
    short i; 
    for(i=0;i<3*ndof;i++){
        get_dof_grad_hat(dofs_grad_hat[i]);
    }    
    //phgExportVTKn(dofs[0]->g, "gradhat.vtk", 3*ndof, dofs_grad_hat);
 
}


/* Old function, no parallel and the logic to get the quad is too complicated*/
//FLOAT
//dgQuadFaceNeighDofDotBas(ELEMENT *e, DOF *u, int face, int N, ELEMENT *neigh_e, 
//DOF *v, int neigh_face, int order)
//{
//    int i, j, e_v[NVert], neigh_v[NVert];
//    FLOAT d, lambda[Dim + 1], neigh_lambda[Dim + 1];
//    FLOAT *dof, *buffer;
//    const FLOAT *bas, *p, *w;
//    QUAD *quad;
//    DOF_TYPE *type;
//
//    assert(!SpecialDofType(u->type));
//    assert(face >= 0 && face <= 3);
//
//    type = (DofIsHP(u) ? u->hp->info->types[u->hp->max_order] : u->type);
//
//    if (order < 0) {
//	    i = DofTypeOrder(v, neigh_e);
//	    j = DofTypeOrder(u, e);
//	    if (i < 0)
//	        i = j;
//	    order = i + j;
//    }
//    quad = phgQuadGetQuad2D(order);
//
//    /*有更好的实现方法吗？*/
//    /*e_v[i]: face号面上的i号顶点在本单元的编号；
//      neigh_v[i]: 该顶点在face面的相邻单元的编号*/
//    for (i=0; i < NVert - 1; i++) {
//        e_v[i] = GetFaceVertex(face, i);
//        for (j=0; j < NVert; j++) {
//            if (neigh_e->verts[j] == e->verts[e_v[i]]) {
//                neigh_v[i] = j;
//                break; 
//            } 
//        }
//    }
//    lambda[face] = 0.;
//    neigh_lambda[neigh_face] = 0.;
//
//    buffer = phgAlloc(sizeof(*buffer));
//    p = quad->points;
//    w = quad->weights;
//
//    d = 0.;
//    for (i = 0; i < quad->npoints; i++) {
//        for (j = 0; j < NVert - 1; j++) {
//        lambda[e_v[j]] = *(p++);
//        neigh_lambda[neigh_v[j]] = lambda[e_v[j]];
//        }
//        
//        dof = phgDofEval(v, neigh_e, neigh_lambda, buffer);
//        bas = (FLOAT *)type->BasFuncs(u, e, N, N + 1, lambda);
//        d += *(bas++) * *(dof++) * *(w++);
//    } 
//    phgFree(buffer);
//
//    return d * phgGeomGetFaceArea(u->g, e, face);
//}


