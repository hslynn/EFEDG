#include "phg.h"

#include <string.h>
#include <math.h>

#define gamma_0 1
#define gamma_1 -1
#define gamma_2 1

const int dims_psi[2] = {4, 4};
const int dims_Pi[2] = {4, 4};
const int dims_Phi[3] = {3, 4, 4};
const int dims_g[2] = {3, 3};
const int dims_Shift[1] = {3};
const int dims_scalar[0];
const int dims_vec[1] = {4};

struct TENSOR
{ 
};

static void 
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) 
{   
    *value = (x+1.)*(x+1.);
}

static void
func_zero(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0.;
}

static void
func_one(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 1.;
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

    if (order < 0) 
	    order = DofTypeOrder(p, e);
	 
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

FLOAT
dgQuadFaceNeighDofDotBas(ELEMENT *e, DOF *u, int face, int N, ELEMENT *neigh_e, 
DOF *v, int neigh_face, int order)
{
    int i, j, e_v[NVert], neigh_v[NVert];
    FLOAT d, lambda[Dim + 1], neigh_lambda[Dim + 1];
    FLOAT *dof, *buffer;
    const FLOAT *bas, *p, *w;
    QUAD *quad;
    DOF_TYPE *type;

    assert(!SpecialDofType(u->type));
    assert(face >= 0 && face <= 3);

    type = (DofIsHP(u) ? u->hp->info->types[u->hp->max_order] : u->type);

    if (order < 0) {
	    i = DofTypeOrder(v, neigh_e);
	    j = DofTypeOrder(u, e);
	if (i < 0)
	    i = j;
	    order = i + j;
    }
    quad = phgQuadGetQuad2D(order);

    /*有更好的实现方法吗？*/
    /*e_v[i]: face号面上的i号顶点在本单元的编号；
      neigh_v[i]: 该顶点在face面的相邻单元的编号*/
    for (i=0; i < NVert - 1; i++) {
        e_v[i] = GetFaceVertex(face, i);
        for (j=0; j < NVert; j++) {
            if (neigh_e->verts[j] == e->verts[e_v[i]]) {
                neigh_v[i] = j;
                break; 
            } 
        }
    }
    lambda[face] = 0.;
    neigh_lambda[neigh_face] = 0.;

    buffer = phgAlloc(sizeof(*buffer));
    p = quad->points;
    w = quad->weights;

    d = 0.;
    for (i = 0; i < quad->npoints; i++) {
        for (j = 0; j < NVert - 1; j++) {
        lambda[e_v[j]] = *(p++);
        neigh_lambda[neigh_v[j]] = lambda[e_v[j]];
        }
        
        dof = phgDofEval(v, neigh_e, neigh_lambda, buffer);
        bas = (FLOAT *)type->BasFuncs(u, e, N, N + 1, lambda);
        d += *(bas++) * *(dof++) * *(w++);
    } 
    phgFree(buffer);

    return d * phgGeomGetFaceArea(u->g, e, face);
}

/*tensor product*/
void
tensorProduct(struct TENSOR A, struct TENSOR B, struct TENSOR C)
{
    int i;
    int num_A = 1, num;

    if(m==0 && n==0){
        *C = phgDofMM(MAT_OP_N, MAT_OP_N, 1, 1, 1, 1., *A, 0, *B, 0., NULL);
    }
    else if(m>0){
        for(i=1; i<m; i++){
           num_A *= *(dims_A + i); 
        }
        
        num = num_A;
        for(i=0; i<n; i++){
           num *= *(dims_B + i); 
        }

        for(i=0; i<dims_A[0]; i++){
            tensorProduct(A + i*num_A, m-1, dims_A + 1, B, n, dims_B, C + i*num);
        }
    }
    else{
        tensorProduct(B, n, dims_B, A, m, dims_A, C);
    }
        
}

void
contract(DOF **A, int s, int *dims, int m, int n, DOF **B)
{   
    int i, j, k, l, num_1 = 1, num_2 = 1, num_3 = 1;
    GRID *g = (*A)->g;
    DOF_TYPE *dof_tp = (*A)->type;
    DOF **ptr_B = B, **ptr_A = A;

    for(i=0; i<m-1; i++)
        num_1 *= dims[i];

    for(i=m; i<n-1; i++)
        num_2 *= dims[i];

    for(i=n; i<s; i++)
        num_3 *= dims[i];

    for(i=0; i<num_1; i++)
        for(j=0; j<num_2; j++)
            for(k=0; k<num_3; k++){  
                *ptr_B = phgDofNew(g, dof_tp, 1, "B", func_zero); 
                for(l=0; l<dims[m-1]; l++){
                    ptr_A = A + i*dims[m-1]*num_2*dims[n-1]*num_3 + l*num_2*dims[m-1]*num_3 + j*dims[m-1]*num_3 + l*num_3 + k;
                    phgDofAXPY(1., *ptr_A, ptr_B); 
                }
                ptr_B++;
            }
}

void
productAndContract()
{
}

void
DofReciprocal(DOF *u_h)
{
    int i;
    int data_len = DofGetDataCount(u_h);
    FLOAT *p = DofData(u_h);
    for(i = 0; i < data_len; i++){
        *p = 1./(*p);
        p++;    
    }
}

void
DofSqrt(DOF *u_h)
{
    int i;
    int data_len = DofGetDataCount(u_h);
    FLOAT *p = DofData(u_h);
    for(i=0; i<data_len; i++){
        *p = sqrt(*p);
        p++;
    }
} 
    
void
computeShift(DOF *psi[4][4], DOF *inverse_metric[3][3], DOF *Shift[3])
{
}

void
computeLapse(DOF *psi[4][4], DOF *space_metric, DOF **Shift, DOF *Lapse)
{
}

void
inverseMetric(DOF *metric[][3], DOF *metric_inverse[][3])
{
    int i, j;
    DOF *inv_determinant;
    DOF *adj[3][3];
    
    /*compute Adj(M)*/ 
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            if(j < i) adj[i][j] = adj[j][i];
            else{            
                adj[i][j] = phgDofMM(MAT_OP_N, MAT_OP_N, 1, 1, 1, 1., metric[(i+1)%3][(j+1)%3], 0, metric[(i+2)%3][(j+2)%3], 0., NULL); 
                phgDofMM(MAT_OP_N, MAT_OP_N, 1, 1, 1, -1., metric[(i+1)%3][(j+2)%3], 0, metric[(i+2)%3][(j+1)%3], 1., &adj[i][j]); 
            }
        }
    }

    /*compute the reciprocal of the determinant of metric matrix*/
    inv_determinant = phgDofMM(MAT_OP_N, MAT_OP_N, 1, 1, 1, 1., metric[0][0], 0, adj[0][0], 0., NULL);  
    phgDofMM(MAT_OP_N, MAT_OP_N, 1, 1, 1, 1., metric[0][1], 0, adj[0][1], 1., &inv_determinant);  
    phgDofMM(MAT_OP_N, MAT_OP_N, 1, 1, 1, 1., metric[0][2], 0, adj[0][2], 1., &inv_determinant);  
    DofReciprocal(inv_determinant); 
   
    /*compute the inverse matrix of the metric using M^-1 = Adj(M)/Det(M).*/ 
    for(i =0; i<3; i++){
        for(j = 0; j<3; j++){
            if(j < i) metric_inverse[i][j] = metric_inverse[j][i];
            else{
                 metric_inverse[i][j] = phgDofMM(MAT_OP_N, MAT_OP_N, 1, 1, 1, 1., adj[i][j], 0, inv_determinant, 0., NULL);
            }
        }
    }

    /*release the temp DOFs*/
    phgDofFree(&inv_determinant);
    for(i=0;i<3;i++)
        for(j=0;j<=i;j++)
            phgDofFree(&adj[i][j]);
}

static void
compute_RHS_psi()
{
     
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


int 
main(int argc, char * argv[])
{
    char *fn = "cube.dat"; 
    int i, j, k;
    //const char *matrix_fn = "mat.m", *var_name = "mat_M";
    GRID *g; 
    DOF_TYPE *dof_tp = DOF_DG2;
    DOF *u_h, *p_0, *p_1, *q_0, *q_1, *w_0, *w_1;
    DOF *psi_00, *psi_01, *psi_02, *psi_03, *psi_11, *psi_12, *psi_13, *psi_22, *psi_23, *psi_33;
    DOF *Pi_00, *Pi_01, *Pi_02, *Pi_03, *Pi_11, *Pi_12, *Pi_13, *Pi_22, *Pi_23, *Pi_33;
    DOF *Phi_100, *Phi_101, *Phi_102, *Phi_103, *Phi_111, *Phi_112, *Phi_113, *Phi_122, *Phi_123, *Phi_133;
    DOF *Phi_200, *Phi_201, *Phi_202, *Phi_203, *Phi_211, *Phi_212, *Phi_213, *Phi_222, *Phi_223, *Phi_233;
    DOF *Phi_300, *Phi_301, *Phi_302, *Phi_303, *Phi_311, *Phi_312, *Phi_313, *Phi_322, *Phi_323, *Phi_333;

    DOF *Lapse;
    DOF *Shift[3];
    SOLVER *solver;
    //MAT *mat_M;

    phgInit(&argc, &argv);	
    
    g = phgNewGrid(-1); 
    phgImport(g, fn, FALSE);

    //phgRefineAllElements(g, 1);
   
    /*creat dofs for all the functions to be solved*/ 
    psi_00 = phgDofNew(g, dof_tp, 1, "psi_00", DofInterpolation);
    psi_01 = phgDofNew(g, dof_tp, 1, "psi_01", DofInterpolation);
    psi_02 = phgDofNew(g, dof_tp, 1, "psi_02", DofInterpolation);
    psi_03 = phgDofNew(g, dof_tp, 1, "psi_03", DofInterpolation);
    psi_11 = phgDofNew(g, dof_tp, 1, "psi_11", DofInterpolation);
    psi_12 = phgDofNew(g, dof_tp, 1, "psi_12", DofInterpolation);
    psi_13 = phgDofNew(g, dof_tp, 1, "psi_13", DofInterpolation);
    psi_22 = phgDofNew(g, dof_tp, 1, "psi_22", DofInterpolation);
    psi_23 = phgDofNew(g, dof_tp, 1, "psi_23", DofInterpolation);
    psi_33 = phgDofNew(g, dof_tp, 1, "psi_33", DofInterpolation);

    DOF *psi_components[4][4] = {{psi_00, psi_01, psi_02, psi_03}, 
                      {psi_01, psi_11, psi_12, psi_13},
                      {psi_02, psi_12, psi_22, psi_23},
                      {psi_03, psi_13, psi_23, psi_33}};

    struct TENSOR psi_tensor = {*psi_components, 2, dims_psi};

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
    
    DOF *Pi_components[4][4] = {{Pi_00, Pi_01, Pi_02, Pi_03}, 
                     {Pi_01, Pi_11, Pi_12, Pi_13},
                     {Pi_02, Pi_12, Pi_22, Pi_23},
                     {Pi_03, Pi_13, Pi_23, Pi_33}};
   
    struct TENSOR Pi_tensor = {*Pi_components, 2, dims_Pi};

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

    DOF *Phi[3][4][4] ={{{Phi_100, Phi_101, Phi_102, Phi_103}, 
                         {Phi_101, Phi_111, Phi_112, Phi_113},
                         {Phi_102, Phi_112, Phi_122, Phi_123},
                         {Phi_103, Phi_113, Phi_123, Phi_133}}, 
                        {{Phi_200, Phi_201, Phi_202, Phi_203},
                         {Phi_201, Phi_211, Phi_212, Phi_213},
                         {Phi_202, Phi_212, Phi_222, Phi_223},
                         {Phi_203, Phi_213, Phi_223, Phi_233}},
                        {{Phi_300, Phi_301, Phi_302, Phi_303},
                         {Phi_301, Phi_311, Phi_312, Phi_313},
                         {Phi_302, Phi_312, Phi_322, Phi_323},
                         {Phi_303, Phi_313, Phi_323, Phi_333}}};

    struct TENSOR Phi_tensor = {*Phi_components, 3, dims_Phi};

    u_h = phgDofNew(g, dof_tp, 1, "u_h", func_u);
    p_0 = phgDofNew(g, dof_tp, 1, "p_0", DofInterpolation); 
    p_1 = phgDofNew(g, dof_tp, 1, "p_1", DofInterpolation); 
    q_0 = phgDofNew(g, dof_tp, 1, "q_0", DofInterpolation); 
    q_1 = phgDofNew(g, dof_tp, 1, "q_1", DofInterpolation); 
    w_0 = phgDofNew(g, dof_tp, 1, "w_0", DofInterpolation); 
    w_1 = phgDofNew(g, dof_tp, 1, "w_1", DofInterpolation); 

      
    /*建立解法器，生成矩阵及右端项*/
    //solver = phgSolverCreate(SOLVER_DEFAULT, p_0, p_1, q_0, q_1, w_0, w_1, NULL);
  
    //build_linear_system(solver, u_h);
    
    /*将解法器中的矩阵输出到文件*/ 
    //mat_M = phgSolverGetMat(solver);
    //phgMatDumpMATLAB(mat_M, var_name, matrix_fn);

    /*解法器求解*/
    //phgSolverSolve(solver, TRUE, p_0, p_1, q_0, q_1, w_0, w_1, NULL);
    //phgSolverDestroy(&solver);

    /*输出结果*/

    //phgDofDump(w_1);

    /*Export VTK*/
    //phgExportVTK(g, "tmp.vtk", u_h, p_0, p_1, q_0, q_1, w_0, w_1, NULL);


    /*release the mem*/
    for(i=0;i<4;i++)
        for(j=0;j<=i;j++){
            phgDofFree(&psi_components[i][j]);
            phgDofFree(&Pi_components[i][j]);
            for(k=0;k<3;k++){
                phgDofFree(&Phi_components[k][i][j]);
            }
        }


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
