#include "phg.h"

#include <string.h>
#include <math.h>


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
