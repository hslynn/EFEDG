#include "grad.c"
#include "Hhat_values.c"

DOF *get_dof_grad_hat(DOF *dof_grad)
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
    return dof_grad;
}

static void
get_dofs_grad_hat(DOF **dofs_grad, DOF **dofs_grad_hat, INT ndof) 
{
   for(INT i=0;i<ndof;i++){
        dofs_grad_hat[i] = get_dof_grad_hat(dofs_grad[i]);
    } 
}

static void
get_dofs_Hhat(DOF **dofs_grad_var, DOF **dofs_g, DOF **dofs_N, DOF **dofs_Hhat)
{
    DOF *dofs_grad_var_hat[50]; 
    get_dofs_grad_hat(dofs_grad_var, dofs_grad_var_hat, 50);
    
    INT i, j, n, data_count = DofGetDataCount(dofs_N[0]);
    FLOAT *p_grad_hat[50], *p_g[6], *p_N[4], *p_Hhat[50];
    FLOAT values_grad_hat[50][6], values_g[6], values_N[4], values_Hhat[50];
    
    for(i=0; i<50; i++){
        p_grad_hat[i] = DofData(dofs_grad_var_hat[i]);
        p_Hhat[i] = DofData(dofs_Hhat[i]);
    }    
    for(i=0; i<6; i++){
        p_g[i] = DofData(dofs_g[i]);
    }
    for(i=0;i<4;i++){
        p_N[i] = DofData(dofs_N[i]);
    }
   
    // 
    for(n=0;n<data_count;n++){
        for(i=0; i<50; i++){
            for(j=0; j<6; j++){
                values_grad_hat[i][j] = *(p_grad_hat[i]++);
            }
        }
        for(i=0; i<6; i++){
            values_g[i] = *(p_g[i]++);
        }
        for(i=0; i<4; i++){
            values_N[i] = *(p_N[i]++);
        }
        
        get_values_Hhat(values_grad_hat, values_g, values_N, values_Hhat);

        for(i=0; i<50; i++){
            *(p_Hhat[i]++) = values_Hhat[i];
        }
    }
        
}
