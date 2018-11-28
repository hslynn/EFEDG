#include "grad.c"
#include "Hhat_values.c"

static void
get_dofs_Hhat(DOF **dofs_var, DOF **dofs_g, DOF **dofs_N, 
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF **dofs_Hhat)
{
    get_dofs_grad_hat(dofs_var, dofs_Psi, 10);
    get_dofs_grad_hat(dofs_var + 10, dofs_Pi, 10);
    get_dofs_grad_hat(dofs_var + 20, dofs_Phi, 30);
    
    INT i, j, n, data_count = DofGetDataCount(dofs_N[0]);
    FLOAT *p_gradPsi[30], *p_gradPi[30], *p_gradPhi[90], *p_g[6], *p_N[4], *p_Hhat[50];
    FLOAT v_gradPsi[30][2], v_gradPi[30][2], v_gradPhi[90][2];
    FLOAT v_g[6], v_N[4], v_Hhat[50];
    
    for(i=0; i<30; i++){
        p_gradPsi[i] = DofData(dofs_gradPsi[i]);
        p_gradPi[i] = DofData(dofs_gradPi[i]);
        p_gradPhi[i] = DofData(dofs_gradPhi[i]);
        p_gradPhi[30+i] = DofData(dofs_gradPhi[30+i]);
        p_gradPhi[60+i] = DofData(dofs_gradPhi[60+i]);
    }    

    for(i=0; i<50; i++){
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
        for(i=0; i<30; i++){
            for(j=0; j<2; j++){
                values_gradPsi[i][j] = *(p_gradPsi[i]++);
                values_gradPi[i][j] = *(p_gradPi[i]++);
                values_gradPhi[i][j] = *(p_gradPhi[i]++);
                values_gradPhi[30+i][j] = *(p_gradPhi[30+i]++);
                values_gradPhi[60+i][j] = *(p_gradPhi[60+i]++);
            }
        }
        for(i=0; i<6; i++){
            values_g[i] = *(p_g[i]++);
        }
        for(i=0; i<4; i++){
            values_N[i] = *(p_N[i]++);
        }
        
        get_values_Hhat(values_gradPsi, values_gradPi, values_gradPhi,  
                values_g, values_N, values_Hhat);

        for(i=0; i<50; i++){
            *(p_Hhat[i]++) = values_Hhat[i];
        }
    }
        
}
