#include "grad.c"
static void
get_values_Hhat(FLOAT values_gradPsi[][2], FLOAT values_gradPi[][2], FLOAT values_gradPhi[][2], FLOAT *values_g, FLOAT *values_N, FLOAT *values_Hhat)
{
    FLOAT list_gradPsi_ave[30], list_gradPi_ave[30], list_gradPhi_ave[90];
    FLOAT list_gradPsi_diff[30], list_gradPi_diff[30], list_gradPhi_diff[90];
    FLOAT tensor_gradPsi_ave[3][4][4], tensor_gradPi_ave[3][4][4], tensor_gradPhi_ave[3][3][4][4];
    FLOAT tensor_gradPsi_diff[3][4][4], tensor_gradPi_diff[3][4][4], tensor_gradPhi_diff[3][3][4][4];
    FLOAT tensor_g[3][3];
    INT i, j, a, b;
    for(i=0;i<30;i++){
        list_gradPsi_ave[i] = values_gradPsi[i][0];
        list_gradPsi_diff[i] = values_gradPsi[i][1];
        list_gradPi_ave[i] = values_gradPi[i][0];
        list_gradPi_diff[i] = values_gradPi[i][1];
    }
    for(i=0;i<90;i++){
        list_gradPhi_ave[i] = values_gradPhi[i][0];
        list_gradPhi_diff[i] = values_gradPhi[i][1];
    }

    for(i=0;i<3;i++){
        list2tensor(list_gradPsi_ave + i*10, tensor_gradPsi_ave[i][0], 4);
        list2tensor(list_gradPi_ave + i*10, tensor_gradPi_ave[i][0], 4);
        list2tensor(list_gradPsi_diff + i*10, tensor_gradPsi_diff[i][0], 4);
        list2tensor(list_gradPi_diff + i*10, tensor_gradPi_diff[i][0], 4);
        for(j=0;i<3;i++){
            list2tensor(list_gradPhi_ave + i*j*10, tensor_gradPhi_ave[i][j][0], 4);
            list2tensor(list_gradPhi_diff + i*j*10, tensor_gradPhi_diff[i][j][0], 4);
        }
    }
    list2tensor(values_g, tensor_g[0], 3);
    
    
    FLOAT *p_Psi = values_Hhat, *p_Pi = values_Hhat + 10; 
    FLOAT *p_xPhi = values_Hhat + 20, *p_yPhi = values_Hhat + 30, *p_zPhi = values_Hhat + 40;
    for(a=0;a<4;a++){
        for(b=a;b<4;b++){
            for(i=0;i<3;i++){
                *p_Psi += -(1+gamma1)*values_N[i+1]*tensor_gradPsi_ave[i][a][b];
                *p_Pi += -values_N[i+1]*tensor_gradPi_ave[i][a][b] 
                        - gamma1*gamma2*values_N[i+1]*tensor_gradPi_ave[i][a][b];
                *p_xPhi += -values_N[i+1]*tensor_gradPhi_ave[i][0][a][b];
                *p_yPhi += -values_N[i+1]*tensor_gradPhi_ave[i][1][a][b];
                *p_zPhi += -values_N[i+1]*tensor_gradPhi_ave[i][2][a][b];
                for(j=0;j<3;j++){
                    *p_Pi += values_N[0]*tensor_g[i][j]*tensor_gradPhi_ave[i][j][a][b];
                }
            }
            *(p_xPhi++) = values_N[0]*tensor_gradPi_ave[0][a][b] - values_N[0]*gamma2*tensor_gradPsi_ave[0][a][b];
            *(p_yPhi++) = values_N[0]*tensor_gradPi_ave[1][a][b] - values_N[0]*gamma2*tensor_gradPsi_ave[1][a][b];
            *(p_zPhi++) = values_N[0]*tensor_gradPi_ave[2][a][b] - values_N[0]*gamma2*tensor_gradPsi_ave[2][a][b];
            p_Psi++;
            p_Pi++;
        
    printf("test1\n\n"); 
        }
    }
       
}

static void
get_dofs_Hhat(DOF **dofs_var, DOF **dofs_g, DOF **dofs_N, 
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF **dofs_Hhat)
{
    get_dofs_grad_hat(dofs_var, dofs_gradPsi, 10);
    get_dofs_grad_hat(dofs_var + 10, dofs_gradPi, 10);
    get_dofs_grad_hat(dofs_var + 20, dofs_gradPhi, 30);

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
                v_gradPsi[i][j] = *(p_gradPsi[i]++);
                v_gradPi[i][j] = *(p_gradPi[i]++);
                v_gradPhi[i][j] = *(p_gradPhi[i]++);
                v_gradPhi[30+i][j] = *(p_gradPhi[30+i]++);
                v_gradPhi[60+i][j] = *(p_gradPhi[60+i]++);
            }
        }
        for(i=0; i<6; i++){
            v_g[i] = *(p_g[i]++);
        }
        for(i=0; i<4; i++){
            v_N[i] = *(p_N[i]++);
        }
         
        get_values_Hhat(v_gradPsi, v_gradPi, v_gradPhi,  
                v_g, v_N, v_Hhat);
        
        printf("test0\n\n"); 
        for(i=0; i<50; i++){
            *(p_Hhat[i]++) = v_Hhat[i];
        }
    }
        
}




