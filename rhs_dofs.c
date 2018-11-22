#include "rhs_values.c"

static void
get_rhs_dof(DOF **dofs_var, DOF **dofs_rhs)
{

    int i, n, data_count;
    FLOAT *p_var[50], *p_rhs[50], values_var[50], values_rhs[50];
    for(int i=0; i < 50; i++){
        p_var[i] = DofData(dofs_var[i]);
        p_rhs[i] = DofData(dofs_rhs[i]); 
    }
    

    //evaluate dofs values at every data point
    data_count = DofGetDataCount(dofs_var[0]);
    for(n=0;n<data_count;n++){
        //compute rhs values at data point
        for(i=0; i < 50; i++){
            values_var[i] = *(p_var[i]);
        }
            
        rhs_values(values_var, values_rhs); 
           
        for(i=0; i<50; i++){
            *p_rhs[i] = values_rhs[i];
            p_rhs[i]++;
            p_var[i]++; 
        }
    }
}

