#include "global_def.h"

void
get_values_Hhat(FLOAT values_gradPsi[][2], FLOAT values_gradPi[][2], FLOAT values_gradPhi[][2], 
        FLOAT *values_g, FLOAT *values_N, FLOAT *values_Hhat);

void
get_dofs_Hhat(DOF **dofs_var, DOF **dofs_bdry_in, DOF **dofs_bdry_out, DOF **dofs_g, DOF **dofs_N, 
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF **dofs_Hhat);
