void
rk2(FLOAT dt, DOF **dofs_var, DOF **dofs_bdry, DOF **dofs_g, DOF **dofs_N, DOF **dofs_H, DOF **dofs_deriH, DOF **dofs_src, DOF **dofs_C,
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF **dofs_Hhat, DOF **dofs_rhs);

void
rk3(FLOAT dt, DOF **dofs_var, DOF **dofs_bdry, DOF **dofs_g, DOF **dofs_N, DOF **dofs_H, DOF **dofs_deriH, DOF **dofs_src, DOF **dofs_C,
        DOF **dofs_gradPsi, DOF **dofs_gradPi, DOF **dofs_gradPhi, DOF **dofs_Hhat, DOF **dofs_rhs);
