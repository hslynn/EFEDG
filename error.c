static void
get_dofs_diff(DOF **dofs_A, DOF **dofs_B, DOF **dofs_diff, INT ndofs)
{
    short i;
    copy_dofs(dofs_B, dofs_diff, "diff", ndofs);
    for(i = 0;i<ndofs;i++){
        phgDofAXPBY(1.0, dofs_A[i], -1.0, dofs_diff + i);
    }
}
