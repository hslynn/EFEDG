static void
get_dofs_diff(DOF **dofs_var, DOF **dofs_sol, DOF **dofs_diff)
{
    int i;
    copy_dofs(dofs_sol, dofs_diff, "diff", NVAR);
    for(i = 0;i<NVAR;i++){
        phgDofAXPBY(1.0, dofs_var[i], -1.0, dofs_diff + i);
    }
}
