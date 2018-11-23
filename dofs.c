static void
create_dofs(GRID *g, DOF_TYPE *type, DOF **dofs_list, char **name_list, INT ndof)
{
    for(INT i=0; i<ndof; i++){
        dofs_list[i] = phgDofNew(g, type, 1, name_list[i], DofInterpolation);
    } 
}
