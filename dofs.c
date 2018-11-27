static void
create_dofs(GRID *g, DOF_TYPE *type, INT dim, DOF **dofs_list, char **name_list, INT ndof)
{
    for(INT i=0; i<ndof; i++){
        dofs_list[i] = phgDofNew(g, type, dim, name_list[i], DofInterpolation);
    } 
}
