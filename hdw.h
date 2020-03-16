void
create_dofs(GRID *g, DOF_TYPE *type, INT dim, DOF **dofs_list, char *name_head, INT ndof);

void
copy_dofs(DOF **dofs_A, DOF **dofs_B, char *name_head, INT ndof);

void
free_dofs(DOF **dofs, INT ndof);

void
list2tensor(FLOAT *list, FLOAT *tensor, int dim);

void
get_dofs_diff(DOF **dofs_A, DOF **dofs_B, DOF **dofs_diff, INT ndofs);
