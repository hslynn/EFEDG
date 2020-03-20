void
phgQuadDofGradBas(ELEMENT *e, DOF *p, DOF *v, int m, int order, FLOAT *values);

void
set_neighbour_data(DOF *u, ELEMENT *e, INT s, DOF *neigh_u, NEIGHBOUR_DATA *nd);

void
reset_neighbour_data(DOF *neigh_u, ELEMENT *e, INT s);

void
build_linear_system(SOLVER *solver, DOF *u_h, NEIGHBOUR_DATA *nd, DOF *dof_bdry_in, DOF **dofs_bdry_out, int coord);

void
dgHJGradDof(DOF *u_h, NEIGHBOUR_DATA *nd, DOF *dof_bdry_in, DOF **dofs_bdry_out, DOF *dof_grad, int coord);

void
get_dofs_grad(DOF **dofs, DOF **dofs_bdry_in, DOF **dofs_bdry_out, DOF **dofs_grad, int ndof);

void
get_dof_grad_hat(DOF *dof_grad);

void
get_dofs_grad_hat(DOF **dofs, DOF **dofs_bdry_in, DOF **dofs_bdry_out, DOF **dofs_grad_hat, INT ndof);
