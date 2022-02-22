#include "phg.h"
#include "global_def.h"

void
filter(DOF **dofs, DOF_TYPE *dg_filter, INT ndof)
{
    FLOAT alpha = 0.1;
    DOF *dof_copy = phgDofNew(dofs[0]->g, dg_filter, dofs[0]->dim, "dof_copy", DofInterpolation);
    for (INT i = 0; i < ndof; i++)
    {
        phgDofCopy(dofs[i], &dof_copy, dg_filter, "dof_copy");
        phgDofAXPBY(alpha, dof_copy, 1.0-alpha, dofs+i); 
    }
    phgDofFree(&dof_copy);
}
