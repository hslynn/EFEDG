/***************************************************************/
/*  			      filter                           */
/***************************************************************/
#include "phg.h"
#include "global_def.h"

void
filter(DOF **dofs_var, DOF_TYPE *dg_type)
{
    FLOAT alpha = 0.1;
    DOF *dof_copy = phgDofNew(dofs_var[0]->g, dg_type, 1, "dof_copy", DofInterpolation);
    for (INT i = 0; i < NVAR; i++)
    {
        phgDofCopy(dofs_var[i], &dof_copy, dg_type, "dof_copy");
        phgDofAXPBY(alpha, dof_copy, 1.0-alpha, dofs_var+i); 
    }
    phgDofFree(&dof_copy);
}
