#include <stdlib.h>
#include <time.h>
#include "global_def.h"

void
set_perturb(DOF **dofs_var, FLOAT amplitude)
{
    INT i, n, np, idx;
    GRID *g = dofs_var[0]->g; 
    ELEMENT *e;
    FLOAT *p_var;
    FLOAT perturb;

    np = dofs_var[0]->type->np_elem;
    srand((unsigned)time(NULL));
    ForAllElements(g, e){ 
        idx = e->index; 
        for(i=0; i<NVAR; i++){
            p_var = DofElementData(dofs_var[i], idx);
            for(n=0;n<np;n++){
                perturb = -amplitude + (rand()/(double) RAND_MAX)*(2*amplitude);
                *p_var += perturb;
                p_var++;
            }
        }    
    }
}

