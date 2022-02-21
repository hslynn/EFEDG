#include "global_def.h"

#define spacetime 5
#if (spacetime == 1)
#include "Sch_Kerr_Schild_ingoing.c"
#elif (spacetime == 2)
#include "Sch_Kerr_Schild_outgoing.c"
#elif (spacetime == 3)
#include "Sch_Harmonic_non_hp.c"
#elif (spacetime == 4)
#include "Sch_Harmonic_hp.c"
#elif (spacetime == 5)
#include "Kerr_Kerr_Schild_ingoing.c"
#endif

INT SPACETIME = spacetime;

void
set_data_var(DOF **dofs_var)
{
    int i;
    DOF_USER_FUNC funcs_Psi[10] = {func_Psi00, func_Psi01, func_Psi02, func_Psi03,
                                                    func_Psi11, func_Psi12, func_Psi13,
                                                                func_Psi22, func_Psi23,
                                                                            func_Psi33};
    DOF_USER_FUNC funcs_Pi[10] = {func_Pi00, func_Pi01, func_Pi02, func_Pi03,
                                             func_Pi11, func_Pi12, func_Pi13,
                                                        func_Pi22, func_Pi23,
                                                                   func_Pi33};
    DOF_USER_FUNC funcs_Phi[30] = {func_Phi100, func_Phi101, func_Phi102, func_Phi103,        
                                               func_Phi111, func_Phi112, func_Phi113,
                                                          func_Phi122, func_Phi123,
                                                                      func_Phi133,
                                   func_Phi200, func_Phi201, func_Phi202, func_Phi203,        
                                               func_Phi211, func_Phi212, func_Phi213,
                                                          func_Phi222, func_Phi223,
                                                                      func_Phi233,
                                   func_Phi300, func_Phi301, func_Phi302, func_Phi303,        
                                               func_Phi311, func_Phi312, func_Phi313,
                                                          func_Phi322, func_Phi323,
                                                                      func_Phi333};
    for(i=0; i<10; i++){
        phgDofSetDataByFunction(dofs_var[i], funcs_Psi[i]);
        phgDofSetDataByFunction(dofs_var[10 + i], funcs_Pi[i]);
        
        phgDofSetDataByFunction(dofs_var[20 + i], funcs_Phi[i]);
        phgDofSetDataByFunction(dofs_var[30 + i], funcs_Phi[10 + i]);
        phgDofSetDataByFunction(dofs_var[40 + i], funcs_Phi[20 + i]);
    }
}

void
set_data_H(DOF **dofs_H)
{
    int i;
    DOF_USER_FUNC funcs_H[4] = {func_H0, func_H1, func_H2, func_H3};
    for(i=0;i<4;i++){
        phgDofSetDataByFunction(dofs_H[i], funcs_H[i]);
    }
}

void
set_data_deriH(DOF **dofs_deriH)
{
    int i;
    DOF_USER_FUNC funcs_deriH[16] = {func_deriH00, func_deriH01, func_deriH02, func_deriH03,
                                     func_deriH10, func_deriH11, func_deriH12, func_deriH13,
                                     func_deriH20, func_deriH21, func_deriH22, func_deriH23,
                                     func_deriH30, func_deriH31, func_deriH32, func_deriH33};
    for(i=0;i<16;i++){
        phgDofSetDataByFunction(dofs_deriH[i], funcs_deriH[i]);
    }
}
