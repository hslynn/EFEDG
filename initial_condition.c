#include "global_def.h"

#if (spacetime == 1)
#include "Sch_Kerr_Schild_ingoing.c"
#elif (spacetime == 2)
#include "Sch_Kerr_Schild_outgoing.c"
#elif (spacetime == 3)
#include "Sch_Harmonic_non_hp.c"
#elif (spacetime == 4)
#include "Sch_Harmonic_hp.c"
#endif

