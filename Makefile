default: EFEDG 
clean:
	/bin/rm -f *.o *.txt *.vtk Einstein_DG *.err *.out
objects = EFEDG_main.o hdw.o auxi_dofs.o grad.o Hhat.o initial_condition.o rhs.o rk.o source_values.o\
		  filter.o runtime.o
EFEDG: $(objects)
	$(LINKER) $(LDFLAGS) -o EFEDG $(objects) $(LIBS)

EFEDG_main.o: global_def.h initial_condition.h hdw.h rk.h rhs.h filter.h runtime.h
auxi_dofs.o: global_def.h source_values.h 
grad.o: global_def.h filter.h
Hhat.o: global_def.h hdw.h grad.h
initial_condition.o: global_def.h Sch_Harmonic_hp.c Sch_Harmonic_non_hp.c Sch_Kerr_Schild_outgoing.c\
	Sch_Kerr_Schild_ingoing.c Minkovski.c
rhs.o: global_def.h Hhat.h auxi_dofs.h hdw.h
rk.o: global_def.h rhs.h hdw.h
source_values.o: global_def.h
filter.o: global_def.h
runtime.o: global_def.h

include /home/hslynn/local/share/phg/Makefile.inc
#mvapich2-2.3b_gcc
#include /soft/phg/phg-0.9.4-mvapich2-20180902/share/phg/Makefile.inc

#mvapich2-2.0.1_icc
#include /soft/phg/phg-0.9.4-mvapich2-20180317/share/phg/Makefile.inc

#openmpi
#include /soft/phg/phg-0.9.4-openmpi-20180317/share/phg/Makefile.inc
