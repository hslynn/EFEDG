default: Einstein_DG 
clean:
	/bin/rm -f *.o *.txt *.vtk Einstein_DG *.err *.out
objects = Einstein_DG.o hdw.o auxi_dofs.o grad.o Hhat.o initial_condition.o rhs.o rk.o source_values.o filter.o runtime.o
Einstein_DG: $(objects)
	${LINKER} ${LDFLAGS} -o Einstein_DG $(objects) ${LIBS}

include /home/hslynn/local/share/phg/Makefile.inc
#mvapich2-2.3b_gcc
#include /soft/phg/phg-0.9.4-mvapich2-20180902/share/phg/Makefile.inc

#mvapich2-2.0.1_icc
#include /soft/phg/phg-0.9.4-mvapich2-20180317/share/phg/Makefile.inc

#openmpi
#include /soft/phg/phg-0.9.4-openmpi-20180317/share/phg/Makefile.inc
