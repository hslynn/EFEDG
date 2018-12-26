default: Einstein_DG 
all: default
clean:
	/bin/rm -f *.o *.txt *.vtk Einstein_DG *.err *.out
Einstein_DG.o: Einstein_DG.c 


#include /home/hslynn/apps/phg/share/phg/Makefile.inc

#mvapich2-2.3b_gcc
#include /soft/phg/phg-0.9.4-mvapich2-20180902/share/phg/Makefile.inc

#mvapich2-2.0.1_icc
#include /soft/phg/phg-0.9.4-mvapich2-20180317/share/phg/Makefile.inc

#openmpi
include /soft/phg/phg-0.9.4-openmpi-20180317/share/phg/Makefile.inc
