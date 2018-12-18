default: Einstein_DG 
all: default
clean:
	/bin/rm -f *.o *.txt *.vtk Einstein_DG
Einstein_DG.o: Einstein_DG.c 


#include /home/hslynn/apps/phg/share/phg/Makefile.inc
include /soft/phg/phg-0.9.4-mvapich2-20180902/share/phg/Makefile.inc
