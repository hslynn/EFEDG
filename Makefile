default: einstein 
all: default
clean:
	/bin/rm -f *.o einstein
u_x: einstein.o 


include /home/hslynn/apps/phg/share/phg/Makefile.inc
