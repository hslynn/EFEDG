default: v1 
all: default
clean:
	/bin/rm -f *.o v1
v1: v1.o 


include /home/hslynn/apps/phg/share/phg/Makefile.inc
