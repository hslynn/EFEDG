default: v2 
all: default
clean:
	/bin/rm -f *.o v2
v2: v2.o 


include /home/hslynn/apps/phg/share/phg/Makefile.inc
