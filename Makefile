default: Einstein_DG 
all: default
clean:
	/bin/rm -f *.o Einstein_DG
Einstein_DG.o: Einstein_DG.c 


include /home/hslynn/Apps/phg/share/phg/Makefile.inc
