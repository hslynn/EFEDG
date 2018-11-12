default: Einstein_DG 
all: default
clean:
	/bin/rm -f *.o Einstein_DG
Einstein_DG: Einstein_DG 


include /home/hslynn/apps/phg/share/phg/Makefile.inc
