
include ../Makefile.in

CFLAGS = -Wall -O2

GL = -framework OpenGL -framework GLUT -I/usr/X11R6/include
H5 = -I$(HDF_HOME)/include -L$(HDF_HOME)/lib -lhdf5

default : pdf-machine polytrope h5view

pdf-machine : pdf-machine.c
	$(CC) $(CFLAGS) -o $@ $^
	cp $@ ../bin

polytrope : polytrope.c
	$(CC) $(CFLAGS) -o $@ $^
	cp $@ ../bin

h5view : h5view.c
	$(CC) $(CFLAGS) -o $@ $^ $(GL) $(H5)
	cp $@ ../bin
