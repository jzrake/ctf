
CFLAGS = -Wall -O2

default : pdf-machine polytrope

pdf-machine : pdf-machine.c
	$(CC) $(CFLAGS) -o $@ $^
	cp $@ ../bin

polytrope : polytrope.c
	$(CC) $(CFLAGS) -o $@ $^
	cp $@ ../bin
