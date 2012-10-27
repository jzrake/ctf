
TSTDIR       ?= ../tests
EXMDIR       ?= ../examples
CC           ?= cc
CFLAGS       ?= -Wall -O3

OBJ = fish.o reconstruct.o fluids.o riemann.o matrix.o
EXE = $(TSTDIR)/testfish $(TSTDIR)/testfluids $(EXMDIR)/euler

default : $(EXE)

%.o : %.c
	$(CC) $(CFLAGS) -o $@ $< -c -std=c99

$(EXMDIR)/euler : euler.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

$(TSTDIR)/testfish : testfish.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

$(TSTDIR)/testfluids : testfluids.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

clean :
	@rm -rf $(EXE) *.o
