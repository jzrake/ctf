
MAKEFILE_IN = Makefile.in
include $(MAKEFILE_IN)

RM ?= rm -f
AR ?= ar
ARSTATIC ?= $(AR) rcu
CFLAGS ?= -Wall
CFLAGS += -std=c99

FISH_A = libfish.a
LUA_I ?= -I$(LUA_HOME)/include

ifeq ($(strip $(USE_FFTW)), 1)
INC += -I$(FFT_HOME)/include
endif

OBJ = fish.o reconstruct.o fluids.o riemann.o matrix.o grav1d.o

default : $(FISH_A) lua-fluids.o lua-fish.o

%.o : %.c
	$(CC) $(CFLAGS) -c $^ $(LUA_I)

grav1d.o : grav1d.c
	$(CC) $(CFLAGS) -c $^ $(LUA_I) $(INC)

$(FISH_A) : $(OBJ)
	$(ARSTATIC) $@ $?

fluidsfuncs.c : fluids.h parse_fluids.py
	python parse_fluids.py

fishfuncs.c : fish.h parse_fish.py
	python parse_fish.py

lua-fluids.o : lua-fluids.c fluidsfuncs.c
	$(CC) $(CFLAGS) -c $< $(LUA_I)

lua-fish.o : lua-fish.c fishfuncs.c
	$(CC) $(CFLAGS) -c $< $(LUA_I)

clean :
	$(RM) $(FISH_A) $(OBJ) fluidsfuncs.c lua-fluids.o fishfuncs.c lua-fish.o
