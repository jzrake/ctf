
MAKEFILE_IN = Makefile.in
include $(MAKEFILE_IN)

RM ?= rm -f
AR ?= ar
ARSTATIC ?= $(AR) rcu
CFLAGS ?= -Wall
CFLAGS += -std=c99

FLUIDS_A = libfluids.a
FISH_A = libfish.a
LUA_I ?= -I$(LUA_HOME)/include

OBJ = fish.o reconstruct.o fluids.o riemann.o matrix.o

default : $(FLUIDS_A) lua-fluids.o $(FISH_A) lua-fish.o

%.o : %.c
	$(CC) $(CFLAGS) -c $^ $(INC)

$(FLUIDS_A) : $(OBJ)
	$(ARSTATIC) $@ $?

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
	$(RM) $(FLUIDS_A) $(OBJ) fluidsfuncs.c lua-fluids.o fishfuncs.c lua-fish.o
