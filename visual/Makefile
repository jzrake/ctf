
MAKEFILE_IN = Makefile.in
include $(MAKEFILE_IN)

RM ?= rm -f
AR ?= ar
ARSTATIC ?= $(AR) rcu

VIS_A = libvisual.a
LUA_I ?= -I$(LUA_HOME)/include

default : $(VIS_A) lic.o

lic.o : lic.cpp
	$(CXX) $(CFLAGS) -c $<

lua-visual.o : lua-visual.c
	$(CC) $(CFLAGS) -c $< $(LUA_I)

$(VIS_A) : lic.o
	$(ARSTATIC) $@ $?

clean :
	$(RM) $(VIS_A) $(OBJ_C) $(OBJ_CPP) lic.o
