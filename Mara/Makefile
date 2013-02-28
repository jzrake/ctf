
MAKEFILE_IN = Makefile.in
include $(MAKEFILE_IN)

RM ?= rm -f
AR ?= ar
ARSTATIC ?= $(AR) rcu
RANLIB ?= ranlib
CFLAGS ?= -Wall

SRC_C   = $(wildcard *.c)
SRC_CPP = $(wildcard *.cpp)
OBJ_C   = $(SRC_C:.c=.o)
OBJ_CPP = $(SRC_CPP:.cpp=.o)

MAR_A = libmara.a
LUA_I ?= -I$(LUA_HOME)/include

ifeq ($(strip $(USE_MPI)), 1)
DEFINES += -DUSE_MPI
endif

default : $(MAR_A) mara.o

%.o : %.c
	$(CC) $(CFLAGS) -c $^ $(DEFINES)

%.o : %.cpp
	$(CXX) $(CFLAGS) -c $^ $(DEFINES)

mara.o : mara.cpp
	$(CXX) $(CFLAGS) -c $< $(LUA_I) $(DEFINES)

$(MAR_A) : $(OBJ_C) $(OBJ_CPP)
	$(ARSTATIC) $@ $?
	$(RANLIB) $@

clean :
	$(RM) $(MAR_A) $(OBJ_C) $(OBJ_CPP) mara.o
