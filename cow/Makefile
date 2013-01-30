
MAKEFILE_IN = Makefile.in
include $(MAKEFILE_IN)

RM ?= rm -f
AR ?= ar
ARSTATIC ?= $(AR) rcu
CFLAGS ?= -Wall
CFLAGS += -std=c99

COW_A = libcow.a
LUA_I ?= -I$(LUA_HOME)/include
HDF_I ?= -I$(HDF_HOME)/include
FFT_I ?= -I$(FFT_HOME)/include

ifeq ($(strip $(USE_FFTW)), 1)
INCLUDE += $(FFT_I)
DEFINES += -DUSE_FFTW
endif

INC = $(HDF_I) $(FFT_I)
OPT = \
	-DCOW_MPI=$(USE_MPI) \
	-DCOW_HDF5=$(USE_HDF5) \
	-DCOW_HDF5_MPI=$(USE_MPIO) \
	-DCOW_FFTW=$(USE_FFTW)

OBJ = cow.o samp.o hist.o io.o fft.o fft_3d.o pack_3d.o remap_3d.o

default : $(COW_A) lua-cow.o

%.o : %.c
	$(CC) $(CFLAGS) -c $^ $(INC) $(OPT)

$(COW_A) : $(OBJ)
	$(ARSTATIC) $@ $?

cowfuncs.c : cow.h parse.py
	python parse.py

lua-cow.o : lua-cow.c cowfuncs.c
	$(CC) $(CFLAGS) -c $< $(LUA_I)

clean :
	$(RM) $(COW_A) $(OBJ) cowfuncs.c lua-cow.o
