
# ------------------------------------------------------------------------------
# CTF build instructions
# ------------------------------------------------------------------------------
# 
# Copy the `Makefile.in.template` file to `Makefile.in`, and modify it for your
# system.
#
# 1. `git submodule update --init`
# 2. `make lua`
# 3. `make`
#
# ------------------------------------------------------------------------------

MAKEFILE_IN = $(PWD)/Makefile.in
include $(MAKEFILE_IN)

CFLAGS ?= -Wall
CURL   ?= curl
UNTAR  ?= tar -xvf
CD     ?= cd
RM     ?= rm -f
OS     ?= generic
LVER   ?= lua-5.2.1

LUA_I = -I$(LUA_HOME)/include
LUA_L = -L$(LUA_HOME)/lib -llua
LUA_A = $(LUA_HOME)/lib/liblua.a

LUA_COW    = cow/lua-cow.o
LUA_MARA   = Mara/mara.o
LUA_MPI    = lua-mpi/lua-mpi.o
LUA_HDF5   = lua-hdf5/lua-hdf5.o
LUA_BUFFER = lua-buffer/lua-buffer.o
LUA_GLUT   = lua-glut

MODULES = $(LUA_BUFFER)
LOCLIBS = $(LUA_A)

ifeq ($(strip $(USE_MPI)), 1)
MODULES += $(LUA_MPI)
DEFINES += -DUSE_MPI
endif

ifeq ($(strip $(USE_HDF5)), 1)
MODULES += $(LUA_HDF5)
DEFINES += -DUSE_HDF5
LIBS += -L$(HDF_HOME)/lib -lz -lhdf5
endif

ifeq ($(strip $(USE_MARA)), 1)
MODULES += $(LUA_MARA)
DEFINES += -DUSE_MARA
LOCLIBS += Mara/libmara.a
endif

ifeq ($(strip $(USE_COW)), 1)
MODULES += $(LUA_COW)
DEFINES += -DUSE_COW
LOCLIBS += cow/libcow.a
endif

ifeq ($(strip $(USE_FFTW)), 1)
LIBS += -L$(FFT_HOME)/lib -lfftw3
endif


default : main

$(LUA_A) :
	$(CD) $(LVER); $(MAKE) $(OS) CC=$(CC); $(MAKE) install INSTALL_TOP=$(PWD)/$(LVER)

cow/libcow.a : $(LUA_A) .FORCE
	$(MAKE) -C cow libcow.a MAKEFILE_IN=$(MAKEFILE_IN)

$(LUA_COW) : $(LUA_A) .FORCE
	$(MAKE) -C cow lua-cow.o MAKEFILE_IN=$(MAKEFILE_IN)

Mara/libmara.a : $(LUA_A) .FORCE
	$(MAKE) -C Mara libmara.a MAKEFILE_IN=$(MAKEFILE_IN)

$(LUA_MARA) : $(LUA_A) .FORCE
	$(MAKE) -C Mara mara.o MAKEFILE_IN=$(MAKEFILE_IN)

$(LUA_MPI) : $(LUA_A) .FORCE
	$(MAKE) -C lua-mpi lua-mpi.o MAKEFILE_IN=$(MAKEFILE_IN)

$(LUA_HDF5) : $(LUA_A) .FORCE
	$(MAKE) -C lua-hdf5 lua-hdf5.o MAKEFILE_IN=$(MAKEFILE_IN)

$(LUA_BUFFER) : $(LUA_A) .FORCE
	$(MAKE) -C lua-buffer lua-buffer.o MAKEFILE_IN=$(MAKEFILE_IN)

main.o : main.c $(LUA_A)
	$(CC) $(CFLAGS) -c -o $@ $< $(LUA_I) $(DEFINES) -DINSTALL_DIR=\"$(PWD)\"

main : main.o $(MODULES) $(LOCLIBS)
	$(CXX) $(CFLAGS) -o $@ $^ $(LIBS)

$(LUA_GLUT) :
	$(MAKE) -C lua-glut DEFS=$(LUA_I)

clean :
	$(MAKE) -C cow clean MAKEFILE_IN=$(MAKEFILE_IN)
	$(MAKE) -C Mara clean MAKEFILE_IN=$(MAKEFILE_IN)
	$(MAKE) -C lua-buffer clean MAKEFILE_IN=$(MAKEFILE_IN)
	$(MAKE) -C lua-mpi clean MAKEFILE_IN=$(MAKEFILE_IN)
	$(MAKE) -C lua-hdf5 clean MAKEFILE_IN=$(MAKEFILE_IN)
	$(MAKE) -C lua-glut clean
	$(MAKE) -C $(LVER) clean
	$(RM) *.o main

show :
	@echo $(MODULES)
	@echo $(LOCLIBS)
	@echo $(DEFINES)

# Also remove local Lua sources
realclean : clean
	$(RM) -r $(LVER)

.PHONY : lua-glut

.FORCE : 
