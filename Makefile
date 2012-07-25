
INSTALL      ?= ..

TSTDIR       ?= ../tests
LIBDIR       ?= $(INSTALL)/lib
BINDIR       ?= $(INSTALL)/bin
INCDIR       ?= $(INSTALL)/include

CC           ?= cc
AR           ?= ar
LDSHARED     ?= $(CC) -arch i386 -dynamiclib -undefined suppress -flat_namespace
ARSTATIC     ?= $(AR) rcu
FPIC         ?= -fPIC
COW_HDF5     ?= 0
COW_MPI      ?= 0
COW_HDF5_MPI ?= 0
COW_FFTW     ?= 0
CFLAGS       ?= -Wall -g -O0
FFTW_INC     ?= 
FFTW_LIB     ?= 
HDF5_INC     ?= 
HDF5_LIB     ?= 

DEFINES = \
	-DCOW_MPI=$(COW_MPI) \
	-DCOW_HDF5=$(COW_HDF5) \
	-DCOW_HDF5_MPI=$(COW_HDF5_MPI) \
	-DCOW_FFTW=$(COW_FFTW)

LIB = $(HDF5_LIB) $(FFTW_LIB)
INC = $(HDF5_INC) $(FFTW_INC)

OBJ = cow.o hist.o io.o samp.o srhdpack.o fft.o fft_3d.o pack_3d.o remap_3d.o
EXE = 	$(BINDIR)/mhdstats \
	$(BINDIR)/srhdhist \
	$(TSTDIR)/testcow \
	$(TSTDIR)/testhist \
	$(TSTDIR)/testfft \
	$(TSTDIR)/testsamp

LIBS = $(LIBDIR)/libcow.so $(LIBDIR)/libcow.a
HEADERS = $(INCDIR)/cow.h $(INCDIR)/srhdpack.h

default : all

exe : $(BINDIR) $(EXE)

lib : $(LIBDIR) $(LIBS)

headers : $(INCDIR) $(HEADERS)

all : exe lib headers

%.o : %.c
	$(CC) $(CFLAGS) -o $@ $< $(DEFINES) $(INC) $(FPIC) -c -std=c99

$(INCDIR)/%.h : %.h $(INCDIR)
	cp $< $@

$(TSTDIR) :
	@mkdir -p $@

$(LIBDIR) :
	@mkdir -p $@

$(BINDIR) :
	@mkdir -p $@

$(INCDIR) :
	@mkdir -p $@

$(LIBDIR)/libcow.so : $(OBJ)
	$(LDSHARED) -o $@ $? $(LIB)

$(LIBDIR)/libcow.a : $(OBJ)
	$(ARSTATIC) $@ $?

$(BINDIR)/mhdstats : mhdstats.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

$(BINDIR)/srhdhist : srhdhist.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

$(TSTDIR)/testcow : testcow.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

$(TSTDIR)/testhist : testhist.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

$(TSTDIR)/testfft : testfft.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

$(TSTDIR)/testsamp : testsamp.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

clean :
	@rm -rf $(EXE) *.o
