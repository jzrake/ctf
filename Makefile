
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
CFLAGS       ?= -Wall -O3

OBJ = fish.o weno.o
EXE = $(BINDIR)/testfish # $(BINDIR)/euler

LIBS = $(LIBDIR)/libfish.so $(LIBDIR)/libfish.a
HEADERS = $(INCDIR)/fish.h

COWLIBDIR = $(HOME)/Work/cow/lib
COWINCDIR = $(HOME)/Work/cow/include
FLUIDSLIBDIR = $(HOME)/Work/fluids/lib
FLUIDSINCDIR = $(HOME)/Work/fluids/include

LIB = -L$(COWLIBDIR) -L$(FLUIDSLIBDIR) -lcow -lfluids
INC = -I$(COWINCDIR) -I$(FLUIDSINCDIR)

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

$(LIBDIR)/libfish.so : $(OBJ)
	$(LDSHARED) -o $@ $? $(LIB)

$(LIBDIR)/libfish.a : $(OBJ)
	$(ARSTATIC) $@ $?

$(BINDIR)/testfish : testfish.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

$(BINDIR)/euler : euler.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

clean :
	@rm -rf $(EXE) *.o
