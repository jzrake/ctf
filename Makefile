
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
CFLAGS       ?= -Wall -g -O0

OBJ = fluids.o matrix.o
EXE = 	$(BINDIR)/testfluids

LIBS = $(LIBDIR)/libfluids.so $(LIBDIR)/libfluids.a
HEADERS = $(INCDIR)/fluids.h

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

$(LIBDIR)/libfluids.so : $(OBJ)
	$(LDSHARED) -o $@ $? $(LIB)

$(LIBDIR)/libfluids.a : $(OBJ)
	$(ARSTATIC) $@ $?

$(BINDIR)/testfluids : testfluids.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

clean :
	@rm -rf $(EXE) *.o
