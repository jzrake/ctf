
MAKEFILE_IN = Makefile.in
include $(MAKEFILE_IN)

RM ?= rm -f
AR ?= ar
ARSTATIC ?= $(AR) rcu
RANLIB ?= ranlib
CFLAGS ?= -Wall

CTF_A = libctf.a
LUA_I ?= -I$(LUA_HOME)/include

OBJ = ctf.o lua-buffer.o

default : $(CTF_A)

%.o : %.c
	$(CC) $(CFLAGS) -c $^ $(INC)

lua-buffer.o : lua-buffer.c
	$(CC) $(CFLAGS) -c $< $(LUA_I)

ctf.o : ctf.c $(LUA_A) $(MAKEFILE_IN)
	$(CC) $(CFLAGS) -c -o $@ $< $(LUA_I) $(DEFINES)

$(CTF_A) : $(OBJ)
	$(ARSTATIC) $@ $?
	$(RANLIB) $@

clean :
	$(RM) ctf.o lua-buffer.o
