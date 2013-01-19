
MAKEFILE_IN = Makefile.in
include $(MAKEFILE_IN)

RM ?= rm -f
CFLAGS ?= -Wall
CFLAGS += -std=c99

LUA_I ?= -I$(LUA_HOME)/include

default : lua-buffer.o

%.o : %.c
	$(CC) $(CFLAGS) -c $^ $(INC)

lua-buffer.o : lua-buffer.c
	$(CC) $(CFLAGS) -c $< $(LUA_I)

clean :
	$(RM) lua-buffer.o
