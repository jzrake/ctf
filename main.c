#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"


#ifndef INSTALL_DIR
#define INSTALL_DIR "."
#endif // INSTALL_DIR


static int traceback(lua_State *L);
int luaopen_buffer(lua_State *L);
int luaopen_hdf5(lua_State *L);
int luaopen_mpi(lua_State *L);
int luaopen_cow(lua_State *L);
int luaopen_Mara(lua_State *L);

#ifndef USE_MPI
int luaopen_mpi(lua_State *L){lua_newtable(L);return 1;}
#endif
#ifndef USE_HDF5
int luaopen_hdf5(lua_State *L){lua_newtable(L);return 1;}
#endif
#ifndef USE_MARA
int luaopen_Mara(lua_State *L){lua_newtable(L);return 1;}
#endif
#ifndef USE_COW
int luaopen_cow(lua_State *L){lua_newtable(L);return 1;}
#endif

int main(int argc, char **argv)
{
  int n;
  lua_State *L = luaL_newstate();
  luaL_openlibs(L);
  luaL_requiref(L, "buffer", luaopen_buffer, 0); lua_pop(L, 1);
  luaL_requiref(L, "HDF5", luaopen_hdf5, 0); lua_pop(L, 1);
  luaL_requiref(L, "MPI", luaopen_mpi, 0); lua_pop(L, 1);
  luaL_requiref(L, "cow", luaopen_cow, 0); lua_pop(L, 1);
  luaL_requiref(L, "Mara", luaopen_Mara, 0); lua_pop(L, 1);


  // Create the global `arg` table
  // ---------------------------------------------------------------------------
  lua_newtable(L);
  for (n=0; n<argc; ++n) {
    lua_pushstring(L, argv[n]);
    lua_rawseti(L, -2, n);
  }
  lua_setglobal(L, "arg");


  // Set the Lua path
  // ---------------------------------------------------------------------------
  lua_getglobal(L, "package");
  lua_pushfstring(L,"%s/?.lua;%s/modules/?.lua;",
		  INSTALL_DIR,
		  INSTALL_DIR);
  lua_setfield(L, -2, "path");
  lua_pop(L, 1);


  // Run the script
  // ---------------------------------------------------------------------------
  lua_pushcfunction(L, traceback);
  if (argc == 1) {
    printf("usage: main script.lua [arg1=val1 arg2=val2]\n");
  }
  else {
    if (luaL_loadfile(L, argv[1])) {
      printf("%s\n", lua_tostring(L, -1));
    }
    else {
      if (lua_pcall(L, 0, 0, -2)) {
	printf("%s\n", lua_tostring(L, -1));
      }
    }
  }
  lua_close(L);
  return 0;
}


int traceback(lua_State *L)
{
  const char *msg = lua_tostring(L, 1);
  if (msg) {
    luaL_traceback(L, L, msg, 1);
  }
  else if (!lua_isnoneornil(L, 1)) {
    if (!luaL_callmeta(L, 1, "__tostring")) {
      lua_pushliteral(L, "(no error message)");
    }
  }
  return 1;
}
