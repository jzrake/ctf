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
int luaopen_fish(lua_State *L);
int luaopen_fluids(lua_State *L);
int luaopen_visual(lua_State *L);


#ifndef USE_MPI
int luaopen_mpi(lua_State *L){lua_newtable(L);return 1;}
#define HAVE_MPI "no"
#else
#define HAVE_MPI "yes"
#endif

#ifndef USE_HDF5
int luaopen_hdf5(lua_State *L){lua_newtable(L);return 1;}
#define HAVE_HDF5 "no"
#else
#define HAVE_HDF5 "yes"
#endif

#ifndef USE_MARA
int luaopen_Mara(lua_State *L){lua_newtable(L);return 1;}
#define HAVE_MARA "no"
#else
#define HAVE_MARA "yes"
#endif

#ifndef USE_COW
int luaopen_cow(lua_State *L){lua_newtable(L);return 1;}
#define HAVE_COW "no"
#else
#define HAVE_COW "yes"
#endif

#ifndef USE_FISH
int luaopen_fish(lua_State *L){lua_newtable(L);return 1;}
int luaopen_fluids(lua_State *L){lua_newtable(L);return 1;}
#define HAVE_FISH "no"
#else
#define HAVE_FISH "yes"
#endif

#ifndef USE_VIS
int luaopen_visual(lua_State *L){lua_newtable(L);return 1;}
#define HAVE_VIS "no"
#else
#define HAVE_VIS "yes"
#endif

#ifndef USE_FFTW
#define HAVE_FFTW "no"
#else
#define HAVE_FFTW "yes"
#endif

#ifndef USE_MPIO
#define HAVE_MPIO "no"
#else
#define HAVE_MPIO "yes"
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
  luaL_requiref(L, "fish", luaopen_fish, 0); lua_pop(L, 1);
  luaL_requiref(L, "fluids", luaopen_fluids, 0); lua_pop(L, 1);
  luaL_requiref(L, "visual", luaopen_visual, 0); lua_pop(L, 1);


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
  lua_pushfstring(L,"%s/?.lua;%s/modules/?.lua;", INSTALL_DIR, INSTALL_DIR);
  lua_setfield(L, -2, "path");
  lua_pop(L, 1);


  // Set the Lua C path
  // ---------------------------------------------------------------------------
  lua_getglobal(L, "package");
  lua_pushfstring(L,"%s/lua-glut/?.so;", INSTALL_DIR);
  lua_setfield(L, -2, "cpath");
  lua_pop(L, 1);


  // Provide some help or run the script
  // ---------------------------------------------------------------------------
  lua_pushcfunction(L, traceback);
  if (argc == 1) {
    printf("**************************************************\n");
    printf("*       Computational Turbulence Framework       *\n");
    printf("*               Jonathan Zrake                   *\n");
    printf("*         New York University 2008-2013          *\n");
    printf("**************************************************\n");
    printf("\nusage: ctf <script.lua> [<options>]\n");
    printf("\nCompiled with support for:\n");
    printf("\tMPI  ... %s\n", HAVE_MPI);
    printf("\tHDF5 ... %s\n", HAVE_HDF5);
    printf("\tFFTW ... %s\n", HAVE_FFTW);
    printf("\tMPIO ... %s\n", HAVE_MPIO);
    printf("\tMARA ... %s\n", HAVE_MARA);
    printf("\tFISH ... %s\n", HAVE_FISH);
    printf("\tVIS  ... %s\n", HAVE_VIS);
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
