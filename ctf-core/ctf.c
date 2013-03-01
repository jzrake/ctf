
#include <string.h>
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)


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
  int args_start = argc;
  int show_splash = 0;
  int no_modules = 0;
  int ret_on_finish_opts = -1;
  int quiet = 0;

  for (n=1; n<argc; ++n) {
    if (argv[n][0] != '-') {
      args_start = n;
      break;
    }
    else if (strcmp(argv[n], "--splash") == 0) {
      show_splash = 1;
      ret_on_finish_opts = 0;
    }
    else if (strcmp(argv[n], "--no-modules") == 0) {
      no_modules = 1;
    }
    else if (strcmp(argv[n], "-q") == 0) {
      quiet = 1;
    }
    else {
      fprintf(stderr, "[ctf]: no such option: %s\n", argv[n]);
      return 2;
    }
  }

  if (show_splash) {
      printf("**************************************************\n");
      printf("*       Computational Turbulence Framework       *\n");
      printf("*               Jonathan Zrake                   *\n");
      printf("*         New York University 2008-2013          *\n");
      printf("**************************************************\n");
      printf("\nCompiled with support for:\n");
      printf("\tMPI    : %s\n", HAVE_MPI);
      printf("\tHDF5   : %s\n", HAVE_HDF5);
      printf("\tFFTW   : %s\n", HAVE_FFTW);
      printf("\tMPIO   : %s\n", HAVE_MPIO);
      printf("\tMARA   : %s\n", HAVE_MARA);
      printf("\tFISH   : %s\n", HAVE_FISH);
      printf("\tVIS    : %s\n", HAVE_VIS);
      printf("\nBuild details:\n");
      printf("\tpath   : %s\n", STR(INSTALL_DIR));
      printf("\tcommit : %s\n", STR(GIT_SHA));
  }

  //  printf("args start on %d, there are %d args\n", args_start, argc);
  if (ret_on_finish_opts != -1) {
    return ret_on_finish_opts;
  }

  lua_State *L = luaL_newstate();
  if (!no_modules) {
    /* create a fully loaded Lua state */
    luaL_openlibs(L);
    luaL_requiref(L, "buffer", luaopen_buffer, 0);
    luaL_requiref(L, "HDF5", luaopen_hdf5, 0);
    luaL_requiref(L, "MPI", luaopen_mpi, 0);
    luaL_requiref(L, "cow", luaopen_cow, 0);
    luaL_requiref(L, "Mara", luaopen_Mara, 0);
    luaL_requiref(L, "fish", luaopen_fish, 0);
    luaL_requiref(L, "fluids", luaopen_fluids, 0);
    luaL_requiref(L, "visual", luaopen_visual, 0);
    lua_pop(L, 8);
  }
  else {
    /* create a minimal Lua state (below is a piece of luaL_openlibs from
       linit.c) */
    luaL_requiref(L, "_G", luaopen_base, 1);
    luaL_requiref(L, LUA_LOADLIBNAME, luaopen_package, 1);
    lua_pop(L, 2);
  }

  // Create the global `arg` table
  // ---------------------------------------------------------------------------
  lua_newtable(L);
  for (n=1; n<argc; ++n) {
    lua_pushstring(L, argv[args_start + n - 1]);
    lua_rawseti(L, -2, n);
  }
  lua_setglobal(L, "arg");


  // Set the Lua path
  // ---------------------------------------------------------------------------
  lua_getglobal(L, "package");
  lua_pushstring(L,
                 STR(INSTALL_DIR)"/?.lua;"
                 STR(INSTALL_DIR)"/modules/?.lua;");
  lua_setfield(L, -2, "path");
  lua_pop(L, 1);


  // Set the Lua C path
  // ---------------------------------------------------------------------------
  lua_getglobal(L, "package");
  lua_pushstring(L, STR(INSTALL_DIR)"/lua-glut/?.so;");
  lua_setfield(L, -2, "cpath");
  lua_pop(L, 1);


  // Provide some help or run the script
  // ---------------------------------------------------------------------------
  lua_pushcfunction(L, traceback);
  if (args_start >= argc) {
    fprintf(stderr, "[ctf]: need a program to run\n");
    return 2;
  }
  else {
    if (luaL_loadfile(L, argv[args_start])) {
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
