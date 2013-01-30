

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

int line_integral_convolution(double *A, int Nx, int Ny, const char *fname);

static int _line_integral_convolution(lua_State *L)
{
  double *A = (double*) lua_touserdata(L, 1);
  int Nx = luaL_checkinteger(L, 2);
  int Ny = luaL_checkinteger(L, 3);
  const char *fname = luaL_checkstring(L, 4);
  unsigned N_have = lua_rawlen(L, 1) / sizeof(double);
  unsigned N_needed = Nx * Ny * 2;
  if (N_have != N_needed) {
    luaL_error(L, "[visual] input buffer must be (Nx, Ny, 2)");
  }

  line_integral_convolution(A, Nx, Ny, fname);
  printf("in the _line_integral_convolution\n");
  return 0;
}

int luaopen_visual(lua_State *L)
{
  luaL_Reg vis_module_funcs[] = {
    {"line_integral_convolution", _line_integral_convolution},
    {NULL, NULL}};

  lua_newtable(L);
  luaL_setfuncs(L, vis_module_funcs, 0);

  return 1;
}
