
#include <stdlib.h>

#define FISH_PRIVATE_DEFS
#include "fish.h"
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

/*
 * -----------------------------------------------------------------------------
 * The [descr|state|riemn]_light function serve two purposes:
 *
 * 1. If no argument is given, they return a NULL pointer having the type of the
 *    corresponding struct. That instance may be used as an alias to an
 *    initialized instance, for example through fluids_getdescr(real_descr,
 *    light_descr). Only the real one should be free'd.
 *
 * 2. If the argument is an instance of the corresponding struct, then a
 *    lightuserdata representing the address of the struct itself (rather than
 *    Lua's full user data containing it) is returned.
 *
 * -----------------------------------------------------------------------------
 */
#define FISH_STRUCT_TYPE(s)						\
  static void luafish_push_fish_##s(lua_State *L, fish_##s *ini)	\
  {                                                                     \
    fish_##s **ud = (fish_##s**) lua_newuserdata(L, sizeof(fish_##s*)); \
      *ud = ini;							\
      luaL_setmetatable(L, "fish::"#s);					\
  }                                                                     \
  static int _fish_##s##_new(lua_State *L)				\
  {                                                                     \
    fish_##s *a = fish_##s##_new();					\
      luafish_push_fish_##s(L, a);					\
	return 1;							\
  }                                                                     \
  static int _fish_##s##_light(lua_State *L)				\
  {                                                                     \
    if (lua_isnoneornil(L, 1)) {					\
      luafish_push_fish_##s(L, NULL);					\
	return 1;							\
    }									\
    else {								\
      fish_##s **ud = (fish_##s**) luaL_checkudata(L, 1, "fish::"#s);	\
	lua_pushlightuserdata(L, *ud);					\
	return 1;							\
    }									\
  }

FISH_STRUCT_TYPE(state)
FISH_STRUCT_TYPE(block)

static int _fish_block_neighbor(lua_State *L)
{
  fish_block *B = *((fish_block**) luaL_checkudata(L, 1, "fish::block"));
  int dim = luaL_checkinteger(L, 2);
  int LR = luaL_checkinteger(L, 3);
  fish_block *B1;
  fish_block_neighbor(B, dim, LR, &B1);
  char *err = fish_block_geterror(B);
  if (err) {
    luaL_error(L, err);
  }
  lua_pushlightuserdata(L, B1);
  return 1;
}

static int _fish_block_getboundaryblock(lua_State *L)
{
  fish_block *B = *((fish_block**) luaL_checkudata(L, 1, "fish::block"));
  int dim = luaL_checkinteger(L, 2);
  int LR = luaL_checkinteger(L, 3);
  fish_block *B1;
  fish_block_getboundaryblock(B, dim, LR, &B1);
  char *err = fish_block_geterror(B);
  if (err) {
    luaL_error(L, err);
  }
  lua_pushlightuserdata(L, B1);
  return 1;
}

static int _fish_block_getboundaryflag(lua_State *L)
{
  fish_block *B = *((fish_block**) luaL_checkudata(L, 1, "fish::block"));
  int dim = luaL_checkinteger(L, 2);
  int LR = luaL_checkinteger(L, 3);
  int flag;
  fish_block_getboundaryflag(B, dim, LR, &flag);
  char *err = fish_block_geterror(B);
  if (err) {
    luaL_error(L, err);
  }
  lua_pushinteger(L, flag);
  return 1;
}

static int _fish_block_map(lua_State *L)
{
  fish_block *B = *((fish_block**) luaL_checkudata(L, 1, "fish::block"));
  luaL_checktype(L, 2, LUA_TFUNCTION);
  int attrflag = luaL_checkinteger(L, 3);

  fluids_state **fluid = B->fluid;
  int Nd = B->rank;
  int Ng = B->guard;
  int Nx = B->size[0];
  int Ny = B->size[1];
  int Nz = B->size[2];
  int Nc = fluids_descr_getncomp(B->descr, attrflag);

  int i0 = Nd >= 1 ? Ng : 0;
  int j0 = Nd >= 2 ? Ng : 0;
  int k0 = Nd >= 3 ? Ng : 0;
  int i1 = i0 + Nx;
  int j1 = j0 + Ny;
  int k1 = k0 + Nz;
  int sx=1,sy=1,sz=1;

  switch (Nd) {
  case 1:
    sx = 1;
    break;
  case 2:
    sx = (Ny + 2*Ng);
    sy = 1;
    break;
  case 3:
    sx = (Ny + 2*Ng) * (Nz + 2*Ng);
    sy = (Nz + 2*Ng);
    sz = 1;
    break;
  }
  double *P = (double*) malloc(Nc * sizeof(double));
  /*
  printf("your rank is %d\n", Nd);
  printf("your starts  are [%d %d %d]\n", i0, j0, k0);
  printf("your stops   are [%d %d %d]\n", i1, j1, k1);
  printf("your sizes   are [%d %d %d]\n", Nx, Ny, Nz);
  printf("your strides are [%d %d %d]\n", sx, sy, sz);
  */

  for (int i=i0; i<i1; ++i) {
    for (int j=j0; j<j1; ++j) {
      for (int k=k0; k<k1; ++k) {

	double x = fish_block_positionatindex(B, 0, i);
	double y = fish_block_positionatindex(B, 1, j);
	double z = fish_block_positionatindex(B, 2, k);

	lua_pushvalue(L, 2);
	lua_pushnumber(L, x);
	lua_pushnumber(L, y);
	lua_pushnumber(L, z);
	lua_call(L, 3, 1);

	if (lua_type(L, -1) != LUA_TTABLE) {
	  luaL_error(L, "function must return a table");
	}

	for (int n=0; n<Nc; ++n) {
	  lua_rawgeti(L, -1, n+1);
	  if (lua_type(L, -1) != LUA_TNUMBER) {
	    luaL_error(L, "all table values must be numbers");
	  }
	  P[n] = lua_tonumber(L, -1);
	  lua_pop(L, 1);
	}

	//	printf("%d %d %d\n", i, j, k );
	//	printf("%f %f %f %f %f\n", P[0], P[1], P[2], P[3], P[4]);

	fluids_state_setattr(fluid[i*sx+j*sy+k*sz], P, attrflag);
      }
    }
  }
  free(P);
  return 0;
}

#include "fishfuncs.c"

int luaopen_fish(lua_State *L)
{
  luaL_Reg fish_aux[] = {
    {"state_new", _fish_state_new},
    {"state_light", _fish_state_light},
    {"block_new", _fish_block_new},
    {"block_light", _fish_block_light},
    {NULL, NULL}};

  luaL_newmetatable(L, "fish::state"); lua_pop(L, 1);
  luaL_newmetatable(L, "fish::block"); lua_pop(L, 1);

  lua_newtable(L);
  luaL_setfuncs(L, fish_aux, 0);
  luaL_setfuncs(L, fish_module_funcs, 0);
  register_constants(L);

  return 1;
}
