

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

static int _fish_block_getneighbor(lua_State *L)
{
  fish_block *B = *((fish_block**) luaL_checkudata(L, 1, "fish::block"));
  int dim = luaL_checkinteger(L, 2);
  int LR = luaL_checkinteger(L, 3);
  fish_block *B1;
  fish_block_getneighbor(B, dim, LR, &B1);
  char *err = fish_block_geterror(B);
  if (err) {
    luaL_error(L, err);
  }
  lua_pushlightuserdata(L, B1);
  return 1;
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
