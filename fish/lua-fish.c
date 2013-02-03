

#include "fish.h"
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"


#define FISH_STRUCT_TYPE(s)						\
  static void luafish_push_fish_##s(lua_State *L, fish_##s *ini)	\
  {                                                                     \
    fish_##s **ud = (fish_##s**) lua_newuserdata(L, sizeof(fish_##s*)); \
      *ud = ini;							\
      luaL_setmetatable(L, "fish::"#s);					\
  }                                                                     \
  static int _fish_##s##_new(lua_State *L)				\
  {                                                                     \
    fish_##s *a = fish_new();						\
      luafish_push_fish_##s(L, a);					\
	return 1;							\
  }                                                                     \
  static int _fish_##s##_light(lua_State *L)				\
  {                                                                     \
    fish_##s **ud = (fish_##s**) luaL_checkudata(L, 1, "fish::"#s);	\
      lua_pushlightuserdata(L, *ud);					\
      return 1;								\
  }

FISH_STRUCT_TYPE(state)
FISH_STRUCT_TYPE(block)

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
