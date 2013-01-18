

#include "fluids.h"
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"


#define FLUIDS_STRUCT_TYPE(s)						\
  static void luafluids_push_fluids_##s(lua_State *L, fluids_##s *ini)	\
  {                                                                     \
    fluids_##s **ud = (fluids_##s**) lua_newuserdata(L, sizeof(fluids_##s*)); \
      *ud = ini;							\
      luaL_setmetatable(L, "fluids::"#s);				\
  }                                                                     \
  static int _fluids_##s##_new(lua_State *L)				\
  {                                                                     \
    fluids_##s *a = fluids_##s##_new();					\
      luafluids_push_fluids_##s(L, a);					\
	return 1;							\
  }                                                                     \
  static int _fluids_##s##_light(lua_State *L)				\
  {                                                                     \
    fluids_##s **ud = (fluids_##s**) luaL_checkudata(L, 1, "fluids::"#s); \
      lua_pushlightuserdata(L, *ud);					\
      return 1;								\
  }

FLUIDS_STRUCT_TYPE(descr)
FLUIDS_STRUCT_TYPE(state)
FLUIDS_STRUCT_TYPE(riemn)

#include "fluidsfuncs.c"

int luaopen_fluids(lua_State *L)
{
  luaL_Reg fluids_aux[] = {
    {"descr_new", _fluids_descr_new},
    {"state_new", _fluids_state_new},
    {"riemn_new", _fluids_riemn_new},
    {"descr_light", _fluids_descr_light},
    {"state_light", _fluids_state_light},
    {"riemn_light", _fluids_riemn_light},
    {NULL, NULL}};

  luaL_newmetatable(L, "fluids::descr"); lua_pop(L, 1);
  luaL_newmetatable(L, "fluids::state"); lua_pop(L, 1);
  luaL_newmetatable(L, "fluids::riemn"); lua_pop(L, 1);

  lua_newtable(L);
  luaL_setfuncs(L, fluids_aux, 0);
  luaL_setfuncs(L, fluids_module_funcs, 0);
  register_constants(L);

  return 1;
}
