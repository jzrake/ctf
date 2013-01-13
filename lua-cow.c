

#include "cow.h"
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"


#define COW_STRUCT_TYPE(s)						\
  static void luacow_push_cow_##s(lua_State *L, cow_##s *ini)		\
  {                                                                     \
    cow_##s **ud = (cow_##s**) lua_newuserdata(L, sizeof(cow_##s*));	\
      *ud = ini;							\
      luaL_setmetatable(L, "cow::"#s);                                  \
  }                                                                     \
  static int _cow_##s##_new(lua_State *L)				\
  {                                                                     \
    cow_##s *a = cow_##s##_new();					\
      luacow_push_cow_##s(L, a);					\
	return 1;							\
  }                                                                     \

COW_STRUCT_TYPE(domain)
COW_STRUCT_TYPE(dfield)
COW_STRUCT_TYPE(histogram)

static int _cow_init(lua_State *L)
{
  int modes = luaL_checkinteger(L, 1);
  cow_init(0, NULL, modes);
  return 0;
}

#include "cowfuncs.c"

int luaopen_cow(lua_State *L)
{
  luaL_Reg cow_aux[] = {
    {"init", _cow_init},
    {"domain_new", _cow_domain_new},
    {"dfield_new", _cow_dfield_new},
    {"histogram_new", _cow_histogram_new},
    {NULL, NULL}};

  luaL_newmetatable(L, "cow::domain"); lua_pop(L, 1);
  luaL_newmetatable(L, "cow::dfield"); lua_pop(L, 1);
  luaL_newmetatable(L, "cow::histogram"); lua_pop(L, 1);

  lua_newtable(L);
  luaL_setfuncs(L, cow_aux, 0);
  luaL_setfuncs(L, cow_module_funcs, 0);
  register_constants(L);
  return 1;
}

