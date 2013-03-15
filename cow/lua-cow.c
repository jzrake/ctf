
#include <string.h>
#include "cow.h"
#include "srhdpack.h"
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
  static int _cow_##s##_light(lua_State *L)				\
  {                                                                     \
    cow_##s **ud = (cow_##s**) luaL_checkudata(L, 1, "cow::"#s);	\
      lua_pushlightuserdata(L, *ud);					\
      return 1;								\
  }

COW_STRUCT_TYPE(domain)
COW_STRUCT_TYPE(dfield)
COW_STRUCT_TYPE(histogram)


static void luacow_push_cow_transform(lua_State *L, cow_transform ini)
{
  cow_transform *op = (cow_transform*) lua_newuserdata(L, sizeof(cow_transform*));
  *op = ini;
  luaL_setmetatable(L, "cow::transform");
}

static int _cow_init(lua_State *L)
{
  int modes = luaL_optinteger(L, 1, 0);
  cow_init(0, NULL, modes);
  return 0;
}

#include "cowfuncs.c"

int _srhdpack_onepointpdfs(lua_State *L)
{
  cow_dfield *f = *((cow_dfield**) luaL_checkudata(L, 1, "cow::dfield"));
  char *which   = (char*) luaL_checkstring(L, 2);
  char *h5fname = (char*) luaL_checkstring(L, 3);
  char *h5gname = (char*) luaL_optstring(L, 4, NULL);
  char *res;
  if (lua_gettop(L) >= 5) {
    double interval[2] = {
      luaL_checknumber(L, 5),
      luaL_checknumber(L, 6) };
    res = srhdpack_onepointpdfs(f, which, h5fname, h5gname, interval);
  }
  else {
    res = srhdpack_onepointpdfs(f, which, h5fname, h5gname, NULL);
  }
  if (res) {
    luaL_error(L, res);
  }
  return 0;
}


int luaopen_cow(lua_State *L)
{
  luaL_Reg cow_aux[] = {
    {"init", _cow_init},
    {"domain_new", _cow_domain_new},
    {"dfield_new", _cow_dfield_new},
    {"histogram_new", _cow_histogram_new},
    {"domain_light", _cow_domain_light},
    {"dfield_light", _cow_dfield_light},
    {"histogram_light", _cow_histogram_light},
    {NULL, NULL}};

  luaL_Reg cow_srhdpack[] = {
    {"onepointpdfs", _srhdpack_onepointpdfs},
    {NULL, NULL}};

  luaL_newmetatable(L, "cow::domain"); lua_pop(L, 1);
  luaL_newmetatable(L, "cow::dfield"); lua_pop(L, 1);
  luaL_newmetatable(L, "cow::histogram"); lua_pop(L, 1);
  luaL_newmetatable(L, "cow::transform"); lua_pop(L, 1);

  lua_newtable(L);
  luaL_setfuncs(L, cow_aux, 0);
  luaL_setfuncs(L, cow_module_funcs, 0);
  register_constants(L);

  lua_newtable(L);
  luaL_setfuncs(L, cow_srhdpack, 0);
  lua_setfield(L, -2, "srhdpack");

  lua_newtable(L);
  luacow_push_cow_transform(L, cow_trans_rot5); lua_setfield(L, -2, "rot5");
  luacow_push_cow_transform(L, cow_trans_div5); lua_setfield(L, -2, "div5");
  luacow_push_cow_transform(L, cow_trans_divcorner); lua_setfield(L, -2, "divcorner");
  luacow_push_cow_transform(L, cow_trans_laplacian); lua_setfield(L, -2, "laplacian");
  luacow_push_cow_transform(L, cow_trans_component); lua_setfield(L, -2, "component");
  luacow_push_cow_transform(L, cow_trans_magnitude); lua_setfield(L, -2, "magnitude");
  luacow_push_cow_transform(L, cow_trans_cross); lua_setfield(L, -2, "cross");
  luacow_push_cow_transform(L, cow_trans_dot3); lua_setfield(L, -2, "dot3");
  lua_setfield(L, -2, "transform");

  return 1;
}

