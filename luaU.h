

#ifndef __LuaUtilities_HEADER__
#define __LuaUtilities_HEADER__

#ifdef __cplusplus
extern "C" {
#endif

#include "lauxlib.h"

void lua_h5_load(lua_State *L);
void lua_mpi_load(lua_State *L);
void lua_measure_load(lua_State *L);
void lua_fft_load(lua_State *L);
void lua_vis_load(lua_State *L);

void    luaU_stack_dump(lua_State *L);
void    luaU_pusharray(lua_State *L, double *A, int N);
void    luaU_pusharray_wshape(lua_State *L, double *A, const int *shape, int Nd);
void    luaU_pusharray_i(lua_State *L, int *A, int N);
void    luaU_pusharray_astable(lua_State *L, double *A, int N);
void    luaU_pusharray_astable_i(lua_State *L, int *A, int N);
double *luaU_checkarray(lua_State *L, int n);
double *luaU_checklarray(lua_State *L, int n, int *N);
int    *luaU_checkarray_i(lua_State *L, int n);
int    *luaU_checklarray_i(lua_State *L, int n, int *N);

#ifdef __cplusplus
}
#endif

#endif // __LuaUtilities_HEADER__
