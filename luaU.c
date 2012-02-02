

#include <stdlib.h>
#include "lunum.h"
#include "luaU.h"

void luaU_pusharray(lua_State *L, double *A, int N)
{
  lunum_pusharray2(L, A, ARRAY_TYPE_DOUBLE, N);
}
void luaU_pusharray_i(lua_State *L, int *A, int N)
{
  lunum_pusharray2(L, A, ARRAY_TYPE_INT, N);
}

void luaU_pusharray_wshape(lua_State *L, double *A, const int *shape, int Nd)
{
  int ntot=1;
  for (int i=0; i<Nd; ++i) ntot *= shape[i];
  lunum_pusharray2(L, A, ARRAY_TYPE_DOUBLE, ntot);
  struct Array *B = lunum_checkarray1(L, -1);
  array_resize(B, shape, Nd);
}

void luaU_pusharray_astable(lua_State *L, double *A, int N)
{
  lunum_pusharray2(L, A, ARRAY_TYPE_DOUBLE, N);
  lunum_astable(L, -1);
  lua_replace(L, -2);
}
void luaU_pusharray_astable_i(lua_State *L, int *A, int N)
{
  lunum_pusharray2(L, A, ARRAY_TYPE_INT, N);
  lunum_astable(L, -1);
  lua_replace(L, -2);
}

double *luaU_checkarray(lua_State *L, int n)
{
  return (double*) lunum_checkarray2(L, n, ARRAY_TYPE_DOUBLE, 0);
}
int *luaU_checkarray_i(lua_State *L, int n)
{
  return (int*) lunum_checkarray2(L, n, ARRAY_TYPE_INT, 0);
}
double *luaU_checklarray(lua_State *L, int n, int *N)
{
  return (double*) lunum_checkarray2(L, n, ARRAY_TYPE_DOUBLE, N);
}
int *luaU_checklarray_i(lua_State *L, int n, int *N)
{
  return (int*) lunum_checkarray2(L, n, ARRAY_TYPE_INT, N);
}
