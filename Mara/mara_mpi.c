

#include <stdlib.h>
#include "config.h"
#include "mara_mpi.h"
//#include "luaU.h"

#if (__MARA_USE_MPI)
#include <mpi.h>
#endif //__MARA_USE_MPI

/*
static int luaC_mpi_dims_create(lua_State *L);
static int luaC_mpi_get_rank(lua_State *L);
static int luaC_mpi_get_size(lua_State *L);
static int luaC_mpi_barrier(lua_State *L);


void lua_mpi_load(lua_State *L)
{
  lua_register(L, "mpi_dims_create", luaC_mpi_dims_create);
  lua_register(L, "mpi_get_rank", luaC_mpi_get_rank);
  lua_register(L, "mpi_get_size", luaC_mpi_get_size);
  lua_register(L, "mpi_barrier", luaC_mpi_barrier);
}
*/

int Mara_mpi_get_rank()
{
#if (__MARA_USE_MPI)
  int rank, run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  if (run_uses_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
  else {
    rank = 0;
  }
  return rank;
#else
  return 0;
#endif //__MARA_USE_MPI
}

int Mara_mpi_get_size()
{
#if (__MARA_USE_MPI)
  int size, run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  if (run_uses_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }
  else {
    size = 1;
  }
  return size;
#else
  return 1;
#endif //__MARA_USE_MPI
}
/*
int luaC_mpi_barrier(lua_State *L)
{
  Mara_mpi_barrier();
  return 0;
}
*/
void Mara_mpi_barrier()
{
#if (__MARA_USE_MPI)
  int run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  if (run_uses_mpi) {
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif //__MARA_USE_MPI
}

int Mara_mpi_int_prod(int myval)
{
#if (__MARA_USE_MPI)
  int val, run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  if (run_uses_mpi) {
    MPI_Allreduce(&myval, &val, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
  }
  else {
    val = myval;
  }
  return val;
#else
  return myval;
#endif //__MARA_USE_MPI
}

int Mara_mpi_int_sum(int myval)
{
#if (__MARA_USE_MPI)
  int val, run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  if (run_uses_mpi) {
    MPI_Allreduce(&myval, &val, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  else {
    val = myval;
  }
  return val;
#else
  return myval;
#endif //__MARA_USE_MPI
}

double Mara_mpi_dbl_min(double myval)
{
#if (__MARA_USE_MPI)
  double val;
  int run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  if (run_uses_mpi) {
    MPI_Allreduce(&myval, &val, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  }
  else {
    val = myval;
  }
  return val;
#else
  return myval;
#endif //__MARA_USE_MPI
}

double Mara_mpi_dbl_max(double myval)
{
#if (__MARA_USE_MPI)
  double val;
  int run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  if (run_uses_mpi) {
    MPI_Allreduce(&myval, &val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
  else {
    val = myval;
  }
  return val;
#else
  return myval;
#endif //__MARA_USE_MPI
}

double Mara_mpi_dbl_sum(double myval)
{
#if (__MARA_USE_MPI)
  double val;
  int run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  if (run_uses_mpi) {
    MPI_Allreduce(&myval, &val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  else {
    val = myval;
  }
  return val;
#else
  return myval;
#endif //__MARA_USE_MPI
}

int Mara_mpi_active()
{
#if (__MARA_USE_MPI)
  int run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  return run_uses_mpi;
#else
  return 0;
#endif
}
/*
int luaC_mpi_dims_create(lua_State *L)
{
  int i;
  const int num_dims = luaL_checkinteger(L, 1);
#if (__MARA_USE_MPI)
  int run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  if (run_uses_mpi) {

    int *dims = (int*) malloc(num_dims*sizeof(int));
    for (i=0; i<num_dims; ++i) dims[i] = 0;

    MPI_Dims_create(Mara_mpi_get_size(), num_dims, dims);
    luaU_pusharray_i(L, dims, num_dims);
    free(dims);

    return 1;
  }
  else {

    int *dims = (int*) malloc(num_dims*sizeof(int));
    for (i=0; i<num_dims; ++i) dims[i] = 1;

    luaU_pusharray_i(L, dims, num_dims);
    free(dims);

    return 1;
  }
#else
  int *dims = (int*) malloc(num_dims*sizeof(int));
  for (i=0; i<num_dims; ++i) dims[i] = 1;

  luaU_pusharray_i(L, dims, num_dims);
  free(dims);

  return 1;
#endif //__MARA_USE_MPI
}


int luaC_mpi_get_rank(lua_State *L)
{
  lua_pushnumber(L, Mara_mpi_get_rank());
  return 1;
}

int luaC_mpi_get_size(lua_State *L)
{
  lua_pushnumber(L, Mara_mpi_get_size());
  return 1;
}
*/


