


#include "config.h"
#include "lualib.h"
#if (__MARA_USE_HDF5)

#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

#include "lauxlib.h"
#include "luaU.h"
#include "lunum.h"


static int luaC_h5_read_string(lua_State *L);
static int luaC_h5_write_string(lua_State *L);
static int luaC_h5_read_numeric_table(lua_State *L);
static int luaC_h5_write_numeric_table(lua_State *L);
static int luaC_h5_read_array(lua_State *L);
static int luaC_h5_write_array(lua_State *L);
static int luaC_h5_open_file(lua_State *L);
static int luaC_h5_open_file(lua_State *L);
static int luaC_h5_close_file(lua_State *L);
static int luaC_h5_open_group(lua_State *L);
static int luaC_h5_close_group(lua_State *L);
static int luaC_h5_get_nsets(lua_State *L);
static int luaC_h5_get_ndims(lua_State *L);


static hid_t PresentFile = 0;
static hid_t PresentGroup = 0;

static lua_State *Lua = NULL;

static herr_t group_to_lua_table(hid_t loc_id, const char *name, void *opdata)
{
  double v;

  hid_t fspc = H5Screate(H5S_SCALAR);
  hid_t dset = H5Dopen(PresentGroup, name, H5P_DEFAULT);
  H5Dread(dset, H5T_NATIVE_DOUBLE, fspc, fspc, H5P_DEFAULT, &v);
  H5Dclose(dset);

  lua_pushstring(Lua, name);
  lua_pushnumber(Lua, v);
  lua_settable(Lua, -3);

  H5Sclose(fspc);
  return 0;
}

static herr_t group_nelem(hid_t loc_id, const char *name, void *opdata)
{
  const int n = lua_tonumber(Lua, -1);
  lua_pop(Lua, 1);
  lua_pushnumber(Lua, n+1);
  return 0;
}

void lua_h5_load(lua_State *L)
{
  lua_register(L, "h5_read_string"        , luaC_h5_read_string);
  lua_register(L, "h5_write_string"       , luaC_h5_write_string);
  lua_register(L, "h5_read_numeric_table" , luaC_h5_read_numeric_table);
  lua_register(L, "h5_write_numeric_table", luaC_h5_write_numeric_table);
  lua_register(L, "h5_read_array"         , luaC_h5_read_array);
  lua_register(L, "h5_write_array"        , luaC_h5_write_array);
  lua_register(L, "h5_open_file"          , luaC_h5_open_file);
  lua_register(L, "h5_close_file"         , luaC_h5_close_file);
  lua_register(L, "h5_open_group"         , luaC_h5_open_group);
  lua_register(L, "h5_close_group"        , luaC_h5_close_group);
  lua_register(L, "h5_get_nsets"          , luaC_h5_get_nsets);
  lua_register(L, "h5_get_ndims"          , luaC_h5_get_ndims);
}

int luaC_h5_read_string(lua_State *L)
{
  const char *dsetnm = luaL_checkstring(L, 1);

  if (!PresentFile) {
    printf("[hdf5] error: no open file.\n");
    return 0;
  }

  hid_t dset = H5Dopen(PresentFile, dsetnm, H5P_DEFAULT);
  hid_t fspc = H5Dget_space(dset);
  hid_t strn = H5Dget_type(dset);
  hsize_t msize = H5Tget_size(strn);

  char *string = (char*) malloc(msize*sizeof(char));

  H5Dread(dset, strn, fspc, fspc, H5P_DEFAULT, string);
  string[msize] = '\0'; // Make sure to null-terminate the string
  lua_pushstring(L, string);

  H5Tclose(strn);
  H5Sclose(fspc);
  H5Dclose(dset);
  free(string);

  return 1;
}

int luaC_h5_write_string(lua_State *L)
{
  const char *dsetnm = luaL_checkstring(L, 1);
  const char *string = luaL_checkstring(L, 2);

  if (!PresentFile) {
    printf("[hdf5] error: no open file.\n");
    return 0;
  }

  hsize_t size = strlen(string);
  hid_t fspc = H5Screate(H5S_SCALAR);
  hid_t strn = H5Tcopy(H5T_C_S1);
  H5Tset_size(strn, size);

  hid_t dset = H5Dcreate(PresentFile, dsetnm, strn, fspc,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, strn, fspc, fspc, H5P_DEFAULT, string);
  H5Dclose(dset);
  H5Tclose(strn);
  H5Sclose(fspc);

  return 0;
}

int luaC_h5_open_file(lua_State *L)
{
  const char *fname = luaL_checkstring(L, 1);
  const char *mode  = luaL_checkstring(L, 2);

  if (PresentFile) {
    H5Fclose(PresentFile);
    PresentFile = 0;
  }

  if (strcmp(mode, "w") == 0) {
    printf("[hdf5] creating file %s\n", fname);
    PresentFile = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }
  else if ((strcmp(mode, "r+") == 0)) {
    printf("[hdf5] opening file %s\n", fname);
    PresentFile = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  }
  else if ((strcmp(mode, "r") == 0)) {
    printf("[hdf5] opening file %s\n", fname);
    PresentFile = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  }
  else {
    printf("[hdf5] error: invalid file access mode '%s'\n", mode);
  }
  return 0;
}

int luaC_h5_close_file(lua_State *L)
{
  if (PresentFile) {
    H5Fclose(PresentFile);
    PresentFile = 0;
    printf("[hdf5] closing file\n");
  }
  else {
    printf("[hdf5] no open file to close\n");
  }
  return 0;
}

int luaC_h5_open_group(lua_State *L)
{
  const char *gname = luaL_checkstring(L, 1);
  const char *mode  = luaL_checkstring(L, 2);

  hid_t grp = 0;

  if (!PresentFile) {
    printf("[hdf5] error: need an open file to open group\n");
  }
  else if (strcmp(mode, "w") == 0) {
    if (H5Lexists(PresentFile, gname, H5P_DEFAULT)) {
      H5Ldelete(PresentFile, gname, H5P_DEFAULT);
    }
    grp = H5Gcreate(PresentFile, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  else if (strcmp(mode, "r+") == 0) {
    if (H5Lexists(PresentFile, gname, H5P_DEFAULT)) {
      grp = H5Gopen(PresentFile, gname, H5P_DEFAULT);
    }
    else {
      grp = H5Gcreate(PresentFile, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
  }
  else {
    printf("[hdf5] error: invalid group access mode '%s'\n", mode);
  }

  lua_pushnumber(L, grp);
  return 1;
}

int luaC_h5_close_group(lua_State *L)
{
  const hid_t grp = luaL_checkinteger(L, 1);
  H5Gclose(grp);
  return 0;
}

int luaC_h5_write_numeric_table(lua_State *L)
{
  const char *gname = luaL_checkstring(L, 1);

  if (!PresentFile) {
    printf("[hdf5] error: no open file.\n");
    return 0;
  }
  hid_t grp = H5Gcreate(PresentFile, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t fspc = H5Screate(H5S_SCALAR);

  lua_pushnil(L);

  while (lua_next(L, 2) != 0) {

    if (lua_type(L, -2) != LUA_TSTRING) {
      printf("[hdf5] error: all keys must be strings.\n");
      return 0;
    }
    if (lua_type(L, -1) != LUA_TNUMBER) {
      printf("[hdf5] error: all values must be numbers.\n");
      return 0;
    }

    const double v = lua_tonumber(L, -1);
    hid_t dset = H5Dcreate(grp, lua_tostring(L, -2), H5T_NATIVE_DOUBLE, fspc,
			   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, fspc, fspc, H5P_DEFAULT, &v);
    H5Dclose(dset);

    lua_pop(L, 1);
  }

  H5Sclose(fspc);
  H5Gclose(grp);
  return 0;
}

int luaC_h5_read_array(lua_State *L)
{
  const char *dsetnm = luaL_checkstring(L, 1);

  if (!PresentFile) {
    printf("[hdf5] error: no open file.\n");
    return 0;
  }

  hid_t dset = H5Dopen(PresentFile, dsetnm, H5P_DEFAULT);
  hid_t fspc = H5Dget_space(dset);
  size_t Nd  = H5Sget_simple_extent_ndims(fspc);
  hsize_t *dims = (hsize_t*) malloc(Nd*sizeof(hsize_t));
  int *dims_int =     (int*) malloc(Nd*sizeof(int));

  H5Sget_simple_extent_dims(fspc, dims, NULL);

  int i, ntot=1;
  for (i=0; i<Nd; ++i) {
    dims_int[i] = dims[i];
    ntot *= dims[i];
  }

  double *data = (double*) malloc(ntot * sizeof(double));
  H5Dread(dset, H5T_NATIVE_DOUBLE, fspc, fspc, H5P_DEFAULT, data);
  luaU_pusharray_wshape(L, data, dims_int, Nd);

  H5Sclose(fspc);
  H5Dclose(dset);
  free(data);

  return 1;
}

int luaC_h5_write_array(lua_State *L)
{
  const char *dsetnm = luaL_checkstring(L, 1);
  struct Array *A = lunum_checkarray1(L, 2);
  int i;

  if (!PresentFile) {
    printf("[hdf5] error: no open file.\n");
    return 0;
  }
  if (A->dtype != ARRAY_TYPE_DOUBLE) {
    luaL_error(L, "[hdf5] only double arrays can be written");
  }

  int ndims = A->ndims;
  hsize_t *sizes = (hsize_t*) malloc(ndims*sizeof(hsize_t));

  for (i=0; i<ndims; ++i) {
    sizes[i] = A->shape[i];
  }
  hid_t fspc = H5Screate_simple(ndims, sizes, H5P_DEFAULT);
  hid_t dset = H5Dcreate(PresentFile, dsetnm, H5T_NATIVE_DOUBLE, fspc,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Dwrite(dset, H5T_NATIVE_DOUBLE, fspc, fspc, H5P_DEFAULT, A->data);
  H5Dclose(dset);
  H5Sclose(fspc);
  free(sizes);

  return 0;
}

int luaC_h5_read_numeric_table(lua_State *L)
{
  const char *gname = luaL_checkstring(L, 1);

  if (!PresentFile) {
    printf("[hdf5] error: no open file.\n");
    return 0;
  }

  Lua = L;
  lua_newtable(L);

  PresentGroup = H5Gopen(PresentFile, gname, H5P_DEFAULT);
  H5Giterate(PresentGroup, ".", NULL, group_to_lua_table, NULL);
  H5Gclose(PresentGroup);

  Lua = NULL;
  return 1;
}

int luaC_h5_get_nsets(lua_State *L)
{
  const char *gname = luaL_checkstring(L, 1);

  if (!PresentFile) {
    printf("[hdf5] error: no open file.\n");
    return 0;
  }

  Lua = L;
  lua_pushnumber(L, 0);

  PresentGroup = H5Gopen(PresentFile, gname, H5P_DEFAULT);
  H5Giterate(PresentGroup, ".", NULL, group_nelem, NULL);
  H5Gclose(PresentGroup);

  Lua = NULL;
  return 1;
}

int luaC_h5_get_ndims(lua_State *L)
{
  const char *dsetnm = luaL_checkstring(L, 1);

  if (!PresentFile) {
    printf("[hdf5] error: no open file.\n");
    return 0;
  }

  hid_t dset = H5Dopen(PresentFile, dsetnm, H5P_DEFAULT);
  hid_t fspc = H5Dget_space(dset);
  size_t Nd  = H5Sget_simple_extent_ndims(fspc);
  hsize_t *dims = (hsize_t*) malloc(Nd*sizeof(hsize_t));
  int *dims_int =     (int*) malloc(Nd*sizeof(int));

  H5Sget_simple_extent_dims(fspc, dims, NULL);

  int i;
  for (i=0; i<Nd; ++i) {
    dims_int[i] = dims[i];
  }
  luaU_pusharray_i(L, dims_int, Nd);

  H5Sclose(fspc);
  H5Dclose(dset);
  free(dims);
  free(dims_int);

  return 1;
}


#else
void lua_h5_load(lua_State *L) { }
#endif // __MARA_USE_HDF5
