
/*------------------------------------------------------------------------------
 * FILE: lua_fft.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 * http://www.cs.sandia.gov/~sjplimp/docs/fft/README.html
 *
 * NOTES:
 *
 * This code makes the parallel FFT and remapping routines of Steve Plimpton at
 * Sandia National Labs available to lua code.
 *
 * - The option to 'SCALED_YES' to divide by N is only regarded for forward
 *   transforms.
 *
 * - The order of indices provided to the FFT's is fast,mid,slow varying. For C
 *   arrays, this means it gets called with Nz, Ny, Nx.
 *
 * - All of the input data is expected to be without ghost zones, in other words
 *   of dimension domain.GetLocalShape().
 *
 *
 *------------------------------------------------------------------------------
 */


#include "config.h"
extern "C" {
#include "lualib.h"
}
#if (__MARA_USE_FFTW && __MARA_USE_MPI)

#include <cstdlib>
#include <mpi.h>
#include <complex>
#include "luaU.h"
#include "hydro.hpp"
#include "histogram.hpp"

#define LUNUM_API_NOCOMPLEX
extern "C" {
#include "lunum.h"
#include "fft_3d.h"
}

#define NBINS 128

// Public functions exported to Lua
// -----------------------------------------------------------------------------
static int luaC_fft_forward(lua_State *L);
static int luaC_fft_reverse(lua_State *L);
static int luaC_fft_helmholtz(lua_State *L);
static int luaC_fft_power_vector_field(lua_State *L);
static int luaC_fft_power_scalar_field(lua_State *L);


void lua_fft_load(lua_State *L)
{
  lua_register(L, "fft_forward"               , luaC_fft_forward);
  lua_register(L, "fft_reverse"               , luaC_fft_reverse);
  lua_register(L, "fft_helmholtz"             , luaC_fft_helmholtz);
  lua_register(L, "fft_power_vector_field"    , luaC_fft_power_vector_field);
  lua_register(L, "fft_power_scalar_field"    , luaC_fft_power_scalar_field);
}


// Private functions used in this module
// -----------------------------------------------------------------------------
#define SCALED_NOT 0
#define SCALED_YES 1

#define PERMUTE_NONE 0
#define FFT_FWD (+1)
#define FFT_REV (-1)

static struct fft_plan_3d *call_fft_plan_3d(int *nbuf);
static double k_at(int i, int j, int k, double *khat);
static double khat_at(int i, int j, int k, double *khat);
static double cnorm(FFT_DATA z);



int luaC_fft_forward(lua_State *L)
{
  int nbuf;
  struct fft_plan_3d *plan = call_fft_plan_3d(&nbuf);

  if (lunum_upcast(L, 1, ARRAY_TYPE_COMPLEX, nbuf)) {
    lua_replace(L, 1);
  }

  struct Array *Fx = lunum_checkarray1(L, 1);
  struct Array Fk = array_new_zeros(nbuf, ARRAY_TYPE_COMPLEX);

  fft_3d((FFT_DATA*)Fx->data, (FFT_DATA*)Fk.data, FFT_FWD, plan);
  fft_3d_destroy_plan(plan);

  lunum_pusharray1(L, &Fk);
  return 1;
}


int luaC_fft_reverse(lua_State *L)
{
  int nbuf;
  struct fft_plan_3d *plan = call_fft_plan_3d(&nbuf);

  if (lunum_upcast(L, 1, ARRAY_TYPE_COMPLEX, nbuf)) {
    lua_replace(L, 1);
  }

  struct Array *Fk = lunum_checkarray1(L, 1);
  struct Array Fx = array_new_zeros(nbuf, ARRAY_TYPE_COMPLEX);

  fft_3d((FFT_DATA*)Fk->data, (FFT_DATA*)Fx.data, FFT_REV, plan);
  fft_3d_destroy_plan(plan);

  lunum_pusharray1(L, &Fx);
  return 1;
}


int luaC_fft_helmholtz(lua_State *L)
// -----------------------------------------------------------------------------
// input : fx, fy, fz                         ... 3 numeric Lua tables
// output: fx_s, fy_s, fz_s, fx_c, fy_c, fz_c ... their Helmholtz decompositions
//
// This function performs 3d a spectral helmholtz on the input vector field
// (fx,fy,fz). It returns the coordinate (not spectral) realization of the
// solenoidal (div-less) and compressive (curl-less) parts as Lua arrays. Note:
// both the compressive 'c' and solenoidal 's' parts of the returned field
// include the zero-mode. That means that f != c + s.
// -----------------------------------------------------------------------------
{
  clock_t start = clock();

  double *fx_in = luaU_checkarray(L, 1);
  double *fy_in = luaU_checkarray(L, 2);
  double *fz_in = luaU_checkarray(L, 3);

  int nbuf;
  struct fft_plan_3d *plan = call_fft_plan_3d(&nbuf);

  FFT_DATA *fx = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *fy = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *fz = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));

  for (int i=0; i<nbuf; ++i) {
    fx[i].re = fx_in[i];
    fx[i].im = 0.0;

    fy[i].re = fy_in[i];
    fy[i].im = 0.0;

    fz[i].re = fz_in[i];
    fz[i].im = 0.0;
  }

  FFT_DATA *gx = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *gy = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *gz = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));

  fft_3d(fx, gx, FFT_FWD, plan);
  fft_3d(fy, gy, FFT_FWD, plan);
  fft_3d(fz, gz, FFT_FWD, plan);

  free(fx);
  free(fy);
  free(fz);

  FFT_DATA *gx_c = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *gy_c = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *gz_c = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));

  FFT_DATA *gx_s = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *gy_s = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *gz_s = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));

  const int Nx = HydroModule::Mara->domain->GetLocalShape()[0];
  const int Ny = HydroModule::Mara->domain->GetLocalShape()[1];
  const int Nz = HydroModule::Mara->domain->GetLocalShape()[2];

  for (int i=0; i<Nx; ++i) {
    for (int j=0; j<Ny; ++j) {
      for (int k=0; k<Nz; ++k) {

        const int m = i*Ny*Nz + j*Nz + k;
        double khat[3];
        khat_at(i, j, k, khat);

        FFT_DATA gdotk;

        gdotk.re = gx[m].re * khat[0] + gy[m].re * khat[1] + gz[m].re * khat[2];
        gdotk.im = gx[m].im * khat[0] + gy[m].im * khat[1] + gz[m].im * khat[2];

        gx_c[m].re = gdotk.re * khat[0];
        gx_c[m].im = gdotk.im * khat[0];

        gy_c[m].re = gdotk.re * khat[1];
        gy_c[m].im = gdotk.im * khat[1];

        gz_c[m].re = gdotk.re * khat[2];
        gz_c[m].im = gdotk.im * khat[2];

        gx_s[m].re = gx[m].re - gx_c[m].re;
        gx_s[m].im = gx[m].im - gx_c[m].im;

        gy_s[m].re = gy[m].re - gy_c[m].re;
        gy_s[m].im = gy[m].im - gy_c[m].im;

        gz_s[m].re = gz[m].re - gz_c[m].re;
        gz_s[m].im = gz[m].im - gz_c[m].im;
      }
    }
  }

  free(gx);
  free(gy);
  free(gz);


  // Putting in the solenoidal projection
  // ---------------------------------------------------------------------------
  FFT_DATA *fx_s = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *fy_s = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *fz_s = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));

  fft_3d(gx_s, fx_s, FFT_REV, plan);
  fft_3d(gy_s, fy_s, FFT_REV, plan);
  fft_3d(gz_s, fz_s, FFT_REV, plan);

  free(gx_s);
  free(gy_s);
  free(gz_s);

  double *Fx_s = (double*) malloc(nbuf*sizeof(double));
  double *Fy_s = (double*) malloc(nbuf*sizeof(double));
  double *Fz_s = (double*) malloc(nbuf*sizeof(double));

  for (int i=0; i<nbuf; ++i) {
    Fx_s[i] = fx_s[i].re;
    Fy_s[i] = fy_s[i].re;
    Fz_s[i] = fz_s[i].re;
  }

  free(fx_s);
  free(fy_s);
  free(fz_s);

  luaU_pusharray(L, Fx_s, nbuf);
  luaU_pusharray(L, Fy_s, nbuf);
  luaU_pusharray(L, Fz_s, nbuf);

  free(Fx_s);
  free(Fy_s);
  free(Fz_s);


  // Putting in the compressive projection
  // ---------------------------------------------------------------------------
  FFT_DATA *fx_c = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *fy_c = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *fz_c = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));

  fft_3d(gx_c, fx_c, FFT_REV, plan);
  fft_3d(gy_c, fy_c, FFT_REV, plan);
  fft_3d(gz_c, fz_c, FFT_REV, plan);

  free(gx_c);
  free(gy_c);
  free(gz_c);

  double *Fx_c = (double*) malloc(nbuf*sizeof(double));
  double *Fy_c = (double*) malloc(nbuf*sizeof(double));
  double *Fz_c = (double*) malloc(nbuf*sizeof(double));

  for (int i=0; i<nbuf; ++i) {
    Fx_c[i] = fx_c[i].re;
    Fy_c[i] = fy_c[i].re;
    Fz_c[i] = fz_c[i].re;
  }

  free(fx_c);
  free(fy_c);
  free(fz_c);

  luaU_pusharray(L, Fx_c, nbuf);
  luaU_pusharray(L, Fy_c, nbuf);
  luaU_pusharray(L, Fz_c, nbuf);

  free(Fx_c);
  free(Fy_c);
  free(Fz_c);


  fft_3d_destroy_plan(plan);

  printf("[fft] helmholtz decomposition took %3.2f seconds\n",
         (double) (clock() - start) / CLOCKS_PER_SEC);

  return 6;
}



int luaC_fft_power_vector_field(lua_State *L)
// -----------------------------------------------------------------------------
// input : fx, fy, fz, hid, name ... vector field f, and the hdf5 target
// output: nothing               ...
//
// This function gets the power spectrum of the 3d vector field given by
// (fx,fy,fz) by binning it spherically. It writes the result to the (assumed to
// be open) hdf5 file or group id given by 'hid'. Note: in the future this
// function could be modified to return the histogram as a Lua array if 'hid' is
// not provided.
// -----------------------------------------------------------------------------
{
  clock_t start = clock();

  double *fx_in = luaU_checkarray(L, 1);
  double *fy_in = luaU_checkarray(L, 2);
  double *fz_in = luaU_checkarray(L, 3);
  int hid = luaL_checkinteger(L, 4);
  const char *name = luaL_checkstring(L, 5);

  int nbuf;
  struct fft_plan_3d *plan = call_fft_plan_3d(&nbuf);

  FFT_DATA *fx = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *fy = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *fz = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));

  for (int i=0; i<nbuf; ++i) {
    fx[i].re = fx_in[i];
    fx[i].im = 0.0;

    fy[i].re = fy_in[i];
    fy[i].im = 0.0;

    fz[i].re = fz_in[i];
    fz[i].im = 0.0;
  }

  FFT_DATA *gx = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *gy = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  FFT_DATA *gz = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));

  fft_3d(fx, gx, FFT_FWD, plan);
  fft_3d(fy, gy, FFT_FWD, plan);
  fft_3d(fz, gz, FFT_FWD, plan);

  free(fx);
  free(fy);
  free(fz);

  const int Nx = HydroModule::Mara->domain->GetLocalShape()[0];
  const int Ny = HydroModule::Mara->domain->GetLocalShape()[1];
  const int Nz = HydroModule::Mara->domain->GetLocalShape()[2];

  const int *N = HydroModule::Mara->domain->GetGlobalShape();
  Histogram1d hist(NBINS, 1.0, 0.5*sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]),
		   Histogram::Logspace);
  hist.binning_mode = Histogram::BinDensity;

  for (int i=0; i<Nx; ++i) {
    for (int j=0; j<Ny; ++j) {
      for (int k=0; k<Nz; ++k) {

        const int m = i*Ny*Nz + j*Nz + k;
        double khat[3];
        khat_at(i, j, k, khat);

	double kvec[3];
	// ---------------------------------------------------------------------
	// Here we are taking the complex norm (absolute value squared) of the
	// vector-valued Fourier amplitude corresponding to the wave-vector, k.
	//
	//                        P(k) = |\vec{f}_\vec{k}|^2
	//
	// ---------------------------------------------------------------------
	const double Kijk = k_at(i, j, k, kvec);
	const double Pijk = cnorm(gx[m]) + cnorm(gy[m]) + cnorm(gz[m]);

	hist.add_sample(Kijk, Pijk);
      }
    }
  }

  free(gx);
  free(gy);
  free(gz);

  hist.nickname = name;
  hist.synchronize();
  hist.dump_hdf5(hid); // Histogram does nothing if hid == 0

  fft_3d_destroy_plan(plan);

  printf("[fft] power_vector_field took %3.2f seconds\n",
         (double) (clock() - start) / CLOCKS_PER_SEC);

  return 0;
}


int luaC_fft_power_scalar_field(lua_State *L)
// -----------------------------------------------------------------------------
// input : f, hid, name ... scalar field f, and the hdf5 target
// output: nothing      ...
//
// This function gets the power spectrum of the 3d scalar field given by 'f' by
// binning it spherically. It writes the result to the (assumed to be open) hdf5
// file or group id given by 'hid'. Note: in the future this function could be
// modified to return the histogram as a Lua array if 'hid' is not provided.
// -----------------------------------------------------------------------------
{
  clock_t start = clock();

  double *f_in = luaU_checkarray(L, 1);
  int hid = luaL_checkinteger(L, 2);
  const char *name = luaL_checkstring(L, 3);

  int nbuf;
  struct fft_plan_3d *plan = call_fft_plan_3d(&nbuf);
  FFT_DATA *f = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));

  for (int i=0; i<nbuf; ++i) {
    f[i].re = f_in[i];
    f[i].im = 0.0;
  }

  FFT_DATA *g = (FFT_DATA*) malloc(nbuf*sizeof(FFT_DATA));
  fft_3d(f, g, FFT_FWD, plan);
  free(f);

  const int Nx = HydroModule::Mara->domain->GetLocalShape()[0];
  const int Ny = HydroModule::Mara->domain->GetLocalShape()[1];
  const int Nz = HydroModule::Mara->domain->GetLocalShape()[2];

  const int *N = HydroModule::Mara->domain->GetGlobalShape();
  Histogram1d hist(NBINS, 1.0, 0.5*sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]),
		   Histogram::Logspace);
  hist.binning_mode = Histogram::BinDensity;

  for (int i=0; i<Nx; ++i) {
    for (int j=0; j<Ny; ++j) {
      for (int k=0; k<Nz; ++k) {

        const int m = i*Ny*Nz + j*Nz + k;
        double khat[3];
        khat_at(i, j, k, khat);

	double kvec[3];
	// ---------------------------------------------------------------------
	// Here we are taking the complex norm (absolute value squared) of the
	// vector-valued Fourier amplitude corresponding to the wave-vector, k.
	//
	//                        P(k) = |\vec{f}_\vec{k}|^2
	//
	// ---------------------------------------------------------------------
	const double Kijk = k_at(i, j, k, kvec);
	const double Pijk = cnorm(g[m]);

	hist.add_sample(Kijk, Pijk);
      }
    }
  }
  free(g);

  hist.nickname = name;
  hist.synchronize();
  hist.dump_hdf5(hid); // Histogram does nothing if hid == 0

  fft_3d_destroy_plan(plan);

  printf("[fft] power_scalar_field took %3.2f seconds\n",
         (double) (clock() - start) / CLOCKS_PER_SEC);

  return 0;
}













struct fft_plan_3d *call_fft_plan_3d(int *nbuf)
{
  const int i0 = HydroModule::Mara->domain->GetGlobalStart()[0];
  const int i1 = HydroModule::Mara->domain->GetLocalShape()[0] + i0-1;

  const int j0 = HydroModule::Mara->domain->GetGlobalStart()[1];
  const int j1 = HydroModule::Mara->domain->GetLocalShape()[1] + j0-1;

  const int k0 = HydroModule::Mara->domain->GetGlobalStart()[2];
  const int k1 = HydroModule::Mara->domain->GetLocalShape()[2] + k0-1;

  const int Nx = HydroModule::Mara->domain->GetGlobalShape()[0];
  const int Ny = HydroModule::Mara->domain->GetGlobalShape()[1];
  const int Nz = HydroModule::Mara->domain->GetGlobalShape()[2];

  return fft_3d_create_plan(MPI_COMM_WORLD,
                            Nz, Ny, Nx,
                            k0,k1, j0,j1, i0,i1,
                            k0,k1, j0,j1, i0,i1,
                            SCALED_YES, PERMUTE_NONE, nbuf);
}

double k_at(int i, int j, int k, double *kvec)
// -----------------------------------------------------------------------------
// Here, we populate the wave vectors on the Fourier lattice. The convention
// used by FFTW is the same as that used by numpy, described at the link
// below. For N odd, the (positive) Nyquist frequency is placed in the middle
// bin.
//
// http://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.fftfreq.html
// -----------------------------------------------------------------------------
{
  i += HydroModule::Mara->domain->GetGlobalStart()[0];
  j += HydroModule::Mara->domain->GetGlobalStart()[1];
  k += HydroModule::Mara->domain->GetGlobalStart()[2];

  const int Nx = HydroModule::Mara->domain->GetGlobalShape()[0];
  const int Ny = HydroModule::Mara->domain->GetGlobalShape()[1];
  const int Nz = HydroModule::Mara->domain->GetGlobalShape()[2];

  kvec[0] = (Nx % 2 == 0) ?
    ((i<  Nx   /2) ? i : i-Nx):  // N even
    ((i<=(Nx-1)/2) ? i : i-Nx);  // N odd

  kvec[1] = (Ny % 2 == 0) ?
    ((j<  Ny   /2) ? j : j-Ny):
    ((j<=(Ny-1)/2) ? j : j-Ny);

  kvec[2] = (Nz % 2 == 0) ?
    ((k<  Nz   /2) ? k : k-Nz):
    ((k<=(Nz-1)/2) ? k : k-Nz);

  return sqrt(kvec[0]*kvec[0] + kvec[1]*kvec[1] + kvec[2]*kvec[2]);
}

double khat_at(int i, int j, int k, double *khat)
{
  const double k0 = k_at(i,j,k,khat);

  if (fabs(k0) > 1e-12) {
    // don't divide by zero
    khat[0] /= k0;
    khat[1] /= k0;
    khat[2] /= k0;
  }

  return k0;
}

double cnorm(FFT_DATA z)
// http://www.cplusplus.com/reference/std/complex/norm
{
  return z.re*z.re + z.im*z.im;
}

#else
void lua_fft_load(lua_State *L) { }
#endif // (__MARA_USE_FFTW && __MARA_USE_MPI)





