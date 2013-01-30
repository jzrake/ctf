
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
 * This code wraps the parallel FFT and remapping routines of Steve Plimpton at
 * Sandia National Labs.
 *
 * - The option to 'SCALED_YES' to divide by N is only regarded for reverse
 *   transforms.
 *
 * - Power spectra do the normalization by hand after the forward FFT.
 *
 * - The order of indices provided to the FFT's is fast,mid,slow varying. For C
 *   arrays, this means it gets called with Nz, Ny, Nx.
 *
 *------------------------------------------------------------------------------
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define COW_PRIVATE_DEFS
#include "cow.h"

#if (COW_FFTW)
#if (COW_MPI)
#include "fft_3d.h"
#else
#include "fftw3.h"
#define FFT_DATA fftw_complex
#endif // COW_MPI
#endif // COW_FFTW
#define MODULE "fft"
#define EFFECTIVELY_ZERO 1e-14
#ifndef M_PI
#define M_PI (4*atan(1.0))
#endif // M_PI


// Private functions used in this module
// -----------------------------------------------------------------------------
#define SCALED_NOT 0
#define SCALED_YES 1

#define PERMUTE_NONE 0
#define FFT_FWD (+1)
#define FFT_REV (-1)

#if (COW_FFTW)
#if (COW_MPI)
static struct fft_plan_3d *call_fft_plan_3d(cow_domain *d, int *nbuf);
#endif // COW_MPI
static double k_at(cow_domain *d, int i, int j, int k, double *khat);
static double k_at_wspc(cow_domain *d, int i, int j, int k, double *khat);
static double khat_at(cow_domain *d, int i, int j, int k, double *khat);
static double cnorm(FFT_DATA z);
static FFT_DATA *_fwd(cow_domain *d, double *fx, int start, int stride);
static double *_rev(cow_domain *d, FFT_DATA *Fk);
#endif // COW_FFTW

void cow_fft_forward(cow_dfield *f, cow_dfield *fkre, cow_dfield *fkim)
{
#if (COW_FFTW)
  int nx = cow_domain_getnumlocalzonesinterior(f->domain, 0);
  int ny = cow_domain_getnumlocalzonesinterior(f->domain, 1);
  int nz = cow_domain_getnumlocalzonesinterior(f->domain, 2);
  int ng = cow_domain_getguard(f->domain);
  int ntot = nx * ny * nz;
  int I0[3] = { ng, ng, ng };
  int I1[3] = { nx + ng, ny + ng, nz + ng };

  double *input = (double*) malloc(ntot * sizeof(double) * f->n_members);
  double *outre = (double*) malloc(ntot * sizeof(double) * f->n_members);
  double *outim = (double*) malloc(ntot * sizeof(double) * f->n_members);
  cow_dfield_extract(f, I0, I1, input);

  for (int m=0; m<f->n_members; ++m) {
    FFT_DATA *g = _fwd(f->domain, input, m, f->n_members); // start, stride
    for (int n=0; n<ntot; ++n) {
      outre[f->n_members*n + m] = g[n][0];
      outim[f->n_members*n + m] = g[n][1];
    }
    free(g);
  }

  cow_dfield_replace(fkre, I0, I1, outre);
  cow_dfield_replace(fkim, I0, I1, outim);
  free(outre);
  free(outim);
  free(input);
#endif // COW_FFTW
}

void cow_fft_reverse(cow_dfield *f, cow_dfield *fkre, cow_dfield *fkim)
{
#if (COW_FFTW)
  int nx = cow_domain_getnumlocalzonesinterior(f->domain, 0);
  int ny = cow_domain_getnumlocalzonesinterior(f->domain, 1);
  int nz = cow_domain_getnumlocalzonesinterior(f->domain, 2);
  int ng = cow_domain_getguard(f->domain);
  int ntot = nx * ny * nz;
  int I0[3] = { ng, ng, ng };
  int I1[3] = { nx + ng, ny + ng, nz + ng };

  double *out = (double*) malloc(ntot * sizeof(double) * f->n_members);
  double *inr = (double*) malloc(ntot * sizeof(double) * f->n_members);
  double *ini = (double*) malloc(ntot * sizeof(double) * f->n_members);
  cow_dfield_extract(fkre, I0, I1, inr);
  cow_dfield_extract(fkim, I0, I1, ini);

  for (int m=0; m<f->n_members; ++m) {
    FFT_DATA *fk = (FFT_DATA*) malloc(ntot * sizeof(FFT_DATA));
    for (int n=0; n<ntot; ++n) {
      fk[n][0] = inr[f->n_members*n + m];
      fk[n][1] = ini[f->n_members*n + m];
    }
    double *fx = _rev(f->domain, fk);
    for (int n=0; n<ntot; ++n) {
      out[f->n_members*n + m] = fx[n];
    }
    free(fx);
    free(fk);
  }

  cow_dfield_replace(f, I0, I1, out);
  free(out);
  free(inr);
  free(ini);
#endif // COW_FFTW
}

void cow_fft_pspecscafield(cow_dfield *f, cow_histogram *hist)
// -----------------------------------------------------------------------------
// This function computes the spherically integrated power spectrum of the
// scalar field represented in `f`. The user needs to supply a half-initialized
// histogram, which has not yet been committed. This function will commit,
// populate, and seal the histogram by doing the FFT's on the vector field
// components. Client code supplies the fields like in the example below, all
// other will be over-written.
//
//  cow_histogram_setnbins(hist, 0, 256);
//  cow_histogram_setspacing(hist, COW_HIST_SPACING_LINEAR); // or LOG
//  cow_histogram_setnickname(hist, "mypspec"); // optional
//
// -----------------------------------------------------------------------------
{
#if (COW_FFTW)
  if (!f->committed) return;
  if (f->n_members != 1) {
    printf("[%s] error: need a 1-component field for %s", MODULE, __FUNCTION__);
    return;
  }

  clock_t start = clock();
  int nx = cow_domain_getnumlocalzonesinterior(f->domain, 0);
  int ny = cow_domain_getnumlocalzonesinterior(f->domain, 1);
  int nz = cow_domain_getnumlocalzonesinterior(f->domain, 2);
  int Nx = cow_domain_getnumglobalzones(f->domain, 0);
  int Ny = cow_domain_getnumglobalzones(f->domain, 1);
  int Nz = cow_domain_getnumglobalzones(f->domain, 2);
  int ng = cow_domain_getguard(f->domain);
  int ntot = nx * ny * nz;
  int I0[3] = { ng, ng, ng };
  int I1[3] = { nx + ng, ny + ng, nz + ng };
  double norm = pow(cow_domain_getnumglobalzones(f->domain, COW_ALL_DIMS), 2.0);
  double *input = (double*) malloc(ntot * sizeof(double));
  cow_dfield_extract(f, I0, I1, input);

  FFT_DATA *gx = _fwd(f->domain, input, 0, 1); // start, stride
  free(input);

  cow_histogram_setlower(hist, 0, 1.0);
  cow_histogram_setupper(hist, 0, 0.5*sqrt(Nx*Nx + Ny*Ny + Nz*Nz));
  cow_histogram_setbinmode(hist, COW_HIST_BINMODE_DENSITY);
  cow_histogram_setdomaincomm(hist, f->domain);
  cow_histogram_commit(hist);
  for (int i=0; i<nx; ++i) {
    for (int j=0; j<ny; ++j) {
      for (int k=0; k<nz; ++k) {
	int m = i*ny*nz + j*nz + k;
	double kvec[3];
        double khat[3];
	khat_at(f->domain, i, j, k, khat);
	// ---------------------------------------------------------------------
	// Here we are taking the complex norm (absolute value squared) of the
	// Fourier amplitude corresponding to the wave-vector, k.
	//
	//                        P(k) = |\vec{f}_\vec{k}|^2
	//
	// ---------------------------------------------------------------------
	double Kijk = k_at(f->domain, i, j, k, kvec);
	double Pijk = cnorm(gx[m]);
	cow_histogram_addsample1(hist, Kijk, Pijk/norm);
      }
    }
  }
  cow_histogram_seal(hist);
  free(gx);
  printf("[%s] %s took %3.2f seconds\n",
	 MODULE, __FUNCTION__, (double) (clock() - start) / CLOCKS_PER_SEC);
#endif // COW_FFTW
}


void cow_fft_pspecvecfield(cow_dfield *f, cow_histogram *hist)
// -----------------------------------------------------------------------------
// This function computes the spherically integrated power spectrum of the
// vector field represented in `f`. The user needs to supply a half-initialized
// histogram, which has not yet been committed. This function will commit,
// populate, and seal the histogram by doing the FFT's on the vector field
// components. The supplies the fields like in the example below, all other will
// be over-written.
//
//  cow_histogram_setnbins(hist, 0, 256);
//  cow_histogram_setspacing(hist, COW_HIST_SPACING_LINEAR); // or LOG
//  cow_histogram_setnickname(hist, "mypspec"); // optional
//
// -----------------------------------------------------------------------------
{
#if (COW_FFTW)
  if (!f->committed) return;
  if (f->n_members != 3) {
    printf("[%s] error: need a 3-component field for %s", MODULE, __FUNCTION__);
    return;
  }

  clock_t start = clock();
  int nx = cow_domain_getnumlocalzonesinterior(f->domain, 0);
  int ny = cow_domain_getnumlocalzonesinterior(f->domain, 1);
  int nz = cow_domain_getnumlocalzonesinterior(f->domain, 2);
  int Nx = cow_domain_getnumglobalzones(f->domain, 0);
  int Ny = cow_domain_getnumglobalzones(f->domain, 1);
  int Nz = cow_domain_getnumglobalzones(f->domain, 2);
  int ng = cow_domain_getguard(f->domain);
  int ntot = nx * ny * nz;
  int I0[3] = { ng, ng, ng };
  int I1[3] = { nx + ng, ny + ng, nz + ng };
  double norm = pow(cow_domain_getnumglobalzones(f->domain, COW_ALL_DIMS), 2.0);
  double *input = (double*) malloc(3 * ntot * sizeof(double));
  cow_dfield_extract(f, I0, I1, input);

  FFT_DATA *gx = _fwd(f->domain, input, 0, 3); // start, stride
  FFT_DATA *gy = _fwd(f->domain, input, 1, 3);
  FFT_DATA *gz = _fwd(f->domain, input, 2, 3);
  free(input);

  cow_histogram_setlower(hist, 0, 1.0);
  cow_histogram_setupper(hist, 0, 0.5*sqrt(Nx*Nx + Ny*Ny + Nz*Nz));
  cow_histogram_setbinmode(hist, COW_HIST_BINMODE_DENSITY);
  cow_histogram_setdomaincomm(hist, f->domain);
  cow_histogram_commit(hist);
  for (int i=0; i<nx; ++i) {
    for (int j=0; j<ny; ++j) {
      for (int k=0; k<nz; ++k) {
	int m = i*ny*nz + j*nz + k;
	double kvec[3];
        double khat[3];
	khat_at(f->domain, i, j, k, khat);
	// ---------------------------------------------------------------------
	// Here we are taking the complex norm (absolute value squared) of the
	// vector-valued Fourier amplitude corresponding to the wave-vector, k.
	//
	//                        P(k) = |\vec{f}_\vec{k}|^2
	//
	// ---------------------------------------------------------------------
	double Kijk = k_at(f->domain, i, j, k, kvec);
	double Pijk = cnorm(gx[m]) + cnorm(gy[m]) + cnorm(gz[m]);
	cow_histogram_addsample1(hist, Kijk, Pijk/norm);
      }
    }
  }
  cow_histogram_seal(hist);
  free(gx);
  free(gy);
  free(gz);
  printf("[%s] %s took %3.2f seconds\n",
	 MODULE, __FUNCTION__, (double) (clock() - start) / CLOCKS_PER_SEC);
#endif // COW_FFTW
}

void cow_fft_helmholtzdecomp(cow_dfield *f, int mode)
{
#if (COW_FFTW)
  if (!f->committed) return;
  if (f->n_members != 3) {
    printf("[%s] error: need a 3-component field for pspecvectorfield", MODULE);
    return;
  }
  clock_t start = clock();
  int nx = cow_domain_getnumlocalzonesinterior(f->domain, 0);
  int ny = cow_domain_getnumlocalzonesinterior(f->domain, 1);
  int nz = cow_domain_getnumlocalzonesinterior(f->domain, 2);
  int ng = cow_domain_getguard(f->domain);
  int ntot = nx * ny * nz;
  int I0[3] = { ng, ng, ng };
  int I1[3] = { nx + ng, ny + ng, nz + ng };

  double *input = (double*) malloc(3 * ntot * sizeof(double));
  cow_dfield_extract(f, I0, I1, input);

  FFT_DATA *gx = _fwd(f->domain, input, 0, 3); // start, stride
  FFT_DATA *gy = _fwd(f->domain, input, 1, 3);
  FFT_DATA *gz = _fwd(f->domain, input, 2, 3);
  free(input);

  FFT_DATA *gx_p = (FFT_DATA*) malloc(ntot * sizeof(FFT_DATA));
  FFT_DATA *gy_p = (FFT_DATA*) malloc(ntot * sizeof(FFT_DATA));
  FFT_DATA *gz_p = (FFT_DATA*) malloc(ntot * sizeof(FFT_DATA));
  for (int i=0; i<nx; ++i) {
    for (int j=0; j<ny; ++j) {
      for (int k=0; k<nz; ++k) {
	int m = i*ny*nz + j*nz + k;
        FFT_DATA gdotk;
        double khat[3];
        khat_at(f->domain, i, j, k, khat);
        gdotk[0] = gx[m][0] * khat[0] + gy[m][0] * khat[1] + gz[m][0] * khat[2];
        gdotk[1] = gx[m][1] * khat[0] + gy[m][1] * khat[1] + gz[m][1] * khat[2];
	switch (mode) {
	case COW_PROJECT_OUT_DIV:
	  gx_p[m][0] = gx[m][0] - gdotk[0] * khat[0];
	  gx_p[m][1] = gx[m][1] - gdotk[1] * khat[0];
	  gy_p[m][0] = gy[m][0] - gdotk[0] * khat[1];
	  gy_p[m][1] = gy[m][1] - gdotk[1] * khat[1];
	  gz_p[m][0] = gz[m][0] - gdotk[0] * khat[2];
	  gz_p[m][1] = gz[m][1] - gdotk[1] * khat[2];
	  break;
	case COW_PROJECT_OUT_CURL:
	  gx_p[m][0] = gdotk[0] * khat[0];
	  gx_p[m][1] = gdotk[1] * khat[0];
	  gy_p[m][0] = gdotk[0] * khat[1];
	  gy_p[m][1] = gdotk[1] * khat[1];
	  gz_p[m][0] = gdotk[0] * khat[2];
	  gz_p[m][1] = gdotk[1] * khat[2];
	  break;
	default: break;
	}
      }
    }
  }
  free(gx);
  free(gy);
  free(gz);
  double *fx_p = _rev(f->domain, gx_p);
  double *fy_p = _rev(f->domain, gy_p);
  double *fz_p = _rev(f->domain, gz_p);
  free(gx_p);
  free(gy_p);
  free(gz_p);

  double *res = (double*) malloc(3 * ntot * sizeof(double));
  for (int i=0; i<ntot; ++i) {
    res[3*i + 0] = fx_p[i];
    res[3*i + 1] = fy_p[i];
    res[3*i + 2] = fz_p[i];
  }
  free(fx_p);
  free(fy_p);
  free(fz_p);

  cow_dfield_replace(f, I0, I1, res);
  cow_dfield_syncguard(f);
  free(res);
  printf("[%s] %s took %3.2f seconds\n",
	 MODULE, __FUNCTION__, (double) (clock() - start) / CLOCKS_PER_SEC);
#endif // COW_FFTW
}


void cow_fft_solvepoisson(cow_dfield *rho, cow_dfield *phi)
{
#if (COW_FFTW)
  if (!rho->committed) return;
  if (!phi->committed) return;

  if (rho->n_members != 1) {
    printf("[%s] error: rho must have 1-component in %s", MODULE, __FUNCTION__);
    return;
  }
  if (phi->n_members != 1) {
    printf("[%s] error: phi must have 1-component in %s", MODULE, __FUNCTION__);
    return;
  }
  if (rho->domain != phi->domain) {
    printf("[%s] error: rho and phi must have the same domain", MODULE);
    return;
  }

  cow_domain *domain = phi->domain;
  int nx = cow_domain_getnumlocalzonesinterior(domain, 0);
  int ny = cow_domain_getnumlocalzonesinterior(domain, 1);
  int nz = cow_domain_getnumlocalzonesinterior(domain, 2);
  int ng = cow_domain_getguard(domain);
  int ntot = nx * ny * nz;
  int I0[3] = { ng, ng, ng };
  int I1[3] = { nx + ng, ny + ng, nz + ng };

  double *rhox = (double*) malloc(ntot * sizeof(double));
  cow_dfield_extract(rho, I0, I1, rhox);

  FFT_DATA *rhok = _fwd(domain, rhox, 0, 1); // start, stride
  FFT_DATA *phik = (FFT_DATA*) malloc(ntot * sizeof(FFT_DATA));

  for (int i=0; i<nx; ++i) {
    for (int j=0; j<ny; ++j) {
      for (int k=0; k<nz; ++k) {
	int m = i*ny*nz + j*nz + k;
        double khat[3];
        double k2 = pow(k_at_wspc(domain, i, j, k, khat), 2.0);
	if (k2 > EFFECTIVELY_ZERO) {
	  phik[m][0] = -rhok[m][0] / k2;
	  phik[m][1] = -rhok[m][1] / k2;
	}
	else {
	  phik[m][0] = 0.0;
	  phik[m][1] = 0.0;
	}
      }
    }
  }

  double *phix = _rev(domain, phik);
  cow_dfield_replace(phi, I0, I1, phix);

  free(rhox);
  free(rhok);
  free(phix);
  free(phik);

#endif // COW_FFTW
}

#if (COW_FFTW)
#if (COW_MPI)
struct fft_plan_3d *call_fft_plan_3d(cow_domain *d, int *nbuf)
{
  const int i0 = cow_domain_getglobalstartindex(d, 0);
  const int i1 = cow_domain_getnumlocalzonesinterior(d, 0) + i0 - 1;
  const int j0 = cow_domain_getglobalstartindex(d, 1);
  const int j1 = cow_domain_getnumlocalzonesinterior(d, 1) + j0 - 1;
  const int k0 = cow_domain_getglobalstartindex(d, 2);
  const int k1 = cow_domain_getnumlocalzonesinterior(d, 2) + k0 - 1;
  const int Nx = cow_domain_getnumglobalzones(d, 0);
  const int Ny = cow_domain_getnumglobalzones(d, 1);
  const int Nz = cow_domain_getnumglobalzones(d, 2);
  return fft_3d_create_plan(d->mpi_cart,
                            Nz, Ny, Nx,
                            k0,k1, j0,j1, i0,i1,
                            k0,k1, j0,j1, i0,i1,
                            SCALED_NOT, PERMUTE_NONE, nbuf);
}
#endif // COW_MPI

FFT_DATA *_fwd(cow_domain *d, double *fx, int start, int stride)
{
  FFT_DATA *Fk = NULL;
  FFT_DATA *Fx = NULL;
  if (cow_mpirunning()) {
#if (COW_MPI)
    int nbuf;
    struct fft_plan_3d *plan = call_fft_plan_3d(d, &nbuf);
    Fx = (FFT_DATA*) malloc(nbuf * sizeof(FFT_DATA));
    Fk = (FFT_DATA*) malloc(nbuf * sizeof(FFT_DATA));
    for (int n=0; n<nbuf; ++n) {
      Fx[n][0] = fx[stride * n + start];
      Fx[n][1] = 0.0;
    }
    fft_3d(Fx, Fk, FFT_FWD, plan);
    free(Fx);
    fft_3d_destroy_plan(plan);
#endif // COW_MPI
  }
  else {
    int nbuf = cow_domain_getnumlocalzonesinterior(d, COW_ALL_DIMS);
    Fx = (FFT_DATA*) malloc(nbuf * sizeof(FFT_DATA));
    Fk = (FFT_DATA*) malloc(nbuf * sizeof(FFT_DATA));
    for (int n=0; n<nbuf; ++n) {
      Fx[n][0] = fx[stride * n + start];
      Fx[n][1] = 0.0;
    }
    int *N = d->L_nint;
    fftw_plan plan = fftw_plan_many_dft(3, N, 1,
					Fx, NULL, 1, 0,
					Fk, NULL, 1, 0,
					FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    free(Fx);
  }
  return Fk;
}
double *_rev(cow_domain *d, FFT_DATA *Fk)
{
  FFT_DATA *Fx = NULL;
  double *fx = NULL;
  if (cow_mpirunning()) {
#if (COW_MPI)
  int nbuf;
  long long ntot = cow_domain_getnumglobalzones(d, COW_ALL_DIMS);
  struct fft_plan_3d *plan = call_fft_plan_3d(d, &nbuf);
  fx = (double*) malloc(nbuf * sizeof(double));
  Fx = (FFT_DATA*) malloc(nbuf * sizeof(FFT_DATA));
  fft_3d(Fk, Fx, FFT_REV, plan);
  for (int n=0; n<nbuf; ++n) {
    fx[n] = Fx[n][0] / ntot;
  }
  free(Fx);
  fft_3d_destroy_plan(plan);
#endif // COW_MPI
  }
  else {
    int nbuf = cow_domain_getnumlocalzonesinterior(d, COW_ALL_DIMS);
    long long ntot = cow_domain_getnumglobalzones(d, COW_ALL_DIMS);
    fx = (double*) malloc(nbuf * sizeof(double));
    Fx = (FFT_DATA*) malloc(nbuf * sizeof(FFT_DATA));
    int *N = d->L_nint;
    fftw_plan plan = fftw_plan_many_dft(3, N, 1,
					Fk, NULL, 1, 0,
					Fx, NULL, 1, 0,
					FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    for (int n=0; n<nbuf; ++n) {
      fx[n] = Fx[n][0] / ntot;
    }
    free(Fx);
    fftw_destroy_plan(plan);
  }
  return fx;
}

double k_at(cow_domain *d, int i, int j, int k, double *kvec)
// -----------------------------------------------------------------------------
// Here, we populate the wave vectors on the Fourier lattice. The convention
// used by FFTW is the same as that used by numpy, described at the link
// below. For N odd, the (positive) Nyquist frequency is placed in the middle
// bin.
//
// http://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.fftfreq.html
// -----------------------------------------------------------------------------
{
  const int Nx = cow_domain_getnumglobalzones(d, 0);
  const int Ny = cow_domain_getnumglobalzones(d, 1);
  const int Nz = cow_domain_getnumglobalzones(d, 2);
  i += cow_domain_getglobalstartindex(d, 0);
  j += cow_domain_getglobalstartindex(d, 1);
  k += cow_domain_getglobalstartindex(d, 2);
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
double k_at_wspc(cow_domain *d, int i, int j, int k, double *kvec)
// -----------------------------------------------------------------------------
// Uses the physical domain lengths to scale the wave-numbers.
// -----------------------------------------------------------------------------
{
  k_at(d, i, j, k, kvec);
  kvec[0] *= 2 * M_PI / (d->glb_upper[0] - d->glb_lower[0]);
  kvec[1] *= 2 * M_PI / (d->glb_upper[1] - d->glb_lower[1]);
  kvec[2] *= 2 * M_PI / (d->glb_upper[2] - d->glb_lower[2]);
  return sqrt(kvec[0]*kvec[0] + kvec[1]*kvec[1] + kvec[2]*kvec[2]);
}
double khat_at(cow_domain *d, int i, int j, int k, double *khat)
{
  const double k0 = k_at(d, i, j, k, khat);
  if (fabs(k0) > EFFECTIVELY_ZERO) { // don't divide by zero
    khat[0] /= k0;
    khat[1] /= k0;
    khat[2] /= k0;
  }
  return k0;
}
double cnorm(FFT_DATA z)
// http://www.cplusplus.com/reference/std/complex/norm
{
  return z[0]*z[0] + z[1]*z[1];
}
#endif // COW_FFTW

