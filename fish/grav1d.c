

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include "fish.h"

static struct fluids_state **fluid;
static fluids_descr *descr;
static fish_state *scheme;
static int TotalZones;
static int NumGhostZones;
static double DomainLength = 1.0;
static double dx;


static void solve_poisson(double *Rho, double *Phi, double *Gph, double *rhobar);
static void timederiv_nogrv(double *L);
static void timederiv_selfg(double *L);
static void set_bc(double *U, int nq);
static void (*timederiv)(double *L);

void fish_grav1d_setscheme(fish_state *S)
{
  scheme = S;
}

void fish_grav1d_init(fluids_descr *descr_, int N)
{
  NumGhostZones = 3;
  TotalZones = N + 2 * NumGhostZones;

  timederiv = 0 ? timederiv_nogrv : timederiv_selfg;
  descr = descr_;

  fluid = (fluids_state**) malloc(TotalZones * sizeof(fluids_state*));
  dx = DomainLength / N;

  for (int n=0; n<TotalZones; ++n) {
    fluid[n] = fluids_state_new();
    fluids_state_setdescr(fluid[n], descr);
  }
}

void fish_grav1d_finalize()
{
  for (int n=0; n<TotalZones; ++n) {
    fluids_state_del(fluid[n]);
  }
  free(fluid);
}

void fish_grav1d_advance(double dt)
{
  double *L  = (double*) malloc(TotalZones * 5 * sizeof(double));
  double *U0 = (double*) malloc(TotalZones * 5 * sizeof(double));
  double U1[5];
  int method;
  fish_getparami(scheme, &method, FISH_TIME_UPDATE);

  for (int m=0; m < 5 * TotalZones; ++m) { U0[m] = 0.0; }
  for (int n=0; n<TotalZones; ++n) {
    fluids_state_derive(fluid[n], &U0[5*n], FLUIDS_CONSERVED);
  }

  switch (method) {
  case FISH_MIDPOINT:
    timederiv(L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
	U1[q] = U0[5*n + q] + L[5*n + q] * dt * 0.5;
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }

    timederiv(L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
        U1[q] = U0[5*n + q] + L[5*n + q] * dt;
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }
  break;

  case FISH_SHUOSHER_RK3:
    /* ******************************* Step 1 ******************************* */
    timederiv(L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
	int m = 5*n + q;
	U1[q] = U0[m] + L[m] * dt;
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }

    /* ******************************* Step 2 ******************************* */
    timederiv(L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
	int m = 5*n + q;
	U1[q] = 3./4*U0[m] + 1./4*U1[q] + 1./4 * dt * L[m];
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }

    /* ******************************* Step 3 ******************************* */
    timederiv(L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
	int m = 5*n + q;
	U1[q] = 1./3*U0[m] + 2./3*U1[q] + 2./3 * dt * L[m];
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }
    break;
  }

  free(L);
  free(U0);
}

double fish_grav1d_maxwavespeed()
{
  double a = 0.0;
  for (int n=0; n<TotalZones; ++n) {
    double A[5];
    fluids_state_derive(fluid[n], A, FLUIDS_EVAL0);
    for (int q=0; q<5; ++q) {
      a = fabs(A[q]) > a ? fabs(A[q]) : a;
    }
  }
  return a;
}

void fish_grav1d_getprim(double *prim, double *grav)
{
  for (int n=0; n<TotalZones; ++n) {
    double P[5];
    double G[4];
    fluids_state_getattr(fluid[n], P, FLUIDS_PRIMITIVE);
    fluids_state_getattr(fluid[n], G, FLUIDS_GRAVITY);
    memcpy(&prim[5*n], P, 5 * sizeof(double));
    memcpy(&grav[4*n], G, 4 * sizeof(double));
  }
}

void fish_grav1d_setprim(double *prim)
{
  for (int n=0; n<TotalZones; ++n) {
    fluids_state_setattr(fluid[n], &prim[5*n], FLUIDS_PRIMITIVE);
  }
}

void fish_grav1d_mapbuffer(double *x, long flag)
{
  int nq = fluids_descr_getncomp(descr, flag);
  for (int n=0; n<TotalZones; ++n) {
    fluids_state_mapbuffer(fluid[n], &x[nq*n], flag);
  }
}

void solve_poisson(double *Rho, double *Phi, double *Gph, double *rhobar)
{
  int Ng = NumGhostZones;
  int N = TotalZones - Ng;

  fftw_complex *Rhox = fftw_alloc_complex(N);
  fftw_complex *Phix = fftw_alloc_complex(N);
  fftw_complex *Gphx = fftw_alloc_complex(N);
  fftw_complex *Rhok = fftw_alloc_complex(N);
  fftw_complex *Phik = fftw_alloc_complex(N);
  fftw_complex *Gphk = fftw_alloc_complex(N);

  fftw_plan fwd = fftw_plan_dft_1d(N, Rhox, Rhok, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan revphi = fftw_plan_dft_1d(N, Phik, Phix, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan revgph = fftw_plan_dft_1d(N, Gphk, Gphx, FFTW_BACKWARD, FFTW_ESTIMATE);

  for (int i=0; i<N; ++i) {
    Rhox[i][0] = Rho[i+Ng];
    Rhox[i][1] = 0.0;
  }
  fftw_execute(fwd);

  *rhobar = Rhok[0][0] / N;
  Phik[0][0] = 0.0;
  Phik[0][1] = 0.0;
  Gphk[0][0] = 0.0;
  Gphk[0][1] = 0.0;

  /*
   * gph[k] = I * k * (phi[k].re + I * phi[k].im)
   */
  for (int i=1; i<N; ++i) {
    double k = 2*M_PI * (i < N/2 ? i : i-N);
    Phik[i][0] = -Rhok[i][0] / (k*k);
    Phik[i][1] = -Rhok[i][1] / (k*k);
    Gphk[i][0] = -Phik[i][1] * k;
    Gphk[i][1] = +Phik[i][0] * k;
  }
  fftw_execute(revphi);
  fftw_execute(revgph);

  for (int i=0; i<N; ++i) {
    Phi[i+Ng] = Phix[i][0] / N;
    Gph[i+Ng] = Gphx[i][0] / N;
  }

  fftw_destroy_plan(fwd);
  fftw_destroy_plan(revphi);
  fftw_destroy_plan(revgph);
  fftw_free(Rhox);
  fftw_free(Phix);
  fftw_free(Gphx);
  fftw_free(Rhok);
  fftw_free(Phik);
  fftw_free(Gphk);
}

void timederiv_nogrv(double *L)
{
  for (int m=0; m < 5 * TotalZones; ++m) {
    L[m] = 0.0;
  }
  fish_timederivative(scheme, fluid, 1, &TotalZones, &dx, L);
  set_bc(L, 5);
}
void timederiv_selfg(double *L)
{
  int N = TotalZones;

  double *Rho = (double*) malloc(N * sizeof(double));
  double *Phi = (double*) malloc(N * sizeof(double));
  double *Gph = (double*) malloc(N * sizeof(double));
  double P[5];
  double rhobar;

  for (int i=0; i<N; ++i) {
    fluids_state_getattr(fluid[i], P, FLUIDS_PRIMITIVE);
    Rho[i] = P[0];
  }
  solve_poisson(Rho, Phi, Gph, &rhobar);
  set_bc(Phi, 1);
  set_bc(Gph, 1);

  for (int i=0; i<N; ++i) {
    double G[4];
    G[0] = Phi[i];
    G[1] = Gph[i];
    G[2] = 0.0;
    G[3] = 0.0;
    fluids_state_setattr(fluid[i], G, FLUIDS_GRAVITY);
  }
  fluids_descr_setrhobar(descr, rhobar);
  free(Rho);
  free(Phi);
  free(Gph);

  for (int m=0; m < 5 * TotalZones; ++m) {
    L[m] = 0.0;
  }
  fish_timederivative(scheme, fluid, 1, &N, &dx, L);

  double source[5];
  for (int i=0; i<TotalZones; ++i) {
    fluids_state_derive(fluid[i], source, FLUIDS_SOURCETERMS);
    for (int q=0; q<5; ++q) {
      L[5*i + q] += source[q];
    }
  }
  set_bc(L, 5);
}

void set_bc(double *U, int nq)
{
  int N = TotalZones;
  int Ng = NumGhostZones;
  int BC;

  fish_getparami(scheme, &BC, FISH_BOUNDARY_CONDITIONS);

  switch (BC) {
  case FISH_PERIODIC:
    for (int i=0; i<Ng; ++i) {
      for (int q=0; q<nq; ++q) {
	U[(    i    )*nq + q] = U[( N - 2*Ng + i)*nq + q];
	U[(N - i - 1)*nq + q] = U[(-1 + 2*Ng - i)*nq + q];
      }
    }
    break;
  case FISH_OUTFLOW:
    for (int i=0; i<Ng; ++i) {
      for (int q=0; q<nq; ++q) {
	U[(    i    )*nq + q] = U[(    Ng    )*nq + q];
	U[(N - i - 1)*nq + q] = U[(N - Ng - 1)*nq + q];
      }
    }
    break;
  default:
    break; // warning!
  }
}
