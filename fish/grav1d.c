

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include "fish.h"

static struct fluids_state *fluid[100];
static fluids_descr *descr;
static double dx = 1.0 / 94;
static double dt = 0.001;

static void solve_poisson(double *Rho, double *Phi, double *Gph, double *rhobar);
static void timederiv(double *L);

void fish_grav1d_init()
{
  descr = fluids_descr_new();
  fluids_descr_setfluid(descr, FLUIDS_NRHYD);
  fluids_descr_setgamma(descr, 1.4);
  fluids_descr_seteos(descr, FLUIDS_EOS_GAMMALAW);
  for (int n=0; n<100; ++n) {
    fluid[n] = fluids_state_new();
    fluids_state_setdescr(fluid[n], descr);
  }
}

void fish_grav1d_finalize()
{
  for (int n=0; n<100; ++n) {
    fluids_state_del(fluid[n]);
  }
  fluids_descr_del(descr);
}

void fish_grav1d_advance()
{
  double L[500], U0[500];

  for (int m=0; m < 5 * 100; ++m) { U0[m] = 0.0; }
  for (int n=0; n<100; ++n) {
    fluids_state_derive(fluid[n], &U0[5*n], FLUIDS_CONSERVED);
  }

  timederiv(L);
  for (int n=0; n<100; ++n) {
    double U[5];
    fluids_state_derive(fluid[n], U, FLUIDS_CONSERVED);
    for (int q=0; q<5; ++q) {
      U[q] += L[5*n + q] * dt * 0.5;
    }
    fluids_state_fromcons(fluid[n], U, FLUIDS_CACHE_DEFAULT);
  }

  timederiv(L);
  for (int n=0; n<100; ++n) {
    double U[5];
    fluids_state_derive(fluid[n], U, FLUIDS_CONSERVED);
    for (int q=0; q<5; ++q) {
      U[q] = U0[5*n + q] + L[5*n + q] * dt;
    }
    fluids_state_fromcons(fluid[n], U, FLUIDS_CACHE_DEFAULT);
  }
}

void fish_grav1d_getprim(double *prim, double *grav)
{
  for (int n=0; n<100; ++n) {
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
  for (int n=0; n<100; ++n) {
    fluids_state_setattr(fluid[n], &prim[5*n], FLUIDS_PRIMITIVE);
  }
}

void fish_grav1d_mapbuffer(double *x, long flag)
{
  for (int n=0; n<100; ++n) {
    fluids_state_mapbuffer(fluid[n], &x[5*n], flag);
  }
}

void solve_poisson(double *Rho, double *Phi, double *Gph, double *rhobar)
{
  int N = 94;
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
    Rhox[i][0] = Rho[i+3];
    Rhox[i][1] = 0.0;
  }
  fftw_execute(fwd);

  *rhobar = Rhok[0][0] / N;
  Phik[0][0] = 0.0;
  Phik[0][1] = 0.0;
  Gphk[0][0] = 0.0;
  Gphk[0][1] = 0.0;

  // gph[k] = I * k * (phi[k].re + I * phi[k].im) =>
  //
  // gph[k].re = -k * phi[k].im
  // gph[k].im = +k * phi[k].re

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
    Phi[i+3] = Phix[i][0] / N;
    Gph[i+3] = Gphx[i][0] / N;
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

void timederiv(double *L)
{
  int N = 100;
  fish_state *S = fish_new();
  fish_setparami(S, FISH_WENO5, FISH_RECONSTRUCTION);
  fish_setparami(S, FLUIDS_RIEMANN_HLLC, FISH_RIEMANN_SOLVER);
  fish_setparami(S, FISH_SPECTRAL, FISH_SOLVER_TYPE);
  fish_setparamd(S, 2.0, FISH_PLM_THETA);

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

  Phi[ 0] = Phi[94]; // set periodic BC's on Phi
  Phi[ 1] = Phi[95];
  Phi[ 2] = Phi[96];
  Phi[97] = Phi[ 3];
  Phi[98] = Phi[ 4];
  Phi[99] = Phi[ 5];

  Gph[ 0] = Gph[94]; // set periodic BC's on Gph
  Gph[ 1] = Gph[95];
  Gph[ 2] = Gph[96];
  Gph[97] = Gph[ 3];
  Gph[98] = Gph[ 4];
  Gph[99] = Gph[ 5];

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

  for (int m=0; m < 5 * 100; ++m) {
    L[m] = 0.0;
  }
  fish_timederivative(S, fluid, 1, &N, &dx, L);

  double source[5];
  for (int i=0; i<100; ++i) {
    fluids_state_derive(fluid[i], source, FLUIDS_SOURCETERMS);
    for (int q=0; q<5; ++q) {
      L[5*i + q] += source[q];
    }
  }
  // periodic BC's
  for (int q=0; q<5; ++q) {
    L[5* 0 + q] = L[5*94 + q];
    L[5* 1 + q] = L[5*95 + q];
    L[5* 2 + q] = L[5*96 + q];
    L[5*97 + q] = L[5* 3 + q];
    L[5*98 + q] = L[5* 4 + q];
    L[5*99 + q] = L[5* 5 + q];
  }

  // outflow BC's
  /*
  for (int q=0; q<5; ++q) {
    L[5* 0 + q] = 0.0;
    L[5* 1 + q] = 0.0;
    L[5* 2 + q] = 0.0;
    L[5*97 + q] = 0.0;
    L[5*98 + q] = 0.0;
    L[5*99 + q] = 0.0;
  }
  */
  fish_del(S);
}
