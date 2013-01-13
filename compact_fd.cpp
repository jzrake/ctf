
#include <cstring>
#include "hydro.hpp"
#include "matrix.h"



#ifndef BUILD_COMPACTFD

void Compact_intercell_flux_sweep(const PhysicalDomain &domain,
				  const FluidEquations &fluid,
				  const double *U, const double *P,
				  const double *F, const double *A,
				  double *Fiph, int Nq, int dim)
{ }
#else


/*------------------------------------------------------------------------------
 *
 * Private inline functions
 *
 */

static void compact_fiph(const double *f, double *Fiph, double s, int Nx, int Ng);

void Compact_intercell_flux_sweep(const PhysicalDomain &domain,
				  const FluidEquations &fluid,
				  const double *U, const double *P,
				  const double *F, const double *A,
				  double *Fiph, int Nq, int dim)
{
  const int Nx = domain.get_N(dim);
  const int Ng = domain.get_Ng();

  // Local memory requirements for WENO flux calculation
  // ---------------------------------------------------------------------------
  double *Fp = new double[Nq*(Nx+2*Ng)]; // Component-wise split fluxes on local cell
  double *Fm = new double[Nq*(Nx+2*Ng)]; // " " Left going
  double *Fiph_p = new double[Nq*(Nx+2*Ng)];
  double *Fiph_m = new double[Nq*(Nx+2*Ng)];
  // ---------------------------------------------------------------------------

  const double ml = *std::max_element(A, A+Nx+2*Ng);

  for (int i=0; i<Nx+2*Ng; ++i) {
    for (int q=0; q<Nq; ++q) {
      const int m = i* Nq       + q; // (i,q)
      const int n = q*(Nx+2*Ng) + i; // (q,i)

      Fp[n] = 0.5*(F[m] + ml*U[m]);
      Fm[n] = 0.5*(F[m] - ml*U[m]);
    }
  }

  for (int q=0; q<Nq; ++q) {
    compact_fiph(Fp+q*(Nx+2*Ng), Fiph_p+q*(Nx+2*Ng), +1, Nx, Ng);
    compact_fiph(Fm+q*(Nx+2*Ng), Fiph_m+q*(Nx+2*Ng), -1, Nx, Ng);
  }

  for (int i=0; i<Nx+2*Ng; ++i) {
    for (int q=0; q<Nq; ++q) {
      const int m = i* Nq       + q; // (i,q)
      const int n = q*(Nx+2*Ng) + i; // (q,i)

      Fiph[m] = Fiph_p[n] + Fiph_m[n];
    }
  }

  // Clean up local memory requirements for WENO flux calculation
  // ---------------------------------------------------------------------------
  delete [] Fp;
  delete [] Fm;

  delete [] Fiph_p;
  delete [] Fiph_m;
  // ---------------------------------------------------------------------------
}


extern "C" {
#include "../CSparse/Include/cs.h"
}

static void compact_fiph(const double *f, double *Fiph, double s, int Nx, int Ng)
{
  double phi[Nx+2*Ng], psi[Nx+2*Ng];
  double b[Nx+2*Ng];
  int i;

  for (i=0; i<Nx+2*Ng-1; ++i) {
    phi[i] = (1.0/3.0) + (1.0/6.0)*s;
    psi[i] = (1.0/3.0) - (1.0/6.0)*s;
  }
  for (i=1; i<Nx+2*Ng-2; ++i) {
    b[i] =
      0.5*(1 + s) * ((1.0/18.0) * f[i-1] + (19.0/18.0) * f[i  ] + (5.0/ 9.0) * f[i+1]) +
      0.5*(1 - s) * ((5.0/ 9.0) * f[i  ] + (19.0/18.0) * f[i+1] + (1.0/18.0) * f[i+2]);
  }
  struct cs_sparse *triplet = cs_spalloc(Nx+1, Nx+1, 3*(Nx+1), 1, 1);

  cs_entry(triplet, Nx, Nx, 1.0);
  cs_entry(triplet,  0, Nx, phi[   Ng-1]);
  cs_entry(triplet, Nx,  0, psi[Nx+Ng-1]);

  for (i=0; i<Nx; ++i) {
    cs_entry(triplet, i, i, 1.0);
    cs_entry(triplet, i, i+1, psi[Ng+i-1]);
    cs_entry(triplet, i+1, i, phi[Ng+i  ]);
  }
  struct cs_sparse *matrix = cs_compress(triplet);

  cs_lusol(0, matrix, b-1+Ng, 1e-12);
  cs_spfree(triplet);
  cs_spfree(matrix);
  memcpy(Fiph, b, (Nx+2*Ng)*sizeof(double));
}


#endif
