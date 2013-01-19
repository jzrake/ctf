

#include <iostream>
#include <cstring>
#include <algorithm>
#include "eulers.hpp"
#include "weno-split.hpp"
#include "riemann_hll.hpp"
#include "matrix.h"
#include "weno.h"

typedef WenoSplit Deriv;

std::valarray<double> Deriv::dUdt(const std::valarray<double> &Uin)
{
  this->prepare_integration();
  Mara->FailureMask = 0;

  Uglb.resize(stride[0]);
  Pglb.resize(stride[0]);
  Lglb.resize(stride[0]);

  Fiph.resize(stride[0]*(ND>=1));
  Giph.resize(stride[0]*(ND>=2));
  Hiph.resize(stride[0]*(ND>=3));

  Uglb = Uin;
  Pglb = Mara->PrimitiveArray;
  int err = ConsToPrim(Uglb, Pglb);

  if (err != 0) {
    printf("c2p failed on %d zones\n", err);
    throw IntermediateFailure();
  }

  switch (ND) {
  case 1: drive_sweeps_1d(); break;
  case 2: drive_sweeps_2d(); break;
  case 3: drive_sweeps_3d(); break;
  }

  return Lglb;
}


void Deriv::intercell_flux_sweep(const double *U, const double *P,
                                 const double *F, const double *A,
                                 double *Fiph, int dim)
{
  const int Nx = Mara->domain->get_N(dim);
  const int Ng = Mara->domain->get_Ng();

  // Local memory requirements for WENO flux calculation
  // ---------------------------------------------------------------------------
  double *lam = new double[NQ];       // Characteristic eigenvalues
  double *Piph = new double[NQ];      // Primitive variables at i+1/2 (averaged)
  double *Uiph = new double[NQ];      // Conserved " " (taken from prim)
  double *Liph = new double[NQ*NQ];   // Left eigenvectors at i+1/2
  double *Riph = new double[NQ*NQ];   // Right eigenvectors at i+1/2
  double *fp = new double[6*NQ];      // Characteristic split fluxes on local stencil
  double *fm = new double[6*NQ];      // " " Left going
  double *Fp = new double[NQ];        // Component-wise split fluxes on local cell
  double *Fm = new double[NQ];        // " " Left going

  double *f = new double[NQ];         // WENO characteristic flux
  double *fpT = new double[NQ*6];     // Transposed split flux to [wave, zone]
  double *fmT = new double[NQ*6];     // " " Left going
  // ---------------------------------------------------------------------------

  for (int i=Ng-1; i<Nx+Ng; ++i) {

    if (fluxsplit_method == FLUXSPLIT_MARQUINA) {

      // Hard-codes NQ=5 for now
      // -----------------------
      double laml[5], lamr[5];
      double Ll[5][5], Rl[5][5];
      double Lr[5][5], Rr[5][5];
      double ul[6][5], ur[6][5];
      double fl[6][5], fr[6][5];
      double fweno_p[5], fweno_m[5];
      double Fweno_p[5], Fweno_m[5];
      double Pl[5], Pr[5];
      double Ul[5], Ur[5];

      for (int q=0; q<NQ; ++q) {
	const int m = i*NQ + q;
	double v[6] = { P[m-2*NQ], P[m-NQ], P[m], P[m+NQ], P[m+2*NQ], P[m+3*NQ] };
	
	switch (GodunovOperator::reconstruct_method) {
	case RECONSTRUCT_PCM:
	  Pl[q] = v[2];
	  Pr[q] = v[3];
	  break;
	case RECONSTRUCT_PLM:
	  Pl[q] = reconstruct(&v[2], PLM_C2R);
	  Pr[q] = reconstruct(&v[3], PLM_C2L);
	  break;
	case RECONSTRUCT_WENO5:
	  Pl[q] = reconstruct(&v[2], WENO5_FV_C2R);
	  Pr[q] = reconstruct(&v[3], WENO5_FV_C2L);
	  break;
	}
      }
      Mara->fluid->PrimToCons(Pl, Ul);
      Mara->fluid->PrimToCons(Pr, Ur);
      Mara->fluid->Eigensystem(Ul, Pl, Ll[0], Rl[0], laml, dim);
      Mara->fluid->Eigensystem(Ur, Pr, Lr[0], Rr[0], lamr, dim);


      for (int j=0; j<6; ++j) {
        matrix_vector_product(Ll[0], &U[(i+j-2)*NQ], ul[j], NQ, NQ);
        matrix_vector_product(Lr[0], &U[(i+j-2)*NQ], ur[j], NQ, NQ);
        matrix_vector_product(Ll[0], &F[(i+j-2)*NQ], fl[j], NQ, NQ);
        matrix_vector_product(Lr[0], &F[(i+j-2)*NQ], fr[j], NQ, NQ);
      }

      for (int q=0; q<NQ; ++q) {
        if (laml[q] > 0.0 && lamr[q] > 0.0) {
          // No sign change, right-going waves only: set fm to zero and fp to f
          for (int j=0; j<6; ++j) {
            fp[j*NQ + q] = fl[j][q];
            fm[j*NQ + q] = 0.0;
          }
        }
        else if (laml[q] < 0.0 && lamr[q] < 0.0) {
          // No sign change, left-going waves only: set fp to zero and fm to f
          for (int j=0; j<6; ++j) {
            fp[j*NQ + q] = 0.0;
            fm[j*NQ + q] = fr[j][q];
          }
        }
        else {
          // There is a sign change in the speed of this characteristic field
          const double a = fabs(laml[q]) > fabs(lamr[q]) ? laml[q] : lamr[q];
          for (int j=0; j<6; ++j) {
            fp[j*NQ + q] = 0.5*(fl[j][q] + fabs(a)*ul[j][q]);
            fm[j*NQ + q] = 0.5*(fr[j][q] - fabs(a)*ur[j][q]);
          }
        }
      }

      for (int q=0; q<NQ; ++q) {
        for (int j=0; j<6; ++j) {
          fpT[q*6 + j] = fp[j*NQ + q];
          fmT[q*6 + j] = fm[j*NQ + q];
        }
      }

      for (int q=0; q<NQ; ++q) {
        fweno_p[q] = reconstruct(fpT+q*6+2, WENO5_FD_C2R);
        fweno_m[q] = reconstruct(fmT+q*6+3, WENO5_FD_C2L);
      }

      matrix_vector_product(Rl[0], fweno_p, Fweno_p, NQ, NQ);
      matrix_vector_product(Rr[0], fweno_m, Fweno_m, NQ, NQ);

      for (int q=0; q<NQ; ++q) {
	Fiph[i*NQ + q] = Fweno_p[q] + Fweno_m[q];
      }
    }

    else if (fluxsplit_method == FLUXSPLIT_LOCAL_LAX_FRIEDRICHS) {

      for (int q=0; q<NQ; ++q) {
        Piph[q] = 0.5*(P[i*NQ + q] + P[(i+1)*NQ + q]);
      }

      Mara->fluid->PrimToCons(Piph, Uiph);
      Mara->fluid->Eigensystem(Uiph, Piph, Liph, Riph, lam, dim);

      // Select the maximum wavespeed on the local stencil
      // -------------------------------------------------------------------------
      const double ml = *std::max_element(A+i-2, A+i+4);

      for (int j=0; j<6; ++j) {
        for (int q=0; q<NQ; ++q) {

          // Local Lax-Friedrichs flux splitting
          // ---------------------------------------------------------------------
          const int m = (i+j-2)*NQ + q;
          Fp[q] = 0.5*(F[m] + ml*U[m]);
          Fm[q] = 0.5*(F[m] - ml*U[m]);
        }
        matrix_vector_product(Liph, Fp, fp+j*NQ, NQ, NQ);
        matrix_vector_product(Liph, Fm, fm+j*NQ, NQ, NQ);
      }

      for (int q=0; q<NQ; ++q) {
        for (int j=0; j<6; ++j) {
          fpT[q*6 + j] = fp[j*NQ + q];
          fmT[q*6 + j] = fm[j*NQ + q];
        }
      }

      for (int q=0; q<NQ; ++q) {
        f[q] =
          reconstruct(fpT+q*6+2, WENO5_FD_C2R) +
          reconstruct(fmT+q*6+3, WENO5_FD_C2L);
      }

      matrix_vector_product(Riph, f, Fiph+i*NQ, NQ, NQ);
    }
  }

  // Clean up local memory requirements for WENO flux calculation
  // ---------------------------------------------------------------------------
  delete [] lam;
  delete [] Piph;
  delete [] Uiph;
  delete [] Liph;
  delete [] Riph;
  delete [] fp;
  delete [] fm;
  delete [] Fp;
  delete [] Fm;

  delete [] f;
  delete [] fpT;
  delete [] fmT;
  // ---------------------------------------------------------------------------
}


void Deriv::drive_single_sweep(const double *Ug, const double *Pg,
                               double *Fiph_g, int dim)
// Ug   .... Pointer to start of conserved variables for this sweep
// Pg   .... "                 " primitive "                      "
// Fiph .... "                 " intercell flux in the dim-direction
// dim  .... Direction along which to take the sweep (x=1, y=2, z=3)
// -----------------------------------------------------------------------------
{
  const int N = Mara->domain->aug_shape()[dim-1];
  const int S = stride[dim];

  double *U = (double*) malloc(N*NQ*sizeof(double)); // Conserved
  double *P = (double*) malloc(N*NQ*sizeof(double)); // Primitive
  double *F = (double*) malloc(N*NQ*sizeof(double)); // Fluxes in along dim-axis
  double *A = (double*) malloc(N*   sizeof(double)); // Max wavespeed for dim

  for (int i=0; i<N; ++i) {
    // In this loop, i is the zone index, not the memory offset
    // -------------------------------------------------------------------------

    memcpy(U+i*NQ, Ug+i*S, NQ*sizeof(double));
    memcpy(P+i*NQ, Pg+i*S, NQ*sizeof(double));

    double ap, am;
    Mara->fluid->FluxAndEigenvalues(U+i*NQ, P+i*NQ, F+i*NQ, &ap, &am, dim);

    A[i] = (fabs(ap)>fabs(am)) ? fabs(ap) : fabs(am);
    if (MaxLambda < A[i]) MaxLambda = A[i];
  }
  double *Fiph_l = (double*) malloc(N*NQ*sizeof(double));
  intercell_flux_sweep(U, P, F, A, Fiph_l, dim);

  // Here we unload the local intercell fluxes, Fiph_l, into the global array,
  // Fiph_g.
  // ---------------------------------------------------------------------------
  for (int i=0; i<N; ++i) {
    memcpy(Fiph_g+i*S, Fiph_l+i*NQ, NQ*sizeof(double));
  }
  free(U); free(P); free(F); free(A); free(Fiph_l);
}


void Deriv::drive_sweeps_1d()
{
  const int Sx = stride[1];

  drive_single_sweep(&Uglb[0], &Pglb[0], &Fiph[0], 1);

  for (int i=Sx; i<stride[0]; ++i) {
    Lglb[i] = -(Fiph[i]-Fiph[i-Sx])/dx;
  }
}
void Deriv::drive_sweeps_2d()
{
  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);
  const int Ng = Mara->domain->get_Ng();
  const int Sx = stride[1];
  const int Sy = stride[2];

  for (int i=0; i<Nx+2*Ng; ++i) {
    drive_single_sweep(&Uglb[i*Sx], &Pglb[i*Sx], &Giph[i*Sx], 2);
  }
  for (int j=0; j<Ny+2*Ng; ++j) {
    drive_single_sweep(&Uglb[j*Sy], &Pglb[j*Sy], &Fiph[j*Sy], 1);
  }

  Mara->fluid->ConstrainedTransport2d(&Fiph[0], &Giph[0], stride);

  for (int i=Sx; i<stride[0]; ++i) {
    Lglb[i] = -(Fiph[i]-Fiph[i-Sx])/dx - (Giph[i]-Giph[i-Sy])/dy;
  }
}
void Deriv::drive_sweeps_3d()
{
  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);
  const int Nz = Mara->domain->get_N(3);
  const int Ng = Mara->domain->get_Ng();
  const int Sx = stride[1];
  const int Sy = stride[2];
  const int Sz = stride[3];

  for (int j=0; j<Ny+2*Ng; ++j) {
    for (int k=0; k<Nz+2*Ng; ++k) {
      const int m = j*Sy + k*Sz;
      drive_single_sweep(&Uglb[m], &Pglb[m], &Fiph[m], 1);
    }
  }
  for (int k=0; k<Nz+2*Ng; ++k) {
    for (int i=0; i<Nx+2*Ng; ++i) {
      const int m = k*Sz + i*Sx;
      drive_single_sweep(&Uglb[m], &Pglb[m], &Giph[m], 2);
    }
  }
  for (int i=0; i<Nx+2*Ng; ++i) {
    for (int j=0; j<Ny+2*Ng; ++j) {
      const int m = i*Sx + j*Sy;
      drive_single_sweep(&Uglb[m], &Pglb[m], &Hiph[m], 3);
    }
  }
  Mara->fluid->ConstrainedTransport3d(&Fiph[0], &Giph[0], &Hiph[0], stride);

  for (int i=Sx; i<stride[0]; ++i) {
    Lglb[i] =
      -(Fiph[i]-Fiph[i-Sx])/dx +
      -(Giph[i]-Giph[i-Sy])/dy +
      -(Hiph[i]-Hiph[i-Sz])/dz;
  }
}

