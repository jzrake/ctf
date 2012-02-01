
#include <iostream>
#include <cstring>
#include "ctu-hancock.hpp"
#include "logging.hpp"
#include "mara_mpi.h"


#define MAXNQ 8 // Used for static array initialization

typedef PlmCtuHancockOperator Deriv;
static double plm_theta = 2.0;
static int DoNotTolerateFailures = 1;

/*------------------------------------------------------------------------------
 *
 * Private inline functions
 *
 */
static inline double sign(double x)
{
  return (x>0)-(x<0);
}
static inline double min3(double a, double b, double c)
{
  const double ab=(a<b)?a:b;
  return (ab<c)?ab:c;
}
static inline double plm_minmod(double ul, double u0, double ur)
{
  const double a = plm_theta * (u0 - ul);
  const double b =     0.5   * (ur - ul);
  const double c = plm_theta * (ur - u0);
  return 0.25*fabs(sign(a)+sign(b))*(sign(a)+sign(c))*min3(fabs(a),fabs(b),fabs(c));
}


void Deriv::SetPlmTheta(double plm)
{
  plm_theta = plm;
}
void Deriv::SetSafetyLevel(int level)
// 0: safest
// 1: ignore c2p and v>c failures by reverting to zone centers
{
  if (level == 0) {
    DoNotTolerateFailures = 1;
  }
  else {
    DoNotTolerateFailures = 0;
  }
}
std::valarray<double> Deriv::dUdt(const std::valarray<double> &Uin)
{
  this->prepare_integration();

  std::valarray<double> U = Uin;
  std::valarray<double> L(U.size());
  std::valarray<double> &P = Mara->PrimitiveArray;

  try {
    this->ConsToPrim(U, P);
  }
  catch (const ConsToPrimFailure &e) {
    throw;
  }

  this->ctu_hancock(&P[0], &U[0], &L[0], TimeStepDt);

  return L;
}
void Deriv::reconstruct_plm(const double *P0, double *Pl, double *Pr, int S)
{
  int i,T=2*S;
  for (i=0; i<NQ; ++i) {
    Pr[i] = P0[S+i] - 0.5 * plm_minmod(P0[ 0+i], P0[S+i], P0[T+i]);
    Pl[i] = P0[0+i] + 0.5 * plm_minmod(P0[-S+i], P0[0+i], P0[S+i]);
  }
}

void Deriv::ctu_hancock(const double *P, const double *U, double *L, double dt)
{
  const int num_dims=ND;
  const int D1=(num_dims>=1);
  const int D2=(num_dims>=2);
  const int D3=(num_dims>=3);
  const int N1=stride[0]*D1;
  const int N2=stride[0]*D2;
  const int N3=stride[0]*D3;

  FluidEquations     &fluid    = *Mara->fluid;
  RiemannSolver      &riemann  = *Mara->riemann;
  BoundaryConditions &boundary = *Mara->boundary;

  std::valarray<double> F(N1), G(N2), H(N3);
  std::valarray<double> Ux(N1), Uy(N2), Uz(N3);
  std::valarray<double> Px(N1), Py(N2), Pz(N3);
  std::valarray<double> dPdx(N1), dPdy(N2), dPdz(N3);

  const int sx=stride[1], sy=stride[2], sz=stride[3];

  Mara->FailureMask = 0;

  /* Step 1
     ---------------------------------------------------------------------------
     Compute the slopes of primitive quantities in each direction.

     Note: These derivatives, computed at the beginning of the time step, are
     used for the reconstructed values in the Hancock and Godunov operators for
     both the predictor and corrector steps.
     ---------------------------------------------------------------------------
  */
  for (int i=sx; i<stride[0]-sx; ++i) {
    if(D1) dPdx[i] = plm_minmod(P[i-sx], P[i], P[i+sx]);
    if(D2) dPdy[i] = plm_minmod(P[i-sy], P[i], P[i+sy]);
    if(D3) dPdz[i] = plm_minmod(P[i-sz], P[i], P[i+sz]);
  }

  /* Step 2
     ---------------------------------------------------------------------------
     Apply the Godunov operator to each face in all directions.

     Notes:

     1) The whole array of intercell fluxes is computed at once, rather than one
     cell at a time. This avoids redundant evaluations of the Godunov operator
     by computing only one solution to the Riemann problem per face.

     2) For the predictor step along a given axis, only the intercell fluxes in
     the transverse directions are used. The Hancock operator is applied along
     the longitudinal direction.
     ---------------------------------------------------------------------------
  */
  for (int i=0; i<stride[0]-sx; i+=NQ) {
    double Pl[MAXNQ], Pr[MAXNQ];

    if (D1) {
      for (int j=0; j<NQ; ++j) {
        Pr[j] = P[i+sx+j] - 0.5*dPdx[i+sx+j];
        Pl[j] = P[i   +j] + 0.5*dPdx[i   +j];
      }
      if (fluid.PrimCheck(Pl) || fluid.PrimCheck(Pr)) {
        // ---------------------------------------------------------------------
        // If the reconstructed primitive state was unphysical,
        // DoNotTolerateFailures=1 will set the fail flag for this zone,
        // triggering an abort for this integration. DoNotTolerateFailures=0
        // indicates that the the zone-centered primitives will be used instead.
        // ---------------------------------------------------------------------
        if (DoNotTolerateFailures) {
          Mara->FailureMask[i/NQ] += 1;
        }
        else {
          riemann.IntercellFlux(&P[i], &P[i+sx], 0, &F[i], 0.0, 1);
        }
      }
      else {
        riemann.IntercellFlux(Pl, Pr, 0, &F[i], 0.0, 1);
      }
    }

    if (D2) {
      for (int j=0; j<NQ; ++j) {
        Pr[j] = P[i+sy+j] - 0.5*dPdy[i+sy+j];
        Pl[j] = P[i   +j] + 0.5*dPdy[i   +j];
      }
      if (fluid.PrimCheck(Pl) || fluid.PrimCheck(Pr)) {
        // ---------------------------------------------------------------------
        // If the reconstructed primitive state was unphysical,
        // DoNotTolerateFailures=1 will set the fail flag for this zone,
        // triggering an abort for this integration. DoNotTolerateFailures=0
        // indicates that the the zone-centered primitives will be used instead.
        // ---------------------------------------------------------------------
        if (DoNotTolerateFailures) {
          Mara->FailureMask[i/NQ] += 1;
        }
        else {
          riemann.IntercellFlux(&P[i], &P[i+sy], 0, &G[i], 0.0, 2);
        }
      }
      else {
        riemann.IntercellFlux(Pl, Pr, 0, &G[i], 0.0, 2);
      }
    }

    if (D3) {
      for (int j=0; j<NQ; ++j) {
        Pr[j] = P[i+sz+j] - 0.5*dPdz[i+sz+j];
        Pl[j] = P[i   +j] + 0.5*dPdz[i   +j];
      }
      if (fluid.PrimCheck(Pl) || fluid.PrimCheck(Pr)) {
        // ---------------------------------------------------------------------
        // If the reconstructed primitive state was unphysical,
        // DoNotTolerateFailures=1 will set the fail flag for this zone,
        // triggering an abort for this integration. DoNotTolerateFailures=0
        // indicates that the the zone-centered primitives will be used instead.
        // ---------------------------------------------------------------------
        if (DoNotTolerateFailures) {
          Mara->FailureMask[i/NQ] += 1;
        }
        else {
          riemann.IntercellFlux(&P[i], &P[i+sz], 0, &H[i], 0.0, 3);
        }
      }
      else {
        riemann.IntercellFlux(Pl, Pr, 0, &H[i], 0.0, 3);
      }
    }
  }
  // ***************************************************************************
  // ABORT POINT
  if (Mara_mpi_int_sum(Mara->FailureMask.sum()))    throw IntermediateFailure();
  // ***************************************************************************


  /* Step 3
     ---------------------------------------------------------------------------
     Apply the Hancock operators and complete the corrector step.

     Notes: This update occurs for a given cell all at once. The Hancock
     operator is evaluated by summing the fluxes on the inner walls of the local
     cell, so there is no danger of redundant calculations.
     ---------------------------------------------------------------------------
  */
  for (int i=sx; i<stride[0]; i+=NQ) {
    double PL[MAXNQ], PR[MAXNQ]; // Capital L/R refers to left and right interior walls of
    double UL[MAXNQ], UR[MAXNQ]; // the local cell, whereas lower-case l/r refer to i_{+-1/2}

    double FL[MAXNQ], FR[MAXNQ];
    double GL[MAXNQ], GR[MAXNQ];
    double HL[MAXNQ], HR[MAXNQ];

    if (D1) {
      for (int j=0; j<NQ; ++j) {
        PL[j] = P[i+j] - 0.5*dPdx[i+j]; // Primitive states on the inner x-facing
        PR[j] = P[i+j] + 0.5*dPdx[i+j]; // walls of the local cell.
      }
      if (fluid.PrimCheck(PL) || fluid.PrimCheck(PR)) {
        if (DoNotTolerateFailures) {
          Mara->FailureMask[i/NQ] += 1;
        }
        else {
          fluid.FluxAndEigenvalues(&U[i], &P[i], FL, 0, 0, 1);
          fluid.FluxAndEigenvalues(&U[i], &P[i], FR, 0, 0, 1);
        }
      }
      else {
        fluid.PrimToCons(PL, UL); // Corresponding conserved quantities and fluxes
        fluid.PrimToCons(PR, UR); // in the x-direction.
        fluid.FluxAndEigenvalues(UL, PL, FL, 0, 0, 1);
        fluid.FluxAndEigenvalues(UR, PR, FR, 0, 0, 1);
      }
    }

    if (D2) {
      for (int j=0; j<NQ; ++j) {
        PL[j] = P[i+j] - 0.5*dPdy[i+j]; // Primitive states on the inner y-facing
        PR[j] = P[i+j] + 0.5*dPdy[i+j]; // walls of the local cell.
      }
      if (fluid.PrimCheck(PL) || fluid.PrimCheck(PR)) {
        if (DoNotTolerateFailures) {
          Mara->FailureMask[i/NQ] += 1;
        }
        else {
          fluid.FluxAndEigenvalues(&U[i], &P[i], GL, 0, 0, 2);
          fluid.FluxAndEigenvalues(&U[i], &P[i], GR, 0, 0, 2);
        }
      }
      else {
        fluid.PrimToCons(PL, UL); // Corresponding conserved quantities and fluxes
        fluid.PrimToCons(PR, UR); // in the y-direction.
        fluid.FluxAndEigenvalues(UL, PL, GL, 0, 0, 2);
        fluid.FluxAndEigenvalues(UR, PR, GR, 0, 0, 2);
      }
    }

    if (D3) {
      for (int j=0; j<NQ; ++j) {
        PL[j] = P[i+j] - 0.5*dPdz[i+j]; // Primitive states on the inner z-facing
        PR[j] = P[i+j] + 0.5*dPdz[i+j]; // walls of the local cell.
      }
      if (fluid.PrimCheck(PL) || fluid.PrimCheck(PR)) {
        if (DoNotTolerateFailures) {
          Mara->FailureMask[i/NQ] += 1;
        }
        else {
          fluid.FluxAndEigenvalues(&U[i], &P[i], HL, 0, 0, 3);
          fluid.FluxAndEigenvalues(&U[i], &P[i], HR, 0, 0, 3);
        }
      }
      else {
        fluid.PrimToCons(PL, UL); // Corresponding conserved quantities and fluxes
        fluid.PrimToCons(PR, UR); // in the z-direction.
        fluid.FluxAndEigenvalues(UL, PL, HL, 0, 0, 3);
        fluid.FluxAndEigenvalues(UR, PR, HR, 0, 0, 3);
      }
    }

    switch (num_dims) {
    case 1:                    /**** 1D ****/
      for (int j=0; j<NQ; ++j) {
        //                   Hancock (normal)
        // ==========================================
        Ux[i+j] = U[i+j] - ((FR[j]-FL[j])/dx)*0.5*dt;
        // ==========================================
      }
      break;
    case 2:                    /**** 2D ****/
      for (int j=0; j<NQ; ++j) {
        //                   Hancock (normal)            Godunov (transverse)
        // ==================================================================
        Ux[i+j] = U[i+j] - ((FR[j]-FL[j])/dx + (G[i+j]-G[i-sy+j])/dy)*0.5*dt;
        Uy[i+j] = U[i+j] - ((GR[j]-GL[j])/dy + (F[i+j]-F[i-sx+j])/dx)*0.5*dt;
        // ==================================================================
      }
      break;
    case 3:                    /**** 3D ****/
      for (int j=0; j<NQ; ++j) {
        //                   Hancock (normal)            Godunov (transverse)
        // ==========================================================================================
        Ux[i+j] = U[i+j] - ((FR[j]-FL[j])/dx + (G[i+j]-G[i-sy+j])/dy + (H[i+j]-H[i-sz+j])/dz)*0.5*dt;
        Uy[i+j] = U[i+j] - ((GR[j]-GL[j])/dy + (H[i+j]-H[i-sz+j])/dz + (F[i+j]-F[i-sx+j])/dx)*0.5*dt;
        Uz[i+j] = U[i+j] - ((HR[j]-HL[j])/dz + (F[i+j]-F[i-sx+j])/dx + (G[i+j]-G[i-sy+j])/dy)*0.5*dt;
        // ==========================================================================================
      }
      break;
    }
  }
  // ***************************************************************************
  // ABORT POINT
  if (Mara_mpi_int_sum(Mara->FailureMask.sum()))    throw IntermediateFailure();
  // ***************************************************************************



  // There is a now predictor conserved state Ux, Uy, Uz at the half-time step:
  // one in each zone for each direction. We need to do a ConsToPrim on each of
  // them, which may fail. If it does fail, then we'll revert to using the input
  // prim and input cons for that zone. Note that the value of P is unchanged by
  // the ConsToPrim call if it fails.
  // ---------------------------------------------------------------------------

  if (D1) boundary.ApplyBoundaries(Ux);
  if (D2) boundary.ApplyBoundaries(Uy);
  if (D3) boundary.ApplyBoundaries(Uz);


  for (int i=0; i<stride[0]; i+=NQ) {
    if (D1) {
      std::memcpy(&Px[i], &P[i], NQ*sizeof(double));
      if (fluid.ConsToPrim(&Ux[i], &Px[i])) {
        // ---------------------------------------------------------------------
        // Unless we are tolerating errors, set the fail flag for this zone.
        // ---------------------------------------------------------------------
        if (DoNotTolerateFailures) {
          Mara->FailureMask[i/NQ] += 1;
        }
        else {
          std::memcpy(&Ux[i], &U[i], NQ*sizeof(double));
        }
      }
    }
    if (D2) {
      std::memcpy(&Py[i], &P[i], NQ*sizeof(double));
      if (fluid.ConsToPrim(&Uy[i], &Py[i])) {
        // ---------------------------------------------------------------------
        // Unless we are tolerating errors, set the fail flag for this zone.
        // ---------------------------------------------------------------------
        if (DoNotTolerateFailures) {
          Mara->FailureMask[i/NQ] += 1;
        }
        else {
          std::memcpy(&Uy[i], &U[i], NQ*sizeof(double));
        }
      }
    }
    if (D3) {
      std::memcpy(&Pz[i], &P[i], NQ*sizeof(double));
      if (fluid.ConsToPrim(&Uz[i], &Pz[i])) {
        // ---------------------------------------------------------------------
        // Unless we are tolerating errors, set the fail flag for this zone.
        // ---------------------------------------------------------------------
        if (DoNotTolerateFailures) {
          Mara->FailureMask[i/NQ] += 1;
        }
        else {
          std::memcpy(&Uz[i], &U[i], NQ*sizeof(double));
        }
      }
    }
  }
  // ***************************************************************************
  // ABORT POINT
  if (Mara_mpi_int_sum(Mara->FailureMask.sum()))    throw IntermediateFailure();
  // ***************************************************************************



  /* Step 4
     ---------------------------------------------------------------------------
     Complete the integration by applying the Godunov fluxes.

     Notes: The final derivative operator is obtained from the Godunov intercell
     fluxes, which in turn are computed using the cell-centered predicted state,
     reconstructed to zone edges using the derivatives from the beginning of the
     time step.
     ---------------------------------------------------------------------------
  */

  for (int i=sx; i<stride[0]-sx; i+=NQ) {
    double Pl[MAXNQ], Pr[MAXNQ];

    if (D1) {
      for (int j=0; j<NQ; ++j) {
        Pr[j] = Px[i+sx+j] - 0.5*dPdx[i+sx+j];
        Pl[j] = Px[i   +j] + 0.5*dPdx[i   +j];
      }
      if (fluid.PrimCheck(Pl) || fluid.PrimCheck(Pr)) {
	if (DoNotTolerateFailures) {
          Mara->FailureMask[i/NQ] += 1;
        }
	else {
	  riemann.IntercellFlux(&P[i], &P[i+sx], 0, &F[i], 0.0, 1);
	}
      }
      else {
        riemann.IntercellFlux(Pl, Pr, 0, &F[i], 0.0, 1);
      }
    }

    if (D2) {
      for (int j=0; j<NQ; ++j) {
        Pr[j] = Py[i+sy+j] - 0.5*dPdy[i+sy+j];
        Pl[j] = Py[i   +j] + 0.5*dPdy[i   +j];
      }
      if (fluid.PrimCheck(Pl) || fluid.PrimCheck(Pr)) {
	if (DoNotTolerateFailures) {
          Mara->FailureMask[i/NQ] += 1;
        }
	else {
	  riemann.IntercellFlux(&P[i], &P[i+sy], 0, &G[i], 0.0, 2);
	}
      }
      else {
        riemann.IntercellFlux(Pl, Pr, 0, &G[i], 0.0, 2);
      }
    }

    if (D3) {
      for (int j=0; j<NQ; ++j) {
        Pr[j] = Pz[i+sz+j] - 0.5*dPdz[i+sz+j];
        Pl[j] = Pz[i   +j] + 0.5*dPdz[i   +j];
      }
      if (fluid.PrimCheck(Pl) || fluid.PrimCheck(Pr)) {
	if (DoNotTolerateFailures) {
          Mara->FailureMask[i/NQ] += 1;
        }
	else {
	  riemann.IntercellFlux(&P[i], &P[i+sz], 0, &H[i], 0.0, 3);
	}
      }
      else {
        riemann.IntercellFlux(Pl, Pr, 0, &H[i], 0.0, 3);
      }
    }
  }

  switch (num_dims) {
    /*---------------------------------- 1D ----------------------------------*/
  case 1:
    // No constrained transport
    for (int i=sx; i<stride[0]; ++i) {
      L[i] = -((F[i]-F[i-sx])/dx);
    }
    break;
    /*---------------------------------- 2D ----------------------------------*/
  case 2:
    fluid.ConstrainedTransport2d(&F[0], &G[0], stride);
    for (int i=sx; i<stride[0]; ++i) {
      L[i] = -((F[i]-F[i-sx])/dx + (G[i]-G[i-sy])/dy);
    }
    break;
    /*---------------------------------- 3D ----------------------------------*/
  case 3:
    fluid.ConstrainedTransport3d(&F[0], &G[0], &H[0], stride);
    for (int i=sx; i<stride[0]; ++i) {
      L[i] = -((F[i]-F[i-sx])/dx + (G[i]-G[i-sy])/dy + (H[i]-H[i-sz])/dz);
    }
    break;
  }
  // ***************************************************************************
  // ABORT POINT
  if (Mara_mpi_int_sum(Mara->FailureMask.sum()))    throw IntermediateFailure();
  // ***************************************************************************
}
