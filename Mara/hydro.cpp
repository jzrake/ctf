
/*------------------------------------------------------------------------------
 * FILE: hydro.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#include <cstdio>
#include <fstream>
#include "mara.hpp"

MaraApplication *HydroModule::Mara;

MaraApplication::MaraApplication()
{
  units    = new PhysicalUnits(1.0, 1.0, 1.0);
  domain   = NULL;
  boundary = NULL;
  fluid    = new AdiabaticIdealEulers;
  eos      = new AdiabaticEos(1.4);
  godunov  = NULL;
  riemann  = NULL;
  advance  = NULL;
  driving  = NULL;
  cooling  = NULL;
  srcterm = NULL;
}

MaraApplication::~MaraApplication()
{
  if (units)    delete units;
  if (domain)   delete domain;
  if (boundary) delete boundary;
  if (fluid)    delete fluid;
  if (eos)      delete eos;
  if (godunov)  delete godunov;
  if (riemann)  delete riemann;
  if (advance)  delete advance;
  if (driving)  delete driving;
  if (cooling)  delete cooling;
  if (srcterm)  delete srcterm;
}



// -----------------------------------------------------------------------------
// RiemannSolver
// -----------------------------------------------------------------------------
double RiemannSolver::GetMaxLambda()
{
  return MaxLambda;
}
void RiemannSolver::ResetMaxLambda()
{
  MaxLambda = 0;
}
double RiemannSolver::MaxLambda;



// -----------------------------------------------------------------------------
// GodunovOperator
// -----------------------------------------------------------------------------
GodunovOperator::ReconstructMethod GodunovOperator::reconstruct_method =
  GodunovOperator::RECONSTRUCT_PLM;
GodunovOperator::FluxSplittingMethod GodunovOperator::fluxsplit_method =
  GodunovOperator::FLUXSPLIT_LOCAL_LAX_FRIEDRICHS;

GodunovOperator::GodunovOperator()
  : emergency_reset(0),
    reset_density(1.0),
    reset_pressure(1.0) { }

void GodunovOperator::prepare_integration()
{
  std::vector<int> N = Mara->domain->aug_shape();
  NQ = Mara->domain->get_Nq();
  ND = N.size();

  N.insert(N.begin(), Mara->domain->get_Ng());
  for (int i=ND; i<=4; ++i) {
    N.push_back(1);
  }

  stride[0] = N[1]*N[2]*N[3]*NQ;
  stride[1] =      N[2]*N[3]*NQ;
  stride[2] =           N[3]*NQ;
  stride[3] =                NQ;

  std::vector<double> x0 = Mara->domain->get_x0();
  std::vector<double> x1 = Mara->domain->get_x1();

  x0.insert(x0.begin(), 0.0);
  x1.insert(x1.begin(), 0.0);

  for (int i=ND; i<=4; ++i) {
    x0.push_back(0.0);
    x1.push_back(1.0);
  }

  dx = (x1[1]-x0[1]) / (N[1]-2*N[0]);
  dy = (x1[2]-x0[2]) / (N[2]-2*N[0]);
  dz = (x1[3]-x0[3]) / (N[3]-2*N[0]);
}

int GodunovOperator::PrimToCons(const std::valarray<double> &P, std::valarray<double> &U)
{
  this->prepare_integration();
  Mara->boundary->SetToReceivePrimitive();
  Mara->boundary->ApplyBoundaries(const_cast<std::valarray<double> &>(P));
  Mara->boundary->SetToReceiveConserved();

  for (int i=0; i<stride[0]; i+=NQ) {
    Mara->fluid->PrimToCons(&P[i], &U[i]);
  }
  return 0;
}

int GodunovOperator::ConsToPrim(const std::valarray<double> &U, std::valarray<double> &P)
{
  this->prepare_integration();
  if (&P != &Mara->PrimitiveArray) P = Mara->PrimitiveArray; // don't copy it to itself

  Mara->boundary->ApplyBoundaries(const_cast<std::valarray<double> &>(U));
  int ttl_error=0;

  for (int i=0; i<stride[0]; i+=NQ) {
    int error = (Mara->fluid->ConsToPrim(&U[i], &P[i]) != 0);
    if (error) {
      int N[3];
      absolute_index_to_3d(i/NQ, N);
      if (1) {
        fprintf(stderr, "[GodunovOperator::ConsToPrim] "
                "recording c2p error at position [%f %f %f]\n"
                "(local zone [%d %d %d])\n",
                Mara->domain->x_at(N[0]),
                Mara->domain->y_at(N[1]),
                Mara->domain->z_at(N[2]), N[0], N[1], N[2]);
      }
    }
    ttl_error += error;
    Mara->FailureMask[i/NQ] = error;
  }

  /* ---------------------------------------------------------------------------
   *                          EMERGENCY STATE REPAIR
   *
   * NOTES:
   *
   * + Zones with the fail flag set will be get their primitive data replaced
   *   with zero velocity, magnetic field of the conserved state, and
   *   user-specified pressure and density.
   *
   * + New primitive states are then converted back to conservered ones at the
   *   end of the step.
   *
   * ---------------------------------------------------------------------------
   */
  if (ttl_error && emergency_reset) {

    std::valarray<double>  P_fixed = P;
    std::valarray<double> &U_fixed = const_cast<std::valarray<double>&>(U);

    for (int i=0; i<stride[0]; i+=NQ) {
      if (Mara->FailureMask[i/NQ]) { /* this is a bad zone */
	P_fixed[i + 0] = this->reset_density;
	P_fixed[i + 1] = this->reset_pressure;
	P_fixed[i + 2] = 0.0; // vx
	P_fixed[i + 3] = 0.0; // vy 
	P_fixed[i + 4] = 0.0; // vz
	if (NQ == 8) { /* MHD */
	  P_fixed[i + 5] = U[i + 5]; // Bx
	  P_fixed[i + 6] = U[i + 6]; // By
	  P_fixed[i + 7] = U[i + 7]; // Bz
	}
      }
    }

    Mara->FailureMask = 0;
    ttl_error = 0;
    P = P_fixed;

    for (int i=0; i<stride[0]; i+=NQ) {
      int perr = Mara->fluid->PrimToCons(&P_fixed[i], &U_fixed[i]);
      if (perr) {
	fprintf(stderr, "[Mara] FATAL: could not recover a primitive state\n");
	exit(1);
      }
    }
  }
  /* ---------------------------------------------------------------------------
   *                       END EMERGENCY STATE REPAIR
   * ---------------------------------------------------------------------------
   */
  return Mara_mpi_int_sum(ttl_error);
}

std::valarray<double> GodunovOperator::
LaxDiffusion(const std::valarray<double> &U, double r)
{
  this->prepare_integration();

  const std::valarray<int> &FM = Mara->FailureMask;
  std::valarray<double> L(U.size());
  Mara->boundary->ApplyBoundaries(const_cast<std::valarray<double>&>(U));

  if (ND == 1) {
    const int Ni=stride[0],sx=stride[1];
    std::valarray<double> Fiph(U.size());

    for (int i=0; i<Ni-sx; i+=NQ) {
      if (FM[i/NQ] || FM[(i+sx)/NQ]) {
        for (int j=0; j<NQ; ++j) {
          Fiph[i+j] = -0.5*r*(U[i+j+sx] - U[i+j])*dx;
        }
      }
    }

    for (int i=sx; i<Ni-sx; ++i) {
      L[i] = -((Fiph[i] - Fiph[i-sx])/dx);
    }
  }
  else if (ND == 2) {
    const int Ni=stride[0],sx=stride[1],sy=stride[2];
    std::valarray<double> Fiph(U.size()), Giph(U.size());

    for (int i=0; i<Ni-sx; i+=NQ) {

      if (FM[i/NQ] || FM[(i+sx)/NQ]) {
        for (int j=0; j<NQ; ++j) {
          Fiph[i+j] = -0.25*r*(U[i+j+sx] - U[i+j])*dx;
        }
      }
      if (FM[i/NQ] || FM[(i+sy)/NQ]) {
        for (int j=0; j<NQ; ++j) {
          Giph[i+j] = -0.25*r*(U[i+j+sy] - U[i+j])*dy;
        }
      }
    }
    Mara->fluid->ConstrainedTransport2d(&Fiph[0], &Giph[0], stride);

    for (int i=sx; i<Ni-sx; ++i) {
      L[i] = -((Fiph[i] - Fiph[i-sx])/dx +
               (Giph[i] - Giph[i-sy])/dy);
    }
  }
  else if (ND == 3) {
    const int Ni=stride[0],sx=stride[1],sy=stride[2],sz=stride[3];
    std::valarray<double> Fiph(U.size()), Giph(U.size()), Hiph(U.size());

    for (int i=0; i<Ni-sx; i+=NQ) {
      if (FM[i/NQ] || FM[(i+sx)/NQ]) {
        for (int j=0; j<NQ; ++j) {
          Fiph[i+j] = -(1.0/6.0)*r*(U[i+j+sx] - U[i+j])*dx;
        }
      }
      if (FM[i/NQ] || FM[(i+sy)/NQ]) {
        for (int j=0; j<NQ; ++j) {
          Giph[i+j] = -(1.0/6.0)*r*(U[i+j+sy] - U[i+j])*dy;
        }
      }
      if (FM[i/NQ] || FM[(i+sz)/NQ]) {
        for (int j=0; j<NQ; ++j) {
          Hiph[i+j] = -(1.0/6.0)*r*(U[i+j+sz] - U[i+j])*dz;
        }
      }
    }
    Mara->fluid->ConstrainedTransport3d(&Fiph[0], &Giph[0], &Hiph[0], stride);

    for (int i=sx; i<Ni-sx; ++i) {
      L[i] = -((Fiph[i] - Fiph[i-sx])/dx +
               (Giph[i] - Giph[i-sy])/dy +
               (Hiph[i] - Hiph[i-sz])/dz);
    }
  }
  return L;
}

void GodunovOperator::
AddFluxSourceTerms(double *Fiph, double *Giph, double *Hiph)
{
  /*
    if (Mara->fluxsrc == NULL) return;
    this->prepare_integration();

    for (int i=0; i<stride[0]; i+=NQ) {
    int N[3];
    absolute_index_to_3d(i/NQ, N);

    double x[3] = {
    Mara->domain->x_at(N[0]),
    Mara->domain->y_at(N[1]),
    Mara->domain->z_at(N[2]) };

    double xx[3] = {x[0], x[1], x[2]};
    double xy[3] = {x[0], x[1], x[2]};
    double xz[3] = {x[0], x[1], x[2]};

    xx[0] += 0.5*dx;
    xy[1] += 0.5*dy;
    xz[2] += 0.5*dz;

    Mara->fluxsrc->AddIntercellFlux(xx, 1, &Fiph[i]);
    Mara->fluxsrc->AddIntercellFlux(xy, 2, &Giph[i]);
    Mara->fluxsrc->AddIntercellFlux(xz, 3, &Hiph[i]);
    }
  */
}

int GodunovOperator::absolute_index_to_3d(const int &m, int ind[3])
{
  const int Nd = HydroModule::Mara->domain->get_Nd();
  const std::vector<int> &N = HydroModule::Mara->domain->aug_shape();

  if (Nd == 1) {
    const int i = m;
    ind[0] = i;
    ind[1] = 0;
    ind[2] = 0;
  }
  else if (Nd == 2) {
    const int i = (m         ) / (N[1]);
    const int j = (m - i*N[1]) / (1);
    ind[0] = i;
    ind[1] = j;
    ind[2] = 0;
  }
  else if (Nd == 3) {
    const int i = (m                       ) / (N[1]*N[2]);
    const int j = (m - i*N[1]*N[2]         ) / (N[2]);
    const int k = (m - i*N[1]*N[2] - j*N[2]) / (1);
    ind[0] = i;
    ind[1] = j;
    ind[2] = k;
  }

  return 0;
}
