
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
  Mara->boundary->ApplyBoundaries(const_cast<std::valarray<double> &>(P));

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
    ttl_error += error;
    Mara->FailureMask[i/NQ] = error;
  }

  return Mara_mpi_int_sum(ttl_error);
}

std::valarray<double> GodunovOperator::LaxDiffusion(const std::valarray<double> &U, double r)
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


