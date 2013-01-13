/*------------------------------------------------------------------------------
 * FILE: euler.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#ifndef __AdiabaticIdealEulers_HEADER__
#define __AdiabaticIdealEulers_HEADER__

#include "hydro.hpp"

class AdiabaticIdealEulers : public FluidEquations
{
public:
  enum { rho, nrg, px, py, pz }; // Conserved
  enum { RHO, pre, vx, vy, vz }; // Primitive

public:
  AdiabaticIdealEulers() { }
  ~AdiabaticIdealEulers() { }
  int ConsToPrim(const double *U, double *P) const;
  int PrimToCons(const double *P, double *U) const;

  void FluxAndEigenvalues(const double *U,
			  const double *P, double *F,
			  double *ap, double *am, int dimension) const;
  void Eigensystem(const double *U, const double *P,
		   double *L, double *R, double *lam, int dim) const;

  int GetNq() const { return 5; }
  std::vector<std::string> GetPrimNames() const;
  std::string PrintPrim(const double *P) const;
  std::string PrintCons(const double *P) const;
} ;

#endif // __AdiabaticIdealEulers_HEADER__
