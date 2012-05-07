/*------------------------------------------------------------------------------
 * FILE: srhd.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#ifndef __AdiabaticIdealSrhd_HEADER__
#define __AdiabaticIdealSrhd_HEADER__

#include "hydro.hpp"

class AdiabaticIdealSrhd : public FluidEquations
{
public:
  enum { ddd, tau, Sx, Sy, Sz }; // Conserved
  enum { rho, pre, vx, vy, vz }; // Primitive
  const double AdiabaticGamma;

public:
  AdiabaticIdealSrhd(double AdiabaticGamma);
  virtual ~AdiabaticIdealSrhd() { }
  virtual int ConsToPrim(const double *U, double *P) const;
  virtual int PrimToCons(const double *P, double *U) const;

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

#endif // __AdiabaticIdealSrhd_HEADER__
