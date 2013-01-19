/*------------------------------------------------------------------------------
 * FILE: advect.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#ifndef __ScalarAdvection_HEADER__
#define __ScalarAdvection_HEADER__

#include "hydro.hpp"

class ScalarAdvection : public FluidEquations
{
private:
  double WaveSpeed;

public:
  ScalarAdvection(double WaveSpeed);
  virtual ~ScalarAdvection() { }
  virtual int ConsToPrim(const double *U, double *P) const;
  virtual int PrimToCons(const double *P, double *U) const;

  void FluxAndEigenvalues(const double *U,
			  const double *P, double *F,
			  double *ap, double *am, int dimension) const;
  void Eigensystem(const double *U, const double *P,
		   double *L, double *R, double *lam, int dim) const;

  int GetNq() const { return 1; }
  std::vector<std::string> GetPrimNames() const;
  std::string PrintPrim(const double *P) const;
  std::string PrintCons(const double *P) const;
} ;

#endif // __ScalarAdvection_HEADER__
