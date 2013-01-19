/*------------------------------------------------------------------------------
 * FILE: rmhd.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#ifndef __AdiabaticIdealRmhd_HEADER__
#define __AdiabaticIdealRmhd_HEADER__

#include "hydro.hpp"

class AdiabaticIdealRmhd : public FluidEquations
{
public:
  enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
  enum { rho, pre, vx, vy, vz };             // Primitive

public:
  AdiabaticIdealRmhd() { }
  ~AdiabaticIdealRmhd() { }
  int ConsToPrim(const double *U, double *P) const;
  int PrimToCons(const double *P, double *U) const;

  void FluxAndEigenvalues(const double *U,
                          const double *P, double *F,
                          double *ap, double *am, int dimension) const;

  void ConstrainedTransport2d(double *Fx, double *Fy,             int stride[4]) const;
  void ConstrainedTransport3d(double *Fx, double *Fy, double *Fz, int stride[4]) const;

  int GetNq() const { return 8; }
  std::vector<std::string> GetPrimNames() const;
  std::string PrintPrim(const double *P) const;
  std::string PrintCons(const double *P) const;

  int PrimCheck(const double *P) const;
  int ConsCheck(const double *U) const;
} ;

#endif // __AdiabaticIdealRmhd_HEADER__
