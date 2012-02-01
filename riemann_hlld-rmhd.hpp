
/*------------------------------------------------------------------------------
 * FILE: riemann_hlld-rmhd.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * PURPOSE: Obtain the intercell flux between two constant states using a 5-wave
 *   approximation.
 *
 * REFERENCES: Mignone, Ugliano, & Bodo (2008) A ﬁve-wave Harten–Lax–van Leer
 *   Riemann solver for relativistic magnetohydrodynamics
 *
 *------------------------------------------------------------------------------
 */

#ifndef __HlldRmhdRiemannSolver_HEADER__
#define __HlldmhdRiemannSolver_HEADER__

#include "rmhd.hpp"

class HlldRmhdRiemannSolver : public RiemannSolver
{
public:
  static int NumberOfSuccess;
  static int NumberOfFailure;

private:
  enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
  enum { rho, pre, vx, vy, vz };             // Primitive

public:
  int IntercellFlux(const double *Pl, const double *Pr, double *U,
                    double *F, double s, int dim);
} ;

#endif // __HlldRmhdRiemannSolver_HEADER__
