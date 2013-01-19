
/*------------------------------------------------------------------------------
 * FILE: riemann_exact-eulers.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION:
 *
 * REFERENCES: E.F Toro (1999)
 *             Riemann Solvers and Numerical Methods for Fluid Dynamics
 *------------------------------------------------------------------------------
 */

#ifndef __ExactRiemannSolverEulers_HEADER__
#define __ExactRiemannSolverEulers_HEADER__

#include <iostream>
#include "nrsolver.hpp"
#include "eulers.hpp"


class ExactEulersRiemannSolver : public RiemannSolver
{
private:
  const AdiabaticIdealEulers &fluid;
  RiemannSolver *BackupRiemannSolver;
  int v1, v2, v3;

  enum { rho, nrg, px, py, pz }; // Conserved
  enum { RHO, pre, vx, vy, vz }; // Primitive

public:
  class NoPressureSolution : public std::exception
  {
  public:
    virtual const char *what() const throw()
    {
      return "Failed to find p* in Riemann solution.";
    }
  } ;
  ExactEulersRiemannSolver(const AdiabaticIdealEulers &fluid)
    : fluid(fluid), BackupRiemannSolver(NULL) { }
  ExactEulersRiemannSolver(const AdiabaticIdealEulers &fluid,
			   RiemannSolver *BackupRiemannSolver)
    : fluid(fluid), BackupRiemannSolver(BackupRiemannSolver) { }
  ~ExactEulersRiemannSolver() { }
  int IntercellFlux(const double *pl, const double *pr, double *U_out,
		    double *F, double s, int dim);
} ;

#endif // __ExactRiemannSolverEulers_HEADER__
