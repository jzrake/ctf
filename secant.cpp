
/*------------------------------------------------------------------------------
 * FILE: secant.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION: Implements the SecantMethodSolver class
 *
 * REFERENCES: http://en.wikipedia.org/wiki/Secant_method
 *
 *
 *------------------------------------------------------------------------------
 */

#include <cstdio>
#include "config.h"
#include "secant.hpp"
#include "nrsolver.hpp"

int SecantMethodSolver::Solve(const EquationSystemBaseClass &eqn, double *x)
{
  const double TOL = 1e-15;
  double &x0 = *x;
  double y0, y1, y2, x1, x2;
  int WorseIteration = 0;

  eqn.Function(&x0, &y0);
  x1 = x0 + 1e-6;

  if (fabs(y0) < TOL) {
    return 0;
  }
  else {
    x1 = x0 + 1e-6;
    eqn.Function(&x1, &y1);
  }

  Iterations = 0;
  while (fabs(y1) > TOL) {

    x2 = x1 - y1 * (x1-x0) / (y1-y0);
    eqn.Function(&x2, &y2);

    x0 = x1;
    y0 = y1;

    x1 = x2;
    y1 = y2;

    WorseIteration += (y1 > y0);
    if (WorseIteration >= 10) {
      throw NewtonRaphesonSolver::NonConvergingSearch();
    }
    if (Mara_isnan_cxx(y2)) {
      throw NewtonRaphesonSolver::NanLocated();
    }

    if (++Iterations > 15) return 1;
  }
  x0 = x1;

  return 0;
}
