
/*------------------------------------------------------------------------------
 * FILE: nrsolver.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION: Implements the NewtonRaphesonSolver class
 *
 * REFERENCES:
 Numerical Recipes, Third Edition: Section 9.6
 *
 *------------------------------------------------------------------------------
 */

#include "config.h"
#include "nrsolver.hpp"
#include "matrix.h"


int NewtonRaphesonSolver::Solve(const EquationSystemBaseClass &eqn, double *x)
/*----------------------------------------------------------------------------
 * This method is based on the first order approximation to the function
 *
 *          f(x) = f(x0) + J(x0).(x-x0) + O(2)
 *
 * where J is the Jacobian matrix, J_ij = df_i/dx_j. Setting the left hand
 * side to zero, and ignoring O(2) we obtain
 *
 *          x = x0 - G(x0).f(x0)
 *
 * where G is the inverse of J. Iterations are carried out this way until the
 * root is reached within some tolerance.
 *----------------------------------------------------------------------------
 */
{
  int Rank = eqn.GetRank();
  int WorseIteration = 0;
  double Epsilon_last;

  double *dx = new double[Rank];
  double *y  = new double[Rank];
  double *J  = new double[Rank*Rank]; // Jacobian
  double *G  = new double[Rank*Rank]; // Jacobian inverse

  Epsilon = 0.0;
  Iterations = 0;
  do {
    Epsilon_last = Epsilon;
    Epsilon = 0.0;

    eqn.Function(x, y);
    eqn.Jacobian(x, J);

    matrix_inverse(J, G, Rank);
    matrix_vector_product(G, y, dx, Rank, Rank);

    Epsilon = eqn.GetErrorUpdateX(x, dx, y);
    WorseIteration += (Epsilon > Epsilon_last);

    if (WorseIteration >= 10 || Iterations == MaxIterations) {
      delete [] dx;
      delete [] y;
      delete [] J;
      delete [] G;
      throw NonConvergingSearch();
    }
    if (Mara_isnan_cxx(Epsilon)) {
      delete [] dx;
      delete [] y;
      delete [] J;
      delete [] G;
      throw NanLocated();
    }
  } while (Epsilon > MinEpsilon &&
	   ++Iterations <= MaxIterations);

  delete [] dx;
  delete [] y;
  delete [] J;
  delete [] G;

  return 0;
}
