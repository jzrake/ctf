
/*------------------------------------------------------------------------------
 * FILE: nrsolver.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION: Contains a class for root-finding using basic Newton-Rapheson
 *
 * REFERENCES:
 Numerical Recipes, Third Edition: Section 9.6
 *
 *------------------------------------------------------------------------------
 */

#ifndef __NewtonRaphesonSolver_HEADER__
#define __NewtonRaphesonSolver_HEADER__

#include <cmath>
#include <exception>
#include "eqnbase.hpp"

class NewtonRaphesonSolver
{
private:
  int Iterations, MaxIterations;
  double Epsilon, MinEpsilon;

public:
  class Failure : public std::exception { } ;
  class SingularJacobian : public Failure
  {
  public:
    static const int ErrorCode = 1;
    virtual const char *what() const throw()
    {
      return "Matrix inversion encountered a singular Jacobian.";
    }
  } ;
  class NonConvergingSearch : public Failure
  {
  public:
    static const int ErrorCode = 2;
    virtual const char *what() const throw()
    {
      return "The root-finder is getting farther from the root.";
    }
  } ;
  class NanLocated : public Failure
  {
  public:
    static const int ErrorCode = 3;
    virtual const char *what() const throw()
    {
      return "The root-finder found a NAN somewhere.";
    }
  } ;

  NewtonRaphesonSolver(int MaxIterations, double MinEpsilon) :
    Iterations(0),
    MaxIterations(MaxIterations),
    Epsilon(0.0),
    MinEpsilon(MinEpsilon)  { }

  int Solve(const EquationSystemBaseClass &eqn, double *x);
  int GetIterations() { return Iterations; }
  double GetEpsilon() { return Epsilon; }
} ;

#endif // __NewtonRaphesonSolver_HEADER__
