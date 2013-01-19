
/*------------------------------------------------------------------------------
 * FILE: secant.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION: Contains a class for root-finding using secant method in 1d
 *
 * REFERENCES: http://en.wikipedia.org/wiki/Secant_method
 *
 *------------------------------------------------------------------------------
 */

#ifndef __SecantMethodSolver_HEADER__
#define __SecantMethodSolver_HEADER__

#include <cmath>
#include <exception>
#include "eqnbase.hpp"

class SecantMethodSolver
{
private:
  int Iterations, MaxIterations;
  double Epsilon, MinEpsilon;

public:
  class Failure : public std::exception { } ;
  SecantMethodSolver(int MaxIterations, double MinEpsilon) :
    Iterations(0),
    MaxIterations(MaxIterations),
    Epsilon(0.0),
    MinEpsilon(MinEpsilon)  { }

  int Solve(const EquationSystemBaseClass &eqn, double *x);
  int GetIterations() { return Iterations; }
  double GetEpsilon() { return Epsilon; }
} ;

#endif // __SecantMethodSolver_HEADER__
