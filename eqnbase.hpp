
/*------------------------------------------------------------------------------
 * FILE: eqnbase.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION: Base class describing a system of equations
 *
 * REFERENCES:
 *------------------------------------------------------------------------------
 */

#ifndef __EquationSystemBaseClass_HEADER__
#define __EquationSystemBaseClass_HEADER__

class EquationSystemBaseClass
{
private:
  const int Rank;

public:
  EquationSystemBaseClass(int Rank) : Rank(Rank) { }
  virtual ~EquationSystemBaseClass() { }

  virtual void Function(const double *x, double *y) const = 0;
  virtual void Jacobian(const double *x, double *J) const = 0;
  virtual double GetErrorUpdateX(double *x, const double *dx, const double *y) const
  {
    double error = 0;
    for (int i=0; i<Rank; ++i) {
      x[i] -= dx[i];
      error += y[i] > 0 ? y[i] : -y[i];
    }
    return error;
  }
  const int &GetRank() const { return Rank; }
} ;

#endif // __EquationSystemBaseClass_HEADER__
