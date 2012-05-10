

#ifndef __MethodOfLinesSplit_HEADER__
#define __MethodOfLinesSplit_HEADER__

#include "hydro.hpp"

class MethodOfLinesSplit : public GodunovOperator
{
private:
  void intercell_flux_sweep(const double *P, double *F, int dim);

  void drive_sweeps_1d(const double *P, double *L);
  void drive_sweeps_2d(const double *P, double *L);
  void drive_sweeps_3d(const double *P, double *L);
  void DriveSweeps(const std::valarray<double> &P, std::valarray<double> &L);
public:
  std::valarray<double> dUdt(const std::valarray<double> &Uin);
} ;

#endif // __MethodOfLinesSplit_HEADER__
