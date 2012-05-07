

#ifndef __PlmMethodOfLinesSplit_HEADER__
#define __PlmMethodOfLinesSplit_HEADER__

#include "hydro.hpp"

class PlmMethodOfLinesSplit : public GodunovOperator
{
private:
  void reconstruct_plm(const double *P0, double *Pl, double *Pr, int S);
  void intercell_flux_sweep(const double *P, double *F, int dim);

  void drive_sweeps_1d(const double *P, double *L);
  void drive_sweeps_2d(const double *P, double *L);
  void drive_sweeps_3d(const double *P, double *L);
  void DriveSweeps(const std::valarray<double> &P, std::valarray<double> &L);

public:
  std::valarray<double> dUdt(const std::valarray<double> &Uin);
} ;

class Weno5RiemannMethodOfLinesSplit : public GodunovOperator
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
/*
class Weno5RamScheme : public GodunovOperator
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
*/

#endif // __PlmMethodOfLinesSplit_HEADER__
