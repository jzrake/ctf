

#ifndef __WenoSplit_HEADER__
#define __WenoSplit_HEADER__

#include <vector>
#include "hydro.hpp"


class WenoSplit : public GodunovOperator, public RiemannSolver
{
private:
  std::valarray<double> Uglb, Pglb, Lglb;
  std::valarray<double> Fiph, Giph, Hiph;

  void intercell_flux_sweep(const double *U, const double *P,
			    const double *F, const double *A,
			    double *Fiph, int dim);
  void drive_single_sweep(const double *Ug, const double *Pg,
			  double *Fiph_g, int dim);
  void drive_sweeps_1d();
  void drive_sweeps_2d();
  void drive_sweeps_3d();

public:
  std::valarray<double> dUdt(const std::valarray<double> &Uin);
  int IntercellFlux(const double *pl, const double *pr, double *U,
		    double *F, double s, int dim) { return 0; }
} ;

#endif // __WenoSplit_HEADER__
