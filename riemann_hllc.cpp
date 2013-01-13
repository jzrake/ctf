

#include "riemann_hllc.hpp"

int HllcRiemannSolver::IntercellFlux(const double *pl, const double *pr, double *U,
                                     double *F, double s, int dim)
{
  FluidEquations &fluid = *Mara->fluid;
  RiemannSolver *riemann;

  if      (typeid(fluid) == typeid(AdiabaticIdealEulers)) {
    riemann = new HllcEulersRiemannSolver;
  }
  else if (typeid(fluid) == typeid(AdiabaticIdealSrhd)) {
    riemann = new HllcSrhdRiemannSolver;
  }
  else if (typeid(fluid) == typeid(AdiabaticIdealRmhd)) {
    riemann = new HllcRmhdRiemannSolver;
  }
  else {
    throw std::bad_typeid();
  }

  int error = riemann->IntercellFlux(pl, pr, U, F, s, dim);
  delete riemann;
  return error;
}
