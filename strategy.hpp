#ifndef __StrategyCollection_HEADER__
#define __StrategyCollection_HEADER__


#include "eulers.hpp"
#include "riemann_hll.hpp"
#include "riemann_hllc.hpp"
#include "riemann_hlld-rmhd.hpp"
#include "riemann_exact-eulers.hpp"
#include "plm-split.hpp"
#include "hybrid-split.hpp"
#include "ctu-hancock.hpp"
#include "runge-kutta.hpp"


class SafestSolutionStrategy : public SolutionStrategy
{
public:
  SafestSolutionStrategy() : SolutionStrategy("Safest") { }
  void BuildSimDriverComponents(const TestProblem &problem,
                                GodunovOperator *&deriv,
                                RiemannSolver *&riemann,
                                RungeKuttaIntegration *&RK) const
  {
    const PhysicalDomain &domain = problem.GetDomain();
    const BoundaryConditions &boundary = problem.GetBoundary();
    const FluidEquations &fluid = problem.GetFluid();

    riemann = new HllRiemannSolver(fluid);
    deriv   = new PlmMethodOfLinesSplit(domain, boundary, fluid, *riemann);
    RK      = new RungeKuttaShuOsherRk3;
  }
} ;
class HllcSolutionStrategy : public SolutionStrategy
{
public:
  HllcSolutionStrategy() : SolutionStrategy("Hllc") { }
  void BuildSimDriverComponents(const TestProblem &problem,
                                GodunovOperator *&deriv,
                                RiemannSolver *&riemann,
                                RungeKuttaIntegration *&RK) const
  {
    const PhysicalDomain &domain = problem.GetDomain();
    const BoundaryConditions &boundary = problem.GetBoundary();
    const FluidEquations &fluid = problem.GetFluid();

    riemann = new HllcRiemannSolver(fluid);
    deriv   = new PlmMethodOfLinesSplit(domain, boundary, fluid, *riemann);
    RK      = new RungeKuttaShuOsherRk3;
  }
} ;
class CtuHllSolutionStrategy : public SolutionStrategy
{
public:
  CtuHllSolutionStrategy() : SolutionStrategy("CtuHll") { }
  void BuildSimDriverComponents(const TestProblem &problem,
                                GodunovOperator *&deriv,
                                RiemannSolver *&riemann,
                                RungeKuttaIntegration *&RK) const
  {
    const PhysicalDomain &domain = problem.GetDomain();
    const BoundaryConditions &boundary = problem.GetBoundary();
    const FluidEquations &fluid = problem.GetFluid();

    riemann = new HllRiemannSolver(fluid);
    deriv   = new PlmCtuHancockOperator(domain, boundary, fluid, *riemann);
    RK      = new RungeKuttaSingleStep;
  }
} ;
class CtuHllcSolutionStrategy : public SolutionStrategy
{
public:
  CtuHllcSolutionStrategy() : SolutionStrategy("CtuHllc") { }
  void BuildSimDriverComponents(const TestProblem &problem,
                                GodunovOperator *&deriv,
                                RiemannSolver *&riemann,
                                RungeKuttaIntegration *&RK) const
  {
    const PhysicalDomain &domain = problem.GetDomain();
    const BoundaryConditions &boundary = problem.GetBoundary();
    const FluidEquations &fluid = problem.GetFluid();

    riemann = new HllcRiemannSolver(fluid);
    deriv   = new PlmCtuHancockOperator(domain, boundary, fluid, *riemann);
    RK      = new RungeKuttaSingleStep;
  }
} ;
class HlldRmhdRk3SolutionStrategy : public SolutionStrategy
{
public:
  HlldRmhdRk3SolutionStrategy() : SolutionStrategy("HlldRmhdRk3") { }
  void BuildSimDriverComponents(const TestProblem &problem,
                                GodunovOperator *&deriv,
                                RiemannSolver *&riemann,
                                RungeKuttaIntegration *&RK) const
  {
    const PhysicalDomain &domain = problem.GetDomain();
    const BoundaryConditions &boundary = problem.GetBoundary();
    const AdiabaticIdealRmhd &fluid =
      dynamic_cast<const AdiabaticIdealRmhd&>(problem.GetFluid());

    riemann = new HlldRmhdRiemannSolver(fluid);
    deriv   = new PlmMethodOfLinesSplit(domain, boundary, fluid, *riemann);
    RK      = new RungeKuttaShuOsherRk3;
  }
} ;
class HlldRmhdCtuSolutionStrategy : public SolutionStrategy
{
public:
  HlldRmhdCtuSolutionStrategy() : SolutionStrategy("HlldRmhdCtu") { }
  void BuildSimDriverComponents(const TestProblem &problem,
                                GodunovOperator *&deriv,
                                RiemannSolver *&riemann,
                                RungeKuttaIntegration *&RK) const
  {
    const PhysicalDomain &domain = problem.GetDomain();
    const BoundaryConditions &boundary = problem.GetBoundary();
    const AdiabaticIdealRmhd &fluid =
      dynamic_cast<const AdiabaticIdealRmhd&>(problem.GetFluid());

    riemann = new HlldRmhdRiemannSolver(fluid);
    deriv   = new PlmCtuHancockOperator(domain, boundary, fluid, *riemann);
    RK      = new RungeKuttaSingleStep;
  }
} ;
class ExactEulerSplitSolutionStrategy : public SolutionStrategy
{
public:
  ExactEulerSplitSolutionStrategy() : SolutionStrategy("ExactEulerSplit") { }
  void BuildSimDriverComponents(const TestProblem &problem,
                                GodunovOperator *&deriv,
                                RiemannSolver *&riemann,
                                RungeKuttaIntegration *&RK) const
  {
    const PhysicalDomain &domain = problem.GetDomain();
    const BoundaryConditions &boundary = problem.GetBoundary();
    const AdiabaticIdealEulers &fluid =
      dynamic_cast<const AdiabaticIdealEulers&>(problem.GetFluid());

    riemann = new ExactEulersRiemannSolver(fluid);
    deriv   = new PlmMethodOfLinesSplit(domain, boundary, fluid, *riemann);
    RK      = new RungeKuttaShuOsherRk3;
  }
} ;
class ExactEulerCtuSolutionStrategy : public SolutionStrategy
{
public:
  ExactEulerCtuSolutionStrategy() : SolutionStrategy("ExactEulerCtu") { }
  void BuildSimDriverComponents(const TestProblem &problem,
                                GodunovOperator *&deriv,
                                RiemannSolver *&riemann,
                                RungeKuttaIntegration *&RK) const
  {
    const PhysicalDomain &domain = problem.GetDomain();
    const BoundaryConditions &boundary = problem.GetBoundary();
    const AdiabaticIdealEulers &fluid =
      dynamic_cast<const AdiabaticIdealEulers&>(problem.GetFluid());

    riemann = new ExactEulersRiemannSolver(fluid);
    deriv   = new PlmCtuHancockOperator(domain, boundary, fluid, *riemann);
    RK      = new RungeKuttaSingleStep;
  }
} ;

class HybridSolutionStrategy : public SolutionStrategy
{
public:
  HybridSolutionStrategy() : SolutionStrategy("Hybrid") { }
  void BuildSimDriverComponents(const TestProblem &problem,
                                GodunovOperator *&deriv,
                                RiemannSolver *&riemann,
                                RungeKuttaIntegration *&RK) const
  {
    const PhysicalDomain &domain = problem.GetDomain();
    const BoundaryConditions &boundary = problem.GetBoundary();
    const FluidEquations &fluid = problem.GetFluid();

    riemann = NULL;
    deriv   = new HybridSplit(domain, boundary, fluid);
    RK      = new RungeKuttaShuOsherRk3;
  }
} ;

#endif // __StrategyCollection_HEADER__
