
/*------------------------------------------------------------------------------
 * FILE: hydro.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#ifndef __HydroBaseInterface_HEADER__
#define __HydroBaseInterface_HEADER__


#include <sstream>
#include <string>
#include <valarray>
#include <vector>
#include <typeinfo>


// -----------------------------------------------------------------------------
// Forward declarations of all classes defined in this header file
// -----------------------------------------------------------------------------
class MaraApplication;
class BoundaryConditions;
class CoolingModule;
class DrivingModule;
class StochasticVectorField;
class EquationOfState;
class FluidEquations;
class PhysicalDomain;
class RiemannSolver;
class GodunovOperator;
class RungeKuttaIntegration;
class PhysicalUnits;
// -----------------------------------------------------------------------------




class MaraApplication
{
private:

public:
  MaraApplication();
  ~MaraApplication();

  template <class T> T &GetDomain() {
    if (domain == NULL) throw std::bad_cast();
    return dynamic_cast<T&>(*domain);
  }

  template <class T> T &GetEos() {
    if (eos == NULL) throw std::bad_cast();
    return dynamic_cast<T&>(*eos);
  }

  template <class T> T &GetFluid() {
    if (fluid == NULL) throw std::bad_cast();
    return dynamic_cast<T&>(*fluid);
  }

  template <class T> T &GetUnits() {
    if (units == NULL) throw std::bad_cast();
    return dynamic_cast<T&>(*units);
  }

  PhysicalUnits         *units;
  PhysicalDomain        *domain;
  BoundaryConditions    *boundary;
  FluidEquations        *fluid;
  EquationOfState       *eos;
  GodunovOperator       *godunov;
  RiemannSolver         *riemann;
  RungeKuttaIntegration *advance;
  CoolingModule         *cooling;
  DrivingModule         *driving;

  std::valarray<double> PrimitiveArray;
  std::valarray<int> FailureMask;
} ;




class HydroModule
// -----------------------------------------------------------------------------
{
public:
  static MaraApplication *Mara;
} ;
class BoundaryConditions : public HydroModule
// -----------------------------------------------------------------------------
{
public:
  virtual ~BoundaryConditions() { }
  virtual void ApplyBoundaries(std::valarray<double> &U) const = 0;
} ;
class EquationOfState : public HydroModule
// -----------------------------------------------------------------------------
{
public:
  virtual ~EquationOfState() { }

  virtual double Pressure      (double D, double T) const = 0;
  virtual double Internal      (double D, double T) const = 0;
  virtual double Entropy       (double D, double T) const = 0;
  virtual double SoundSpeed2Nr (double D, double T) const = 0;
  virtual double SoundSpeed2Sr (double D, double T) const = 0;
  virtual double Temperature_u (double D, double u) const = 0;
  virtual double Temperature_p (double D, double p) const = 0;
  virtual double TemperatureMeV(double D, double p) const = 0;
  virtual double TemperatureArb(double D, double T_MeV) const = 0;

  virtual double Derivatives_u(double D, double T, double *J) const = 0;
  virtual double Derivatives_p(double D, double T, double *J) const = 0;

  virtual double DensLower() const { return 0.0; } // returns lower/upper bound on D in code units
  virtual double DensUpper() const { return 0.0; }

  virtual double TempLower() const { return 0.0; } // returns lower/upper bound on T in code units
  virtual double TempUpper() const { return 0.0; }
} ;
class FluidEquations : public HydroModule
// -----------------------------------------------------------------------------
{
public:
  virtual ~FluidEquations() { }
  virtual int PrimToCons(const double *P, double *U) const = 0;
  virtual int ConsToPrim(const double *U, double *P) const = 0;
  virtual void FluxAndEigenvalues(const double *U,
                                  const double *P, double *F,
                                  double *ap, double *am, int dim) const = 0;
  virtual void Eigensystem(const double *U, const double *P,
                           double *L, double *R, double *lam, int dim) const { }
  virtual void ConstrainedTransport2d(double *Fx, double *Fy,
                                      int stride[4]) const { }
  virtual void ConstrainedTransport3d(double *Fx, double *Fy, double *Fz,
                                      int stride[4]) const { }
  virtual std::vector<std::string> GetPrimNames() const = 0;
  virtual int GetNq() const = 0;
  virtual std::string PrintPrim(const double *P) const { return ""; }
  virtual std::string PrintCons(const double *P) const { return ""; }
  virtual int PrimCheck(const double *P) const { return 0; }
  virtual int ConsCheck(const double *U) const { return 0; }
} ;
class PhysicalDomain : public HydroModule
// -----------------------------------------------------------------------------
{
public:
  struct SubdomainSpecs
  {
    int *A_nint;
    int *L_ntot;
    int *L_strt;
    int *G_ntot;
    int *G_strt;
    int  n_dims;
    int  n_prim;
  } ;
protected:
  SubdomainSpecs Specs;
public:
  const SubdomainSpecs &GetSpecs() const { return Specs; }
  virtual ~PhysicalDomain() { }
  virtual double get_dx(int d) const = 0;
  virtual double get_min_dx()  const = 0;
  virtual int get_N(int d)     const = 0; // not including guard
  virtual int get_Ng()         const = 0; // guard zones
  virtual int get_Nq()         const = 0; // quantities per zone
  virtual int get_Nd()         const = 0; // number of dimensions
  virtual double x_at(int i)   const = 0;
  virtual double y_at(int j)   const = 0;
  virtual double z_at(int k)   const = 0;
  virtual std::vector<double> get_x0() const = 0;
  virtual std::vector<double> get_x1() const = 0;
  virtual std::vector<int> aug_shape() const = 0; // including guard
  virtual void Synchronize(std::valarray<double> &U) const = 0;
  virtual int SubgridRank() const = 0;
  virtual int SubgridSize() const = 0;
  virtual int GetNumberOfZones() const = 0; // including guard
  virtual int GetGlobalNumberOfZones() const = 0; // over all subgrids
  virtual int GetSubgridIndex(int i) const = 0;
  virtual int GetSubgridSizes(int i) const = 0;

  virtual const double *GetGlobalX0() const = 0;
  virtual const double *GetGlobalX1() const = 0;

  virtual const int *GetGlobalShape() const = 0;
  virtual const int *GetGlobalStart() const = 0;
  virtual const int *GetLocalShape() const = 0; // not including guard

  virtual int SubgridAtPosition(const double *r) const = 0;
  virtual int IndexAtPosition(const double *r, int d) const = 0; // d: 0,1,2 for x,y,z
} ;
class RiemannSolver : public HydroModule
// -----------------------------------------------------------------------------
{
protected:
  static double MaxLambda;
public:
  virtual ~RiemannSolver() { }
  static double GetMaxLambda();
  static void ResetMaxLambda();
  virtual int IntercellFlux(const double *pl, const double *pr, double *U,
                            double *F, double s, int dim) = 0;
} ;
class GodunovOperator : public HydroModule
// -----------------------------------------------------------------------------
{
protected:
  double dx,dy,dz;
  int stride[4],NQ,ND;

public:
  class ConsToPrimFailure : public std::exception
  {
  public:
    virtual const char *what() const throw()
    {
      return "Failed to invert vector of conserved quantities.";
    }
  } ;
  class IntermediateFailure : public std::exception
  {
  public:
    virtual const char *what() const throw()
    {
      return "The integration failed on an intermediate step.";
    }
  } ;
  virtual ~GodunovOperator() { }
  virtual std::valarray<double> dUdt(const std::valarray<double> &Uin) = 0;
  virtual std::valarray<double> LaxDiffusion(const std::valarray<double> &U, double r);
  virtual int PrimToCons(const std::valarray<double> &P, std::valarray<double> &U);
  virtual int ConsToPrim(const std::valarray<double> &U, std::valarray<double> &P);
  virtual void SetTimeStepDt(double dt) { };
  virtual void SetPlmTheta(double plm) { }
  virtual void SetSafetyLevel(int level) { }
  const int *GetStride() const { return stride; }
  int GetNq() const { return NQ; }
  int GetNd() const { return ND; }

protected:
  void prepare_integration();
} ;
class RungeKuttaIntegration : public HydroModule
// -----------------------------------------------------------------------------
{
public:
  virtual ~RungeKuttaIntegration() { }
  virtual void AdvanceState(std::valarray<double> &U, double dt) const = 0;
} ;
class StochasticVectorField
{
public:
  virtual ~StochasticVectorField() { }
  virtual void AdvanceField(double dt) = 0;
  virtual std::vector<double> SampleField(double x, double y, double z) const = 0;
  virtual void Serialize(std::ostream &stream) const = 0;
} ;
class DrivingModule : public HydroModule
// -----------------------------------------------------------------------------
{
public:
  virtual ~DrivingModule() { }
  virtual void Drive(std::valarray<double> &P, double dt) = 0;
  virtual StochasticVectorField *GetField() = 0;
  virtual void ResampleField() = 0;
  virtual double AveragePowerInFields() const = 0;
} ;
class CoolingModule : public HydroModule
// -----------------------------------------------------------------------------
{
public:
  virtual ~CoolingModule() { }
  virtual void Cool(std::valarray<double> &P, double dt) = 0;
  virtual double EnergyRemoved() = 0;
} ;
class PhysicalUnits
// -----------------------------------------------------------------------------
{
public:
  PhysicalUnits(double c, double g, double s);
  ~PhysicalUnits() { }

  double Gram() const;
  double Centimeter() const;
  double Second() const;

  double Erg() const;
  double MeV() const;
  double Gauss() const;
  double Velocity() const;
  double GramsPerCubicCentimeter() const;
  double ErgPerCubicCentimeter() const;
  double MeVPerCubicCentimeter() const;
  double MeVPerCubicFemtometer() const;
  double BoltzmannConstant() const;
  double LightSpeed() const;
  double ProtonMass() const;

  void PrintUnits();

private:
  const double Length;
  const double Mass;
  const double Time;
} ;



#endif // __HydroBaseInterface_HEADER__
