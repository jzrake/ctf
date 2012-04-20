
#ifndef __EquationOfStates_HEADER__
#define __EquationOfStates_HEADER__

#include "hydro.hpp"



// -----------------------------------------------------------------------------
// Physics constants (gram, centimeter, second, kelvin)
// -----------------------------------------------------------------------------
#define LIGHT_SPEED         2.99792458000e+10
#define PROTON_MASS         1.67262177774e-24
#define BOLTZMANN_CONSTANT  1.38064881300e-16
#define MEV_TO_ERG          1.60217648740e-06
#define FM3_TO_CM3          1.00000000000e-39
#define EFFECTIVELY_ZERO    1.00000000000e-16




class PhysicalUnitsFactory_cgs
{
public:
  PhysicalUnits BuildPhysicalUnits();
} ;
class PhysicalUnitsFactory_MeV_gm_s
{
public:
  PhysicalUnits BuildPhysicalUnits();
} ;
class PhysicalUnitsFactory_MeV_c_ms
{
public:
  PhysicalUnits BuildPhysicalUnits();
} ;
class PhysicalUnitsFactory_erg_c_ms
{
public:
  PhysicalUnits BuildPhysicalUnits();
} ;
class PhysicalUnitsFactory_neutron_star
{
public:
  PhysicalUnits BuildPhysicalUnits();
} ;



class AdiabaticEos : public EquationOfState
{
public:
  AdiabaticEos(double Gamma);
  ~AdiabaticEos() { }

  double Pressure      (double D, double T) const;
  double Internal      (double D, double T) const;
  double Entropy       (double D, double T) const;
  double SoundSpeed2Nr (double D, double T) const;
  double SoundSpeed2Sr (double D, double T) const;
  double Temperature_u (double D, double p) const;
  double Temperature_p (double D, double p) const;
  double TemperatureMeV(double D, double p) const;
  double TemperatureArb(double D, double T_MeV) const;

  double Derivatives_u(double D, double T, double *J) const;
  double Derivatives_p(double D, double T, double *J) const;

  const double Gamma;
} ;

class ThermalBarotropicEos : public EquationOfState
{
private:
  const double GammaT;
  const double GammaB;
  const double Kappa;

public:
  ThermalBarotropicEos(double GammaT, double GammaB, double Kappa);
  ~ThermalBarotropicEos() { }

  double Pressure      (double D, double T) const;
  double Internal      (double D, double T) const;
  double Entropy       (double D, double T) const;
  double SoundSpeed2Nr (double D, double T) const;
  double SoundSpeed2Sr (double D, double T) const;
  double Temperature_u (double D, double p) const;
  double Temperature_p (double D, double p) const;
  double TemperatureMeV(double D, double p) const;
  double TemperatureArb(double D, double T_MeV) const;

  double Derivatives_u(double D, double T, double *J) const;
  double Derivatives_p(double D, double T, double *J) const;
} ;



struct TabulatedEos
{
  std::vector<double> logD_values;
  std::vector<double> logT_values;
  std::vector<double> EOS_p;
  std::vector<double> EOS_s;
  std::vector<double> EOS_u;
} ;

class ShenTabulatedNuclearEos : public EquationOfState
{

public:
  static bool verbose;
  static TabulatedEos LoadTable(const char *fname, double YpExtract,
                                const double *TempRange, const double *DensRange);

  ShenTabulatedNuclearEos(const TabulatedEos &tab);
  ~ShenTabulatedNuclearEos() { }

  // Public interface
  // ---------------------------------------------------------------------------

  double Pressure      (double D, double T) const;
  double Internal      (double D, double T) const;
  double Entropy       (double D, double T) const;
  double SoundSpeed2Nr (double D, double T) const;
  double SoundSpeed2Sr (double D, double T) const;
  double Temperature_u (double D, double u) const;
  double Temperature_p (double D, double p) const;
  double TemperatureMeV(double D, double p) const;
  double TemperatureArb(double D, double T_MeV) const;

  double Derivatives_u(double D, double T, double *J) const;
  double Derivatives_p(double D, double T, double *J) const;

  double DensLower() const; // returns lower/upper bound on D in code units
  double DensUpper() const;

  double TempLower() const; // returns lower/upper bound on T in code units
  double TempUpper() const;



  // Unit self-tests
  // ---------------------------------------------------------------------------
  void self_test_derivatives();
  void self_test_inversion();
  void self_test_interpolation();
  void self_test_soundspeed();



  // Exceptions
  // ---------------------------------------------------------------------------
  class UnableToLoadTable : public std::exception {
  public: virtual const char *what() const throw() {
    return "The ASCII table for Shen nuclear EOS could not be loaded."; } } ;

  class SampledOutOfRangeDensity : public std::exception {
  public: virtual const char *what() const throw() {
    return "Sampled the density out of the table's domain."; } } ;

  class SampledOutOfRangeTemperature : public std::exception {
  public: virtual const char *what() const throw() {
    return "Sampled the temperature out of the table's domain."; } } ;

  class SampledNegativePressure : public std::exception {
  public: virtual const char *what() const throw() {
    return "Table contains negative pressure entries for Yp>=0.17. Got one."; } } ;

  class FailedInverseLookupTemperature : public std::exception {
  public: virtual const char *what() const throw() {
    return "Could not find a T for input density and total internal energy."; } } ;

private:

  // Private interface
  // ---------------------------------------------------------------------------

  enum VariableIndex { ilogD, inB, iYp, iF, iEint, iS, iA, iZ, iMN,
                       iXn, iXp, iXa, iXA, ip, iun, iup, iML, iXL } ;


  //  void load_from_table();
  void tabulate_derivatives();
  int find_upper_index_D(double logD) const;
  int find_upper_index_T(double logT) const;
  double sample_EOS(const std::vector<double> &EOS,
                    double logD, double logT, double *J=NULL) const;
  double approx_EOS(const std::vector<double> &EOS,
                    double logD, double logT, double *J=NULL) const;
  double tabled_EOS(const std::vector<double> &EOS,
                    double logD, double logT, double *J=NULL) const;
  double inverse_lookup_T(const std::vector<double> &EOS,
                          double logD, double logF) const;


  // Private member data
  // ---------------------------------------------------------------------------
  std::vector<double> logD_values; // in code units
  std::vector<double> logT_values;


  // All values are in log10
  // ---------------------------------------------------------------------------
  std::vector<double> EOS_p; // gas pressure    (MeV/fm^3)
  std::vector<double> EOS_s; // entropy per particle (kB)
  std::vector<double> EOS_u; // energy density  (MeV/fm^3) not including rest mass

  std::vector<double> EOS_cs2; // sound speed squared
} ;




class GenericTabulatedEos : public EquationOfState
{

public:
  static bool verbose;

  GenericTabulatedEos(std::vector<double> &D_values,
		      std::vector<double> &T_values,
		      std::vector<double> &p,
                      std::vector<double> &u,
                      std::vector<double> &c);
  ~GenericTabulatedEos() { }

  // Public interface
  // ---------------------------------------------------------------------------
  double Pressure      (double D, double T) const;
  double Internal      (double D, double T) const;
  double Entropy       (double D, double T) const;
  double SoundSpeed2Nr (double D, double T) const;
  double SoundSpeed2Sr (double D, double T) const;
  double Temperature_u (double D, double u) const;
  double Temperature_p (double D, double p) const;
  double TemperatureMeV(double D, double p) const;
  double TemperatureArb(double D, double T_MeV) const;

  double Derivatives_u(double D, double T, double *J) const;
  double Derivatives_p(double D, double T, double *J) const;

  class SampledOutOfRangeDensity : public std::exception {
  public: virtual const char *what() const throw() {
    return "Sampled the density out of the table's domain."; } } ;

  class SampledOutOfRangeTemperature : public std::exception {
  public: virtual const char *what() const throw() {
    return "Sampled the temperature out of the table's domain."; } } ;

private:

  // Private interface
  // ---------------------------------------------------------------------------
  int find_upper_index_D(double D) const;
  int find_upper_index_T(double T) const;
  double sample_EOS(const std::vector<double> &EOS,
                    double D, double T, double *J=NULL) const;
  double approx_EOS(const std::vector<double> &EOS,
                    double D, double T, double *J=NULL) const;
  double tabled_EOS(const std::vector<double> &EOS,
                    double D, double T, double *J=NULL) const;
  double inverse_lookup_T(const std::vector<double> &EOS,
                          double D, double F) const;


  // Private member data
  // ---------------------------------------------------------------------------
  const std::vector<double> D_values; // in code units
  const std::vector<double> T_values;

  // EOS variables
  // ---------------------------------------------------------------------------
  const std::vector<double> EOS_p; // gas pressure    (MeV/fm^3)
  const std::vector<double> EOS_u; // energy density  (MeV/fm^3) no rest mass
  const std::vector<double> EOS_c; // sound speed     (units of light-speed)
} ;

#endif // __EquationOfStates_HEADER__
