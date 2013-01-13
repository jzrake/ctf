
/*------------------------------------------------------------------------------
 * FILE: shen.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * DESCRIPTION:
 *
 * This code provides a class which reads the 3d tabulated nuclear equation of
 * state developed by Shen (2011). It selects a single, global value of the
 * proton fraction Yp, and reads a 2d slice of the table in logD := log10(rho)
 * and logT := log10(T).
 *
 * REFERENCES:
 *
 * http://user.numazu-ct.ac.jp/~sumi/eos/table2/guide_EOS3.pdf
 *
 * CAUTION:
 *
 * (1) For densities greater than 10^14 gm/cm^3 the pressure at constant density
 *     can become multi-valued in T, making the inverse lookup impossible.
 *
 * (2) The sound speed can become imaginary for high densities and low
 *     temperatures.
 *
 * CONVENTIONS:
 *
 * The unit conversion facility works by providing member functions such as
 * units.GramsPerCubicCentimeter(). If the variable D_phys contains the density
 * in gm/cm^3 units, then D_code := D_phys * units.GramsPerCubicCentimeter().
 *
 * All public member functions receive the density (D) in code units defined
 * through the PhysicalUnits class instance, and the log10 of temperature (T) in
 * MeV. They also return the pressure (p), internal energy (u) and sound speed
 * (squared, cs2) in code units, but the log10 of temperature in MeV. Private
 * member functions all receive log10 of the quantities in the physical units
 * used by Shen's table, which are MeV for temperature, gm/cm^3 for density, and
 * MeV/fm^3 for pressure and internal energy density.
 *
 *------------------------------------------------------------------------------
 */


#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include "eos.hpp"


// Alias to the global Mara application units instance
#define units (*Mara->units)

// Set this to 0.5 or something to test off-table lookups (but expect warnings)
#define OFF_TABLE_VAL (0)

// Enable to use direct lookup on D,T (assumes equal spacing)
#define NOSEARCH_DT

bool ShenTabulatedNuclearEos::verbose = false;




double ShenTabulatedNuclearEos::DensLower() const
{
  return pow(10.0, logD_values.front()) * units.GramsPerCubicCentimeter();
}
double ShenTabulatedNuclearEos::DensUpper() const
{
  return pow(10.0, logD_values.back()) * units.GramsPerCubicCentimeter();
}
double ShenTabulatedNuclearEos::TempLower() const
{
  return logT_values.front();
}
double ShenTabulatedNuclearEos::TempUpper() const
{
  return logT_values.back();
}



ShenTabulatedNuclearEos::ShenTabulatedNuclearEos(const TabulatedEos &tab)
  : logD_values(tab.logD_values),
    logT_values(tab.logT_values),
    EOS_p(tab.EOS_p),
    EOS_s(tab.EOS_s),
    EOS_u(tab.EOS_u)
{
  this->tabulate_derivatives();
}

void ShenTabulatedNuclearEos::tabulate_derivatives()
{
  const size_t ND = logD_values.size();
  const size_t NT = logT_values.size();

  EOS_cs2.resize(ND*NT);

  for (size_t i=0; i<ND; ++i) {
    for (size_t j=0; j<NT; ++j) {

      const int im1 = (i !=    0) ? i-1 : 0;
      const int ip1 = (i != ND-1) ? i+1 : ND-1;

      const int jm1 = (j !=    0) ? j-1 : 0;
      const int jp1 = (j != NT-1) ? j+1 : NT-1;

      const double dlogD = logD_values[ip1] - logD_values[im1];
      const double dlogT = logT_values[jp1] - logT_values[jm1];

      double Jp[2], Js[2];

      Jp[0] = (EOS_p[ip1 + j*ND] - EOS_p[im1 + j*ND]) / dlogD;
      Js[0] = (EOS_s[ip1 + j*ND] - EOS_s[im1 + j*ND]) / dlogD;

      Jp[1] = (EOS_p[i + jp1*ND] - EOS_p[i + jm1*ND]) / dlogT;
      Js[1] = (EOS_s[i + jp1*ND] - EOS_s[i + jm1*ND]) / dlogT;

      const double c2 = LIGHT_SPEED*LIGHT_SPEED;
      const double f  = MEV_TO_ERG / FM3_TO_CM3;

      const double p = pow(10.0, EOS_p[i + j*ND])*f; // erg/cm^3
      const double u = pow(10.0, EOS_u[i + j*ND])*f; // erg/cm^3
      const double D = pow(10.0, logD_values[i]);    // gm/cm^3

      const double Dh = D + u/c2 + p/c2;
      const double GammaEff = (Jp[0]*Js[1] - Jp[1]*Js[0])/Js[1];
      const double cs2 = GammaEff * p / Dh;

      EOS_cs2[i + j*ND] = cs2/c2; // (cm/s)^2
    }
  }
}
#ifdef NOSEARCH_DT
int ShenTabulatedNuclearEos::find_upper_index_D(double logD) const
{
  int n0=0, n1=logD_values.size()-1;

  if (logD < logD_values[n0]) {
    throw SampledOutOfRangeDensity();
  }
  if (logD > logD_values[n1]) {
    throw SampledOutOfRangeDensity();
  }

  return 1.0 + (logD - logD_values[0]) / (logD_values[1] - logD_values[0]);
}

int ShenTabulatedNuclearEos::find_upper_index_T(double logT) const
{
  int n0=0, n1=logT_values.size()-1;

  if (logT < logT_values[n0]) {
    throw SampledOutOfRangeTemperature();
  }
  if (logT > logT_values[n1]) {
    throw SampledOutOfRangeTemperature();
  }

  return 1.0 + (logT - logT_values[0]) / (logT_values[1] - logT_values[0]);
}
#else

int ShenTabulatedNuclearEos::find_upper_index_D(double logD) const
{
  int n0=0, n1=logD_values.size()-1;

  if (logD < logD_values[n0]) {
    throw SampledOutOfRangeDensity();
  }
  if (logD > logD_values[n1]) {
    throw SampledOutOfRangeDensity();
  }

  while (n1 - n0 > 1) {

    const double logD_mid = logD_values[(n0+n1)/2];

    if (logD > logD_mid) n0 = (n0+n1)/2;
    else                 n1 = (n0+n1)/2;
  }

  return n1;
}

int ShenTabulatedNuclearEos::find_upper_index_T(double logT) const
{
  int n0=0, n1=logT_values.size()-1;

  if (logT < logT_values[n0]) {
    throw SampledOutOfRangeTemperature();
  }
  if (logT > logT_values[n1]) {
    throw SampledOutOfRangeTemperature();
  }

  while (n1 - n0 > 1) {

    const double logT_mid = logT_values[(n0+n1)/2];

    if (logT > logT_mid) n0 = (n0+n1)/2;
    else                 n1 = (n0+n1)/2;
  }

  return n1;
}
#endif

double ShenTabulatedNuclearEos::sample_EOS(const std::vector<double> &EOS,
                                           double logD, double logT, double *J) const
{
  try {
    return this->tabled_EOS(EOS, logD, logT, J);
  }
  catch (const SampledOutOfRangeDensity &e) {
    printf("[shen] warning: density out, using approximate. logD=%e\n", logD);
    return this->approx_EOS(EOS, logD, logT, J);
  }
  catch (const SampledOutOfRangeTemperature &e) {
    printf("[shen] warning: temperature out, using approximate logT=%e\n", logT);
    return this->approx_EOS(EOS, logD, logT, J);
  }
}
double ShenTabulatedNuclearEos::tabled_EOS(const std::vector<double> &EOS,
                                           double logD, double logT, double *J) const
{
  // Receives one of the lookup tables for pressure, temperature, or internal
  // energy, and performs a bilinear interpolation on the nearest 4 samples in
  // order to construct the needed EOS variable. Sampling is done in log10-space
  // for both the density and temperature.
  // ---------------------------------------------------------------------------

  const int Di = find_upper_index_D(logD);
  const int Tj = find_upper_index_T(logT);

  const int ND = logD_values.size();

  const double dlogD = logD_values[Di] - logD_values[Di-1];
  const double dlogT = logT_values[Tj] - logT_values[Tj-1];

  // http://en.wikipedia.org/wiki/Bilinear_interpolation
  // ---------------------------------------------------------------------------
  const double f00 = EOS[(Di-1) + (Tj-1)*ND];
  const double f01 = EOS[(Di-1) + (Tj-0)*ND];
  const double f10 = EOS[(Di-0) + (Tj-1)*ND];
  const double f11 = EOS[(Di-0) + (Tj-0)*ND];

  const double x = (logD - logD_values[Di-1]) / dlogD;
  const double y = (logT - logT_values[Tj-1]) / dlogT;

  const double b1 = f00;
  const double b2 = f10 - f00;
  const double b3 = f01 - f00;
  const double b4 = f00 - f10 - f01 + f11;

  // Fill in the derivative of the EOS variable with respect to the log of rho
  // and T, if requested by the user.
  // ---------------------------------------------------------------------------
  if (J != NULL) {
    const double dfdlogD = (b2 + b4*y)/dlogD; // (df/dx)*(dx/dlogD)
    const double dfdlogT = (b3 + b4*x)/dlogT; // (df/dy)*(dy/dlogT)
    J[0] = dfdlogD;
    J[1] = dfdlogT;
  }

  return b1 + b2*x + b3*y + b4*x*y;
}

double ShenTabulatedNuclearEos::approx_EOS(const std::vector<double> &EOS,
                                           double logD, double logT, double *J) const
{
  const int ND = logD_values.size();
  const int NT = logT_values.size();

  const double dlogD = logD_values[ND-1] - logD_values[0];
  const double dlogT = logT_values[NT-1] - logT_values[0];

  const double dfdlogD = (EOS[ND-1 + (NT/2)*ND] - EOS[0 + (NT/2)*ND]) / dlogD;
  const double dfdlogT = (EOS[ND/2 + (NT-1)*ND] - EOS[ND/2 +   0*ND]) / dlogT;

  if (J != NULL) {
    J[0] = dfdlogD;
    J[1] = dfdlogT;
  }

  const double fc = EOS[ND/2 + (NT/2)*ND];
  const double xc = logD_values[ND/2];
  const double yc = logT_values[NT/2];

  return fc + (logD-xc)*dfdlogD + (logT-yc)*dfdlogT;
}

double ShenTabulatedNuclearEos::inverse_lookup_T(const std::vector<double> &EOS,
                                                 double logD, double logF) const
{
  // First do a bisection to identify the nearest tabulated temperature for the
  // density and variable (F := { p,s,u }) requested.
  // ---------------------------------------------------------------------------
  int n0=0, n1=logT_values.size()-1;

  while (n1 - n0 > 1) {

    const double logF_mid = sample_EOS(EOS, logD, logT_values[(n0+n1)/2]);

    if (logF > logF_mid) n0 = (n0+n1)/2;
    else                 n1 = (n0+n1)/2;
  }

  // The variable 'off_table' below indicates that F could not be bracketed; the
  // corresponding temperature is off the table. When this happens, we use an
  // safe/approximate evalutation instead of a table lookup.
  // ---------------------------------------------------------------------------
  const double logF0 = sample_EOS(EOS, logD, logT_values[n0]);
  const double logF1 = sample_EOS(EOS, logD, logT_values[n1]);

  const double logD0 = logD_values.front();
  const double logD1 = logD_values.back();

  const int off_table =
    (logF < logF0 || logF1 < logF) +
    (logD < logD0 || logD1 < logD);

  double logT = 0.5*(logT_values[n0] + logT_values[n1]);

  // Now refine the guess using a single Newton-Rapheson iteration. The use of
  // bilinear interpolation guarentees that a single iteration gets the root.
  // ---------------------------------------------------------------------------
  double J[2];
  double f = (off_table ?
              approx_EOS(EOS, logD, logT, J) :
              tabled_EOS(EOS, logD, logT, J)) - logF;
  double g = fabs(J[1]) > EFFECTIVELY_ZERO ? J[1] : EFFECTIVELY_ZERO;

  logT -= f/g;

  if (off_table) {
    printf("[shen] warning: inverse lookup used approximate. logT=%e\n", logT);
  }

  return logT;
}


double ShenTabulatedNuclearEos::Pressure(double D, double logT) const
{
  D /= units.GramsPerCubicCentimeter();
  const double p = pow(10.0, sample_EOS(EOS_p, log10(D), logT));
  return p * units.MeVPerCubicFemtometer();
}
double ShenTabulatedNuclearEos::Internal(double D, double logT) const
{
  D /= units.GramsPerCubicCentimeter();
  const double u = pow(10.0, sample_EOS(EOS_u, log10(D), logT));
  return u * units.MeVPerCubicFemtometer();
}
double ShenTabulatedNuclearEos::Entropy(double D, double logT) const
// -----------------------------------------------------------------------------
// Returns the entropy per baryon
// -----------------------------------------------------------------------------
{
  D /= units.GramsPerCubicCentimeter();
  const double s = pow(10.0, sample_EOS(EOS_s, log10(D), logT));
  return s * units.BoltzmannConstant();
}

double ShenTabulatedNuclearEos::SoundSpeed2Nr(double D, double logT) const
{
  return this->SoundSpeed2Sr(D, logT);
}
double ShenTabulatedNuclearEos::SoundSpeed2Sr(double D, double logT) const
{
  D /= units.GramsPerCubicCentimeter();
  const double cs2 = sample_EOS(EOS_cs2, log10(D), logT); // in units of light speed
  return cs2 * pow(units.LightSpeed(), 2.0); // in code units
}

double ShenTabulatedNuclearEos::Temperature_u(double D, double u) const
{
  const double logD = log10(D/units.GramsPerCubicCentimeter());
  const double logu = log10(u/units.MeVPerCubicFemtometer());

  return this->inverse_lookup_T(EOS_u, logD, logu);
}
double ShenTabulatedNuclearEos::Temperature_p(double D, double p) const
{
  const double logD = log10(D/units.GramsPerCubicCentimeter());
  const double logp = log10(p/units.MeVPerCubicFemtometer());

  return this->inverse_lookup_T(EOS_p, logD, logp);
}
double ShenTabulatedNuclearEos::TemperatureMeV(double D, double p) const
{
  const double logD = log10(D/units.GramsPerCubicCentimeter());
  const double logp = log10(p/units.MeVPerCubicFemtometer());
  const double logT = this->inverse_lookup_T(EOS_p, logD, logp);
  return pow(10.0, logT);
}
double ShenTabulatedNuclearEos::TemperatureArb(double D, double T_MeV) const
{
  return log10(T_MeV);
}

double ShenTabulatedNuclearEos::Derivatives_u(double D, double logT, double *J) const
{
  const double logD = log10(D/units.GramsPerCubicCentimeter());

  // J will first hold dlogu / dlogx = (x/u) (du/dx)
  // ---------------------------------------------------------------------------
  const double u = pow(10.0, sample_EOS(EOS_u, logD, logT, J));

  J[0] *= u/D * units.MeVPerCubicFemtometer();
  J[1] *= u   * units.MeVPerCubicFemtometer() * log(10);

  return u * units.MeVPerCubicFemtometer();
}
double ShenTabulatedNuclearEos::Derivatives_p(double D, double logT, double *J) const
{
  const double logD = log10(D/units.GramsPerCubicCentimeter());

  // J will first hold dlogp / dlogx = (x/u) (du/dx)
  // ---------------------------------------------------------------------------
  const double p = pow(10.0, sample_EOS(EOS_p, logD, logT, J));

  J[0] *= p/D * units.MeVPerCubicFemtometer();
  J[1] *= p   * units.MeVPerCubicFemtometer() * log(10);

  return p * units.MeVPerCubicFemtometer();
}










void ShenTabulatedNuclearEos::self_test_derivatives()
{

  int old_verbose = this->verbose;
  this->verbose = false;

  const double logD0 = log10(this->DensLower()) - OFF_TABLE_VAL * 0.5;
  const double logD1 = log10(this->DensUpper()) + OFF_TABLE_VAL * 0.5;

  const double logT0 = this->TempLower() - OFF_TABLE_VAL * 0.5;
  const double logT1 = this->TempUpper() + OFF_TABLE_VAL * 0.5;

  for (int n=0; n<10; ++n) {

    const double logD = logD0 + (logD1 - logD0) * rand() / RAND_MAX;
    const double logT = logT0 + (logT1 - logT0) * rand() / RAND_MAX;

    const double D = pow(10.0, logD);

    const double dx = 0.01 * D;
    const double dy = 0.01 * logT;

    double Ju[2], Jp[2]; // sampled analytically
    double Gu[2], Gp[2]; // sampled numerically

    this->Derivatives_u(D, logT, Ju);
    this->Derivatives_p(D, logT, Jp);

    Gu[0] = (this->Internal(D+dx, logT) - this->Internal(D-dx, logT)) / (2*dx);
    Gu[1] = (this->Internal(D, logT+dy) - this->Internal(D, logT-dy)) / (2*dy);

    Gp[0] = (this->Pressure(D+dx, logT) - this->Pressure(D-dx, logT)) / (2*dx);
    Gp[1] = (this->Pressure(D, logT+dy) - this->Pressure(D, logT-dy)) / (2*dy);

    printf("Ju[0] ?= Gu[0] : %e ?= %e    Ju[1] ?= Gu[1] : %e ?= %e\n", Ju[0], Gu[0], Ju[1], Gu[1]);
    printf("Jp[0] ?= Gp[0] : %e ?= %e    Jp[1] ?= Gp[1] : %e ?= %e\n", Jp[0], Gp[0], Jp[1], Gp[1]);
  }

  this->verbose = old_verbose;
}


void ShenTabulatedNuclearEos::self_test_inversion()
{
  int old_verbose = this->verbose;
  this->verbose = false;

  const double logD0 = log10(this->DensLower());
  const double logD1 = log10(this->DensUpper());

  const double logT0 = this->TempLower();
  const double logT1 = this->TempUpper();

  for (int n=0; n<10; ++n) {

    const double logD = logD0 + (logD1 - logD0) * rand() / RAND_MAX;
    const double logT = logT0 + (logT1 - logT0) * rand() / RAND_MAX;

    const double D = pow(10.0, logD);
    const double T = pow(10.0, logT);

    printf("attempting to invert from D=%le, T=%le... ", D, T);

    const double p = this->Pressure(D, logT);
    const double newlogT = this->Temperature_p(D, p);
    printf("got back T=%le\n", pow(10.0, newlogT));
  }

  for (int n=0; n<10; ++n) {

    const double logD = logD0 + (logD1 - logD0) * rand() / RAND_MAX;
    const double logT = logT0 + (logT1 - logT0) * rand() / RAND_MAX;

    const double D = pow(10.0, logD);
    const double T = pow(10.0, logT);

    printf("attempting to invert from D=%le, T=%le... ", D, T);

    const double u = this->Internal(D, logT);
    const double newlogT = this->Temperature_u(D, u);
    printf("got back T=%le\n", pow(10.0, newlogT));
  }

  if (false) {
    const double logD = 0.5*(logD0 + logD1);
    const double logT = 0.5*(logT0 + logT1);
    const double D = pow(10.0, logD);

    printf("off-table inverse lookup in p, u...\n");

    const double p = this->Pressure(D, logT) * 1e20;
    const double u = this->Internal(D, logT) * 1e20;

    const double logTp = this->Temperature_p(D, p);
    const double logTu = this->Temperature_u(D, u);

    const double newp = this->Pressure(D, logTp);
    const double newu = this->Internal(D, logTu);

    printf("got back p=%e (expect %e)\n", newp, p);
    printf("got back u=%e (expect %e)\n", newu, u);
  }

  if (false) {
    const double logD =      logD1 + 2.0;
    const double logT = 0.5*(logT0 + logT1);
    const double D = pow(10.0, logD);

    printf("off-table inverse lookup in D...\n");

    const double p = this->Pressure(D, logT);
    const double u = this->Internal(D, logT);

    const double logTp = this->Temperature_p(D, p);
    const double logTu = this->Temperature_u(D, u);

    const double newp = this->Pressure(D, logTp);
    const double newu = this->Internal(D, logTu);

    printf("got back p=%e (expect %e)\n", newp, p);
    printf("got back u=%e (expect %e)\n", newu, u);
  }

  if (false) {
    printf("This test fails because the extrapolation being used presently "
           "makes the pressure multi-valued.\n");

    printf("off-table inverse lookup in D... limits are:\n");
    printf("%e %e\n", logD0, logD1);
    printf("%e %e\n", logT0, logT1);

    const double D = 100.0;
    const double p =  10.0;

    const double logT = this->Temperature_p(D, p);
    const double newp = this->Pressure(D, logT);

    printf("got back p=%e (expect %e)\n", newp, p);
  }

  this->verbose = old_verbose;
}

void ShenTabulatedNuclearEos::self_test_interpolation()
{
  FILE *outp = fopen("shen_p.txt", "w");
  FILE *outs = fopen("shen_s.txt", "w");
  FILE *outu = fopen("shen_u.txt", "w");

  const int nsampT = 64;
  const int nsampD = 64;

  const double logD0 = log10(this->DensLower()) - OFF_TABLE_VAL * 0.5;
  const double logD1 = log10(this->DensUpper()) + OFF_TABLE_VAL * 0.5;

  const double logT0 = this->TempLower() - OFF_TABLE_VAL * 0.5;
  const double logT1 = this->TempUpper() + OFF_TABLE_VAL * 0.5;

  const double dlogD = (logD1 - logD0) / nsampD;
  const double dlogT = (logT1 - logT0) / nsampT;

  for (double logD=logD0 + 1e-12; logD < logD1; logD += dlogD) {
    for (double logT=logT0 + 1e-12; logT < logT1; logT += dlogT) {

      const double D = pow(10.0, logD);
      const double T = pow(10.0, logT);

      fprintf(outp, "%le %le %le\n", D, T, this->Pressure(D, logT));
      fprintf(outs, "%le %le %le\n", D, T, this->Entropy(D, logT));
      fprintf(outu, "%le %le %le\n", D, T, this->Internal(D, logT));
    }
  }

  fclose(outp);
  fclose(outs);
  fclose(outu);
}

void ShenTabulatedNuclearEos::self_test_soundspeed()
{
  int old_verbose = this->verbose;
  this->verbose = false;

  FILE *out_cs = fopen("shen_cs.txt", "w");

  const int nsampD = 64;
  const int nsampT = 64;

  const double logD0 = log10(this->DensLower());
  const double logD1 = log10(this->DensUpper());

  const double logT0 = this->TempLower();
  const double logT1 = this->TempUpper();

  const double dlogD = (logD1 - logD0) / nsampD;
  const double dlogT = (logT1 - logT0) / nsampT;

  for (double logD=logD0 + 1e-12; logD < logD1; logD += dlogD) {
    for (double logT=logT0 + 1e-12; logT < logT1; logT += dlogT) {

      const double D = pow(10.0, logD);
      const double T = pow(10.0, logT);

      const double cs2 = this->SoundSpeed2Sr(D, logT);
      fprintf(out_cs, "%le %le %le\n", D, T, sqrt(cs2));

    }
  }

  fclose(out_cs);
  this->verbose = old_verbose;
}



TabulatedEos ShenTabulatedNuclearEos::LoadTable(const char *fname, double YpExtract,
                                                const double *DensRange_,
                                                const double *TempRange_)
{
  if (verbose) {
    printf("loading EOS table %s from disk...\n", fname);
  }

  FILE *shentab = fopen(fname, "r");
  TabulatedEos tab;

  if (shentab == NULL) {
    throw UnableToLoadTable();
  }

  double DensRange[2];
  double TempRange[2];

  if (DensRange_ == NULL) {
    DensRange[0] = 1e-01; // g/cm^3;
    DensRange[1] = 1e+16;
  }
  else {
    std::memcpy(DensRange, DensRange_, 2*sizeof(double));
  }

  if (TempRange_ == NULL) {
    TempRange[0] =   0.1; // MeV
    TempRange[1] = 200.0;
  }
  else {
    std::memcpy(TempRange, TempRange_, 2*sizeof(double));
  }


  char line[1024];
  int n = 0;
  double logT, T;

  while (!feof(shentab)) {

    fgets(line, sizeof(line), shentab);

    if (strncmp(line, " cccccccccccc", 13) == 0) {

      fgets(line, sizeof(line), shentab); //  Log10(Temp)   Temp
      fgets(line, sizeof(line), shentab); // -1.000000E+00  1.000000E-01

      sscanf(line, "%le %le\n", &logT, &T);

      if ((TempRange[0] < T && T < TempRange[1])) {
        tab.logT_values.push_back(logT);
      }
    }
    else if (strlen(line) == 3) {
      // pass
    }
    else {
      double d[18];

      sscanf(line,
             "%le %le %le %le %le %le %le %le %le %le %le %le %le %le "
             "%le %le %le %le\n",
             d+ 0, d+ 1, d+ 2, d+ 3, d+ 4, d+ 5, d+ 6, d+ 7, d +8, d +9,
             d+10, d+11, d+12, d+13, d+14, d+15, d+16, d+17);

      const double nB    = d[inB];          // baryon number density
      const double Yp    = d[iYp];          // proton fraction
      const double logD  = d[ilogD];        // log_10 of density
      const double D     = pow(10.0, logD); // density
      const double S     = d[iS];           // entropy per baryon
      const double E     = d[iEint];        // total energy per baryon (u = E * nB)
      const double p     = d[ip];           // gas pressure

      if ( (fabs(Yp - YpExtract) < 1e-14) &&
           (TempRange[0] < T && T < TempRange[1]) &&
           (DensRange[0] < D && D < DensRange[1]) ) {

        // If this was the first pass through the densities, then populate the
        // list of logD values.
        // ---------------------------------------------------------------------
        if (tab.logT_values.size() == 1) {
          tab.logD_values.push_back(logD);
        }

        tab.EOS_p.push_back(log10(p));
        tab.EOS_s.push_back(log10(S)); // entropy per baryon
        tab.EOS_u.push_back(log10(E * nB));

        ++n;
      }
    }

    // The 2d-arrays tabulated by the EOS are stored in stored in 'Fortran'
    // order with respect to the ordering EOS(D, T). For example,
    //
    //  Pressure(Di, Tj) := Pressure[i + j*ND];
    // -------------------------------------------------------------------------
  }
  if (verbose) {
    printf("added %d = (%ld x %ld) entries to the table\n",
           n, tab.logD_values.size(), tab.logT_values.size());
  }
  fclose(shentab);

  return tab;
}
