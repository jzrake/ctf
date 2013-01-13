
/*------------------------------------------------------------------------------
 * FILE: tabulated_eos.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * DESCRIPTION:
 *
 * Builds a class which reads from a tabulated equation of state in the
 * independent variables density (D) and temperature (T). The user-upplied data
 * are arrays representing the pressure (p), internal energy density (u), and
 * sound speed (c). Each of these arrays must be data buffers of dimension
 * (NDxNT) in standard C (row-major) ordering. The arrays D_values and T_values
 * are equally-spaced arrays of density and temperature, and are of length ND
 * and NT respectively. All units, with the exception of temperature, are
 * assumed to already be converted into code units. Temperature is in units of
 * MeV, although this may easily be changed.
 *
 *------------------------------------------------------------------------------
 */


#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include "eos.hpp"



bool GenericTabulatedEos::verbose = false;



GenericTabulatedEos::GenericTabulatedEos(std::vector<double> &D_values,
					 std::vector<double> &T_values,
					 std::vector<double> &p,
					 std::vector<double> &u,
					 std::vector<double> &c)
  : D_values(D_values),
    T_values(T_values),
    EOS_p(p),
    EOS_u(u),
    EOS_c(c)
{
  printf("[eos] building new GenericTabulatedEos\n");
  printf("[eos] density points: %ld\n", D_values.size());
  printf("[eos] temperature points: %ld\n", T_values.size());

  size_t npts = D_values.size() * T_values.size();

  if (npts != EOS_p.size()) {
    printf("[eos] warning: pressure array has %ld points, expected %ld\n",
	   EOS_p.size(), npts);
  }
  if (npts != EOS_u.size()) {
    printf("[eos] warning: internal energy array has %ld points, expected %ld\n",
	   EOS_u.size(), npts);
  }
  if (npts != EOS_c.size()) {
    printf("[eos] warning: sound speed array has %ld points, expected %ld\n",
	   EOS_c.size(), npts);
  }
}

int GenericTabulatedEos::find_upper_index_D(double D) const
{
  int n0=0, n1=D_values.size()-1;

  if (D < D_values[n0]) {
    throw SampledOutOfRangeDensity();
  }
  if (D > D_values[n1]) {
    throw SampledOutOfRangeDensity();
  }

  return 1.0 + (D - D_values[0]) / (D_values[1] - D_values[0]);
}

int GenericTabulatedEos::find_upper_index_T(double T) const
{
  int n0=0, n1=T_values.size()-1;

  if (T < T_values[n0]) {
    throw SampledOutOfRangeTemperature();
  }
  if (T > T_values[n1]) {
    throw SampledOutOfRangeTemperature();
  }

  return 1.0 + (T - T_values[0]) / (T_values[1] - T_values[0]);
}

double GenericTabulatedEos::sample_EOS(const std::vector<double> &EOS,
                                       double D, double T, double *J) const
{
  try {
    return this->tabled_EOS(EOS, D, T, J);
  }
  catch (const SampledOutOfRangeDensity &e) {
    if (verbose) printf("[eos] warning: density out, using approximate. D=%e\n", D);
    return this->approx_EOS(EOS, D, T, J);
  }
  catch (const SampledOutOfRangeTemperature &e) {
    if (verbose) printf("[eos] warning: temperature out, using approximate T=%e\n", T);
    return this->approx_EOS(EOS, D, T, J);
  }
}
double GenericTabulatedEos::tabled_EOS(const std::vector<double> &EOS,
                                       double D, double T, double *J) const
// -----------------------------------------------------------------------------
// Receives one of the lookup tables for pressure, temperature, or internal
// energy, and performs a bilinear interpolation on the nearest 4 samples in
// order to construct the needed EOS variable.
// -----------------------------------------------------------------------------
{
  const int Di = find_upper_index_D(D);
  const int Tj = find_upper_index_T(T);
  const int NT = T_values.size();
  const double dD = D_values[Di] - D_values[Di-1];
  const double dT = T_values[Tj] - T_values[Tj-1];

  // http://en.wikipedia.org/wiki/Bilinear_interpolation
  // ---------------------------------------------------------------------------
  const double f00 = EOS[(Di-1)*NT + (Tj-1)];
  const double f01 = EOS[(Di-1)*NT + (Tj-0)];
  const double f10 = EOS[(Di-0)*NT + (Tj-1)];
  const double f11 = EOS[(Di-0)*NT + (Tj-0)];

  const double x = (D - D_values[Di-1]) / dD;
  const double y = (T - T_values[Tj-1]) / dT;

  const double b1 = f00;
  const double b2 = f10 - f00;
  const double b3 = f01 - f00;
  const double b4 = f00 - f10 - f01 + f11;

  // Fill in the derivative of the EOS variable with respect to D and T, if
  // requested by the user.
  // ---------------------------------------------------------------------------
  if (J != NULL) {
    const double dfdD = (b2 + b4*y) / dD; // (df/dx)*(dx/dD)
    const double dfdT = (b3 + b4*x) / dT; // (df/dy)*(dy/dT)
    J[0] = dfdD;
    J[1] = dfdT;
  }

  return b1 + b2*x + b3*y + b4*x*y;
}

double GenericTabulatedEos::approx_EOS(const std::vector<double> &EOS,
                                       double D, double T, double *J) const
{
  const int ND = D_values.size();
  const int NT = T_values.size();

  const double dD = D_values[ND-1] - D_values[0];
  const double dT = T_values[NT-1] - T_values[0];

  const double dfdD = (EOS[(ND-1)*NT + NT/2] - EOS[(0   )*NT + NT/2]) / dD;
  const double dfdT = (EOS[(ND/2)*NT + NT-1] - EOS[(ND/2)*NT +    0]) / dT;

  if (J != NULL) {
    J[0] = dfdD;
    J[1] = dfdT;
  }

  const double fc = EOS[(ND/2)*NT + NT/2];
  const double xc = D_values[ND/2];
  const double yc = T_values[NT/2];

  return fc + (D-xc) * dfdD + (T-yc) * dfdT;
}

double GenericTabulatedEos::inverse_lookup_T(const std::vector<double> &EOS,
                                             double D, double F) const
{
  // First do a bisection to identify the nearest tabulated temperature for the
  // density and variable (F := { p,u }) requested.
  // ---------------------------------------------------------------------------
  int n0=0, n1=T_values.size()-1;

  while (n1 - n0 > 1) {

    const double F_mid = sample_EOS(EOS, D, T_values[(n0+n1)/2]);

    if (F > F_mid) n0 = (n0+n1)/2;
    else           n1 = (n0+n1)/2;
  }

  // The variable 'off_table' below indicates that F could not be bracketed; the
  // corresponding temperature is off the table. When this happens, we use an
  // safe/approximate evalutation instead of a table lookup.
  // ---------------------------------------------------------------------------
  const double F0 = sample_EOS(EOS, D, T_values[n0]);
  const double F1 = sample_EOS(EOS, D, T_values[n1]);

  const double D0 = D_values.front();
  const double D1 = D_values.back();

  const int off_table =
    (F < F0 || F1 < F) +
    (D < D0 || D1 < D);

  double T = 0.5*(T_values[n0] + T_values[n1]);

  // Now refine the guess using a single Newton-Rapheson iteration. The use of
  // bilinear interpolation guarantees that a single iteration gets the root.
  // ---------------------------------------------------------------------------
  double J[2];
  double f = (off_table ?
              approx_EOS(EOS, D, T, J) :
              tabled_EOS(EOS, D, T, J)) - F;
  double g = fabs(J[1]) > EFFECTIVELY_ZERO ? J[1] : EFFECTIVELY_ZERO;

  T -= f/g;

  if (off_table && verbose) {
    printf("[eos] warning: inverse lookup on (D,F) = (%f,%f) used approximate. "
    	   "T=%e\n", D, F, T);
  }

  return T;
}


double GenericTabulatedEos::Pressure(double D, double T) const
{
  return sample_EOS(EOS_p, D, T);
}
double GenericTabulatedEos::Internal(double D, double T) const
{
  return sample_EOS(EOS_u, D, T);
}
double GenericTabulatedEos::Entropy(double D, double T) const // per baryon
{
  return 0.0; // not implemented
}
double GenericTabulatedEos::Derivatives_u(double D, double T, double *J) const
{
  return sample_EOS(EOS_u, D, T, J);
}
double GenericTabulatedEos::Derivatives_p(double D, double T, double *J) const
{
  return sample_EOS(EOS_p, D, T, J);
}



double GenericTabulatedEos::SoundSpeed2Nr(double D, double T) const
{
  return this->SoundSpeed2Sr(D, T);
}
double GenericTabulatedEos::SoundSpeed2Sr(double D, double T) const
{
  return pow(sample_EOS(EOS_c, D, T), 2.0); // in units of light speed
}
double GenericTabulatedEos::Temperature_u(double D, double u) const
{
  return this->inverse_lookup_T(EOS_u, D, u);
}
double GenericTabulatedEos::Temperature_p(double D, double p) const
{
  return this->inverse_lookup_T(EOS_p, D, p);
}
double GenericTabulatedEos::TemperatureMeV(double D, double p) const
{
  return this->inverse_lookup_T(EOS_p, D, p);
}
double GenericTabulatedEos::TemperatureArb(double D, double T_MeV) const
{
  return T_MeV;
}
