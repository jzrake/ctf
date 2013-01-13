
/*------------------------------------------------------------------------------
 * FILE: ou-field.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION: Contains a class for generating vector fields from an OU-process
 *
 * REFERENCES:
 *
 *   http://en.wikipedia.org/wiki/Ornstein–Uhlenbeck_process
 *   Wolfram Schmidt, Wolfgang Hillebrandt, Jens C. Niemeyerb (2005)
 *
 *------------------------------------------------------------------------------
 */

#ifndef __StochasticVectorField_HEADER__
#define __StochasticVectorField_HEADER__

#include <valarray>
#include <vector>
#include <complex>
#include <fstream>
#include "hydro.hpp"
#include "ornuhl.hpp"

class StochasticVectorField2d : public StochasticVectorField
{
private:
  int k1, numk;
  double P0, zeta;
  std::vector<OrnsteinUhlenbeckProcess> OuProcesses;
  std::valarray<std::complex<double> > Fx, Fy;
  std::valarray<std::complex<double> > Kx, Ky;

public:
  StochasticVectorField2d(double P0, double zeta, int k1,
			  int seed, double (*PS)(double)=NULL);
  StochasticVectorField2d(std::istream &stream);

  void AdvanceField(double dt);
  std::vector<double> SampleField(double x, double y, double z) const;
  void Serialize(std::ostream &stream) const;

  const std::valarray<std::complex<double> > &GetFx() { return Fx; };
  const std::valarray<std::complex<double> > &GetFy() { return Fy; };

  const std::valarray<std::complex<double> > &GetKx() { return Kx; };
  const std::valarray<std::complex<double> > &GetKy() { return Ky; };
} ;


class StochasticVectorField3d : public StochasticVectorField
{
private:
  int k1, numk;
  double P0, zeta;
  std::vector<OrnsteinUhlenbeckProcess> OuProcesses;
  std::valarray<std::complex<double> > Fx, Fy, Fz;
  std::valarray<std::complex<double> > Kx, Ky, Kz;

public:
  StochasticVectorField3d(double P0, double zeta, int k1, int seed,
			  double (*PS)(double)=NULL);
  StochasticVectorField3d(std::istream &stream);

  void AdvanceField(double dt);
  std::vector<double> SampleField(double x, double y, double z) const;
  void Serialize(std::ostream &stream) const;

  const std::valarray<std::complex<double> > &GetFx() { return Fx; };
  const std::valarray<std::complex<double> > &GetFy() { return Fy; };
  const std::valarray<std::complex<double> > &GetFz() { return Fz; };

  const std::valarray<std::complex<double> > &GetKx() { return Kx; };
  const std::valarray<std::complex<double> > &GetKy() { return Ky; };
  const std::valarray<std::complex<double> > &GetKz() { return Kz; };
} ;

#endif // __StochasticVectorField_HEADER__
