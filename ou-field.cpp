
/*------------------------------------------------------------------------------
 * FILE: ou-field.cpp
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


#include "ou-field.hpp"


// Constructors
// -----------------------------------------------------------------------------
StochasticVectorField2d::StochasticVectorField2d(double P0, double zeta, int k1,
						 int seed, double (*PS)(double))
  : k1(k1), numk((2*k1+1)*(2*k1+1)),
    P0(P0), zeta(zeta),
    Fx(numk), Fy(numk),
    Kx(numk), Ky(numk)
{
  const int si = (2*k1+1);
  const int sj = 1;
  double totpower = 0.0;

  for (int i=-k1; i<=k1; ++i) {
    for (int j=-k1; j<=k1; ++j) {

      const int m = (i+k1)*si + (j+k1)*sj;
      Kx[m] = 2*M_PI*i;
      Ky[m] = 2*M_PI*j;

      const double K[2] = { Kx[m].real(), Ky[m].real() };
      const double k = sqrt(K[0]*K[0] + K[1]*K[1]);

      totpower += PS ? PS(k) : 1.0;
    }
  }

  for (int m=0; m<numk; ++m) {
    const double K[3] = { Kx[m].real(), Ky[m].real() };
    const double k = sqrt(K[0]*K[0] + K[1]*K[1]);
    const double Pk = PS ? PS(k) : 1.0;

    OrnsteinUhlenbeckProcess process(1.0, sqrt(2*P0*Pk/totpower), seed+m);
    OuProcesses.push_back(process);
  }
}

StochasticVectorField3d::StochasticVectorField3d(double P0, double zeta, int k1, int seed,
                                                 double (*PS)(double))
  : k1(k1), numk((2*k1+1)*(2*k1+1)*(2*k1+1)),
    P0(P0), zeta(zeta),
    Fx(numk), Fy(numk), Fz(numk),
    Kx(numk), Ky(numk), Kz(numk)
{
  const int si = (2*k1+1)*(2*k1+1);
  const int sj = (2*k1+1);
  const int sk = 1;
  double totpower = 0.0;

  for (int i=-k1; i<=k1; ++i) {
    for (int j=-k1; j<=k1; ++j) {
      for (int k=-k1; k<=k1; ++k) {

        const int m = (i+k1)*si + (j+k1)*sj + (k+k1)*sk;
        Kx[m] = 2*M_PI*i;
        Ky[m] = 2*M_PI*j;
        Kz[m] = 2*M_PI*k;

        const double K[3] = { Kx[m].real(), Ky[m].real(), Kz[m].real() };
        const double k = sqrt(K[0]*K[0] + K[1]*K[1] + K[2]*K[2]);

        totpower += PS ? PS(k) : 1.0;
      }
    }
  }

  for (int m=0; m<numk; ++m) {
    const double K[3] = { Kx[m].real(), Ky[m].real(), Kz[m].real() };
    const double k = sqrt(K[0]*K[0] + K[1]*K[1] + K[2]*K[2]);
    const double Pk = PS ? PS(k) : 1.0;

    OrnsteinUhlenbeckProcess process(1.0, sqrt(2*P0*Pk/totpower), seed+m);
    OuProcesses.push_back(process);
  }
}

StochasticVectorField2d::StochasticVectorField2d(std::istream &stream)
{
  stream.read((char*)&k1, sizeof(int));
  numk = (2*k1+1)*(2*k1+1);

  Kx.resize(numk);
  Ky.resize(numk);

  Fx.resize(numk);
  Fy.resize(numk);

  OrnsteinUhlenbeckProcess dummy(0,0,0);
  OuProcesses.assign(numk, dummy);

  size_t bytes = numk*sizeof(std::complex<double>);

  stream.read((char*)&P0, sizeof(double));
  stream.read((char*)&zeta, sizeof(double));

  stream.read((char*)&Kx[0], bytes);
  stream.read((char*)&Ky[0], bytes);

  stream.read((char*)&Fx[0], bytes);
  stream.read((char*)&Fy[0], bytes);

  stream.read((char*)&OuProcesses[0],
              sizeof(OrnsteinUhlenbeckProcess)*OuProcesses.size());
}

StochasticVectorField3d::StochasticVectorField3d(std::istream &stream)
{
  stream.read((char*)&k1, sizeof(int));
  numk = (2*k1+1)*(2*k1+1)*(2*k1+1);

  Kx.resize(numk);
  Ky.resize(numk);
  Kz.resize(numk);

  Fx.resize(numk);
  Fy.resize(numk);
  Fz.resize(numk);

  OrnsteinUhlenbeckProcess dummy(0,0,0);
  OuProcesses.assign(numk, dummy);

  size_t bytes = numk*sizeof(std::complex<double>);

  stream.read((char*)&P0, sizeof(double));
  stream.read((char*)&zeta, sizeof(double));

  stream.read((char*)&Kx[0], bytes);
  stream.read((char*)&Ky[0], bytes);
  stream.read((char*)&Kz[0], bytes);

  stream.read((char*)&Fx[0], bytes);
  stream.read((char*)&Fy[0], bytes);
  stream.read((char*)&Fz[0], bytes);

  stream.read((char*)&OuProcesses[0],
              sizeof(OrnsteinUhlenbeckProcess)*OuProcesses.size());
}

// AdvanceField
// -----------------------------------------------------------------------------
void StochasticVectorField2d::AdvanceField(double dt)
{
  for (int m=0; m<numk; ++m) {
    std::complex<double> dW[2], dV[2];

    const double K[2] = { Kx[m].real(), Ky[m].real() };
    const double kdotk = K[0]*K[0] + K[1]*K[1];

    if (kdotk < 1e-8) continue; // don't advance the zero mode

    dW[0] = OuProcesses[m].Deviate(Fx[m], dt);
    dW[1] = OuProcesses[m].Deviate(Fy[m], dt);

    for (int p=0; p<2; ++p) {
      dV[p] = 0.0;
      for (int q=0; q<2; ++q) {
        const double P = (1 - 2*zeta)*K[p]*K[q]/kdotk + zeta*(p==q);
        dV[p] += P * dW[q];
      }
    }

    Fx[m] += dV[0];
    Fy[m] += dV[1];
  }
}
void StochasticVectorField3d::AdvanceField(double dt)
{
  for (int m=0; m<numk; ++m) {
    std::complex<double> dW[3], dV[3];

    const double K[3] = { Kx[m].real(), Ky[m].real(), Kz[m].real() };
    const double kdotk = K[0]*K[0] + K[1]*K[1] + K[2]*K[2];

    if (kdotk < 1e-8) continue; // don't advance the zero mode

    dW[0] = OuProcesses[m].Deviate(Fx[m], dt);
    dW[1] = OuProcesses[m].Deviate(Fy[m], dt);
    dW[2] = OuProcesses[m].Deviate(Fz[m], dt);

    for (int p=0; p<3; ++p) {
      dV[p] = 0.0;
      for (int q=0; q<3; ++q) {
        const double P = (1 - 2*zeta)*K[p]*K[q]/kdotk + zeta*(p==q);
        dV[p] += P * dW[q] / sqrt(1 - 1.25*zeta + 3*zeta*zeta);
      }
    }

    Fx[m] += dV[0];
    Fy[m] += dV[1];
    Fz[m] += dV[2];
  }
}

// SampleField
// -----------------------------------------------------------------------------
std::vector<double> StochasticVectorField2d::SampleField(double x, double y, double z) const
{
  const std::complex<double> X(0,x);
  const std::complex<double> Y(0,y);

  std::vector<double> F(2);
  std::valarray<std::complex<double> > expfac = exp(Kx*X + Ky*Y);

  F[0] = (Fx * expfac).sum().real();
  F[1] = (Fy * expfac).sum().real();

  return F;
}
std::vector<double> StochasticVectorField3d::SampleField(double x, double y, double z) const
{
  const std::complex<double> X(0,x);
  const std::complex<double> Y(0,y);
  const std::complex<double> Z(0,z);

  std::vector<double> F(3);
  std::valarray<std::complex<double> > expfac = exp(Kx*X + Ky*Y + Kz*Z);

  F[0] = (Fx * expfac).sum().real();
  F[1] = (Fy * expfac).sum().real();
  F[2] = (Fz * expfac).sum().real();

  return F;
}



// Serialize field and random number states to bit array.
// -----------------------------------------------------------------------------
void StochasticVectorField2d::Serialize(std::ostream &stream) const
{
  size_t bytes = numk*sizeof(std::complex<double>);

  stream.write((char*)&k1, sizeof(int));
  stream.write((char*)&P0, sizeof(double));
  stream.write((char*)&zeta, sizeof(double));

  stream.write((char*)&Kx[0], bytes);
  stream.write((char*)&Ky[0], bytes);

  stream.write((char*)&Fx[0], bytes);
  stream.write((char*)&Fy[0], bytes);

  stream.write((char*)&OuProcesses[0],
               sizeof(OrnsteinUhlenbeckProcess)*OuProcesses.size());
}
void StochasticVectorField3d::Serialize(std::ostream &stream) const
{
  size_t bytes = numk*sizeof(std::complex<double>);

  stream.write((char*)&k1, sizeof(int));
  stream.write((char*)&P0, sizeof(double));
  stream.write((char*)&zeta, sizeof(double));

  stream.write((char*)&Kx[0], bytes);
  stream.write((char*)&Ky[0], bytes);
  stream.write((char*)&Kz[0], bytes);

  stream.write((char*)&Fx[0], bytes);
  stream.write((char*)&Fy[0], bytes);
  stream.write((char*)&Fz[0], bytes);

  stream.write((char*)&OuProcesses[0],
               sizeof(OrnsteinUhlenbeckProcess)*OuProcesses.size());
}
