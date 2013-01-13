

#include "eos.hpp"
static const double Kelvin_to_MeV = 8.621738e-11;


AdiabaticEos::AdiabaticEos(double Gamma) : Gamma(Gamma) { }

double AdiabaticEos::Pressure(double rho, double T) const
{
  return T;
}
double AdiabaticEos::Internal(double rho, double T) const
{
  return T / (Gamma - 1);
}
double AdiabaticEos::Entropy(double D, double T) const
{
  return 0.0; // Neither needed nor implemented.
}
double AdiabaticEos::SoundSpeed2Nr(double rho, double T) const
{
  return Gamma*T / rho;
}
double AdiabaticEos::SoundSpeed2Sr(double rho, double T) const
{
  const double rhoh = rho + Internal(rho, T) + Pressure(rho, T);
  return Gamma*T / rhoh;
}

double AdiabaticEos::Temperature_u(double rho, double u) const
{
  return u*(Gamma - 1);
}
double AdiabaticEos::Temperature_p(double rho, double P) const
{
  return P;
}
double AdiabaticEos::TemperatureMeV(double rho, double P) const
{
  P   /= Mara->units->ErgPerCubicCentimeter();
  rho /= Mara->units->GramsPerCubicCentimeter();

  const double mp = Mara->units->ProtonMass();
  const double kb = Mara->units->BoltzmannConstant();
  const double T_Kelvin = (P/rho) * (mp/kb);

  return T_Kelvin * Kelvin_to_MeV;
}
double AdiabaticEos::TemperatureArb(double rho, double T_MeV) const
{
  rho /= Mara->units->GramsPerCubicCentimeter();

  const double mp = Mara->units->ProtonMass();
  const double kb = Mara->units->BoltzmannConstant();

  const double T_Kelvin = T_MeV / Kelvin_to_MeV;
  const double P = (kb*T_Kelvin) * (rho/mp);
  return P * Mara->units->ErgPerCubicCentimeter();
}

double AdiabaticEos::Derivatives_u(double D, double T, double *J) const
{
  J[0] = 0.0;
  J[1] = 1.0 / (Gamma - 1);

  return this->Internal(D, T);
}
double AdiabaticEos::Derivatives_p(double D, double T, double *J) const
{
  J[0] = 0.0;
  J[1] = 1.0;

  return this->Pressure(D, T);
}



ThermalBarotropicEos::ThermalBarotropicEos(double GammaT, double GammaB, double Kappa)
  : GammaT(GammaT), GammaB(GammaB), Kappa(Kappa) { }

double ThermalBarotropicEos::Pressure(double rho, double T) const
{
  return T;
}
double ThermalBarotropicEos::Internal(double rho, double T) const
{
  // returns rho * (eb + et)
  // ---------------------------------------------------------------------------
  const double Pb = Kappa*pow(rho, GammaB);
  const double Pt = T - Pb;
  return Pt / (GammaT - 1) + Pb / (GammaB - 1);
}
double ThermalBarotropicEos::Entropy(double D, double T) const
{
  return 0.0; // Neither needed nor implemented.
}
double ThermalBarotropicEos::SoundSpeed2Nr(double rho, double T) const
{
  const double Pb = Kappa*pow(rho, GammaB);
  const double Pt = T - Pb;
  return (Pt*GammaT + Pb*GammaB) / rho;
}
double ThermalBarotropicEos::SoundSpeed2Sr(double rho, double T) const
{
  const double Pb = Kappa*pow(rho, GammaB);
  const double Pt = T - Pb;
  const double rhoh = rho + Internal(rho, T) + Pressure(rho, T);
  return (Pt*GammaT + Pb*GammaB) / rhoh;
}

double ThermalBarotropicEos::Temperature_u(double rho, double u) const
{
  const double Pb = Kappa*pow(rho, GammaB);
  const double ub = Pb / (GammaB - 1);
  return (u - ub)*(GammaT - 1) + Pb;
}
double ThermalBarotropicEos::Temperature_p(double rho, double P) const
{
  return P;
}
double ThermalBarotropicEos::TemperatureMeV(double rho, double P) const
{
  P   /= Mara->units->ErgPerCubicCentimeter();
  rho /= Mara->units->GramsPerCubicCentimeter();

  const double mp = Mara->units->ProtonMass();
  const double kb = Mara->units->BoltzmannConstant();
  const double T_Kelvin = (P/rho) * (mp/kb);
  return T_Kelvin * Kelvin_to_MeV;
}
double ThermalBarotropicEos::TemperatureArb(double rho, double T_MeV) const
{
  rho /= Mara->units->GramsPerCubicCentimeter();

  const double mp = Mara->units->ProtonMass();
  const double kb = Mara->units->BoltzmannConstant();

  const double T_Kelvin = T_MeV / Kelvin_to_MeV;
  const double P = (kb*T_Kelvin) * (rho/mp);
  return P * Mara->units->ErgPerCubicCentimeter();
}

double ThermalBarotropicEos::Derivatives_u(double D, double T, double *J) const
{
  const double gB = GammaB;
  const double gT = GammaT;
  const double kp = Kappa;

  J[0] = -(gB*(gB - gT) / ((gB - 1) * (gT - 1))) * kp*pow(D, gB - 1);
  J[1] = 1.0 / (gT - 1);

  return this->Internal(D, T);
}
double ThermalBarotropicEos::Derivatives_p(double D, double T, double *J) const
{
  J[0] = 0.0;
  J[1] = 1.0;

  return this->Pressure(D, T);
}
