

#include <cstdio>
#include "eos.hpp"


PhysicalUnits::PhysicalUnits(double c, double g, double s)
  : Length(c), Mass(g), Time(s)
{

}
void PhysicalUnits::PrintUnits()
{
  printf("\n");
  printf("%24s = %-14.3e %s\n", "[Length]", Length, "cm");
  printf("%24s = %-14.3e %s\n", "[Mass]", Mass, "gm");
  printf("%24s = %-14.3e %s\n", "[Time]", Time, "s");
  printf("%24s = %-14.3e %s\n", "[Velocity]", Length/Time, "cm/s");
  printf("%24s = %-14.3e %s\n", "[Energy]", 1.0/Erg(), "erg");
  printf("%24s = %-14.3e %s\n", "[B-field]", 1.0/Gauss(), "Gauss");
  printf("%24s = %-14.3e %s\n", "MeV", MeV(), "[Energy]");
  printf("%24s = %-14.3e %s\n", "LightSpeed", LightSpeed(), "[Velocity]");
  printf("%24s = %-14.3e %s\n", "GramsPerCubicCentimeter", GramsPerCubicCentimeter(), "[Mass]/[Length]^3");
}





PhysicalUnits PhysicalUnitsFactory_cgs::
BuildPhysicalUnits()
{
  return PhysicalUnits(1.0, 1.0, 1.0);
}

PhysicalUnits PhysicalUnitsFactory_MeV_gm_s::
BuildPhysicalUnits()
{
  const double E      = MEV_TO_ERG; // erg
  const double Mass   = 1.0;        // gm
  const double Length = 1.0;        // cm
  const double V      = sqrt(E/Mass);
  const double Time   = Length / V;

  return PhysicalUnits(Length, Mass, Time);
}

PhysicalUnits PhysicalUnitsFactory_MeV_c_ms::
BuildPhysicalUnits()
{
  const double E = MEV_TO_ERG;  // erg
  const double V = LIGHT_SPEED; // cm/s
  const double T = 1e-3;        // s

  const double Length = V*T;
  const double Mass   = E / (V*V);
  const double Time   = T;

  return PhysicalUnits(Length, Mass, Time);
}

PhysicalUnits PhysicalUnitsFactory_erg_c_ms::
BuildPhysicalUnits()
{
  const double E = 1.0;         // erg
  const double V = LIGHT_SPEED; // cm/s
  const double T = 1e-3;        // s

  const double Length = V*T;
  const double Mass   = E / (V*V);
  const double Time   = T;

  return PhysicalUnits(Length, Mass, Time);
}

PhysicalUnits PhysicalUnitsFactory_neutron_star::
BuildPhysicalUnits()
{
  const double Density = 1e13;        // gm/cm^3
  const double V       = LIGHT_SPEED; // cm/s
  const double Length  = 1e2;         // cm
  const double Mass    = Density * pow(Length, 3.0);
  const double Time    = Length / V;

  return PhysicalUnits(Length, Mass, Time);
}



double PhysicalUnits::Gram() const
{
  return 1.0 / Mass;
}
double PhysicalUnits::Centimeter() const
{
  return 1.0 / Length;
}
double PhysicalUnits::Second() const
{
  return 1.0 / Time;
}
double PhysicalUnits::Velocity() const
{
  return this->Centimeter() / this->Second();
}
double PhysicalUnits::Erg() const
{
  return this->Gram() * pow(this->Centimeter() / this->Second(), 2.0);
}
double PhysicalUnits::MeV() const
{
  return MEV_TO_ERG * this->Erg();
}
double PhysicalUnits::Gauss() const
{
  return sqrt(this->Erg()/pow(this->Centimeter(), 3.0)/(4*M_PI));
}



double PhysicalUnits::GramsPerCubicCentimeter() const
{
  return this->Gram() / pow(this->Centimeter(), 3.0);
}
double PhysicalUnits::ErgPerCubicCentimeter() const
{
  return this->MeVPerCubicFemtometer() * 1e-33;
}
double PhysicalUnits::MeVPerCubicCentimeter() const
{
  return this->MeV() / pow(this->Centimeter(), 3);
}
double PhysicalUnits::MeVPerCubicFemtometer() const
{
  return this->MeV() / pow(1e-13*this->Centimeter(), 3);
}
double PhysicalUnits::BoltzmannConstant() const
{
  return BOLTZMANN_CONSTANT * this->Erg(); // per Kelvin := 1
}
double PhysicalUnits::LightSpeed() const
{
  return LIGHT_SPEED * this->Centimeter() / this->Second();
}
double PhysicalUnits::ProtonMass() const
{
  return PROTON_MASS * this->Gram();
}

