

#include "cooling.hpp"
#include "valman.hpp"
#include "eulers.hpp"
#include "srhd.hpp"
#include "rmhd.hpp"



enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive


CoolingModuleT4::CoolingModuleT4(double Tref, double t0)
  : Tref(Tref), t0(t0), energy_removed(0.0)
{

}

void CoolingModuleT4::Cool(std::valarray<double> &P, double dt)
{
  const int Nq = Mara->domain->get_Nq();
  const double dx = Mara->domain->get_dx(1);
  const double dy = Mara->domain->get_dx(2);
  const double dz = Mara->domain->get_dx(3);
  const double dV = dx*dy*dz;

  double *U0 = new double[Nq];
  double *U1 = new double[Nq];

  for (size_t i=0; i<P.size()/Nq; ++i) {

    std::slice m(i*Nq, Nq, 1);
    std::valarray<double> P0 = P[m];

    const double T0 = Mara->eos->TemperatureMeV(P0[rho], P0[pre]);
    const double T1 = T0 - dt * (Tref/t0) * pow(T0/Tref, 4);

    Mara->fluid->PrimToCons(&P0[0], &U0[0]);
    P0[pre] = Mara->eos->Pressure(P0[rho], Mara->eos->TemperatureArb(P0[rho], T1));
    Mara->fluid->PrimToCons(&P0[0], &U1[0]);

    energy_removed += (U0[tau] - U1[tau])*dV;
    P[m] = P0;
  }

  delete [] U0;
  delete [] U1;
}

double CoolingModuleT4::EnergyRemoved()
{
  const double ret = energy_removed;
  energy_removed = 0.0;
  return ret;
}



CoolingModuleE4::CoolingModuleE4(double eref, double t0)
  : eref(eref), t0(t0), energy_removed(0.0)
{

}

void CoolingModuleE4::Cool(std::valarray<double> &P, double dt)
{
  const int Nq = Mara->domain->get_Nq();
  const double dx = Mara->domain->get_dx(1);
  const double dy = Mara->domain->get_dx(2);
  const double dz = Mara->domain->get_dx(3);
  const double dV = dx*dy*dz;

  double *U0 = new double[Nq];
  double *U1 = new double[Nq];

  for (size_t i=0; i<P.size()/Nq; ++i) {

    std::slice m(i*Nq, Nq, 1);
    std::valarray<double> P0 = P[m];

    const double T0 = Mara->eos->Temperature_p(P0[rho], P0[pre]);
    const double e0 = Mara->eos->Internal(P0[rho], T0) / P0[rho];

    const double e1 = e0 - dt * (eref/t0) * pow(e0/eref, 4);
    const double T1 = Mara->eos->Temperature_u(P0[rho], P0[rho]*e1);

    Mara->fluid->PrimToCons(&P0[0], &U0[0]);
    P0[pre] = Mara->eos->Pressure(P0[rho], T1);
    Mara->fluid->PrimToCons(&P0[0], &U1[0]);

    energy_removed += (U0[tau] - U1[tau])*dV;
    P[m] = P0;
  }

  delete [] U0;
  delete [] U1;
}

double CoolingModuleE4::EnergyRemoved()
{
  const double ret = energy_removed;
  energy_removed = 0.0;
  return ret;
}


