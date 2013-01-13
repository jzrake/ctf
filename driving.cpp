
#include <iostream>
#include <typeinfo>
#include "driving.hpp"
#include "valman.hpp"
#include "eulers.hpp"
#include "srhd.hpp"
#include "rmhd.hpp"
#include "eos.hpp"
#include "mara_mpi.h"


DrivingProcedure::DrivingProcedure(StochasticVectorField *field)
  : num_dims(Mara->domain->get_Nd()),
    field(field),
    Fx(Mara->domain->GetNumberOfZones()),
    Fy(Mara->domain->GetNumberOfZones()),
    Fz(Mara->domain->GetNumberOfZones() * (num_dims == 3))
{
  CoolingRatio = 0.0;
  this->ResampleField();
}
DrivingProcedure::~DrivingProcedure()
{
  delete field;
}
StochasticVectorField *DrivingProcedure::GetField()
{
  return field;
}
double DrivingProcedure::AveragePowerInFields() const
{
  if (num_dims == 2) {
    return (Fx*Fx + Fy*Fy).sum() / Fx.size();
  }
  else {
    return (Fx*Fx + Fy*Fy + Fz*Fz).sum() / Fx.size();
  }
}
double DrivingProcedure::EnergyInjectionRate() const
{
  return InjectionRate;
}
void DrivingProcedure::ResampleField()
{
  const std::vector<int> Nx = Mara->domain->aug_shape();
  ValarrayIndexer N(Nx);

  if (num_dims == 2) {
    for (int i=0; i<Nx[0]; ++i) {
      for (int j=0; j<Nx[1]; ++j) {

        const double c = 1.0;//LightSpeed;
        const double x = Mara->domain->x_at(i);
        const double y = Mara->domain->y_at(j);

        std::vector<double> F0 = field->SampleField(x*c, y*c, 0.0);

        Fx[ N(i,j) ] = F0[0];
        Fy[ N(i,j) ] = F0[1];
      }
    }
  }
  else if (num_dims == 3) {

    for (int i=0; i<Nx[0]; ++i) {
      for (int j=0; j<Nx[1]; ++j) {
        for (int k=0; k<Nx[2]; ++k) {

          const double c = 1.0;//LightSpeed;
          const double x = Mara->domain->x_at(i);
          const double y = Mara->domain->y_at(j);
          const double z = Mara->domain->z_at(k);

          std::vector<double> F0 = field->SampleField(x*c, y*c, z*c);

          Fx[ N(i,j,k) ] = F0[0];
          Fy[ N(i,j,k) ] = F0[1];
          Fz[ N(i,j,k) ] = F0[2];
        }
      }
    }
  }
}

void DrivingProcedure::Drive(std::valarray<double> &P, double dt)
{
  if (typeid(*Mara->fluid) == typeid(AdiabaticIdealEulers)) {
    Drive_eulers(P, dt);
  }
  else if (typeid(*Mara->fluid) == typeid(AdiabaticIdealSrhd)) {
    Drive_srhd(P, dt);
  }
  else if (typeid(*Mara->fluid) == typeid(AdiabaticIdealRmhd)) {
    Drive_rmhd(P, dt);
  }
}

void DrivingProcedure::Drive_eulers(std::valarray<double> &P, double dt)
{
  typedef AdiabaticIdealEulers Eul;
  const int Nq = Mara->domain->get_Nq();

  for (size_t i=0; i<Fx.size(); ++i) {

    std::slice m(i*Nq, Nq, 1);
    std::valarray<double> P0 = P[m];

    P0[Eul::vx] += Fx[i] * dt;
    P0[Eul::vy] += Fy[i] * dt;
    if (num_dims == 3) {
      P0[Eul::vz] += Fz[i] * dt;
    }

    P[m] = P0;
  }
}

void DrivingProcedure::Drive_srhd(std::valarray<double> &P, double dt)
{
  typedef AdiabaticIdealRmhd Srhd;
  const int Nq = Mara->domain->get_Nq();

  for (size_t i=0; i<Fx.size(); ++i) {

    std::slice m(i*Nq, Nq, 1);
    std::valarray<double> P0 = P[m];

    double uf[4] = { 0.0, P0[Srhd::vx], P0[Srhd::vy], P0[Srhd::vz] };
    double v2 = uf[1]*uf[1] + uf[2]*uf[2] + uf[3]*uf[3];
    double W0 = 1.0 / sqrt(1.0 - v2);

    // To conserve particle number we to keep the quantity rho*W == D fixed,
    // i.e. rho0*W0 = rho1*W1
    // -------------------------------------------------------------------------
    uf[1] *= W0; uf[1] += Fx[i] * (dt/W0);// / LightSpeed; // uf is now a 4-velocity
    uf[2] *= W0; uf[2] += Fy[i] * (dt/W0);// / LightSpeed;
    if (num_dims == 3) {
      uf[3] *= W0; uf[3] += Fz[i] * (dt/W0);// / LightSpeed;
    }

    double u2 = uf[1]*uf[1] + uf[2]*uf[2] + uf[3]*uf[3];
    double W1 = sqrt(1.0 + u2);

    P0[Srhd::vx] = uf[1]/W1;
    P0[Srhd::vy] = uf[2]/W1;
    P0[Srhd::vz] = uf[3]/W1;
    P0[Srhd::rho] *= W0/W1; // fix the particle number

    P[m] = P0;
  }
}

void DrivingProcedure::Drive_rmhd(std::valarray<double> &P, double dt)
{
  typedef AdiabaticIdealRmhd Rmhd;
  const int Nq = Mara->domain->get_Nq();

  for (size_t i=0; i<Fx.size(); ++i) {

    std::slice m(i*Nq, Nq, 1);
    std::valarray<double> P0 = P[m];

    double uf[4] = { 0.0, P0[Rmhd::vx], P0[Rmhd::vy], P0[Rmhd::vz] };
    double v2 = uf[1]*uf[1] + uf[2]*uf[2] + uf[3]*uf[3];
    double W0 = 1.0 / sqrt(1.0 - v2);

    // To conserve particle number we to keep the quantity rho*W == D fixed,
    // i.e. rho0*W0 = rho1*W1
    // -------------------------------------------------------------------------
    uf[1] *= W0; uf[1] += Fx[i] * (dt/W0); // uf is now a 4-velocity
    uf[2] *= W0; uf[2] += Fy[i] * (dt/W0);
    if (num_dims == 3) {
      uf[3] *= W0; uf[3] += Fz[i] * (dt/W0);
    }

    double u2 = uf[1]*uf[1] + uf[2]*uf[2] + uf[3]*uf[3];
    double W1 = sqrt(1.0 + u2);

    P0[Rmhd::vx] = uf[1]/W1;
    P0[Rmhd::vy] = uf[2]/W1;
    P0[Rmhd::vz] = uf[3]/W1;
    P0[Rmhd::rho] *= W0/W1; // fix the particle number

    P[m] = P0;
  }
}
