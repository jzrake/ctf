

#include "luaU.h"
#include "mara_mpi.h"
#include "valman.hpp"
#include "eulers.hpp"
#include "srhd.hpp"
#include "rmhd.hpp"


#define TTL_ZONES (HydroModule::Mara->domain->GetGlobalNumberOfZones())


enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive


static int luaC_mean_velocity(lua_State *L);
static int luaC_mean_cons(lua_State *L);
static int luaC_mean_prim(lua_State *L);
static int luaC_mean_energies(lua_State *L);
static int luaC_mean_max_temperature(lua_State *L);
static int luaC_mean_max_sonic_mach(lua_State *L);
static int luaC_mean_min_alfvenic_mach(lua_State *L);
static int luaC_mean_max_magnetic_field(lua_State *L);
static int luaC_max_lorentz_factor(lua_State *L);
static int luaC_mean_max_divB(lua_State *L);

static inline int absolute_index_is_ghost(const int &m);


void lua_measure_load(lua_State *L)
{
  lua_register(L, "measure_mean_velocity"          , luaC_mean_velocity);
  lua_register(L, "measure_mean_cons"              , luaC_mean_cons);
  lua_register(L, "measure_mean_prim"              , luaC_mean_prim);
  lua_register(L, "measure_mean_energies"          , luaC_mean_energies);
  lua_register(L, "measure_mean_max_temperature"   , luaC_mean_max_temperature);
  lua_register(L, "measure_mean_max_sonic_mach"    , luaC_mean_max_sonic_mach);
  lua_register(L, "measure_mean_min_alfvenic_mach" , luaC_mean_min_alfvenic_mach);
  lua_register(L, "measure_mean_max_magnetic_field", luaC_mean_max_magnetic_field);
  lua_register(L, "measure_max_lorentz_factor"     , luaC_max_lorentz_factor);
  lua_register(L, "measure_mean_max_divB"          , luaC_mean_max_divB);
}



int luaC_mean_velocity(lua_State *L)
{
  const int Nq                   = HydroModule::Mara->domain->get_Nq();
  const std::valarray<double> &P = HydroModule::Mara->PrimitiveArray;

  double Vrms[4]={0,0,0,0};

  for (size_t m=0; m<P.size(); m+=Nq) {
    if (absolute_index_is_ghost(m/Nq)) continue;

    const std::slice M(m, Nq, 1);
    const std::valarray<double> P0 = P[M];
    const double v2 = P0[vx]*P0[vx] + P0[vy]*P0[vy] + P0[vz]*P0[vz];

    Vrms[0] += 1.0 / sqrt(1.0 - v2);
    Vrms[1] += fabs(P0[vx]);
    Vrms[2] += fabs(P0[vy]);
    Vrms[3] += fabs(P0[vz]);
  }

  for (int i=0; i<4; ++i) {
    Vrms[i] = Mara_mpi_dbl_sum(Vrms[i]) / TTL_ZONES;
  }

  luaU_pusharray(L, Vrms, 4);
  return 1;
}

int luaC_mean_cons(lua_State *L)
{
  const FluidEquations &fluid    = *HydroModule::Mara->fluid;
  const int Nq                   =  HydroModule::Mara->domain->get_Nq();
  const std::valarray<double> &P =  HydroModule::Mara->PrimitiveArray;

  double *U0   = new double[Nq];
  double *Uavg = new double[Nq];
  for (int i=0; i<Nq; ++i) Uavg[i] = 0.0;

  for (size_t m=0; m<P.size(); m+=Nq) {
    if (absolute_index_is_ghost(m/Nq)) continue;

    const std::slice M(m, Nq, 1);
    const std::valarray<double> P0 = P[M];

    fluid.PrimToCons(&P0[0], &U0[0]);

    for (int i=0; i<Nq; ++i) {
      Uavg[i] += U0[i];
    }
  }

  for (int i=0; i<Nq; ++i) {
    Uavg[i] = Mara_mpi_dbl_sum(Uavg[i]) / TTL_ZONES;
  }

  luaU_pusharray(L, Uavg, Nq);
  delete [] U0;
  delete [] Uavg;
  return 1;
}

int luaC_mean_prim(lua_State *L)
{
  const int Nq                   = HydroModule::Mara->domain->get_Nq();
  const std::valarray<double> &P = HydroModule::Mara->PrimitiveArray;

  double *Pavg = new double[Nq];
  for (int i=0; i<Nq; ++i) Pavg[i] = 0.0;

  for (size_t m=0; m<P.size(); m+=Nq) {
    if (absolute_index_is_ghost(m/Nq)) continue;

    const std::slice M(m, Nq, 1);
    const std::valarray<double> P0 = P[M];

    for (int i=0; i<Nq; ++i) {
      Pavg[i] += P0[i];
    }
  }

  for (int i=0; i<Nq; ++i) {
    Pavg[i] = Mara_mpi_dbl_sum(Pavg[i]) / TTL_ZONES;
  }

  luaU_pusharray(L, Pavg, Nq);
  delete [] Pavg;
  return 1;
}

int luaC_mean_energies(lua_State *L)
{
  const FluidEquations &fluid    = *HydroModule::Mara->fluid;
  const EquationOfState &eos     = *HydroModule::Mara->eos;
  const int Nq                   =  HydroModule::Mara->domain->get_Nq();
  const std::valarray<double> &P =  HydroModule::Mara->PrimitiveArray;

  double kinetic  = 0.0;
  double internal = 0.0;
  double magnetic = 0.0;
  double total    = 0.0;

  for (size_t m=0; m<P.size(); m+=Nq) {
    if (absolute_index_is_ghost(m/Nq)) continue;

    const std::slice M(m, Nq, 1);
    const std::valarray<double> P0 = P[M];
    double *U0 = new double[Nq];

    fluid.PrimToCons(&P0[0], &U0[0]);

    const double v2 = P0[vx]*P0[vx] + P0[vy]*P0[vy] + P0[vz]*P0[vz];
    const double T0 = eos.Temperature_p(P0[rho], P0[pre]);
    const double u0 = eos.Internal(P0[rho], T0);

    if (typeid(fluid) == typeid(AdiabaticIdealSrhd) ||
        typeid(fluid) == typeid(AdiabaticIdealRmhd)) {

      const double rhoh = P0[rho] + u0 + P0[pre];
      const double W0   = 1.0 / sqrt(1.0 - v2);
      const double e_f  = rhoh * W0*W0 - P0[pre] - U0[ddd];
      const double e_m  = U0[tau] - e_f;

      magnetic += e_m;
      kinetic  += P0[rho] * W0 * (W0-1);
      internal += W0 * u0;
      total    += U0[tau];
    }
    else if (typeid(fluid) == typeid(AdiabaticIdealEulers)) {

      magnetic += 0.0;
      kinetic  += 0.5 * P0[rho] * v2;
      internal += u0;
      total    += U0[tau];
    }
    delete [] U0;
  }

  magnetic = Mara_mpi_dbl_sum(magnetic) / TTL_ZONES;
  kinetic  = Mara_mpi_dbl_sum(kinetic ) / TTL_ZONES;
  internal = Mara_mpi_dbl_sum(internal) / TTL_ZONES;
  total    = Mara_mpi_dbl_sum(total   ) / TTL_ZONES;

  lua_newtable(L);
  lua_pushstring(L, "magnetic"); lua_pushnumber(L, magnetic); lua_settable(L, -3);
  lua_pushstring(L, "kinetic");  lua_pushnumber(L, kinetic);  lua_settable(L, -3);
  lua_pushstring(L, "internal"); lua_pushnumber(L, internal); lua_settable(L, -3);
  lua_pushstring(L, "total");    lua_pushnumber(L, total);    lua_settable(L, -3);

  return 1;
}

int luaC_mean_max_temperature(lua_State *L)
{
  const EquationOfState &eos     = *HydroModule::Mara->eos;
  const int Nq                   =  HydroModule::Mara->domain->get_Nq();
  const std::valarray<double> &P =  HydroModule::Mara->PrimitiveArray;

  double ttl=0.0, max=0.0;

  for (size_t m=0; m<P.size(); m+=Nq) {
    if (absolute_index_is_ghost(m/Nq)) continue;

    const std::slice M(m, Nq, 1);
    const std::valarray<double> P0 = P[M];
    const double T0 = eos.TemperatureMeV(P0[rho], P0[pre]);
    ttl += T0;
    if (T0 > max) max = T0;
  }

  lua_pushnumber(L, Mara_mpi_dbl_sum(ttl) / TTL_ZONES);
  lua_pushnumber(L, Mara_mpi_dbl_max(max));
  return 2;
}

int luaC_mean_max_sonic_mach(lua_State *L)
{
  const FluidEquations &fluid    = *HydroModule::Mara->fluid;
  const EquationOfState &eos     = *HydroModule::Mara->eos;
  const int Nq                   =  HydroModule::Mara->domain->get_Nq();
  const std::valarray<double> &P =  HydroModule::Mara->PrimitiveArray;

  double ttl=0.0, max=0.0;

  for (size_t m=0; m<P.size(); m+=Nq) {
    if (absolute_index_is_ghost(m/Nq)) continue;

    const std::slice M(m, Nq, 1);
    const std::valarray<double> P0 = P[M];
    const double v2 = P0[vx]*P0[vx] + P0[vy]*P0[vy] + P0[vz]*P0[vz];
    const double T0 = eos.Temperature_p(P0[rho], P0[pre]);

    if (typeid(fluid) == typeid(AdiabaticIdealSrhd) ||
        typeid(fluid) == typeid(AdiabaticIdealRmhd)) {

      const double cs2 = eos.SoundSpeed2Sr(P0[rho], T0);
      const double M0 = sqrt((v2/(1-v2)) / (cs2/(1-cs2)));
      ttl += M0;
      if (M0 > max) max = M0;
    }
    else if (typeid(fluid) == typeid(AdiabaticIdealEulers)) {

      const double cs2 = eos.SoundSpeed2Nr(P0[rho], T0);
      const double M0 = sqrt(v2/cs2);
      ttl += M0;
      if (M0 > max) max = M0;
    }
  }

  lua_pushnumber(L, Mara_mpi_dbl_sum(ttl) / TTL_ZONES);
  lua_pushnumber(L, Mara_mpi_dbl_max(max));
  return 2;
}

int luaC_mean_min_alfvenic_mach(lua_State *L)
{
  const FluidEquations &fluid    = *HydroModule::Mara->fluid;
  const EquationOfState &eos     = *HydroModule::Mara->eos;
  const int Nq                   =  HydroModule::Mara->domain->get_Nq();
  const std::valarray<double> &P =  HydroModule::Mara->PrimitiveArray;

  if (typeid(fluid) != typeid(AdiabaticIdealRmhd)) {
    lua_pushnumber(L, 0.0);
    return 1;
  }

  double ttl=0.0, min=1e20;

  for (size_t m=0; m<P.size(); m+=Nq) {
    if (absolute_index_is_ghost(m/Nq)) continue;

    const std::slice M(m, Nq, 1);
    const std::valarray<double> P0 = P[M];

    const double v2 = P0[vx]*P0[vx] + P0[vy]*P0[vy] + P0[vz]*P0[vz];
    const double B2 = P0[Bx]*P0[Bx] + P0[By]*P0[By] + P0[Bz]*P0[Bz];
    const double T0 = eos.Temperature_p(P0[rho], P0[pre]);
    const double u0 = eos.Internal(P0[rho], T0);
    const double rhoh = P0[rho] + u0 + P0[pre];
    const double va2 = B2 / (B2 + rhoh);
    const double M0 = sqrt((v2/(1-v2)) / (va2/(1-va2)));
    ttl += M0;
    if (M0 < min) min = M0;
  }

  lua_pushnumber(L, Mara_mpi_dbl_sum(ttl) / TTL_ZONES);
  lua_pushnumber(L, Mara_mpi_dbl_min(min));
  return 2;
}

int luaC_mean_max_magnetic_field(lua_State *L)
{
  const FluidEquations &fluid    = *HydroModule::Mara->fluid;
  const int Nq                   =  HydroModule::Mara->domain->get_Nq();
  const std::valarray<double> &P =  HydroModule::Mara->PrimitiveArray;

  if (typeid(fluid) != typeid(AdiabaticIdealRmhd)) {
    lua_pushnumber(L, 0.0);
    return 1;
  }

  double ttl=0.0, max=0.0;

  for (size_t m=0; m<P.size(); m+=Nq) {
    if (absolute_index_is_ghost(m/Nq)) continue;

    const std::slice M(m, Nq, 1);
    const std::valarray<double> P0 = P[M];
    const double B0 = sqrt(P0[Bx]*P0[Bx] + P0[By]*P0[By] + P0[Bz]*P0[Bz]);
    ttl += B0;
    if (B0 > max) max = B0;
  }
  const PhysicalUnits *units = HydroModule::Mara->units;

  lua_pushnumber(L, (1.0/units->Gauss())*Mara_mpi_dbl_sum(ttl) / TTL_ZONES);
  lua_pushnumber(L, (1.0/units->Gauss())*Mara_mpi_dbl_max(max));
  return 2;
}

int luaC_max_lorentz_factor(lua_State *L)
{
  const FluidEquations &fluid    = *HydroModule::Mara->fluid;
  const int Nq                   =  HydroModule::Mara->domain->get_Nq();
  const std::valarray<double> &P =  HydroModule::Mara->PrimitiveArray;

  if (typeid(fluid) == typeid(AdiabaticIdealEulers)) {
    lua_pushnumber(L, 1.0);
    return 1;
  }

  double max=0.0;

  for (size_t m=0; m<P.size(); m+=Nq) {
    if (absolute_index_is_ghost(m/Nq)) continue;

    const std::slice M(m, Nq, 1);
    const std::valarray<double> P0 = P[M];
    const double v2 = P0[vx]*P0[vx] + P0[vy]*P0[vy] + P0[vz]*P0[vz];
    const double W0 = 1.0 / sqrt(1.0 - v2);
    if (W0 > max) max = W0;
  }

  lua_pushnumber(L, Mara_mpi_dbl_max(max));
  return 1;
}






int luaC_mean_max_divB(lua_State *L)
{
  const FluidEquations &fluid = *HydroModule::Mara->fluid;
  if (typeid(fluid) != typeid(AdiabaticIdealRmhd)) {
    lua_pushnumber(L, 0.0);
    lua_pushnumber(L, 0.0);
    return 2;
  }

  const std::valarray<double> &P = HydroModule::Mara->PrimitiveArray;
  const int Ng = HydroModule::Mara->domain->get_Ng();
  const int Nq = HydroModule::Mara->domain->get_Nq();
  const int Nd = HydroModule::Mara->domain->get_Nd();
  const int Nz = HydroModule::Mara->domain->GetNumberOfZones();
  const std::vector<int> N = HydroModule::Mara->domain->aug_shape();

  const std::valarray<double> fx = P[std::slice(Bx, Nz, Nq)];
  const std::valarray<double> fy = P[std::slice(By, Nz, Nq)];
  const std::valarray<double> fz = P[std::slice(Bz, Nz, Nq)];

  ValarrayIndexer M(N);
  std::valarray<double> div(Nz);

  if (Nd == 2) {
    for (int i=Ng; i<N[0]-Ng-1; ++i) {
      for (int j=Ng; j<N[1]-Ng-1; ++j) {

        div[ M(i,j) ] = (fx[ M(i+1,j  ) ] + fx[ M(i+1,j+1) ] -
                         fx[ M(i  ,j  ) ] - fx[ M(i  ,j+1) ]) / 2.0
          +             (fy[ M(i  ,j+1) ] + fy[ M(i+1,j+1) ] -
                         fy[ M(i  ,j  ) ] - fy[ M(i+1,j  ) ]) / 2.0;
      }
    }
  }
  else if (Nd == 3) {
    for (int i=Ng; i<N[0]-Ng-1; ++i) {
      for (int j=Ng; j<N[1]-Ng-1; ++j) {
        for (int k=Ng; k<N[2]-Ng-1; ++k) {

          div[ M(i,j,k) ] = ((fx[ M(i+1,j  ,k  ) ] + fx[ M(i+1,j+1,k  ) ] +
                              fx[ M(i+1,j  ,k+1) ] + fx[ M(i+1,j+1,k+1) ]) -
                             (fx[ M(i  ,j  ,k  ) ] + fx[ M(i  ,j+1,k  ) ] +
                              fx[ M(i  ,j  ,k+1) ] + fx[ M(i  ,j+1,k+1) ])) / 4.0
            +               ((fy[ M(i  ,j+1,k  ) ] + fy[ M(i  ,j+1,k+1) ] +
                              fy[ M(i+1,j+1,k  ) ] + fy[ M(i+1,j+1,k+1) ]) -
                             (fy[ M(i  ,j  ,k  ) ] + fy[ M(i  ,j  ,k+1) ] +
                              fy[ M(i+1,j  ,k  ) ] + fy[ M(i+1,j  ,k+1) ])) / 4.0
            +               ((fz[ M(i  ,j  ,k+1) ] + fz[ M(i+1,j  ,k+1) ] +
                              fz[ M(i  ,j+1,k+1) ] + fz[ M(i+1,j+1,k+1) ]) -
                             (fz[ M(i  ,j  ,k  ) ] + fz[ M(i+1,j  ,k  ) ] +
                              fz[ M(i  ,j+1,k  ) ] + fz[ M(i+1,j+1,k  ) ])) / 4.0;
        }
      }
    }
  }
  lua_pushnumber(L, Mara_mpi_dbl_sum(div.sum()/div.size()) / TTL_ZONES);
  lua_pushnumber(L, Mara_mpi_dbl_max(div.max()));

  return 2;
}





int absolute_index_is_ghost(const int &m)
{
  const int Nd = HydroModule::Mara->domain->get_Nd();
  const int Ng = HydroModule::Mara->domain->get_Ng();
  const std::vector<int> &N = HydroModule::Mara->domain->aug_shape();

  if (Nd == 1) {
    const int i = m;
    if (i < Ng || i >= N[0]-Ng) return 1;
  }
  else if (Nd == 2) {
    const int i = (m         ) / (N[1]);
    const int j = (m - i*N[1]) / (1);
    if (i < Ng || i >= N[0]-Ng ||
        j < Ng || j >= N[1]-Ng) return 1;
  }
  else if (Nd == 3) {
    const int i = (m                       ) / (N[1]*N[2]);
    const int j = (m - i*N[1]*N[2]         ) / (N[2]);
    const int k = (m - i*N[1]*N[2] - j*N[2]) / (1);
    if (i < Ng || i >= N[0]-Ng ||
        j < Ng || j >= N[1]-Ng ||
        k < Ng || k >= N[2]-Ng) return 1;
  }

  return 0;
}

