
/*
 * -----------------------------------------------------------------------------
 *                      M A R A    H Y D R O    C O D E
 * -----------------------------------------------------------------------------
 *        __  __    __    ____    __      _   _  _  _  ____  ____  _____
 *       (  \/  )  /__\  (  _ \  /__\    ( )_( )( \/ )(  _ \(  _ \(  _  )
 *        )    (  /(__)\  )   / /(__)\    ) _ (  \  /  )(_) ))   / )(_)(
 *       (_/\/\_)(__)(__)(_)\_)(__)(__)  (_) (_) (__) (____/(_)\_)(_____)
 *                            ___  _____  ____  ____
 *                           / __)(  _  )(  _ \( ___)
 *                          ( (__  )(_)(  )(_) ))__)
 *                           \___)(_____)(____/(____)
 *
 * -----------------------------------------------------------------------------
 *
 * FILE: mara.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 *
 *------------------------------------------------------------------------------
 */

#include "config.h"
#if (__MARA_USE_MPI)
#include <mpi.h>
#endif

#include <iostream>
#include <map>
extern "C" {
#include "lualib.h"
#include "lauxlib.h"
}
#include "mara.hpp"



static void    luaU_pusharray(lua_State *L, double *A, int N);
static void    luaU_pusharray_i(lua_State *L, int *A, int N);
static double *luaU_checkarray(lua_State *L, int n);
static double *luaU_checklarray(lua_State *L, int n, int *N);
static int    *luaU_checklarray_i(lua_State *L, int n, int *N);

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



extern "C"
{
  static int luaC_Mara_start(lua_State *L);
  static int luaC_Mara_close(lua_State *L);
  static int luaC_Mara_version(lua_State *L);
  static int luaC_Mara_show(lua_State *L);
  static int luaC_Mara_advance(lua_State *L);
  static int luaC_Mara_diffuse(lua_State *L);
  static int luaC_Mara_init_prim(lua_State *L);
  static int luaC_Mara_set_gravity(lua_State *L);
  static int luaC_Mara_prim_at_point(lua_State *L);
  static int luaC_Mara_get_timestep(lua_State *L);

  static int luaC_set_domain(lua_State *L);
  static int luaC_set_boundary(lua_State *L);
  static int luaC_set_fluid(lua_State *L);
  static int luaC_set_eos(lua_State *L);
  static int luaC_set_units(lua_State *L);
  static int luaC_set_riemann(lua_State *L);
  static int luaC_set_godunov(lua_State *L);
  static int luaC_set_advance(lua_State *L);
  static int luaC_set_driving(lua_State *L);
  static int luaC_set_cooling(lua_State *L);
  static int luaC_set_primitive(lua_State *L);
  static int luaC_cooling_rate(lua_State *L);
  static int luaC_config_solver(lua_State *L);

  static int luaC_new_ou_field(lua_State *L);
  static int luaC_load_shen(lua_State *L);
  static int luaC_test_shen(lua_State *L);
  static int luaC_test_rmhd_c2p(lua_State *L);
  static int luaC_test_sampling(lua_State *L);
  static int luaC_test_sampling_many(lua_State *L);

  static int luaC_fluid_PrimToCons(lua_State *L);
  static int luaC_fluid_ConsToPrim(lua_State *L);
  static int luaC_fluid_Eigensystem(lua_State *L);
  static int luaC_fluid_FluxFunction(lua_State *L);
  static int luaC_fluid_GetPrimNames(lua_State *L);

  //  static int luaC_boundary_ApplyBoundaries(lua_State *L);

  static int luaC_driving_Advance(lua_State *L);
  static int luaC_driving_Resample(lua_State *L);
  static int luaC_driving_Serialize(lua_State *L);
  static int luaC_driving_Power(lua_State *L);

  static int luaC_eos_TemperatureMeV(lua_State *L);
  static int luaC_eos_TemperatureArb(lua_State *L);
  static int luaC_eos_Temperature_p(lua_State *L);
  static int luaC_eos_Internal(lua_State *L);
  static int luaC_eos_Pressure(lua_State *L);
  static int luaC_eos_SoundSpeed2Nr(lua_State *L);
  static int luaC_eos_SoundSpeed2Sr(lua_State *L);
  static int luaC_eos_DensUpper(lua_State *L);
  static int luaC_eos_DensLower(lua_State *L);
  static int luaC_eos_TempUpper(lua_State *L);
  static int luaC_eos_TempLower(lua_State *L);

  static int luaC_units_Print(lua_State *L);
  static int luaC_units_Gauss(lua_State *L);
  static int luaC_units_Velocity(lua_State *L);
  static int luaC_units_MeV(lua_State *L);
  static int luaC_units_MeVPerCubicFemtometer(lua_State *L);
  static int luaC_units_ErgPerCubicCentimeter(lua_State *L);
  static int luaC_units_GramsPerCubicCentimeter(lua_State *L);

  static int luaC_write_ppm(lua_State *L);

  int luaopen_Mara(lua_State *L);
}

static MaraApplication *Mara = NULL;

int luaopen_Mara(lua_State *L)
{
  luaL_Reg Mara_module[] = {
    {"start"        , luaC_Mara_start},
    {"close"        , luaC_Mara_close},
    {"version"      , luaC_Mara_version},
    {"show"         , luaC_Mara_show},
    {"advance"      , luaC_Mara_advance},
    {"diffuse"      , luaC_Mara_diffuse},
    {"init_prim"    , luaC_Mara_init_prim},
    {"set_gravity"  , luaC_Mara_set_gravity},
    {"prim_at_point", luaC_Mara_prim_at_point},
    {"get_timestep" , luaC_Mara_get_timestep},

    {"set_domain"   , luaC_set_domain},
    {"set_boundary" , luaC_set_boundary},
    {"set_fluid"    , luaC_set_fluid},
    {"set_eos"      , luaC_set_eos},
    {"set_units"    , luaC_set_units},
    {"set_godunov"  , luaC_set_godunov},
    {"set_riemann"  , luaC_set_riemann},
    {"set_advance"  , luaC_set_advance},
    {"set_driving"  , luaC_set_driving},
    {"set_cooling"  , luaC_set_cooling},
    {"set_primitive", luaC_set_primitive},
    {"config_solver", luaC_config_solver},

    {"cooling_rate" , luaC_cooling_rate},
    {"new_ou_field" , luaC_new_ou_field},
    {"load_shen"    , luaC_load_shen},
    {"test_shen"    , luaC_test_shen},
    {"test_rmhd_c2p", luaC_test_rmhd_c2p},
    {"test_sampling", luaC_test_sampling},
    {"test_sampling_many", luaC_test_sampling_many},
    {NULL, NULL}};

  luaL_Reg Mara_measure[] = {
    {"mean_cons"              , luaC_mean_cons},
    {"mean_prim"              , luaC_mean_prim},
    {"mean_energies"          , luaC_mean_energies},
    {"mean_velocity"          , luaC_mean_velocity},
    {"mean_max_temperature"   , luaC_mean_max_temperature},
    {"mean_max_sonic_mach"    , luaC_mean_max_sonic_mach},
    {"mean_min_alfvenic_mach" , luaC_mean_min_alfvenic_mach},
    {"mean_max_magnetic_field", luaC_mean_max_magnetic_field},
    {"max_lorentz_factor"     , luaC_max_lorentz_factor},
    {"mean_max_divB"          , luaC_mean_max_divB},
    {NULL, NULL}};

  luaL_Reg Mara_fluid[] = {
    {"PrimToCons", luaC_fluid_PrimToCons},
    {"ConsToPrim", luaC_fluid_ConsToPrim},
    {"Eigensystem", luaC_fluid_Eigensystem},
    {"FluxFunction", luaC_fluid_FluxFunction},
    {"GetPrimNames", luaC_fluid_GetPrimNames},
    {NULL, NULL}};

  luaL_Reg Mara_boundary[] = {
    //    {"ApplyBoundaries", luaC_boundary_ApplyBoundaries},
    {NULL, NULL}};

  luaL_Reg Mara_eos[] = {
    {"TemperatureMeV", luaC_eos_TemperatureMeV},
    {"TemperatureArb", luaC_eos_TemperatureArb},
    {"Temperature_p", luaC_eos_Temperature_p},
    {"Internal", luaC_eos_Internal},
    {"Pressure", luaC_eos_Pressure},
    {"SoundSpeed2Nr", luaC_eos_SoundSpeed2Nr},
    {"SoundSpeed2Sr", luaC_eos_SoundSpeed2Sr},
    {"TempUpper", luaC_eos_TempUpper},
    {"TempLower", luaC_eos_TempLower},
    {"DensUpper", luaC_eos_DensUpper},
    {"DensLower", luaC_eos_DensLower},
    {NULL, NULL}};

  luaL_Reg Mara_units[] = {
    {"Print", luaC_units_Print},
    {"Gauss", luaC_units_Gauss},
    {"Velocity", luaC_units_Velocity},
    {"MeV", luaC_units_MeV},
    {"MeVPerCubicFemtometer", luaC_units_MeVPerCubicFemtometer},
    {"ErgPerCubicCentimeter", luaC_units_ErgPerCubicCentimeter},
    {"GramsPerCubicCentimeter", luaC_units_GramsPerCubicCentimeter},
    {NULL, NULL}};

  luaL_Reg Mara_driving[] = {
    {"Advance", luaC_driving_Advance},
    {"Resample", luaC_driving_Resample},
    {"Serialize", luaC_driving_Serialize},
    {"Power", luaC_driving_Power},
    {NULL, NULL}};

  luaL_Reg Mara_image[] = {
    {"write_ppm", luaC_write_ppm},
    {NULL, NULL}};


  lua_newtable(L); // Mara
  luaL_setfuncs(L, Mara_module, 0);

  lua_newtable(L); // Mara.measure
  luaL_setfuncs(L, Mara_measure, 0);
  lua_setfield(L, -2, "measure");

  lua_newtable(L); // Mara.fluid
  luaL_setfuncs(L, Mara_fluid, 0);
  lua_setfield(L, -2, "fluid");

  lua_newtable(L); // Mara.boundary
  luaL_setfuncs(L, Mara_boundary, 0);
  lua_setfield(L, -2, "boundary");

  lua_newtable(L); // Mara.eos
  luaL_setfuncs(L, Mara_eos, 0);
  lua_setfield(L, -2, "eos");

  lua_newtable(L); // Mara.units
  luaL_setfuncs(L, Mara_units, 0);
  lua_setfield(L, -2, "units");

  lua_newtable(L); // Mara.driving
  luaL_setfuncs(L, Mara_driving, 0);
  lua_setfield(L, -2, "driving");

  lua_newtable(L); // Mara.image
  luaL_setfuncs(L, Mara_image, 0);
  lua_setfield(L, -2, "image");

  return 1;
}


int luaC_Mara_start(lua_State *L)
{
  printf("\n");
  printf("\t***********************************\n");
  printf("\tMara Astrophysical gasdynamics code\n");
  printf("\t(C) Jonathan Zrake, NYU CCPP\n");
  printf("\tVersion %s\n", __MARA_BASE_VERSION);
  printf("\t***********************************\n\n");

  if (Mara != NULL) {
    delete Mara;
  }

  Mara = new MaraApplication;
  HydroModule::Mara = Mara;
  return 0;
}
int luaC_Mara_close(lua_State *L)
{
  if (Mara != NULL) {
    delete Mara;
    Mara = NULL;
  }
  return 0;
}


int luaC_Mara_advance(lua_State *L)
{
  if (Mara->domain == NULL) {
    luaL_error(L, "[Mara] need a domain to run this, use set_domain");
  }

  clock_t start = clock();
  double *UserPrim = (double*) lua_touserdata(L, 1);
  double dt = luaL_checknumber(L, 2);
  unsigned N = lua_rawlen(L, 1) / sizeof(double); // total number of zones
  unsigned N_needed = Mara->domain->GetNumberOfZones() * Mara->domain->get_Nq();

  if (N != N_needed) {
    luaL_error(L, "[Mara] primitive has the wrong size for the domain");
  }

  std::valarray<double> P(N);
  std::valarray<double> U(N);
  int errors;

  memcpy(&P[0], UserPrim, N * sizeof(double));

  RiemannSolver::ResetMaxLambda();
  Mara->FailureMask.resize(Mara->domain->GetNumberOfZones());
  Mara->PrimitiveArray.resize(P.size());
  Mara->godunov->PrimToCons(P, U);
  try {
    Mara->advance->AdvanceState(U, dt);
    errors = Mara->godunov->ConsToPrim(U, P);
  }
  catch (const GodunovOperator::IntermediateFailure &e) {
    errors = Mara_mpi_int_sum(Mara->FailureMask.sum());
  }
  if (errors == 0) {
    if (Mara->driving) Mara->driving->Drive(P, dt);
    if (Mara->cooling) Mara->cooling->Cool(P, dt);
    memcpy(UserPrim, &P[0], N * sizeof(double));
  }

  const double sec = (double) (clock() - start) / CLOCKS_PER_SEC;
  lua_pushnumber(L, 1e-3*Mara->domain->GetNumberOfZones()/sec);
  lua_pushnumber(L, errors);

  return 2;
}

int luaC_Mara_diffuse(lua_State *L)
{
  double *UserPrim = (double*) lua_touserdata(L, 1);
  double r = luaL_checknumber(L, 2);
  unsigned N = lua_rawlen(L, 1) / sizeof(double); // total number of zones
  unsigned N_needed = Mara->domain->GetNumberOfZones() * Mara->domain->get_Nq();

  if (N != N_needed) {
    luaL_error(L, "[Mara] primitive has the wrong size for the domain");
  }

  std::valarray<double> P(N);
  std::valarray<double> U(N);
  memcpy(&P[0], UserPrim, N * sizeof(double));

  Mara->godunov->PrimToCons(P, U);
  U += Mara->godunov->LaxDiffusion(U, r);
  Mara->godunov->ConsToPrim(U, P);
  memcpy(UserPrim, &P[0], N * sizeof(double));

  return 0;
}

int luaC_Mara_get_timestep(lua_State *L)
{
  const double CFL = luaL_checknumber(L, 1);
  const double dt = Mara_mpi_dbl_min(CFL * Mara->domain->get_min_dx() /
                                     RiemannSolver::GetMaxLambda());
  lua_pushnumber(L, dt);
  return 1;
}

int luaC_Mara_version(lua_State *L)
{
  char str[256];
  sprintf(str, "%s", __MARA_BASE_VERSION);
  lua_pushstring(L, str);
  return 1;
}

std::string Demangle(const std::string &mname);
int luaC_Mara_show(lua_State *L)
{
  std::cout << std::endl;
  std::cout << "\tunits: "
            << (Mara->units ? Demangle(typeid(*Mara->units).name()) : "NULL")
            << std::endl;
  std::cout << "\tdomain: "
            << (Mara->domain ? Demangle(typeid(*Mara->domain).name()) : "NULL")
            << std::endl;
  std::cout << "\tboundary: "
            << (Mara->boundary ? Demangle(typeid(*Mara->boundary).name()) : "NULL")
            << std::endl;
  std::cout << "\tfluid: "
            << (Mara->fluid ? Demangle(typeid(*Mara->fluid).name()) : "NULL")
            << std::endl;
  std::cout << "\teos: "
            << (Mara->eos ? Demangle(typeid(*Mara->eos).name()) : "NULL")
            << std::endl;
  std::cout << "\tgodunov: "
            << (Mara->godunov ? Demangle(typeid(*Mara->godunov).name()) : "NULL")
            << std::endl;
  std::cout << "\triemann: "
            << (Mara->riemann ? Demangle(typeid(*Mara->riemann).name()) : "NULL")
            << std::endl;
  std::cout << "\tadvance: "
            << (Mara->advance ? Demangle(typeid(*Mara->advance).name()) : "NULL")
            << std::endl;
  std::cout << "\tdriving: "
            << (Mara->driving ? Demangle(typeid(*Mara->driving).name()) : "NULL")
            << std::endl;
  std::cout << "\tcooling: "
            << (Mara->cooling ? Demangle(typeid(*Mara->cooling).name()) : "NULL")
            << std::endl;
  std::cout << std::endl;

  return 0;
}


EquationOfState *BuildAdiabaticEos(lua_State *L)
{
  const double Gamma = luaL_checknumber(L, 2);
  return new AdiabaticEos(Gamma);
}

EquationOfState *BuildThermalBarotropicEos(lua_State *L)
{
  const double GammaT = luaL_checknumber(L, 2);
  const double GammaB = luaL_checknumber(L, 3);
  const double Kappa  = luaL_checknumber(L, 4);
  return new ThermalBarotropicEos(GammaT, GammaB, Kappa);
}

EquationOfState *BuildShenTabulatedNuclearEos(lua_State *L)
{
  int N;
  TabulatedEos tab;

  lua_pushstring(L, "logD_values");
  lua_gettable(L, -2);
  double *logD_values = luaU_checklarray(L, -1, &N);
  lua_pop(L, 1);
  tab.logD_values.assign(logD_values, logD_values+N);
  free(logD_values);

  lua_pushstring(L, "logT_values");
  lua_gettable(L, -2);
  double *logT_values = luaU_checklarray(L, -1, &N);
  lua_pop(L, 1);
  tab.logT_values.assign(logT_values, logT_values+N);
  free(logT_values);

  lua_pushstring(L, "EOS_p");
  lua_gettable(L, -2);
  double *EOS_p = luaU_checklarray(L, -1, &N);
  lua_pop(L, 1);
  tab.EOS_p.assign(EOS_p, EOS_p+N);
  free(EOS_p);

  lua_pushstring(L, "EOS_s");
  lua_gettable(L, -2);
  double *EOS_s = luaU_checklarray(L, -1, &N);
  lua_pop(L, 1);
  tab.EOS_s.assign(EOS_s, EOS_s+N);
  free(EOS_s);

  lua_pushstring(L, "EOS_u");
  lua_gettable(L, -2);
  double *EOS_u = luaU_checklarray(L, -1, &N);
  lua_pop(L, 1);
  tab.EOS_u.assign(EOS_u, EOS_u+N);
  free(EOS_u);

  return new ShenTabulatedNuclearEos(tab);
}

EquationOfState *BuildGenericTabulatedEos(lua_State *L)
{
  int N;
  double *tmp;

  lua_getfield(L, -1, "D");
  tmp = luaU_checklarray(L, -1, &N);
  std::vector<double> D_values(tmp, tmp + N);
  lua_pop(L, 1);
  free(tmp);

  lua_getfield(L, -1, "T");
  tmp = luaU_checklarray(L, -1, &N);
  std::vector<double> T_values(tmp, tmp + N);
  lua_pop(L, 1);
  free(tmp);

  lua_getfield(L, -1, "p");
  tmp = luaU_checklarray(L, -1, &N);
  std::vector<double> EOS_p(tmp, tmp + N);
  lua_pop(L, 1);
  free(tmp);

  lua_getfield(L, -1, "u");
  tmp = luaU_checklarray(L, -1, &N);
  std::vector<double> EOS_u(tmp, tmp + N);
  lua_pop(L, 1);
  free(tmp);

  lua_getfield(L, -1, "c");
  tmp = luaU_checklarray(L, -1, &N);
  std::vector<double> EOS_c(tmp, tmp + N);
  lua_pop(L, 1);
  free(tmp);

  return new GenericTabulatedEos(D_values, T_values, EOS_p, EOS_u, EOS_c);
}


int luaC_Mara_set_gravity(lua_State *L)
{
  if (Mara->domain == NULL) {
    luaL_error(L, "[Mara] need a domain to run this, use set_domain");
  }
  if (lua_type(L, 1) != LUA_TFUNCTION) {
    luaL_error(L, "[Mara] argument must be a function");
  }

  const PhysicalDomain *domain = Mara->domain;
  const int Ng = domain->get_Ng();
  const std::vector<int> Ninter(domain->GetLocalShape(),
                                domain->GetLocalShape()+domain->get_Nd());
  int ngrav = 4;
  ValarrayManager M(domain->aug_shape(), ngrav);
  Mara->GravityArray.resize(domain->GetNumberOfZones() * ngrav);

  switch (domain->get_Nd()) {
  case 1:

    for (int i=0; i<Ninter[0]+2*Ng; ++i) {

      const double x = domain->x_at(i);
      const double y = 0.0;
      const double z = 0.0;

      lua_pushvalue(L, 1);
      lua_pushnumber(L, x);
      lua_pushnumber(L, y);
      lua_pushnumber(L, z);
      lua_call(L, 3, 1);

      double *P0 = luaU_checkarray(L, 2);
      std::valarray<double> P(P0, ngrav);
      free(P0);

      Mara->GravityArray[ M(i) ] = P;
      lua_pop(L, 1);
    }
    break;

  case 2:
    for (int i=0; i<Ninter[0]+2*Ng; ++i) {
      for (int j=0; j<Ninter[1]+2*Ng; ++j) {

        const double x = domain->x_at(i);
        const double y = domain->y_at(j);
        const double z = 0.0;

	lua_pushvalue(L, 1);
	lua_pushnumber(L, x);
	lua_pushnumber(L, y);
	lua_pushnumber(L, z);
	lua_call(L, 3, 1);

	double *P0 = luaU_checkarray(L, 2);
	std::valarray<double> P(P0, ngrav);
	free(P0);

	Mara->GravityArray[ M(i,j) ] = P;
	lua_pop(L, 1);
      }
    }
    break;

  case 3:
    for (int i=0; i<Ninter[0]+2*Ng; ++i) {
      for (int j=0; j<Ninter[1]+2*Ng; ++j) {
        for (int k=0; k<Ninter[2]+2*Ng; ++k) {

          const double x = domain->x_at(i);
          const double y = domain->y_at(j);
          const double z = domain->z_at(k);

	  lua_pushvalue(L, 1);
	  lua_pushnumber(L, x);
	  lua_pushnumber(L, y);
	  lua_pushnumber(L, z);
	  lua_call(L, 3, 1);

	  double *P0 = luaU_checkarray(L, 2);
	  std::valarray<double> P(P0, ngrav);
	  free(P0);
	  
	  Mara->GravityArray[ M(i,j,k) ] = P;
	  lua_pop(L, 1);
        }
      }
    }
    break;
  }

  FILE *fpre = fopen("fpre.dat", "w");
  FILE *fphi = fopen("fphi.dat", "w");
  FILE *gphx = fopen("gphx.dat", "w");
  FILE *gphy = fopen("gphy.dat", "w");

  ValarrayManager PM(domain->aug_shape(), 5);

  for (int i=0; i<Ninter[0]+2*Ng; ++i) {
    for (int j=0; j<Ninter[1]+2*Ng; ++j) {
      std::valarray<double> P = Mara->PrimitiveArray[ PM(i,j) ];
      std::valarray<double> G = Mara->GravityArray[ M(i,j) ];

      fprintf(fpre, "%14.12f ", P[1]);
      fprintf(fphi, "%14.12f ", G[0]);
      fprintf(gphx, "%14.12f ", G[1]);
      fprintf(gphy, "%14.12f ", G[2]);

    }
    fprintf(fpre, "\n");
    fprintf(fphi, "\n");
    fprintf(gphx, "\n");
    fprintf(gphy, "\n");
  }

  fclose(fpre);
  fclose(fphi);
  fclose(gphx);
  fclose(gphy);
  return 0;
}


int luaC_Mara_init_prim(lua_State *L)
{
  if (Mara->domain == NULL) {
    luaL_error(L, "[Mara] need a domain to run this, use set_domain");
  }
  if (lua_type(L, 2) != LUA_TFUNCTION) {
    luaL_error(L, "[Mara] argument must be a function");
  }

  const PhysicalDomain *domain = Mara->domain;
  double *UserPrim = (double*) lua_touserdata(L, 1);
  unsigned N_have = lua_rawlen(L, 1) / sizeof(double); // total number of zones
  unsigned N_needed = Mara->domain->GetNumberOfZones() * Mara->domain->get_Nq();

  if (N_have != N_needed) {
    luaL_error(L, "[Mara] primitive has the wrong size for the domain");
  }

  const int Ng = domain->get_Ng();
  const std::vector<int> Ninter(domain->GetLocalShape(),
                                domain->GetLocalShape()+domain->get_Nd());
  ValarrayIndexer N(Ninter);
  ValarrayManager M(domain->aug_shape(), domain->get_Nq());
  Mara->PrimitiveArray.resize(domain->GetNumberOfZones() * domain->get_Nq());

  switch (domain->get_Nd()) {
  case 1:

    for (int i=0; i<Ninter[0]; ++i) {

      const double x = domain->x_at(i+Ng);
      const double y = 0.0;
      const double z = 0.0;

      lua_pushvalue(L, 2);
      lua_pushnumber(L, x);
      lua_pushnumber(L, y);
      lua_pushnumber(L, z);
      lua_call(L, 3, 1);

      double *P0 = luaU_checkarray(L, 3);
      std::valarray<double> P(P0, domain->get_Nq());
      free(P0);

      Mara->PrimitiveArray[ M(i+Ng) ] = P;
      lua_pop(L, 1);
    }
    break;

  case 2:
    for (int i=0; i<Ninter[0]; ++i) {
      for (int j=0; j<Ninter[1]; ++j) {

        const double x = domain->x_at(i+Ng);
        const double y = domain->y_at(j+Ng);
        const double z = 0.0;

	lua_pushvalue(L, 2);
	lua_pushnumber(L, x);
	lua_pushnumber(L, y);
	lua_pushnumber(L, z);
	lua_call(L, 3, 1);

	double *P0 = luaU_checkarray(L, 3);
	std::valarray<double> P(P0, domain->get_Nq());
	free(P0);

	Mara->PrimitiveArray[ M(i+Ng,j+Ng) ] = P;
	lua_pop(L, 1);
      }
    }
    break;

  case 3:
    for (int i=0; i<Ninter[0]; ++i) {
      for (int j=0; j<Ninter[1]; ++j) {
        for (int k=0; k<Ninter[2]; ++k) {

          const double x = domain->x_at(i+Ng);
          const double y = domain->y_at(j+Ng);
          const double z = domain->z_at(k+Ng);

	  lua_pushvalue(L, 2);
	  lua_pushnumber(L, x);
	  lua_pushnumber(L, y);
	  lua_pushnumber(L, z);
	  lua_call(L, 3, 1);

	  double *P0 = luaU_checkarray(L, 3);
	  std::valarray<double> P(P0, domain->get_Nq());
	  free(P0);
	  
	  Mara->PrimitiveArray[ M(i+Ng,j+Ng,k+Ng) ] = P;
	  lua_pop(L, 1);
        }
      }
    }
    break;
  }
  memcpy(UserPrim, &Mara->PrimitiveArray[0], N_have * sizeof(double));
  return 0;
}

int luaC_Mara_prim_at_point(lua_State *L)
{
  double *r1 = luaU_checkarray(L, 1);
  int Nq = Mara->domain->get_Nq();
  double *P1 = (double*) malloc(Nq * sizeof(double));
  Mara_prim_at_point(r1, P1);
  luaU_pusharray(L, P1, Nq);
  free(r1);
  free(P1);
  return 1;
}

int luaC_test_sampling(lua_State *L)
{
  if (Mara->domain == NULL) {
    luaL_error(L, "[Mara] need a domain to run this, use set_domain");
  }
  if (Mara->domain->get_Nd() != 3) {
    luaL_error(L, "[Mara] need a 3d domain to run this");
  }

  const int numsamp = luaL_checkinteger(L, 1);
  const int Nq = Mara->domain->get_Nq();

  const double *gx0 = Mara->domain->GetGlobalX0();
  const double *gx1 = Mara->domain->GetGlobalX1();

  double *P1 = new double[Nq];
  RandomNumberStream rand;

  const clock_t start = clock();

  for (int i=0; i<numsamp; ++i) {
    double r1[3] = { rand.RandomDouble(gx0[0], gx1[0]),
                     rand.RandomDouble(gx0[1], gx1[1]),
                     rand.RandomDouble(gx0[2], gx1[2]) };
    Mara_prim_at_point(r1, P1);
  }

  const double trun = (double) (clock() - start) / CLOCKS_PER_SEC;
  lua_pushnumber(L, trun);

  delete [] P1;
  return 1;
}

int luaC_test_sampling_many(lua_State *L)
{
  if (Mara->domain == NULL) {
    luaL_error(L, "[Mara] need a domain to run this, use set_domain");
  }
  if (Mara->domain->get_Nd() != 3) {
    luaL_error(L, "[Mara] need a 3d domain to run this");
  }

  const int numsamp = luaL_checkinteger(L, 1);
  const int Nq = Mara->domain->get_Nq();

  const double *gx0 = Mara->domain->GetGlobalX0();
  const double *gx1 = Mara->domain->GetGlobalX1();

  RandomNumberStream rand;
  double *Rinpt = new double[ 3*numsamp];
  double *Rlist = new double[ 3*numsamp];
  double *Plist = new double[Nq*numsamp];

  for (int i=0; i<numsamp; ++i) {
    Rinpt[3*i + 0] = rand.RandomDouble(gx0[0], gx1[0]);
    Rinpt[3*i + 1] = rand.RandomDouble(gx0[1], gx1[1]);
    Rinpt[3*i + 2] = rand.RandomDouble(gx0[2], gx1[2]);
  }

  const clock_t start = clock();

  Mara_prim_at_point_many(Rinpt, Rlist, Plist, numsamp);

  for (int i=0; i<numsamp; ++i) {
    //    printf("(%f %f %f)\n", Rinpt[3*i+0], Rinpt[3*i+1], Rinpt[3*i+2]);
    //    printf("(%f %f %f)\n", Rlist[3*i+0], Rlist[3*i+1], Rlist[3*i+2]);
    //    std::cout << Mara->fluid->PrintPrim(&Plist[Nq*i]) << std::endl;
  }

  const double trun = (double) (clock() - start) / CLOCKS_PER_SEC;
  lua_pushnumber(L, trun);

  delete [] Rinpt;
  delete [] Rlist;
  delete [] Plist;

  return 1;
}

int luaC_set_domain(lua_State *L)
{
  int Nd;
  double *x0 = luaU_checkarray(L, 1);
  double *x1 = luaU_checkarray(L, 2);

  int *N = luaU_checklarray_i(L, 3, &Nd);
  int Nq = luaL_checkinteger(L, 4);
  int Ng = luaL_checkinteger(L, 5);
  void *cart_comm = lua_touserdata(L, 6); // optional, OK if NULL

  PhysicalDomain *new_f = NULL;

  if (Mara_mpi_active()) {
    new_f = new DecomposedCartesianDomain(x0, x1, N, Nd, Nq, Ng, cart_comm);
  }
  else {
    new_f = new SimpleCartesianDomain(x0, x1, N, Nd, Nq, Ng);
  }

  if (new_f) {
    if (Mara->domain) delete Mara->domain;
    Mara->domain = new_f;
  }

  free(x0);
  free(x1);
  free(N);
  return 0;
}

int luaC_write_ppm(lua_State *L)
{
  if (Mara->domain == NULL) {
    luaL_error(L, "[Mara] need a domain to run this, use set_domain");
  }
  if (Mara->domain->get_Nd() != 2) {
    luaL_error(L, "[Mara] need a 2d domain to write ppm images");
  }

  const int narg = lua_gettop(L);
  const char *fname = luaL_checkstring(L, 1);
  double *data = luaU_checkarray(L, 2);
  int cmap = narg == 3 ? luaL_checkinteger(L, 3) : 0;
  double *range = narg == 4 ? luaU_checkarray(L, 4) : NULL;

  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);

  Mara_image_write_ppm(fname, data, cmap, Nx, Ny, range);
  free(data);
  free(range);
  return 0;
}

int luaC_set_fluid(lua_State *L)
{
  const char *key = luaL_checkstring(L, 1);
  FluidEquations *new_f = NULL;

  if (strcmp("euler", key) == 0) {
    new_f = new AdiabaticIdealEulers;
  }
  else if (strcmp("srhd", key) == 0) {
    new_f = new AdiabaticIdealSrhd;
  }
  else if (strcmp("rmhd", key) == 0) {
    new_f = new AdiabaticIdealRmhd;
  }

  if (new_f) {
    if (Mara->fluid) delete Mara->fluid;
    Mara->fluid = new_f;
  }
  else {
    luaL_error(L, "[Mara] no such fluid '%s'", key);
  }

  return 0;
}


int luaC_set_eos(lua_State *L)
{
  const char *key = luaL_checkstring(L, 1);
  EquationOfState *new_f = NULL;

  if (strcmp("gamma-law", key) == 0) {
    new_f = BuildAdiabaticEos(L);
  }
  else if (strcmp("gamma-law + polytrope", key) == 0) {
    new_f = BuildThermalBarotropicEos(L);
  }
  else if (strcmp("shen", key) == 0) {
    new_f = BuildShenTabulatedNuclearEos(L);
  }
  else if (strcmp("tabulated", key) == 0) {
    new_f = BuildGenericTabulatedEos(L);
  }

  if (new_f) {
    if (Mara->eos) delete Mara->eos;
    Mara->eos = new_f;
  }
  else {
    luaL_error(L, "[Mara] no such eos '%s'", key);
  }

  return 0;
}



int luaC_set_units(lua_State *L)
{
  const double Mass   = luaL_checknumber(L, 1);
  const double Length = luaL_checknumber(L, 2);
  const double Time   = luaL_checknumber(L, 3);

  PhysicalUnits *new_f = new PhysicalUnits(Mass, Length, Time);

  if (new_f) {
    if (Mara->units) delete Mara->units;
    Mara->units = new_f;
  }
  return 0;
}


int luaC_set_boundary(lua_State *L)
{
  const char *key = luaL_checkstring(L, 1);
  BoundaryConditions *new_f = NULL;

  if (strcmp(key, "periodic") == 0) {
    new_f = new PeriodicBoundary;
  }
  else if (strcmp(key, "outflow") == 0) {
    new_f = new OutflowBoundary;
  }
  else if (strcmp(key, "perxouty") == 0) {
    new_f = new PeriodicXOutflowY2d;
  }
  else if (strcmp(key, "reflect2d") == 0) {
    int revx = 2;//luaL_checkinteger(L, 2);
    int revy = 3;//luaL_checkinteger(L, 3);
    new_f = new ReflectingBoundary2d(revx, revy);
  }

  if (new_f) {
    if (Mara->boundary) delete Mara->boundary;
    Mara->boundary = new_f;
  }
  else {
    luaL_error(L, "[Mara] no such boundary conditions '%s'", key);
  }

  return 0;
}

int luaC_set_riemann(lua_State *L)
{
  const char *key = luaL_checkstring(L, 1);
  RiemannSolver *new_f = NULL;

  if (strcmp(key, "hll") == 0) {
    new_f = new HllRiemannSolver;
  }
  else if (strcmp(key, "hllc") == 0) {
    new_f = new HllcRiemannSolver;
  }
  else if (strcmp(key, "hlld") == 0) {
    new_f = new HlldRmhdRiemannSolver;
  }

  if (new_f) {
    if (Mara->riemann) delete Mara->riemann;
    Mara->riemann = new_f;
  }
  else {
    luaL_error(L, "[Mara] no such riemann solver '%s'", key);
  }

  return 0;
}

int luaC_set_godunov(lua_State *L)
{
  const char *key = luaL_checkstring(L, 1);
  GodunovOperator *new_f = NULL;

  if (strcmp(key, "plm-split") == 0) {
    new_f = new MethodOfLinesSplit;
  }
  else if (strcmp(key, "plm-muscl") == 0) {
    new_f = new PlmCtuHancockOperator;
  }
  else if (strcmp(key, "weno-split") == 0) {
    new_f = new WenoSplit;
  }

  if (new_f) {
    if (Mara->godunov) delete Mara->godunov;
    Mara->godunov = new_f;
  }
  else {
    luaL_error(L, "[Mara] no such integration scheme '%s'", key);
  }

  return 0;
}

int luaC_set_plm_theta(lua_State *L)
{
  if (Mara->godunov == NULL) {
    luaL_error(L, "[Mara] need a godunov operator for this");
  }
  double theta = luaL_checknumber(L, 1);
  Mara->godunov->SetPlmTheta(theta);
  return 0;
}

int luaC_config_solver(lua_State *L)
// -----------------------------------------------------------------------------
// Configures the GodunovOperator::reconstruct_method flag as well as internal
// parameters specific to the reconstruction library. The input is a table which
// may contain any of the following keys:
//
// fsplit (string) : one of [llf, marq]       ... flux splitting mode
// extrap (string) : one of [pcm, plm, weno5] ... reconstruction type
// theta  (number) : must be [0,2]            ... theta value for PLM/minmod
// IS     (string) : one of [js96, b08, sz10] ... smoothness indicator
// sz10A  (number) : should be in [0,100]     ... used by sz10 (see weno.c)
//
// A second positional argument, quiet (bool) may be provided.
// -----------------------------------------------------------------------------
{
  typedef std::map<std::string, GodunovOperator::FluxSplittingMethod> FSmap;
  typedef std::map<std::string, GodunovOperator::ReconstructMethod> RMmap;
  typedef std::map<std::string, SmoothnessIndicator> ISmap;
  luaL_checktype(L, 1, LUA_TTABLE);

  int quiet = 0;
  if (lua_gettop(L) == 2) {
    quiet = lua_toboolean(L, 2);
  }

  FSmap FSmodes;
  FSmodes["llf"] = GodunovOperator::FLUXSPLIT_LOCAL_LAX_FRIEDRICHS;
  FSmodes["marq"] = GodunovOperator::FLUXSPLIT_MARQUINA;

  RMmap RMmodes;
  RMmodes["pcm"] = GodunovOperator::RECONSTRUCT_PCM;
  RMmodes["plm"] = GodunovOperator::RECONSTRUCT_PLM;
  RMmodes["weno5"] = GodunovOperator::RECONSTRUCT_WENO5;

  ISmap ISmodes;
  ISmodes["js96"] = OriginalJiangShu96;
  ISmodes["b08"] = ImprovedBorges08;
  ISmodes["sz10"] = ImprovedShenZha10;

  lua_getfield(L, 1, "fsplit");
  if (lua_isstring(L, -1)) {
    const char *key = lua_tostring(L, -1);
    FSmap::iterator it = FSmodes.find(key);
    if (it != FSmodes.end()) {
      if (!quiet) printf("[config] setting fsplit=%s\n", it->first.c_str());
      GodunovOperator::fluxsplit_method = it->second;
    }
    else {
      luaL_error(L, "[Mara] no such fsplit '%s'", key);
    }
  }
  lua_pop(L, 1);

  lua_getfield(L, 1, "extrap");
  if (lua_isstring(L, -1)) {
    const char *key = lua_tostring(L, -1);
    RMmap::iterator it = RMmodes.find(key);
    if (it != RMmodes.end()) {
      if (!quiet) printf("[config] setting extrap=%s\n", it->first.c_str());
      GodunovOperator::reconstruct_method = it->second;
    }
    else {
      luaL_error(L, "[Mara] no such extrap '%s'", key);
    }
  }
  lua_pop(L, 1);

  lua_getfield(L, 1, "theta");
  if (lua_isnumber(L, -1)) {
    const double theta = lua_tonumber(L, -1);
    if (!quiet) printf("[Mara] setting theta=%f\n", theta);
    reconstruct_set_plm_theta(theta);
    if (Mara->godunov) {
      Mara->godunov->SetPlmTheta(theta);
    }
  }
  lua_pop(L, 1);

  lua_getfield(L, 1, "IS");
  if (lua_isstring(L, -1)) {
    const char *key = lua_tostring(L, -1);
    ISmap::iterator it = ISmodes.find(key);
    if (it != ISmodes.end()) {
      if (!quiet) printf("[Mara] setting IS=%s\n", it->first.c_str());
      reconstruct_set_smoothness_indicator(it->second);
    }
    else {
      luaL_error(L, "no such IS: %s", key);
    }
  }
  lua_pop(L, 1);

  lua_getfield(L, 1, "sz10A");
  if (lua_isnumber(L, -1)) {
    const double A = lua_tonumber(L, -1);
    if (!quiet) printf("[Mara] setting sz10A=%f\n", A);
    reconstruct_set_shenzha10_A(A);
  }
  lua_pop(L, 1);

  return 0;
}

int luaC_set_advance(lua_State *L)
{
  const char *key = luaL_checkstring(L, 1);
  RungeKuttaIntegration *new_f = NULL;

  if (strcmp(key, "single") == 0) {
    new_f = new RungeKuttaSingleStep;
  }
  else if (strcmp(key, "rk2") == 0) {
    new_f = new RungeKuttaRk2Tvd;
  }
  else if (strcmp(key, "rk3") == 0) {
    new_f = new RungeKuttaShuOsherRk3;
  }
  else if (strcmp(key, "rk4") == 0) {
    new_f = new RungeKuttaClassicRk4;
  }

  if (new_f) {
    if (Mara->advance) delete Mara->advance;
    Mara->advance = new_f;
  }
  else {
    luaL_error(L, "[Mara] no such advance option '%s'", key);
  }

  return 0;
}

int luaC_set_driving(lua_State *L)
{
  if (Mara->domain == NULL) {
    luaL_error(L, "[Mara] need a domain to run this, use set_domain");
  }

  size_t len;
  DrivingModule *new_f = NULL;
  const char *buf = luaL_checklstring(L, 1, &len);

  std::stringstream stream;
  stream.write(buf, sizeof(char)*len);

  if (Mara->domain->get_Nd() == 2) {
    StochasticVectorField *field = new StochasticVectorField2d(stream);
    new_f = new DrivingProcedure(field);
  }
  else if (Mara->domain->get_Nd() == 3) {
    StochasticVectorField *field = new StochasticVectorField3d(stream);
    new_f = new DrivingProcedure(field);
  }

  if (new_f) {
    if (Mara->driving) delete Mara->driving;
    Mara->driving = new_f;
  }
  return 0;
}

int luaC_set_cooling(lua_State *L)
{
  const char *key   = luaL_checkstring(L, 1);
  const double eref = luaL_checknumber(L, 2);
  const double t0   = luaL_checknumber(L, 3);

  CoolingModule *new_f = NULL;

  if (strcmp("T4", key) == 0) {
    new_f = new CoolingModuleT4(eref, t0);
  }
  else if (strcmp("E4", key) == 0) {
    new_f = new CoolingModuleE4(eref, t0);
  }

  if (new_f) {
    if (Mara->cooling) delete Mara->cooling;
    Mara->cooling = new_f;
  }
  return 0;
}

int luaC_set_primitive(lua_State *L)
{
  double *UserPrim = (double*) lua_touserdata(L, 1);
  unsigned N = lua_rawlen(L, 1) / sizeof(double); // total number of zones
  unsigned N_needed = Mara->domain->GetNumberOfZones() * Mara->domain->get_Nq();
  if (N != N_needed) {
    luaL_error(L, "[Mara] primitive has the wrong size for the domain");
  }
  Mara->PrimitiveArray.resize(N);
  memcpy(&Mara->PrimitiveArray[0], UserPrim, N * sizeof(double));
  return 0;
}

int luaC_cooling_rate(lua_State *L)
{
  const double dt = luaL_checknumber(L, 1);
  double res;
  if (Mara->cooling == NULL) {
    res = 0.0;
  }
  else {
    res = Mara->cooling->EnergyRemoved() / dt;
  }
  lua_pushnumber(L, res);
  return 1;
}

int luaC_load_shen(lua_State *L)
{
  const int narg = lua_gettop(L);
  const char *fname = luaL_checkstring(L, 1);
  double  YpExtract = luaL_checknumber(L, 2);
  double *DensRange = narg == 2 ? NULL : luaU_checkarray(L, 3);
  double *TempRange = narg == 2 ? NULL : luaU_checkarray(L, 4);

  ShenTabulatedNuclearEos::verbose = 1;
  TabulatedEos tab =
    ShenTabulatedNuclearEos::LoadTable(fname, YpExtract, DensRange, TempRange);

  free(DensRange);
  free(TempRange);

  lua_newtable(L);

  lua_pushstring(L, "logD_values");
  luaU_pusharray(L, &tab.logD_values[0], tab.logD_values.size());
  lua_settable(L, -3);
  lua_pushvalue(L, -1);

  lua_pushstring(L, "logT_values");
  luaU_pusharray(L, &tab.logT_values[0], tab.logT_values.size());
  lua_settable(L, -3);
  lua_pushvalue(L, -1);

  lua_pushstring(L, "EOS_p");
  luaU_pusharray(L, &tab.EOS_p[0], tab.EOS_p.size());
  lua_settable(L, -3);
  lua_pushvalue(L, -1);

  lua_pushstring(L, "EOS_s");
  luaU_pusharray(L, &tab.EOS_s[0], tab.EOS_s.size());
  lua_settable(L, -3);
  lua_pushvalue(L, -1);

  lua_pushstring(L, "EOS_u");
  luaU_pusharray(L, &tab.EOS_u[0], tab.EOS_u.size());
  lua_settable(L, -3);
  lua_pushvalue(L, -1);

  return 1;
}

int luaC_test_shen(lua_State *L)
{
  if (Mara->eos == NULL) {
    luaL_error(L, "set shen as the eos before testing it");
  }
  if (typeid(*Mara->eos) != typeid(ShenTabulatedNuclearEos)) {
    luaL_error(L, "set shen as the eos before testing it");
  }

  Mara->GetEos<ShenTabulatedNuclearEos>().self_test_derivatives();
  Mara->GetEos<ShenTabulatedNuclearEos>().self_test_inversion();
  Mara->GetEos<ShenTabulatedNuclearEos>().self_test_interpolation();
  Mara->GetEos<ShenTabulatedNuclearEos>().self_test_soundspeed();

  return 0;
}

static double StoneProfile(double k)
{
  const double KvecPeak = 2.0;
  return pow(k,6) * exp(-8*k/(2*M_PI*KvecPeak));
}
int luaC_new_ou_field(lua_State *L)
{
  const int Nd = luaL_checkinteger(L, 1);
  const double P0 = luaL_checknumber(L, 2);
  const double zeta = luaL_checknumber(L, 3);
  const int k1 = luaL_checkinteger(L, 4);
  const int seed = luaL_checkinteger(L, 5);

  std::stringstream stream;
  if (Nd == 2) {
    printf("building new 2d driving field\n");
    StochasticVectorField2d field(P0, zeta, k1, seed, StoneProfile);
    field.Serialize(stream);
  }
  else if (Nd == 3) {
    printf("building new 3d driving field\n");
    StochasticVectorField3d field(P0, zeta, k1, seed, StoneProfile);
    field.Serialize(stream);
  }
  std::string d = stream.str();

  lua_pushlstring(L, d.c_str(), d.length());
  return 1;
}

int luaC_fluid_PrimToCons(lua_State *L)
{
  if (Mara->fluid == NULL) {
    luaL_error(L, "[Mara] need a fluid to run this, use set_fluid");
  }
  else {
    int Nq = Mara->fluid->GetNq();
    double *P = luaU_checkarray(L, 1);
    double *U = (double*) malloc(Nq*sizeof(double));
    Mara->fluid->PrimToCons(P, U);
    luaU_pusharray(L, U, Nq);
    free(P);
    free(U);
  }
  return 1;
}
int luaC_fluid_ConsToPrim(lua_State *L)
{
  if (Mara->fluid == NULL) {
    luaL_error(L, "[Mara] need a fluid to run this, use set_fluid");
  }
  else {
    int Nq = Mara->fluid->GetNq();
    double *U = luaU_checkarray(L, 1);
    double *P = (double*) malloc(Nq*sizeof(double));
    for (int q=0; q<Nq; ++q) {
      P[q] = 0.0; // send a bad, but deteriministic guess state
    }
    int err = Mara->fluid->ConsToPrim(U, P);
    luaU_pusharray(L, P, Nq);
    free(U);
    free(P);
    if (err) {
      luaL_error(L, "[Mara] ConsToPrim failed");
    }
  }
  return 1;
}

int luaC_fluid_Eigensystem(lua_State *L)
{
  if (Mara->fluid == NULL) {
    luaL_error(L, "[Mara] need a fluid to run this, use set_fluid");
  }

  double *P = luaU_checkarray(L, 1);
  int dim = luaL_checkinteger(L, 2);

  int Nq = Mara->fluid->GetNq();
  double *U  = (double*) malloc(   Nq*sizeof(double));
  double *Lv = (double*) malloc(Nq*Nq*sizeof(double));
  double *Rv = (double*) malloc(Nq*Nq*sizeof(double));
  double *lm = (double*) malloc(   Nq*sizeof(double));

  Mara->fluid->PrimToCons(P, U);
  Mara->fluid->Eigensystem(U, P, Lv, Rv, lm, dim);

  lua_newtable(L);
  for (int i=0; i<Nq; ++i) {
    lua_pushnumber(L, i+1);
    luaU_pusharray(L, Lv+i*Nq, Nq);
    lua_settable(L, -3);
  }

  lua_newtable(L);
  for (int i=0; i<Nq; ++i) {
    lua_pushnumber(L, i+1);
    luaU_pusharray(L, Rv+i*Nq, Nq);
    lua_settable(L, -3);
  }
  luaU_pusharray(L, lm, Nq);

  free(P);
  free(U);
  free(Lv);
  free(Rv);
  free(lm);

  return 3;
}

int luaC_fluid_FluxFunction(lua_State *L)
{
  if (Mara->fluid == NULL) {
    luaL_error(L, "[Mara] need a fluid to run this, use set_fluid");
  }

  double *P = luaU_checkarray(L, 1);
  int dim = luaL_checkinteger(L, 2);

  int Nq = Mara->fluid->GetNq();
  double *U = (double*) malloc(Nq*sizeof(double));
  double *F = (double*) malloc(Nq*sizeof(double));

  double ap, am;
  Mara->fluid->PrimToCons(P, U);
  Mara->fluid->FluxAndEigenvalues(U, P, F, &ap, &am, dim);
  luaU_pusharray(L, F, Nq);
  lua_pushnumber(L, ap);
  lua_pushnumber(L, am);

  free(P);
  free(U);
  free(F);
  
  return 3;
}

int luaC_fluid_GetPrimNames(lua_State *L)
{
  if (Mara->fluid == NULL) {
    luaL_error(L, "[Mara] need a fluid to run this, use set_fluid");
  }
  std::vector<std::string> pnames = Mara->fluid->GetPrimNames();
  lua_newtable(L);
  for (unsigned int n=0; n<pnames.size(); ++n) {
    lua_pushstring(L, pnames[n].c_str());
    lua_rawseti(L, -2, n+1);
  }
  return 1;
}
/*
int luaC_boundary_ApplyBoundaries(lua_State *L)
{
  if (Mara->boundary == NULL) {
    luaL_error(L, "[Mara] need a boundary to run this, use set_boundary");
  }
  else {
    Mara->boundary->ApplyBoundaries(Mara->PrimitiveArray);
  }
  return 0;
}
*/
int luaC_eos_TemperatureMeV(lua_State *L)
{
  const double D = luaL_checknumber(L, 1);
  const double p = luaL_checknumber(L, 2);

  if (Mara->eos == NULL) {
    luaL_error(L, "[Mara] need an eos to run this, use set_eos");
  }
  else {
    const double T = Mara->eos->TemperatureMeV(D, p);
    lua_pushnumber(L, T);
  }
  return 1;
}

int luaC_eos_TemperatureArb(lua_State *L)
{
  const double D = luaL_checknumber(L, 1);
  const double T_MeV = luaL_checknumber(L, 2);

  if (Mara->eos == NULL) {
    luaL_error(L, "[Mara] need an eos to run this, use set_eos");
  }
  else {
    const double T = Mara->eos->TemperatureArb(D, T_MeV);
    lua_pushnumber(L, T);
  }
return 1;
}

int luaC_eos_Temperature_p(lua_State *L)
{
  const double D = luaL_checknumber(L, 1);
  const double p = luaL_checknumber(L, 2);

  if (Mara->eos == NULL) {
    luaL_error(L, "[Mara] need an eos to run this, use set_eos");
    return 0;
  }
  else {
    const double T = Mara->eos->Temperature_p(D, p);
    lua_pushnumber(L, T);
    return 1;
  }
}
int luaC_eos_Internal(lua_State *L)
{
  const double D = luaL_checknumber(L, 1);
  const double T = luaL_checknumber(L, 2);

  if (Mara->eos == NULL) {
    luaL_error(L, "[Mara] need an eos to run this, use set_eos");
    return 0;
  }
  else {
    const double u = Mara->eos->Internal(D, T);
    lua_pushnumber(L, u);
    return 1;
  }
}
int luaC_eos_Pressure(lua_State *L)
{
  const double D = luaL_checknumber(L, 1);
  const double T = luaL_checknumber(L, 2);

  if (Mara->eos == NULL) {
    luaL_error(L, "[Mara] need an eos to run this, use set_eos");
    return 0;
  }
  else {
    const double p = Mara->eos->Pressure(D, T);
    lua_pushnumber(L, p);
    return 1;
  }
}
int luaC_eos_SoundSpeed2Nr(lua_State *L)
{
  const double D = luaL_checknumber(L, 1);
  const double T = luaL_checknumber(L, 2);

  if (Mara->eos == NULL) {
    luaL_error(L, "[Mara] need an eos to run this, use set_eos");
    return 0;
  }
  else {
    const double p = Mara->eos->SoundSpeed2Nr(D, T);
    lua_pushnumber(L, p);
    return 1;
  }
}
int luaC_eos_SoundSpeed2Sr(lua_State *L)
{
  const double D = luaL_checknumber(L, 1);
  const double T = luaL_checknumber(L, 2);

  if (Mara->eos == NULL) {
    luaL_error(L, "[Mara] need an eos to run this, use set_eos");
    return 0;
  }
  else {
    const double p = Mara->eos->SoundSpeed2Sr(D, T);
    lua_pushnumber(L, p);
    return 1;
  }
}
int luaC_eos_DensUpper(lua_State *L)
{
  if (Mara->eos == NULL) {
    luaL_error(L, "[Mara] need an eos to run this, use set_eos");
    return 0;
  }
  else {
    lua_pushnumber(L, Mara->eos->DensUpper());
    return 1;
  }
}
int luaC_eos_DensLower(lua_State *L)
{
  if (Mara->eos == NULL) {
    luaL_error(L, "[Mara] need an eos to run this, use set_eos");
    return 0;
  }
  else {
    lua_pushnumber(L, Mara->eos->DensLower());
    return 1;
  }
}
int luaC_eos_TempUpper(lua_State *L)
{
  if (Mara->eos == NULL) {
    luaL_error(L, "[Mara] need an eos to run this, use set_eos");
    return 0;
  }
  else {
    lua_pushnumber(L, Mara->eos->TempUpper());
    return 1;
  }
}
int luaC_eos_TempLower(lua_State *L)
{
  if (Mara->eos == NULL) {
    luaL_error(L, "[Mara] need an eos to run this, use set_eos");
    return 0;
  }
  else {
    lua_pushnumber(L, Mara->eos->TempLower());
    return 1;
  }
}

int luaC_units_Print(lua_State *L)
{
  if (Mara->units == NULL) {
    luaL_error(L, "[Mara] need a units system to run this, use set_units");
    return 0;
  }
  else {
    Mara->units->PrintUnits();
    return 0;
  }
}

int luaC_units_Gauss(lua_State *L)
{
  if (Mara->units == NULL) {
    luaL_error(L, "[Mara] need a units system to run this, use set_units");
    return 0;
  }
  else {
    lua_pushnumber(L, Mara->units->Gauss());
    return 1;
  }
}

int luaC_units_Velocity(lua_State *L)
{
  if (Mara->units == NULL) {
    luaL_error(L, "[Mara] need a units system to run this, use set_units");
    return 0;
  }
  else {
    lua_pushnumber(L, Mara->units->Velocity());
    return 1;
  }
}

int luaC_units_MeV(lua_State *L)
{
  if (Mara->units == NULL) {
    luaL_error(L, "[Mara] need a units system to run this, use set_units");
    return 0;
  }
  else {
    lua_pushnumber(L, Mara->units->MeV());
    return 1;
  }
}

int luaC_units_MeVPerCubicFemtometer(lua_State *L)
{
  if (Mara->units == NULL) {
    luaL_error(L, "[Mara] need a units system to run this, use set_units");
    return 0;
  }
  else {
    lua_pushnumber(L, Mara->units->MeVPerCubicFemtometer());
    return 1;
  }
}

int luaC_units_ErgPerCubicCentimeter(lua_State *L)
{
  if (Mara->units == NULL) {
    luaL_error(L, "[Mara] need a units system to run this, use set_units");
    return 0;
  }
  else {
    lua_pushnumber(L, Mara->units->ErgPerCubicCentimeter());
    return 1;
  }
}

int luaC_units_GramsPerCubicCentimeter(lua_State *L)
{
  if (Mara->units == NULL) {
    luaL_error(L, "[Mara] need a units system to run this, use set_units");
    return 0;
  }
  else {
    lua_pushnumber(L, Mara->units->GramsPerCubicCentimeter());
    return 1;
  }
}


int luaC_driving_Advance(lua_State *L)
{
  const double dt = luaL_checknumber(L, 1);
  if (Mara->driving) {
    Mara->driving->GetField()->AdvanceField(dt);
  }
  return 0;
}
int luaC_driving_Resample(lua_State *L)
{
  if (Mara->driving) {
    Mara->driving->ResampleField();
  }
  return 0;
}
int luaC_driving_Serialize(lua_State *L)
{
  if (!Mara->driving) {
    return 0;
  }

  std::stringstream stream;
  Mara->driving->GetField()->Serialize(stream);
  std::string d = stream.str();
  lua_pushlstring(L, d.c_str(), d.length());
  return 1;
}

int luaC_driving_Power(lua_State *L)
{
  if (!Mara->driving) {
    lua_pushnumber(L, 0.0);
  }
  else {
    lua_pushnumber(L, Mara->driving->AveragePowerInFields());
  }
  return 1;
}

static double get_L2_error(const double *P, const double *Q, int Nq)
{
  double L2 = 0.0;
  for (int i=0; i<Nq; ++i) {
    L2 += (P[i]-Q[i])*(P[i]-Q[i]);
  }
  return sqrt(L2);
}

int luaC_test_rmhd_c2p(lua_State *L)
{
  double *P_ = luaU_checkarray(L, 1);

  double P[8], U[8], Q[8];
  memcpy(P, P_, 8*sizeof(double));
  free(P_);

  Mara->fluid->PrimToCons(P, U);

  rmhd_c2p_eos_set_eos(Mara->eos);
  rmhd_c2p_eos_new_state(U);
  rmhd_c2p_eos_estimate_from_cons();

  int codes[4], iters[4];
  double error[4], times[4];
  clock_t start;

  start = clock();
  codes[0] = rmhd_c2p_eos_solve_noble2dzt(Q);
  iters[0] = rmhd_c2p_eos_get_iterations();
  error[0] = get_L2_error(P, Q, 8);
  times[0] = (double) (clock() - start) / CLOCKS_PER_SEC;

  start = clock();
  codes[1] = rmhd_c2p_eos_solve_duffell3d(Q);
  iters[1] = rmhd_c2p_eos_get_iterations();
  error[1] = get_L2_error(P, Q, 8);
  times[1] = (double) (clock() - start) / CLOCKS_PER_SEC;

  if (typeid(*Mara->eos) != typeid(AdiabaticEos)) {
    luaU_pusharray_i(L, codes, 2);
    luaU_pusharray_i(L, iters, 2);
    luaU_pusharray  (L, error, 2);
    luaU_pusharray  (L, times, 2);
    return 4;
  }

  rmhd_c2p_set_gamma(Mara->GetEos<AdiabaticEos>().Gamma);
  rmhd_c2p_new_state(U);
  rmhd_c2p_estimate_from_cons();

  start = clock();
  codes[2] = rmhd_c2p_solve_noble1dw(Q);
  iters[2] = rmhd_c2p_get_iterations();
  error[2] = get_L2_error(P, Q, 8);
  times[2] = (double) (clock() - start) / CLOCKS_PER_SEC;

  start = clock();
  codes[3] = rmhd_c2p_solve_anton2dzw(Q);
  iters[3] = rmhd_c2p_get_iterations();
  error[3] = get_L2_error(P, Q, 8);
  times[3] = (double) (clock() - start) / CLOCKS_PER_SEC;

  luaU_pusharray_i(L, codes, 4);
  luaU_pusharray_i(L, iters, 4);
  luaU_pusharray  (L, error, 4);
  luaU_pusharray  (L, times, 4);

  return 4;
}


void luaU_pusharray(lua_State *L, double *A, int N)
{
  lua_newtable(L);
  for (int i=0; i<N; ++i) {
    lua_pushnumber(L, A[i]);
    lua_rawseti(L, -2, i+1);
  }
}
void luaU_pusharray_i(lua_State *L, int *A, int N)
{
  lua_newtable(L);
  for (int i=0; i<N; ++i) {
    lua_pushnumber(L, A[i]);
    lua_rawseti(L, -2, i+1);
  }
}
double *luaU_checkarray(lua_State *L, int n)
{
  n = lua_absindex(L, n);
  if (lua_type(L, n) != LUA_TTABLE) {
    luaL_error(L, "[Mara] expected a table as argument %d", n);
  }
  int N = lua_rawlen(L, n);
  double *ret = (double*) malloc(N * sizeof(double));
  for (int i=0; i<N; ++i) {
    lua_rawgeti(L, n, i+1);
    ret[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
  return ret;
}
double *luaU_checklarray(lua_State *L, int n, int *N)
{
  n = lua_absindex(L, n);
  if (lua_type(L, n) != LUA_TTABLE) {
    luaL_error(L, "[Mara] expected a table as argument %d", n);
  }
  *N = lua_rawlen(L, n);
  double *ret = (double*) malloc(*N * sizeof(double));
  for (int i=0; i<*N; ++i) {
    lua_rawgeti(L, n, i+1);
    ret[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
  return ret;
}
int *luaU_checklarray_i(lua_State *L, int n, int *N)
{
  n = lua_absindex(L, n);
  if (lua_type(L, n) != LUA_TTABLE) {
    luaL_error(L, "[Mara] expected a table as argument %d", n);
  }
  *N = lua_rawlen(L, n);
  int *ret = (int*) malloc(*N * sizeof(int));
  for (int i=0; i<*N; ++i) {
    lua_rawgeti(L, n, i+1);
    ret[i] = lua_tointeger(L, -1);
    lua_pop(L, 1);
  }
  return ret;
}






enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive
#define TTL_ZONES (HydroModule::Mara->domain->GetGlobalNumberOfZones())


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





#ifdef __GNUC__
#include <cxxabi.h>
std::string Demangle(const std::string &mname)
{
  int     status;
  char   *realname;

  realname = abi::__cxa_demangle(mname.c_str(), 0, 0, &status);
  std::string retname(realname);
  free(realname);

  return retname;
}
#else
std::string Demangle(const std::string &mname)
{
  return mname;
}
#endif // __GNUC__
