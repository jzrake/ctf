
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
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#include "config.h"
#if (__MARA_USE_MPI)
#include <mpi.h>
#endif

#include <iostream>
#include <map>
#include <readline/readline.h>
#include <readline/history.h>
extern "C" {
#include "lualib.h"
#include "lauxlib.h"
}
#include "luaU.h"
#include "mara.hpp"




extern "C"
// -----------------------------------------------------------------------------
// These functions constitute most of the Lua interface. Other functions are
// registered to the Lua environment though the fluid, eos, and units Lua tables
// in this file, while more are added lua_*_load functions declared in luaU.h,
// where * may be mpi, h5, fft, etc.
// -----------------------------------------------------------------------------
{
  static int luaC_mara_version(lua_State *L);
  static int luaC_print_mara(lua_State *L);

  static int luaC_advance(lua_State *L);
  static int luaC_diffuse(lua_State *L);
  static int luaC_init_prim(lua_State *L);
  static int luaC_read_prim(lua_State *L);
  static int luaC_write_prim(lua_State *L);
  static int luaC_get_prim(lua_State *L);
  static int luaC_prim_at_point(lua_State *L);
  static int luaC_get_timestep(lua_State *L);
  static int luaC_streamline(lua_State *L);

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

  static int luaC_new_ou_field(lua_State *L);
  static int luaC_load_shen(lua_State *L);
  static int luaC_test_shen(lua_State *L);
  static int luaC_test_rmhd_c2p(lua_State *L);
  static int luaC_test_sampling(lua_State *L);
  static int luaC_test_sampling_many(lua_State *L);

  static int luaC_fluid_PrimToCons(lua_State *L);
  static int luaC_fluid_ConsToPrim(lua_State *L);
  static int luaC_fluid_Eigensystem(lua_State *L);

  static int luaC_boundary_ApplyBoundaries(lua_State *L);

  static int luaC_driving_Advance(lua_State *L);
  static int luaC_driving_Resample(lua_State *L);
  static int luaC_driving_Serialize(lua_State *L);
  static int luaC_driving_Power(lua_State *L);

  static int luaC_eos_TemperatureMeV(lua_State *L);
  static int luaC_eos_Temperature_p(lua_State *L);
  static int luaC_eos_Internal(lua_State *L);
  static int luaC_eos_Pressure(lua_State *L);
  static int luaC_eos_DensUpper(lua_State *L);
  static int luaC_eos_DensLower(lua_State *L);
  static int luaC_eos_TempUpper(lua_State *L);
  static int luaC_eos_TempLower(lua_State *L);

  static int luaC_units_Print(lua_State *L);
  static int luaC_units_Gauss(lua_State *L);
  static int luaC_units_Velocity(lua_State *L);
  static int luaC_units_MeVPerCubicFemtometer(lua_State *L);
  static int luaC_units_ErgPerCubicCentimeter(lua_State *L);
  static int luaC_units_GramsPerCubicCentimeter(lua_State *L);

  static int luaC_write_ppm(lua_State *L);

  int luaopen_lunum(lua_State *L);
}

static void mara_prim_io(lua_State *L, char mode);
static MaraApplication *Mara;


// Variables whose value may be set by command line flags on startup
// -----------------------------------------------------------------------------
static int PrintVersion          = 0;
static int InteractiveMode       = 0;
static char CommandlineFluid[32] = "euler";
static char LuaProgramName[128]  = "";
// -----------------------------------------------------------------------------



int main(int argc, char **argv)
{
  // ---------------------------------------------------------------------------
  // Collect the entries from the command line
  // ---------------------------------------------------------------------------
  HydroModule::Mara = Mara = new MaraApplication;
  OptionParser parser;

  parser.Register(&InteractiveMode, Opt::Bool, 'i', "interactive",
                  "enter interactive mode after executing 'script'");
  parser.Register(&PrintVersion, Opt::Bool, 'v', "version",
                  "print version number and exit");
  parser.Register(&CommandlineFluid, Opt::String, 'f', "fluid",
                  "startup value for fluid equations");
  std::vector<std::string> args = parser.ParseOptions(argc, argv);


#if (__MARA_USE_MPI)
  // ---------------------------------------------------------------------------
  // Configuring MPI usage
  // ---------------------------------------------------------------------------
  if (InteractiveMode) {
    printf("-i: disabling MPI for interactive mode\n");
  }
  else {
    MPI_Init(&argc, &argv);
    if (Mara_mpi_get_rank() != 0) {
      freopen("/dev/null", "w", stdout);
    }
  }
#endif // ----------------------------------------------------------------------



  printf("\n");
  printf("Mara: Astrophysical gasdynamics code\n");
  printf("(C) Jonathan Zrake, NYU CCPP\n");
  printf("version %s.%s\n\n", __MARA_BASE_VERSION, __MARA_HG_CHANGESET);

  if (parser.HelpMessage) {
    puts("usage: mara [options] program.lua [arg1 [...]] [opt1=val1 [...]]");
    parser.PrintUsageMessage();
    return 0;
  }
  if (PrintVersion) {
    return 0;
  }


  // ---------------------------------------------------------------------------
  // Configure the Lua environment: open Lua stdlib, the Mara libraries, and
  // include the paths to Lua modules in the __MARA_INSTALL_DIR.
  // ---------------------------------------------------------------------------
  lua_State *L = luaL_newstate();
  luaL_openlibs(L);
  luaopen_lunum(L);

  lua_h5_load(L);
  lua_mpi_load(L);
  lua_measure_load(L);
  lua_fft_load(L);
  lua_vis_load(L);

  lua_getglobal(L, "package");
  lua_pushfstring(L, "./?.lua;%s/conf/?.lua;%s/lib/lua/5.2/?.lua",
                  __MARA_INSTALL_DIR,  __MARA_INSTALL_DIR);
  lua_setfield(L, 1, "path");
  lua_remove(L, 1);



  // ---------------------------------------------------------------------------
  // Collect command line arguments. If the first argument after
  // 
  // $> mara [options] 
  //
  // ends in the extension '.lua', then regeister it as the LuaProgramName. All
  // subsequent command line entries are dropped into the 'cmdline.opts' table
  // if they contains an equals sign, and the 'cmdline.args' table otherwise.
  // ---------------------------------------------------------------------------
  lua_newtable(L); // 'cmdline' stack index := 1

  lua_pushstring(L, "opts");
  lua_newtable(L); // 'cmdline.opts' stack index := 3
  for (size_t i=0; i<args.size(); ++i) {

    const char *a = args[i].c_str();
    const size_t p = strcspn(a, "=");

    if (p != strlen(a)) {
      lua_pushlstring(L, a, p);
      lua_pushstring(L, a+p+1);
      lua_settable(L, 3);
    }
  }
  lua_settable(L, 1);


  lua_pushstring(L, "args");
  lua_newtable(L); // 'cmdline.args' stack index := 3
  if (args.size() > 0) {
    if (args[0].compare(args[0].length()-4, 4, ".lua") == 0) {
      strcpy(LuaProgramName, args[0].c_str());
      args.erase(args.begin());
    }
  }
  for (size_t i=0; i<args.size(); ++i) {
      
    const char *a = args[i].c_str();
    const size_t p = strcspn(a, "=");

    if (p == strlen(a)) {
      lua_pushnumber(L, lua_rawlen(L, 3) + 1);
      lua_pushstring(L, a);
      lua_settable(L, 3);
    }
  }
  lua_settable(L, 1);

  lua_setglobal(L, "cmdline");
  // ---------------------------------------------------------------------------



  lua_register(L, "mara_version", luaC_mara_version);
  lua_register(L, "print_mara", luaC_print_mara);
  lua_register(L, "advance", luaC_advance);
  lua_register(L, "diffuse", luaC_diffuse);

  lua_register(L, "init_prim"    , luaC_init_prim);
  lua_register(L, "read_prim"    , luaC_read_prim);
  lua_register(L, "write_prim"   , luaC_write_prim);
  lua_register(L, "write_ppm"    , luaC_write_ppm);
  lua_register(L, "prim_at_point", luaC_prim_at_point);
  lua_register(L, "get_prim"     , luaC_get_prim);
  lua_register(L, "get_timestep" , luaC_get_timestep);
  lua_register(L, "streamline"   , luaC_streamline);

  lua_register(L, "set_domain"   , luaC_set_domain);
  lua_register(L, "set_boundary" , luaC_set_boundary);
  lua_register(L, "set_fluid"    , luaC_set_fluid);
  lua_register(L, "set_eos"      , luaC_set_eos);
  lua_register(L, "set_units"    , luaC_set_units);
  lua_register(L, "set_godunov"  , luaC_set_godunov);
  lua_register(L, "set_riemann"  , luaC_set_riemann);
  lua_register(L, "set_advance"  , luaC_set_advance);
  lua_register(L, "set_driving"  , luaC_set_driving);
  lua_register(L, "set_cooling"  , luaC_set_cooling);

  lua_register(L, "new_ou_field" , luaC_new_ou_field);
  lua_register(L, "load_shen"    , luaC_load_shen);
  lua_register(L, "test_shen"    , luaC_test_shen);
  lua_register(L, "test_rmhd_c2p", luaC_test_rmhd_c2p);
  lua_register(L, "test_sampling", luaC_test_sampling);
  lua_register(L, "test_sampling_many", luaC_test_sampling_many);


  // Expose the fluid interface
  // ---------------------------------------------------------------------------
  lua_newtable(L);

  lua_pushstring(L, "PrimToCons");
  lua_pushcfunction(L, luaC_fluid_PrimToCons);
  lua_settable(L, 1);

  lua_pushstring(L, "ConsToPrim");
  lua_pushcfunction(L, luaC_fluid_ConsToPrim);
  lua_settable(L, 1);

  lua_pushstring(L, "Eigensystem");
  lua_pushcfunction(L, luaC_fluid_Eigensystem);
  lua_settable(L, 1);

  lua_setglobal(L, "fluid");


  // Expose the boundary conditions interface
  // ---------------------------------------------------------------------------
  lua_newtable(L);

  lua_pushstring(L, "ApplyBoundaries");
  lua_pushcfunction(L, luaC_boundary_ApplyBoundaries);
  lua_settable(L, 1);

  lua_setglobal(L, "boundary");


  // Expose the eos interface
  // ---------------------------------------------------------------------------
  lua_newtable(L);

  lua_pushstring(L, "TemperatureMeV");
  lua_pushcfunction(L, luaC_eos_TemperatureMeV);
  lua_settable(L, 1);

  lua_pushstring(L, "Temperature_p");
  lua_pushcfunction(L, luaC_eos_Temperature_p);
  lua_settable(L, 1);

  lua_pushstring(L, "Internal");
  lua_pushcfunction(L, luaC_eos_Internal);
  lua_settable(L, 1);

  lua_pushstring(L, "Pressure");
  lua_pushcfunction(L, luaC_eos_Pressure);
  lua_settable(L, 1);

  lua_pushstring(L, "DensUpper");
  lua_pushcfunction(L, luaC_eos_DensUpper);
  lua_settable(L, 1);

  lua_pushstring(L, "DensLower");
  lua_pushcfunction(L, luaC_eos_DensLower);
  lua_settable(L, 1);

  lua_pushstring(L, "TempUpper");
  lua_pushcfunction(L, luaC_eos_TempUpper);
  lua_settable(L, 1);

  lua_pushstring(L, "TempLower");
  lua_pushcfunction(L, luaC_eos_TempLower);
  lua_settable(L, 1);

  lua_setglobal(L, "eos");


  // Expose the units interface
  // ---------------------------------------------------------------------------
  lua_newtable(L);

  lua_pushstring(L, "Print");
  lua_pushcfunction(L, luaC_units_Print);
  lua_settable(L, 1);

  lua_pushstring(L, "Gauss");
  lua_pushcfunction(L, luaC_units_Gauss);
  lua_settable(L, 1);

  lua_pushstring(L, "Velocity");
  lua_pushcfunction(L, luaC_units_Velocity);
  lua_settable(L, 1);

  lua_pushstring(L, "MeVPerCubicFemtometer");
  lua_pushcfunction(L, luaC_units_MeVPerCubicFemtometer);
  lua_settable(L, 1);

  lua_pushstring(L, "ErgPerCubicCentimeter");
  lua_pushcfunction(L, luaC_units_ErgPerCubicCentimeter);
  lua_settable(L, 1);

  lua_pushstring(L, "GramsPerCubicCentimeter");
  lua_pushcfunction(L, luaC_units_GramsPerCubicCentimeter);
  lua_settable(L, 1);

  lua_setglobal(L, "units");



  // Expose the driving interface
  // ---------------------------------------------------------------------------
  lua_newtable(L);

  lua_pushstring(L, "Advance");
  lua_pushcfunction(L, luaC_driving_Advance);
  lua_settable(L, 1);

  lua_pushstring(L, "Resample");
  lua_pushcfunction(L, luaC_driving_Resample);
  lua_settable(L, 1);

  lua_pushstring(L, "Serialize");
  lua_pushcfunction(L, luaC_driving_Serialize);
  lua_settable(L, 1);

  lua_pushstring(L, "Power");
  lua_pushcfunction(L, luaC_driving_Power);
  lua_settable(L, 1);

  lua_setglobal(L, "driving");


  // Set the default fluid equations
  // ---------------------------------------------------------------------------
  lua_getglobal(L, "set_fluid");
  lua_pushstring(L, CommandlineFluid);
  lua_call(L, 1, 0);

  int listen = 1;

  if (strlen(LuaProgramName) != 0) {
    if (luaL_dofile(L, LuaProgramName)) {
      printf("%s\n", lua_tostring(L, -1));
    }
    listen = InteractiveMode;
  }


  // Configure the readline library and process stdin
  // ---------------------------------------------------------------------------
  rl_bind_key('\t', rl_complete);
  char *input;
  int numhist = 1;

  while (listen) {

    char prompt[256];
    sprintf(prompt, "Mara [%d]: ", numhist);
    input = readline(prompt);

    if (!input) {
      printf("\n\n");
      break;
    }

    if (strlen(input) > 0) {

      if (luaL_dostring(L, input)) {
        fprintf(stderr, "%s\n", lua_tostring(L, -1));
      }

      ++numhist;
      add_history(input);
    }
    printf("\n");
  }


  // Clean up and close libraries
  // ---------------------------------------------------------------------------
  lua_close(L);
  delete Mara;

#if (__MARA_USE_MPI)
  if (Mara_mpi_active()) MPI_Finalize();
#endif

  return 0;
  // ---------------------------------------------------------------------------
}



int luaC_advance(lua_State *L)
{
  const clock_t start = clock();
  const double dt = luaL_checknumber(L, 1);

  std::valarray<double> P = Mara->PrimitiveArray;
  std::valarray<double> U(P.size());
  int errors;

  Mara->riemann->ResetMaxLambda();
  Mara->FailureMask.resize(Mara->domain->GetNumberOfZones());
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
    Mara->PrimitiveArray = P;
  }

  const double sec = (double) (clock() - start) / CLOCKS_PER_SEC;

  lua_pushnumber(L, 1e-3*Mara->domain->GetNumberOfZones()/sec);
  lua_pushnumber(L, errors);

  return 2;
}

int luaC_diffuse(lua_State *L)
{
  const double r = luaL_checknumber(L, 1);

  std::valarray<double> &P = Mara->PrimitiveArray;
  std::valarray<double> U(P.size());

  Mara->godunov->PrimToCons(P, U);
  U += Mara->godunov->LaxDiffusion(U, r);
  Mara->godunov->ConsToPrim(U, P);

  return 0;
}

int luaC_get_timestep(lua_State *L)
{
  const double CFL = luaL_checknumber(L, 1);
  const double dt = Mara_mpi_dbl_min(CFL*Mara->domain->get_min_dx() /
                                     RiemannSolver::GetMaxLambda());
  lua_pushnumber(L, dt);

  return 1;
}

int luaC_streamline(lua_State *L)
{
  const double *r0 = luaU_checkarray (L, 1);
  const double  s  = luaL_checknumber(L, 2);
  const double ds  = luaL_checknumber(L, 3);
  const char *type = luaL_checkstring(L, 4);


  if (strcmp(type, "velocity") == 0) {
    std::vector<double> strm =
      Mara_streamline_velocity(r0, s, ds, Mara_streamline_scalars_velocity);

    int shape[2] = { strm.size()/4, 4 };
    luaU_pusharray_wshape(L, &strm[0], shape, 2);
  }
  else if (strcmp(type, "magnetic") == 0) {
    std::vector<double> strm =
      Mara_streamline_magnetic(r0, s, ds, Mara_streamline_scalars_magnetic);

    int shape[2] = { strm.size()/4, 4 };
    luaU_pusharray_wshape(L, &strm[0], shape, 2);
  }
  else {
    luaL_error(L, "[mara] please choose either 'velocity' or 'magnetic'");
  }

  return 1;
}

int luaC_mara_version(lua_State *L)
{
  char str[256];
  sprintf(str, "%s.%s", __MARA_BASE_VERSION, __MARA_HG_CHANGESET);
  lua_pushstring(L, str);
  return 1;
}

std::string Demangle(const std::string &mname);
int luaC_print_mara(lua_State *L)
{
  std::cout << std::endl;
  std::cout << "units: "
            << (Mara->units ? Demangle(typeid(*Mara->units).name()) : "NULL")
            << std::endl;


  std::cout << "domain: "
            << (Mara->domain ? Demangle(typeid(*Mara->domain).name()) : "NULL")
            << std::endl;


  std::cout << "boundary: "
            << (Mara->boundary ? Demangle(typeid(*Mara->boundary).name()) : "NULL")
            << std::endl;


  std::cout << "fluid: "
            << (Mara->fluid ? Demangle(typeid(*Mara->fluid).name()) : "NULL")
            << std::endl;


  std::cout << "eos: "
            << (Mara->eos ? Demangle(typeid(*Mara->eos).name()) : "NULL")
            << std::endl;


  std::cout << "godunov: "
            << (Mara->godunov ? Demangle(typeid(*Mara->godunov).name()) : "NULL")
            << std::endl;


  std::cout << "riemann: "
            << (Mara->riemann ? Demangle(typeid(*Mara->riemann).name()) : "NULL")
            << std::endl;


  std::cout << "advance: "
            << (Mara->advance ? Demangle(typeid(*Mara->advance).name()) : "NULL")
            << std::endl;


  std::cout << "driving: "
            << (Mara->driving ? Demangle(typeid(*Mara->driving).name()) : "NULL")
            << std::endl;


  std::cout << "cooling: "
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

  lua_pushstring(L, "logT_values");
  lua_gettable(L, -2);
  double *logT_values = luaU_checklarray(L, -1, &N);
  lua_pop(L, 1);
  tab.logT_values.assign(logT_values, logT_values+N);

  lua_pushstring(L, "EOS_p");
  lua_gettable(L, -2);
  double *EOS_p = luaU_checklarray(L, -1, &N);
  lua_pop(L, 1);
  tab.EOS_p.assign(EOS_p, EOS_p+N);

  lua_pushstring(L, "EOS_s");
  lua_gettable(L, -2);
  double *EOS_s = luaU_checklarray(L, -1, &N);
  lua_pop(L, 1);
  tab.EOS_s.assign(EOS_s, EOS_s+N);

  lua_pushstring(L, "EOS_u");
  lua_gettable(L, -2);
  double *EOS_u = luaU_checklarray(L, -1, &N);
  lua_pop(L, 1);
  tab.EOS_u.assign(EOS_u, EOS_u+N);

  return new ShenTabulatedNuclearEos(tab);
}

GodunovOperator *BuildPlmCtuHancockOperator(lua_State *L)
{
  const double plm_theta = luaL_checknumber(L, 2);
  const int safety = luaL_checkinteger(L, 3);

  PlmCtuHancockOperator *new_f = new PlmCtuHancockOperator;
  new_f->SetPlmTheta(plm_theta);
  new_f->SetSafetyLevel(safety);

  return new_f;
}

GodunovOperator *BuildPlmMethodOfLinesSplit(lua_State *L)
{
  const double plm_theta = luaL_checknumber(L, 2);
  const int safety = luaL_checkinteger(L, 3);

  PlmMethodOfLinesSplit *new_f = new PlmMethodOfLinesSplit;
  new_f->SetPlmTheta(plm_theta);
  new_f->SetSafetyLevel(safety);

  return new_f;
}




int luaC_init_prim(lua_State *L)
{
  const PhysicalDomain *domain = Mara->domain;
  double **data = NULL;

  if (domain == NULL) {
    luaL_error(L, "need a domain to run this, use set_domain");
  }
  if (lua_type(L, 1) == LUA_TFUNCTION) {
    // Handled by if (data == NULL) cases below
  }
  else if (lua_type(L, 1) == LUA_TTABLE) {
    std::vector<std::string> pnames = Mara->fluid->GetPrimNames();
    Mara->PrimitiveArray.resize(domain->GetNumberOfZones() * domain->get_Nq());
    data = (double**) malloc(pnames.size()*sizeof(double*));

    for (size_t q=0; q<pnames.size(); ++q) {
      lua_pushstring(L, pnames[q].c_str());
      lua_gettable(L, 1);
      data[q] = luaU_checkarray(L, 2);
      lua_pop(L, 1);
    }
  }
  else {
    luaL_error(L, "argument must be function or table");
  }

  Mara->PrimitiveArray.resize(domain->GetNumberOfZones() * domain->get_Nq());

  const int Ng = domain->get_Ng();
  const std::vector<int> Ninter(domain->GetLocalShape(),
                                domain->GetLocalShape()+domain->get_Nd());
  ValarrayIndexer N(Ninter);
  ValarrayManager M(domain->aug_shape(), domain->get_Nq());

  switch (domain->get_Nd()) {
  case 1:

    for (int i=0; i<Ninter[0]; ++i) {

      const double x = domain->x_at(i+Ng);
      const double y = 0.0;
      const double z = 0.0;

      if (data == NULL) {
        lua_pushvalue(L, 1);
        lua_pushnumber(L, x);
        lua_pushnumber(L, y);
        lua_pushnumber(L, z);
        lua_call(L, 3, 1);

        double *P0 = luaU_checkarray(L, 2);
        std::valarray<double> P(P0, domain->get_Nq());

        Mara->PrimitiveArray[ M(i+Ng) ] = P;
        lua_pop(L, 1);
      }
      else {
        std::valarray<double> P(domain->get_Nq());
        for (size_t q=0; q<P.size(); ++q) {
          P[q] = data[q][N(i)];
        }
        Mara->PrimitiveArray[ M(i+Ng) ] = P;
      }
    }
    break;

  case 2:
    for (int i=0; i<Ninter[0]; ++i) {
      for (int j=0; j<Ninter[1]; ++j) {

        const double x = domain->x_at(i+Ng);
        const double y = domain->y_at(j+Ng);
        const double z = 0.0;

        if (data == NULL) {
          lua_pushvalue(L, 1);
          lua_pushnumber(L, x);
          lua_pushnumber(L, y);
          lua_pushnumber(L, z);
          lua_call(L, 3, 1);

          double *P0 = luaU_checkarray(L, 2);
          std::valarray<double> P(P0, domain->get_Nq());

          Mara->PrimitiveArray[ M(i+Ng,j+Ng) ] = P;
          lua_pop(L, 1);
        }
        else {
          std::valarray<double> P(domain->get_Nq());
          for (size_t q=0; q<P.size(); ++q) {
            P[q] = data[q][N(i,j)];
          }
          Mara->PrimitiveArray[ M(i+Ng,j+Ng) ] = P;
        }
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

          if (data == NULL) {
            lua_pushvalue(L, 1);
            lua_pushnumber(L, x);
            lua_pushnumber(L, y);
            lua_pushnumber(L, z);
            lua_call(L, 3, 1);

            double *P0 = luaU_checkarray(L, 2);
            std::valarray<double> P(P0, domain->get_Nq());

            Mara->PrimitiveArray[ M(i+Ng,j+Ng,k+Ng) ] = P;
            lua_pop(L, 1);
          }
          else {
            std::valarray<double> P(domain->get_Nq());
            for (size_t q=0; q<P.size(); ++q) {
              P[q] = data[q][N(i,j,k)];
            }
            Mara->PrimitiveArray[ M(i+Ng,j+Ng,k+Ng) ] = P;
          }
        }
      }
    }
    break;
  }

  if (data != NULL) free(data);
  return 0;
}

int luaC_prim_at_point(lua_State *L)
{
  const double *r1 = luaU_checkarray(L, 1);
  const int Nq = Mara->domain->get_Nq();
  double *P1 = new double[Nq];
  Mara_prim_at_point(r1, P1);
  luaU_pusharray(L, P1, Nq);
  delete [] P1;
  return 1;
}

int luaC_test_sampling(lua_State *L)
{
  if (Mara->domain == NULL) {
    luaL_error(L, "[mara] error: need a domain to run this, use set_domain\n");
  }
  if (Mara->domain->get_Nd() != 3) {
    luaL_error(L, "[mara] error: need a 3d domain to run this\n");
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
    //    printf("(%f %f %f) ", r1[0], r1[1], r1[2]);
    //    std::cout << Mara->fluid->PrintPrim(P1) << std::endl;
  }

  const double trun = (double) (clock() - start) / CLOCKS_PER_SEC;
  lua_pushnumber(L, trun);

  delete [] P1;
  return 1;
}

int luaC_test_sampling_many(lua_State *L)
{
  if (Mara->domain == NULL) {
    luaL_error(L, "[mara] error: need a domain to run this, use set_domain\n");
  }
  if (Mara->domain->get_Nd() != 3) {
    luaL_error(L, "[mara] error: need a 3d domain to run this\n");
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
    //    printf("(%f %f %f) ", Rlist[3*i+0], Rlist[3*i+1], Rlist[3*i+2]);
    //    std::cout << Mara->fluid->PrintPrim(&Plist[Nq*i]) << std::endl;
  }

  const double trun = (double) (clock() - start) / CLOCKS_PER_SEC;
  lua_pushnumber(L, trun);

  delete [] Rinpt;
  delete [] Rlist;
  delete [] Plist;

  return 1;
}


int luaC_get_prim(lua_State *L)
{
  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  std::vector<std::string> pnames = Mara->fluid->GetPrimNames();

  lua_newtable(L);

  const int Nd = Mara->domain->get_Nd();
  const int Ng = Mara->domain->get_Ng();
  const int Nx = Mara->domain->get_N(1);
  const int Ny = Nd >= 2 ? Mara->domain->get_N(2) : 1;
  const int Nz = Nd >= 3 ? Mara->domain->get_N(3) : 1;

  for (size_t n=0; n<pnames.size(); ++n) {

    std::vector<double> var;

    switch (Nd) {
    case 1:
      for (int i=Ng; i<Nx+Ng; ++i) {
        std::valarray<double> P0 = Mara->PrimitiveArray[ M(i) ];
        var.push_back(P0[n]);
      }
      break;
    case 2:
      for (int i=Ng; i<Nx+Ng; ++i) {
        for (int j=Ng; j<Ny+Ng; ++j) {
          std::valarray<double> P0 = Mara->PrimitiveArray[ M(i,j) ];
          var.push_back(P0[n]);
        }
      }
      break;
    case 3:
      for (int i=Ng; i<Nx+Ng; ++i) {
        for (int j=Ng; j<Ny+Ng; ++j) {
          for (int k=Ng; k<Nz+Ng; ++k) {
            std::valarray<double> P0 = Mara->PrimitiveArray[ M(i,j,k) ];
            var.push_back(P0[n]);
          }
        }
      }
      break;
    }
    lua_pushstring(L, pnames[n].c_str());
    luaU_pusharray_wshape(L, &var[0], Mara->domain->GetLocalShape(), Nd);
    lua_settable(L, 1);
  }

  return 1;
}




int luaC_read_prim(lua_State *L)
{
  const PhysicalDomain *domain = Mara->domain;

  if (domain == NULL) {
    luaL_error(L, "need a domain to run this, use set_domain");
  }

  Mara->PrimitiveArray.resize(domain->GetNumberOfZones() * domain->get_Nq());
  clock_t start = clock();
  mara_prim_io(L, 'r');
  lua_pushnumber(L, (double) (clock() - start) / CLOCKS_PER_SEC);
  return 1;
}
int luaC_write_prim(lua_State *L) {
  clock_t start = clock();
  mara_prim_io(L, 'w');
  lua_pushnumber(L, (double) (clock() - start) / CLOCKS_PER_SEC);
  return 1;
}
void mara_prim_io(lua_State *L, char mode)
{
  const int narg = lua_gettop(L);
  const char *fname = luaL_checkstring(L, 1);


  // Input/output function modes are given string aliases to be used from lua.
  // ---------------------------------------------------------------------------
  std::map<std::string, enum MaraIoFunction> io_modes;
  io_modes["BINARY"] = MARA_IO_FUNC_BINARY;
  io_modes["H5MPI" ] = MARA_IO_FUNC_H5MPI;
  io_modes["H5SER" ] = MARA_IO_FUNC_H5SER;


  // We will load these from the input lua table if it was provided (narg == 2),
  // otherwise these are the default values.
  // ---------------------------------------------------------------------------
  MaraIoFunction output_function = MARA_IO_FUNC_BINARY;
  MaraIoFunction  input_function = MARA_IO_FUNC_BINARY;
  int disk_align_threshold = 0;
  int stripe_size_mb = 0;
  int enable_chunking = 0;
  int enable_alignment = 0;


  // If a table with additional options was provided as input, execute this.
  // ---------------------------------------------------------------------------
  if (narg == 2) {

    lua_pushstring(L, "input_function");
    lua_gettable(L, 2);
    input_function = io_modes[lua_tostring(L, -1)];
    lua_pop(L, 1);

    lua_pushstring(L, "output_function");
    lua_gettable(L, 2);
    output_function = io_modes[lua_tostring(L, -1)];
    lua_pop(L, 1);

    lua_pushstring(L, "disk_align_threshold");
    lua_gettable(L, 2);
    disk_align_threshold = lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_pushstring(L, "stripe_size_mb");
    lua_gettable(L, 2);
    stripe_size_mb = lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_pushstring(L, "enable_chunking");
    lua_gettable(L, 2);
    enable_chunking = lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_pushstring(L, "enable_alignment");
    lua_gettable(L, 2);
    enable_alignment = lua_tointeger(L, -1);
    lua_pop(L, 1);
  }


  // Gather the necessary information from Mara application to pass to Mara_io.
  // ---------------------------------------------------------------------------
  const PhysicalDomain::SubdomainSpecs d = Mara->domain->GetSpecs();
  const std::vector<std::string> PNames = Mara->fluid->GetPrimNames();
  const char **pnames = (const char**) malloc(PNames.size()*sizeof(char*));

  for (size_t i=0; i<PNames.size(); ++i) {
    pnames[i] = PNames[i].c_str();
  }


  // Open, run, and close the io library.
  // ---------------------------------------------------------------------------
  Mara_io_init(0, d.n_dims, d.n_prim, d.A_nint, d.L_ntot, d.L_strt, d.G_ntot, d.G_strt);
  Mara_io_set_logfile(stdout);
  Mara_io_set_output_function(output_function);
  Mara_io_set_input_function(input_function);
  Mara_io_set_disk_block_size(1024*1024*stripe_size_mb);
  Mara_io_set_disk_align_threshold(disk_align_threshold);
  Mara_io_set_chunk_size(d.A_nint);
  Mara_io_set_enable_chunking(enable_chunking);
  Mara_io_set_enable_alignment(enable_alignment);

  if (mode == 'r') {
    Mara_io_read_prim(fname, pnames, &Mara->PrimitiveArray[0]);
  }
  else if (mode == 'w') {
    Mara_io_write_prim(fname, pnames, &Mara->PrimitiveArray[0]);
  }
  Mara_io_free();

  free(pnames);
}



int luaC_set_domain(lua_State *L)
{
  int Nd;
  double *x0 = luaU_checkarray(L, 1);
  double *x1 = luaU_checkarray(L, 2);

  int *N = luaU_checklarray_i(L, 3, &Nd);
  int Nq = luaL_checknumber(L, 4);
  int Ng = luaL_checknumber(L, 5);

  PhysicalDomain *new_f = NULL;

  if (Mara_mpi_active()) {
    new_f = new DecomposedCartesianDomain(x0, x1, N, Nd, Nq, Ng);
  }
  else {
    new_f = new SimpleCartesianDomain(x0, x1, N, Nd, Nq, Ng);
  }

  if (new_f) {
    if (Mara->domain) delete Mara->domain;
    Mara->domain = new_f;
  }

  return 0;
}
int luaC_write_ppm(lua_State *L)
{
  const int narg = lua_gettop(L);
  const char *fname = luaL_checkstring(L, 1);
  double *data = luaU_checkarray(L, 2);
  int cmap = narg == 3 ? luaL_checkinteger(L, 3) : 0;
  double *range = narg == 4 ? luaU_checkarray(L, 4) : NULL;

  if (Mara->domain == NULL) {
    luaL_error(L, "need a domain to run this, use set_domain");
  }
  if (Mara->domain->get_Nd() != 2) {
    luaL_error(L, "need a 2d domain to write ppm images");
  }

  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);

  Mara_image_write_ppm(fname, data, cmap, Nx, Ny, range);

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
    new_f = new AdiabaticIdealSrhd(1.4);
  }
  else if (strcmp("rmhd", key) == 0) {
    new_f = new AdiabaticIdealRmhd;
  }
  else {
    luaL_error(L, "no such fluid: %s", key);
  }

  if (new_f) {
    if (Mara->fluid) delete Mara->fluid;
    Mara->fluid = new_f;
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
  else {
    luaL_error(L, "no such eos: %s", key);
  }

  if (new_f) {
    if (Mara->eos) delete Mara->eos;
    Mara->eos = new_f;
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
  else if (strcmp(key, "reflect2d") == 0) {
    int revx = luaL_checkinteger(L, 2);
    int revy = luaL_checkinteger(L, 3);
    new_f = new ReflectingBoundary2d(revx, revy);
  }

  if (new_f) {
    if (Mara->boundary) delete Mara->boundary;
    Mara->boundary = new_f;
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

  return 0;
}

int luaC_set_godunov(lua_State *L)
{
  const char *key = luaL_checkstring(L, 1);
  GodunovOperator *new_f = NULL;

  if (strcmp(key, "plm-split") == 0) {
    new_f = BuildPlmMethodOfLinesSplit(L);
  }
  else if (strcmp(key, "plm-muscl") == 0) {
    new_f = BuildPlmCtuHancockOperator(L);
  }
  else if (strcmp(key, "weno-split") == 0) {
    new_f = new WenoSplit;
  }
  else {
    luaL_error(L, "no such integration scheme: %s", key);
  }

  if (new_f) {
    if (Mara->godunov) delete Mara->godunov;
    Mara->godunov = new_f;
  }

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

  if (new_f) {
    if (Mara->advance) delete Mara->advance;
    Mara->advance = new_f;
  }

  return 0;
}

int luaC_set_driving(lua_State *L)
{
  size_t len;
  DrivingModule *new_f = NULL;
  const char *buf = luaL_checklstring(L, 1, &len);

  if (Mara->domain == NULL) {
    luaL_error(L, "need a domain to run this, use set_domain");
  }

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
  const double k1 = luaL_checknumber(L, 4);
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
  double *P = luaU_checkarray(L, 1);

  if (Mara->fluid == NULL) {
    luaL_error(L, "need a fluid to run this, use set_fluid");
  }
  else {
    int Nq = Mara->fluid->GetNq();
    double *U = (double*) malloc(Nq*sizeof(double));
    Mara->fluid->PrimToCons(P, U);
    luaU_pusharray(L, U, Nq);
    free(U);
  }

  return 1;
}
int luaC_fluid_ConsToPrim(lua_State *L)
{
  double *U = luaU_checkarray(L, 1);

  if (Mara->fluid == NULL) {
    luaL_error(L, "need a fluid to run this, use set_fluid");
  }

  int Nq = Mara->fluid->GetNq();
  double *P = (double*) malloc(Nq*sizeof(double));
  Mara->fluid->ConsToPrim(U, P);
  luaU_pusharray(L, P, Nq);
  free(P);

  return 1;
}

int luaC_fluid_Eigensystem(lua_State *L)
{
  double *P = luaU_checkarray(L, 1);
  int dim = luaL_checkinteger(L, 2);

  if (Mara->fluid == NULL) {
    luaL_error(L, "need a fluid to run this, use set_fluid");
  }

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
    luaU_pusharray_astable(L, Lv+i*Nq, Nq);
    lua_settable(L, -3);
  }

  lua_newtable(L);
  for (int i=0; i<Nq; ++i) {
    lua_pushnumber(L, i+1);
    luaU_pusharray_astable(L, Rv+i*Nq, Nq);
    lua_settable(L, -3);
  }
  luaU_pusharray(L, lm, Nq);

  free(U);
  free(Lv);
  free(Rv);
  free(lm);

  return 3;
}

int luaC_boundary_ApplyBoundaries(lua_State *L)
{
  if (Mara->boundary == NULL) {
    printf("[mara] error: need a boundary conditions to run this, use set_boundary.\n");
  }
  else {
    Mara->boundary->ApplyBoundaries(Mara->PrimitiveArray);
  }
  return 0;
}

int luaC_eos_TemperatureMeV(lua_State *L)
{
  const double D = luaL_checknumber(L, 1);
  const double p = luaL_checknumber(L, 2);

  if (Mara->eos == NULL) {
    luaL_error(L, "need an eos to run this, use set_eos");
    return 0;
  }
  else {
    const double T = Mara->eos->TemperatureMeV(D, p);
    lua_pushnumber(L, T);
    return 1;
  }
}
int luaC_eos_Temperature_p(lua_State *L)
{
  const double D = luaL_checknumber(L, 1);
  const double p = luaL_checknumber(L, 2);

  if (Mara->eos == NULL) {
    luaL_error(L, "need an eos to run this, use set_eos");
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
    luaL_error(L, "need an eos to run this, use set_eos");
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
    luaL_error(L, "need an eos to run this, use set_eos");
    return 0;
  }
  else {
    const double p = Mara->eos->Pressure(D, T);
    lua_pushnumber(L, p);
    return 1;
  }
}
int luaC_eos_DensUpper(lua_State *L)
{
  if (Mara->eos == NULL) {
    luaL_error(L, "need an eos to run this, use set_eos");
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
    luaL_error(L, "need an eos to run this, use set_eos");
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
    luaL_error(L, "need an eos to run this, use set_eos");
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
    luaL_error(L, "need an eos to run this, use set_eos");
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
    luaL_error(L, "need a units system to run this, use set_units");
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
    luaL_error(L, "need a units system to run this, use set_units");
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
    luaL_error(L, "need a units system to run this, use set_units");
    return 0;
  }
  else {
    lua_pushnumber(L, Mara->units->Velocity());
    return 1;
  }
}

int luaC_units_MeVPerCubicFemtometer(lua_State *L)
{
  if (Mara->units == NULL) {
    luaL_error(L, "need a units system to run this, use set_units");
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
    luaL_error(L, "need a units system to run this, use set_units");
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
    luaL_error(L, "need a units system to run this, use set_units");
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
  double *P = luaU_checkarray(L, 1);

  double U[8], Q[8];

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

