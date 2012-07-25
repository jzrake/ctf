
#ifndef FLUIDS_HEADER_INCLUDED
#define FLUIDS_HEADER_INCLUDED

#define FLUIDS_LOCATION          (2<<0)
#define FLUIDS_PASSIVE           (2<<1)
#define FLUIDS_CONSERVED         (2<<2)
#define FLUIDS_PRIMITIVE         (2<<3)
#define FLUIDS_FLUX0             (2<<4)
#define FLUIDS_FLUX1             (2<<5)
#define FLUIDS_FLUX2             (2<<6)
#define FLUIDS_MAGNETIC          (2<<7)
#define FLUIDS_FOURVELOCITY      (2<<8)
#define FLUIDS_EIGENVALUES       (2<<9)
#define FLUIDS_EIGENVECTORS      (2<<10)
#define FLUIDS_SOUNDSPEEDSQUARED (2<<11)
#define FLUIDS_TEMPERATURE       (2<<12)
#define FLUIDS_SPECIFICENTHALPY  (2<<13)
#define FLUIDS_SPECIFICINTERNAL  (2<<14)
#define FLUIDS_FLUXALL           (FLUIDS_FLUX0|FLUIDS_FLUX1|FLUIDS_FLUX2)
#define FLUIDS_FLAGSALL          ((2<<31) - 1)

#define FLUIDS_SCALAR_ADVECTION  -41
#define FLUIDS_SCALAR_BURGERS    -42
#define FLUIDS_SHALLOW_WATER     -43
#define FLUIDS_NRHYD             -44
#define FLUIDS_SRHYD             -45
#define FLUIDS_URHYD             -46
#define FLUIDS_GRHYD             -47
#define FLUIDS_NRMHD             -48
#define FLUIDS_SRMHD             -49
#define FLUIDS_GRMHD             -50

#define FLUIDS_EOS_GAMMALAW      -51
#define FLUIDS_EOS_TABULATED     -52

#define FLUIDS_COORD_CARTESIAN   -53
#define FLUIDS_COORD_SPHERICAL   -54
#define FLUIDS_COORD_CYLINDRICAL -55

#define FLUIDS_ERROR_BADARG      -66
#define FLUIDS_ERROR_BADREQUEST  -67

struct fluid_state;
typedef struct fluid_state fluid_state;

fluid_state *fluids_new(void);
int fluids_del(fluid_state *S);
int fluids_setfluid(fluid_state *S, int fluid);
int fluids_seteos(fluid_state *S, int eos);
int fluids_setcoordsystem(fluid_state *S, int coordsystem);
int fluids_setnpassive(fluid_state *S, int n);
int fluids_c2p(fluid_state *S);
int fluids_p2c(fluid_state *S);
int fluids_getlocation(fluid_state *S, double *x);
int fluids_setlocation(fluid_state *S, double *x);
int fluids_getpassive(fluid_state *S, double *x);
int fluids_setpassive(fluid_state *S, double *x);
int fluids_getconserved(fluid_state *S, double *x);
int fluids_setconserved(fluid_state *S, double *x);
int fluids_getprimitive(fluid_state *S, double *x);
int fluids_setprimitive(fluid_state *S, double *x);
int fluids_getflux0(fluid_state *S, double *x);
int fluids_setflux0(fluid_state *S, double *x);
int fluids_getflux1(fluid_state *S, double *x);
int fluids_setflux1(fluid_state *S, double *x);
int fluids_getflux2(fluid_state *S, double *x);
int fluids_setflux2(fluid_state *S, double *x);
int fluids_getmagnetic(fluid_state *S, double *x);
int fluids_setmagnetic(fluid_state *S, double *x);
int fluids_getfourvelocity(fluid_state *S, double *x);
int fluids_setfourvelocity(fluid_state *S, double *x);
int fluids_geteigenvalues(fluid_state *S, double *x);
int fluids_seteigenvalues(fluid_state *S, double *x);
int fluids_getleigenvectors(fluid_state *S, double *x);
int fluids_setleigenvectors(fluid_state *S, double *x);
int fluids_getreigenvectors(fluid_state *S, double *x);
int fluids_setreigenvectors(fluid_state *S, double *x);


int fluids_update(fluid_state *S, long flags);

#ifdef FLUIDS_PRIVATE_DEFS
struct fluid_state {
  int fluid;
  int eos;
  int coordsystem;
  int nwaves;
  int npassive;
  double *location;
  double *passive;
  double *conserved;
  double *primitive;
  double *flux0;
  double *flux1;
  double *flux2;
  double *magnetic;
  double *fourvelocity;
  double *eigenvalues;
  double *leigenvectors;
  double *reigenvectors;
  double soundspeedsquared;
  double temperature;
  double specificenthalpy;
  double specificinternal;
} ;
#endif // FLUIDS_PRIVATE_DEFS
#endif // FLUIDS_HEADER_INCLUDED
