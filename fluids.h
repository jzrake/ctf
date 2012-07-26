
#ifndef FLUIDS_HEADER_INCLUDED
#define FLUIDS_HEADER_INCLUDED

#define FLUIDS_LOCATION          (1<<0)
#define FLUIDS_PASSIVE           (1<<1)
#define FLUIDS_CONSERVED         (1<<2)
#define FLUIDS_PRIMITIVE         (1<<3)
#define FLUIDS_MAGNETIC          (1<<4)
#define FLUIDS_FOURVELOCITY      (1<<5)
#define FLUIDS_FLUX0             (1<<6)
#define FLUIDS_FLUX1             (1<<7)
#define FLUIDS_FLUX2             (1<<8)
#define FLUIDS_EIGENVALUES0      (1<<9)
#define FLUIDS_EIGENVALUES1      (1<<10)
#define FLUIDS_EIGENVALUES2      (1<<11)
#define FLUIDS_LEIGENVECTORS0    (1<<12)
#define FLUIDS_LEIGENVECTORS1    (1<<13)
#define FLUIDS_LEIGENVECTORS2    (1<<14)
#define FLUIDS_REIGENVECTORS0    (1<<15)
#define FLUIDS_REIGENVECTORS1    (1<<16)
#define FLUIDS_REIGENVECTORS2    (1<<17)
#define FLUIDS_SOUNDSPEEDSQUARED (1<<18)
#define FLUIDS_TEMPERATURE       (1<<19)
#define FLUIDS_SPECIFICENTHALPY  (1<<20)
#define FLUIDS_SPECIFICINTERNAL  (1<<21)
#define FLUIDS_GAMMALAWINDEX     (1<<22)
#define FLUIDS_FLUXALL           (FLUIDS_FLUX0|	\
				  FLUIDS_FLUX1|	\
				  FLUIDS_FLUX2)
#define FLUIDS_EIGENVALUESALL    (FLUIDS_EIGENVALUES0|	\
				  FLUIDS_EIGENVALUES1|	\
				  FLUIDS_EIGENVALUES2)
#define FLUIDS_LEIGENVECTORSALL  (FLUIDS_LEIGENVECTORS0|	\
				  FLUIDS_LEIGENVECTORS1|	\
				  FLUIDS_LEIGENVECTORS2)
#define FLUIDS_REIGENVECTORSALL  (FLUIDS_REIGENVECTORS0|	\
				  FLUIDS_REIGENVECTORS1|	\
				  FLUIDS_REIGENVECTORS2)
#define FLUIDS_FLAGSALL          ((1<<30) - 1)

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

#ifdef FLUIDS_INDEX_VARS
enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive
#endif // FLUIDS_INDEX_VARS

struct fluid_state;
typedef struct fluid_state fluid_state;

fluid_state *fluids_new(void);
int fluids_del(fluid_state *S);
int fluids_update(fluid_state *S, long flags);
int fluids_c2p(fluid_state *S);
int fluids_p2c(fluid_state *S);
int fluids_setfluid(fluid_state *S, int fluid);
int fluids_seteos(fluid_state *S, int eos);
int fluids_setcoordsystem(fluid_state *S, int coordsystem);
int fluids_setnpassive(fluid_state *S, int n);
int fluids_getattrib(fluid_state *S, double *x, long flag);
int fluids_setattrib(fluid_state *S, double *x, long flag);

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
  double *magnetic;
  double *fourvelocity;
  double *flux[3];
  double *eigenvalues[3];
  double *leigenvectors[3];
  double *reigenvectors[3];
  double soundspeedsquared;
  double temperature;
  double specificenthalpy;
  double specificinternal;
  double gammalawindex;
} ;
#endif // FLUIDS_PRIVATE_DEFS
#endif // FLUIDS_HEADER_INCLUDED
