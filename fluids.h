
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
#define FLUIDS_EVAL0             (1<<9)
#define FLUIDS_EVAL1             (1<<10)
#define FLUIDS_EVAL2             (1<<11)
#define FLUIDS_LEVECS0           (1<<12)
#define FLUIDS_LEVECS1           (1<<13)
#define FLUIDS_LEVECS2           (1<<14)
#define FLUIDS_REVECS0           (1<<15)
#define FLUIDS_REVECS1           (1<<16)
#define FLUIDS_REVECS2           (1<<17)
#define FLUIDS_JACOBIAN0         (1<<18)
#define FLUIDS_JACOBIAN1         (1<<19)
#define FLUIDS_JACOBIAN2         (1<<20)
#define FLUIDS_SOUNDSPEEDSQUARED (1<<21)
#define FLUIDS_TEMPERATURE       (1<<22)
#define FLUIDS_SPECIFICENTHALPY  (1<<23)
#define FLUIDS_SPECIFICINTERNAL  (1<<24)
#define FLUIDS_GAMMALAWINDEX     (1<<25)
#define FLUIDS_FLAGSALL          ((1<<30) - 1)

#define FLUIDS_FLUXALL           (FLUIDS_FLUX0|FLUIDS_FLUX1|FLUIDS_FLUX2)
#define FLUIDS_EVALSALL          (FLUIDS_EVAL0|FLUIDS_EVAL1|FLUIDS_EVAL2)
#define FLUIDS_LEVECSALL         (FLUIDS_LEVECS0|FLUIDS_LEVECS1|FLUIDS_LEVECS2)
#define FLUIDS_REVECSALL         (FLUIDS_REVECS0|FLUIDS_REVECS1|FLUIDS_REVECS2)
#define FLUIDS_JACOBIANALL       (FLUIDS_JACOBIAN0|FLUIDS_JACOBIAN1|FLUIDS_JACOBIAN2)

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
#define FLUIDS_ERROR_RIEMANN     -68

#define FLUIDS_RIEMANN_HLL       -69
#define FLUIDS_RIEMANN_HLLC      -70
#define FLUIDS_RIEMANN_EXACT     -71


#ifdef FLUIDS_INDEX_VARS
enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive
#endif // FLUIDS_INDEX_VARS

struct fluid_state;
struct fluid_riemann;
typedef struct fluid_state fluid_state;
typedef struct fluid_riemann fluid_riemann;

fluid_state *fluids_new(void);
int fluids_del(fluid_state *S);
int fluids_update(fluid_state *S, long flags);
int fluids_resetcache(fluid_state *S);
int fluids_getlastupdate(fluid_state *S, long *flags);
int fluids_alloc(fluid_state *S, long flags);
int fluids_c2p(fluid_state *S);
int fluids_p2c(fluid_state *S);
int fluids_setfluid(fluid_state *S, int fluid);
int fluids_seteos(fluid_state *S, int eos);
int fluids_setcoordsystem(fluid_state *S, int coordsystem);
int fluids_setnpassive(fluid_state *S, int n);
int fluids_getattrib(fluid_state *S, double *x, long flag);
int fluids_setattrib(fluid_state *S, double *x, long flag);

fluid_riemann *fluids_riemann_new(void);
int fluids_riemann_del(fluid_riemann *R);
int fluids_riemann_setstateL(fluid_riemann *R, fluid_state *S);
int fluids_riemann_setstateR(fluid_riemann *R, fluid_state *S);
int fluids_riemann_setdim(fluid_riemann *R, int dim);
int fluids_riemann_execute(fluid_riemann *R);
int fluids_riemann_sample(fluid_riemann *R, fluid_state *S, double s);
int fluids_riemann_setsolver(fluid_riemann *R, int solver);

#ifdef FLUIDS_PRIVATE_DEFS
struct fluid_state {
  int fluid;
  int eos;
  int coordsystem;
  int nwaves;
  int npassive;
  long needsupdateflags;
  long lastupdatedflags;
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
  double *jacobian[3];
  double soundspeedsquared;
  double temperature;
  double specificenthalpy;
  double specificinternal;
  double gammalawindex;
} ;
#endif // FLUIDS_PRIVATE_DEFS

#endif // FLUIDS_HEADER_INCLUDED
