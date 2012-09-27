
#ifndef FLUIDS_HEADER_INCLUDED
#define FLUIDS_HEADER_INCLUDED

#define FLUIDS_PRIMITIVE         (1<<1)
#define FLUIDS_PASSIVE           (1<<2)
#define FLUIDS_GRAVITY           (1<<3)
#define FLUIDS_MAGNETIC          (1<<4)
#define FLUIDS_LOCATION          (1<<5)
#define FLUIDS_CONSERVED         (1<<6)
#define FLUIDS_FOURVELOCITY      (1<<7)
#define FLUIDS_FLUX0             (1<<8)
#define FLUIDS_FLUX1             (1<<9)
#define FLUIDS_FLUX2             (1<<10)
#define FLUIDS_EVAL0             (1<<11)
#define FLUIDS_EVAL1             (1<<12)
#define FLUIDS_EVAL2             (1<<13)
#define FLUIDS_LEVECS0           (1<<14)
#define FLUIDS_LEVECS1           (1<<15)
#define FLUIDS_LEVECS2           (1<<16)
#define FLUIDS_REVECS0           (1<<17)
#define FLUIDS_REVECS1           (1<<18)
#define FLUIDS_REVECS2           (1<<19)
#define FLUIDS_JACOBIAN0         (1<<20)
#define FLUIDS_JACOBIAN1         (1<<21)
#define FLUIDS_JACOBIAN2         (1<<22)
#define FLUIDS_SOUNDSPEEDSQUARED (1<<23)
#define FLUIDS_TEMPERATURE       (1<<24)
#define FLUIDS_SPECIFICENTHALPY  (1<<25)
#define FLUIDS_SPECIFICINTERNAL  (1<<26)
#define FLUIDS_FLAGSALL          ((1<<30) - 1)

#define FLUIDS_FLUXALL           (FLUIDS_FLUX0|FLUIDS_FLUX1|FLUIDS_FLUX2)
#define FLUIDS_EVALSALL          (FLUIDS_EVAL0|FLUIDS_EVAL1|FLUIDS_EVAL2)
#define FLUIDS_LEVECSALL         (FLUIDS_LEVECS0|FLUIDS_LEVECS1|FLUIDS_LEVECS2)
#define FLUIDS_REVECSALL         (FLUIDS_REVECS0|FLUIDS_REVECS1|FLUIDS_REVECS2)
#define FLUIDS_JACOBIANALL       (FLUIDS_JACOBIAN0|FLUIDS_JACOBIAN1|FLUIDS_JACOBIAN2)

#define FLUIDS_SCADV             -41 // Scalar advection
#define FLUIDS_SCBRG             -42 // Burgers equation
#define FLUIDS_SHWAT             -43 // Shallow water equations
#define FLUIDS_NRHYD             -44 // Euler equations
#define FLUIDS_SRHYD             -45 // Special relativistic
#define FLUIDS_URHYD             -46 // Ultra relativistic
#define FLUIDS_GRHYD             -47 // General relativistic
#define FLUIDS_NRMHD             -48 // Magnetohydrodynamic (MHD)
#define FLUIDS_SRMHD             -49 // Special relativistic MHD
#define FLUIDS_GRMHD             -50 // General relativistic MHD

#define FLUIDS_EOS_GAMMALAW      -51
#define FLUIDS_EOS_TABULATED     -52

#define FLUIDS_COORD_CARTESIAN   -53
#define FLUIDS_COORD_SPHERICAL   -54
#define FLUIDS_COORD_CYLINDRICAL -55

#define FLUIDS_ERROR_BADARG      -66
#define FLUIDS_ERROR_BADREQUEST  -67
#define FLUIDS_ERROR_RIEMANN     -68
#define FLUIDS_ERROR_INCOMPLETE  -69

#define FLUIDS_RIEMANN_HLL       -70
#define FLUIDS_RIEMANN_HLLC      -71
#define FLUIDS_RIEMANN_EXACT     -72

#define FLUIDS_CACHE_NOTOUCH     -73
#define FLUIDS_CACHE_RESET       -74
#define FLUIDS_CACHE_ERASE       -75

#ifdef FLUIDS_INDEX_VARS
enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive
enum { phi, phidot, gradphi };             // Gravity
#endif // FLUIDS_INDEX_VARS


typedef struct fluids_descr fluids_descr;
typedef struct fluids_cache fluids_cache;
typedef struct fluids_state fluids_state;
typedef struct fluids_riemn fluids_riemn;


/* fluids_descr member functions */
fluids_descr *fluids_descr_new(void);
int fluids_descr_del(fluids_descr *D);
int fluids_descr_getfluid(fluids_descr *D, int *fluid);
int fluids_descr_setfluid(fluids_descr *D, int fluid);
int fluids_descr_geteos(fluids_descr *D, int *eos);
int fluids_descr_seteos(fluids_descr *D, int eos);
int fluids_descr_getcoordsystem(fluids_descr *D, int *coordsystem);
int fluids_descr_setcoordsystem(fluids_descr *D, int coordsystem);
int fluids_descr_getgamma(fluids_descr *D, double *gam);
int fluids_descr_setgamma(fluids_descr *D, double gam);
int fluids_descr_getncomp(fluids_descr *D, long flag);

/* fluids_state member functions */
fluids_state *fluids_state_new(void);
int fluids_state_del(fluids_state *S);
int fluids_state_getdescr(fluids_state *S, fluids_descr **D);
int fluids_state_setdescr(fluids_state *S, fluids_descr *D);
int fluids_state_resetcache(fluids_state *S);
int fluids_state_erasecache(fluids_state *S);
int fluids_state_getattr(fluids_state *S, double *x, long flag);
int fluids_state_setattr(fluids_state *S, double *x, long flag);
int fluids_state_fromcons(fluids_state *S, double *U, int cachebehavior);
int fluids_state_derive(fluids_state *S, double *x, int flag);

/* fluids_riemn member functions */
fluids_riemn *fluids_riemn_new(void);
int fluids_riemn_del(fluids_riemn *R);
int fluids_riemn_setstateL(fluids_riemn *R, fluids_state *S);
int fluids_riemn_setstateR(fluids_riemn *R, fluids_state *S);
int fluids_riemn_setdim(fluids_riemn *R, int dim);
int fluids_riemn_execute(fluids_riemn *R);
int fluids_riemn_sample(fluids_riemn *R, fluids_state *S, double s);
int fluids_riemn_getsolver(fluids_riemn *R, int *solver);
int fluids_riemn_setsolver(fluids_riemn *R, int solver);


#ifdef FLUIDS_PRIVATE_DEFS

struct fluids_cache {
  double *conserved;
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
  long needsupdateflags;
  fluids_state *state;
} ;

struct fluids_descr {
  int fluid;
  int eos;
  int coordsystem;
  int nprimitive;
  int npassive;
  int ngravity;
  int nmagnetic;
  int nlocation;
  long cacheflags;
  double gammalawindex;
} ;

struct fluids_state {
  double *primitive;
  double *passive;
  double *gravity;
  double *magnetic;
  double *location;
  fluids_cache *cache;
  fluids_descr *descr;
} ;

/* http://en.wikipedia.org/wiki/Bitwise_operation#NOT */
#define BITWISENOT(x) (-(x) - 1)
#define MAPBUF 2
#define ALLOC 1
#define DEALLOC 0

#endif // FLUIDS_PRIVATE_DEFS

#endif // FLUIDS_HEADER_INCLUDED
