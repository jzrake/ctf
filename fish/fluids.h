
#ifndef FLUIDS_HEADER_INCLUDED
#define FLUIDS_HEADER_INCLUDED

enum {
  FLUIDS_PRIMITIVE         =  1<<1,
  FLUIDS_PASSIVE           =  1<<2,
  FLUIDS_GRAVITY           =  1<<3,
  FLUIDS_MAGNETIC          =  1<<4,
  FLUIDS_LOCATION          =  1<<5,
  FLUIDS_USERFLAG          =  1<<6,
  FLUIDS_CONSERVED         =  1<<7,
  FLUIDS_SOURCETERMS       =  1<<8,
  FLUIDS_FOURVELOCITY      =  1<<9,
  FLUIDS_FLUX0             =  1<<10,
  FLUIDS_FLUX1             =  1<<11,
  FLUIDS_FLUX2             =  1<<12,
  FLUIDS_EVAL0             =  1<<13,
  FLUIDS_EVAL1             =  1<<14,
  FLUIDS_EVAL2             =  1<<15,
  FLUIDS_LEVECS0           =  1<<16,
  FLUIDS_LEVECS1           =  1<<17,
  FLUIDS_LEVECS2           =  1<<18,
  FLUIDS_REVECS0           =  1<<19,
  FLUIDS_REVECS1           =  1<<20,
  FLUIDS_REVECS2           =  1<<21,
  FLUIDS_JACOBIAN0         =  1<<22,
  FLUIDS_JACOBIAN1         =  1<<23,
  FLUIDS_JACOBIAN2         =  1<<24,
  FLUIDS_SOUNDSPEEDSQUARED =  1<<25,
  FLUIDS_TEMPERATURE       =  1<<26,
  FLUIDS_SPECIFICENTHALPY  =  1<<27,
  FLUIDS_SPECIFICINTERNAL  =  1<<28,
  FLUIDS_FLAGSALL          = (1<<30) - 1,
  FLUIDS_FLUXALL     = FLUIDS_FLUX0|FLUIDS_FLUX1|FLUIDS_FLUX2,
  FLUIDS_EVALSALL    = FLUIDS_EVAL0|FLUIDS_EVAL1|FLUIDS_EVAL2,
  FLUIDS_LEVECSALL   = FLUIDS_LEVECS0|FLUIDS_LEVECS1|FLUIDS_LEVECS2,
  FLUIDS_REVECSALL   = FLUIDS_REVECS0|FLUIDS_REVECS1|FLUIDS_REVECS2,
  FLUIDS_JACOBIANALL = FLUIDS_JACOBIAN0|FLUIDS_JACOBIAN1|FLUIDS_JACOBIAN2,
} ;

enum {
  FLUIDS_SUCCESS,
  FLUIDS_ERROR_BADARG,
  FLUIDS_ERROR_BADREQUEST,
  FLUIDS_ERROR_RIEMANN,
  FLUIDS_ERROR_INCOMPLETE,
  FLUIDS_ERROR_NEGATIVE_DENSITY_CONS,
  FLUIDS_ERROR_NEGATIVE_DENSITY_PRIM,
  FLUIDS_ERROR_NEGATIVE_ENERGY,
  FLUIDS_ERROR_NEGATIVE_PRESSURE,
  FLUIDS_ERROR_SUPERLUMINAL,
  FLUIDS_ERROR_C2P_MAXITER,
  FLUIDS_ERROR_NOT_IMPLEMENTED,

  FLUIDS_COORD_CARTESIAN,
  FLUIDS_COORD_SPHERICAL,
  FLUIDS_COORD_CYLINDRICAL,

  FLUIDS_SCADV, // Scalar advection
  FLUIDS_SCBRG, // Burgers equation
  FLUIDS_SHWAT, // Shallow water equations
  FLUIDS_NRHYD, // Euler equations
  FLUIDS_GRAVS, // Gravitating Euler equations (with source terms on p and E)
  FLUIDS_GRAVP, // " "                         (no source term on p)
  FLUIDS_GRAVE, // " "                         (no source terms at all)
  FLUIDS_SRHYD, // Special relativistic
  FLUIDS_URHYD, // Ultra relativistic
  FLUIDS_GRHYD, // General relativistic
  FLUIDS_NRMHD, // Magnetohydrodynamic (MHD)
  FLUIDS_SRMHD, // Special relativistic MHD
  FLUIDS_GRMHD, // General relativistic MHD

  FLUIDS_EOS_GAMMALAW,
  FLUIDS_EOS_TABULATED,

  FLUIDS_RIEMANN_HLL,
  FLUIDS_RIEMANN_HLLC,
  FLUIDS_RIEMANN_EXACT,

  FLUIDS_CACHE_DEFAULT,
  FLUIDS_CACHE_NOTOUCH,
  FLUIDS_CACHE_CREATE,
  FLUIDS_CACHE_STEAL,
  FLUIDS_CACHE_RESET,
  FLUIDS_CACHE_ERASE,
} ;

#ifdef FLUIDS_INDEX_VARS
enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive
enum { phi=0, gph=1, phd=4, gpd=5 };       // Gravity
#endif // FLUIDS_INDEX_VARS


typedef struct fluids_descr fluids_descr;
typedef struct fluids_cache fluids_cache;
typedef struct fluids_state fluids_state;
typedef struct fluids_svect fluids_svect;
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
int fluids_descr_getrhobar(fluids_descr *D, double *rhobar);
int fluids_descr_setrhobar(fluids_descr *D, double rhobar);
int fluids_descr_getncomp(fluids_descr *D, long flag);


/* fluids_state member functions */
fluids_state *fluids_state_new(void);
int fluids_state_del(fluids_state *S);
int fluids_state_copy(fluids_state *S0, fluids_state *S1);
int fluids_state_getdescr(fluids_state *S, fluids_descr **D);
int fluids_state_setdescr(fluids_state *S, fluids_descr *D);
int fluids_state_getattr(fluids_state *S, double *x, long flag);
int fluids_state_setattr(fluids_state *S, double *x, long flag);
int fluids_state_fromcons(fluids_state *S, double *U, int cache);
int fluids_state_derive(fluids_state *S, double *x, long flags);
int fluids_state_getcached(fluids_state *S, double *x, long flag);
int fluids_state_mapbuffer(fluids_state *S, double *buffer, long flag);
int fluids_state_mapbufferuserflag(fluids_state *S, int *buffer);
int fluids_state_getuserflag(fluids_state *S, int *x);
int fluids_state_setuserflag(fluids_state *S, int *x);
int fluids_state_cache(fluids_state *S, int operation);


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
  double *sourceterms;
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
  double rhobar;
  fluids_cache *cache;
} ;

struct fluids_state {
  double *primitive;
  double *passive;
  double *gravity;
  double *magnetic;
  double *location;
  int *userflag;
  char ownscache;
  long ownsbufferflags;
  fluids_cache *cache;
  fluids_descr *descr;
} ;

#define MAPBUF 2
#define ALLOC 1
#define DEALLOC 0

#endif // FLUIDS_PRIVATE_DEFS
#endif // FLUIDS_HEADER_INCLUDED
