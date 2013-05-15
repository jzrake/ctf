

#ifndef COW_HEADER_INCLUDED
#define COW_HEADER_INCLUDED
#include <stdlib.h>

#ifdef COW_PRIVATE_DEFS
#if (COW_MPI)
#include <mpi.h>
#endif // COW_MPI
#if (COW_HDF5)
#include <hdf5.h>
#endif // COW_HDF5
#endif // COW_PRIVATE_DEFS


#define COW_NOREOPEN_STDOUT   (1<<0)
#define COW_DISABLE_MPI       (1<<1)
#define COW_HASNAN            (1<<2)
#define COW_HASINF            (1<<3)

#define COW_ALL_DIMS             -41
#define COW_HIST_SPACING_LINEAR  -42
#define COW_HIST_SPACING_LOG     -43
#define COW_HIST_BINMODE_COUNTS  -44 // traditional histogram
#define COW_HIST_BINMODE_DENSITY -45 // divides by bin width
#define COW_HIST_BINMODE_AVERAGE -46
#define COW_PROJECT_OUT_DIV      -47 // used for Helmholtz decomposition
#define COW_PROJECT_OUT_CURL     -48
#define COW_SAMPLE_NEAREST       -49 // sample the nearest zone center
#define COW_SAMPLE_LINEAR        -50 // use (uni/bi/tri) linear interp
#define COW_SAMPLE_ERROR_OUT     -51 // out-of-bounds sample request
#define COW_SAMPLE_ERROR_WRONGD  -52 // wrong number of dims on sample coords

// -----------------------------------------------------------------------------
//
// These prototypes constitute the C.O.W. interface
//
// -----------------------------------------------------------------------------
struct cow_domain; // forward declarations (for opaque data structure)
struct cow_dfield;
typedef struct cow_domain cow_domain;
typedef struct cow_dfield cow_dfield;
typedef struct cow_histogram cow_histogram;
typedef void (*cow_transform)(double *result, double **args, int **strides,
			      void *udata);

void cow_init(int argc, char **argv, int modes);
void cow_finalize(void);
int cow_mpirunning(void);

cow_domain *cow_domain_new(void);
void cow_domain_commit(cow_domain *d);
void cow_domain_del(cow_domain *d);
void cow_domain_setsize(cow_domain *d, int dim, int size);
void cow_domain_setndim(cow_domain *d, int ndim);
void cow_domain_setguard(cow_domain *d, int guard);
void cow_domain_setprocsizes(cow_domain *d, int dim, int size);
void cow_domain_setcollective(cow_domain *d, int mode);
void cow_domain_setchunk(cow_domain *d, int mode);
void cow_domain_setalign(cow_domain *d, int alignthreshold, int diskblocksize);
void cow_domain_readsize(cow_domain *d, char *fname, char *dname);
int cow_domain_getndim(cow_domain *d);
int cow_domain_getsize(cow_domain *d, int dim);
int cow_domain_getguard(cow_domain *d);
long long cow_domain_getnumlocalzonesincguard(cow_domain *d, int dim);
long long cow_domain_getnumlocalzonesinterior(cow_domain *d, int dim);
long long cow_domain_getnumglobalzones(cow_domain *d, int dim);
int cow_domain_getglobalstartindex(cow_domain *d, int dim);
double cow_domain_getlowercoord(cow_domain *d, int dim); // local physical extent
double cow_domain_getuppercoord(cow_domain *d, int dim);
double cow_domain_getgridspacing(cow_domain *d, int dim);
int cow_domain_getcartrank(cow_domain *d);
int cow_domain_getcartsize(cow_domain *d);
int cow_domain_subgridatposition(cow_domain *d, double x, double y, double z);
int cow_domain_indexatposition(cow_domain *d, int dim, double x);
double cow_domain_positionatindex(cow_domain *d, int dim, int index);
void cow_domain_barrier(cow_domain *d);
int cow_domain_intprod(cow_domain *d, int myval);
int cow_domain_intsum(cow_domain *d, int myval);
double cow_domain_dblmin(cow_domain *d, double myval);
double cow_domain_dblmax(cow_domain *d, double myval);
double cow_domain_dblsum(cow_domain *d, double myval);
void cow_domain_getcomm(cow_domain *d, void *comm);

cow_dfield *cow_dfield_new(void);
cow_dfield *cow_dfield_dup(cow_dfield *f);
void cow_dfield_commit(cow_dfield *f);
void cow_dfield_del(cow_dfield *f);
void cow_dfield_setdomain(cow_dfield *f, cow_domain *d);
void cow_dfield_addmember(cow_dfield *f, char *name);
void cow_dfield_setname(cow_dfield *f, char *name);
void cow_dfield_extract(cow_dfield *f, int *I0, int *I1, void *out);
void cow_dfield_replace(cow_dfield *f, int *I0, int *I1, void *out);
void cow_dfield_loop(cow_dfield *f, cow_transform op, void *udata);
void cow_dfield_settransform(cow_dfield *f, cow_transform op);
void cow_dfield_clearargs(cow_dfield *f);
void cow_dfield_pusharg(cow_dfield *f, cow_dfield *arg);
void cow_dfield_setuserdata(cow_dfield *f, void *userdata);
void cow_dfield_setiparam(cow_dfield *f, int p);
void cow_dfield_setfparam(cow_dfield *f, double p);
void cow_dfield_transformexecute(cow_dfield *f);
char *cow_dfield_iteratemembers(cow_dfield *f);
char *cow_dfield_nextmember(cow_dfield *f);
char *cow_dfield_getname(cow_dfield *f);
cow_domain *cow_dfield_getdomain(cow_dfield *f);
int cow_dfield_getstride(cow_dfield *f, int dim);
int cow_dfield_getnmembers(cow_dfield *f);
size_t cow_dfield_getdatabytes(cow_dfield *f);
void cow_dfield_setdatabuffer(cow_dfield *f, void *buffer);
void cow_dfield_sampleglobalind(cow_dfield *f, int i, int j, int k, double **x, int *n0);
int cow_dfield_setsamplecoords(cow_dfield *f, double *x, int n0, int n1);
void cow_dfield_getsamplecoords(cow_dfield *f, double **x, int *n0, int *n1);
void cow_dfield_getsampleresult(cow_dfield *f, double **x, int *n0, int *n1);
void cow_dfield_setsamplemode(cow_dfield *f, int mode);
void cow_dfield_sampleexecute(cow_dfield *f);
int cow_dfield_getownsdata(cow_dfield *f);
void *cow_dfield_getdatabuffer(cow_dfield *f);
void cow_dfield_syncguard(cow_dfield *f);
void cow_dfield_reduce(cow_dfield *f, double *x);
int cow_dfield_write(cow_dfield *f, char *fname);
int cow_dfield_read(cow_dfield *f, char *fname);


cow_histogram *cow_histogram_new(void);
void cow_histogram_commit(cow_histogram *h);
void cow_histogram_del(cow_histogram *h);
void cow_histogram_setbinmode(cow_histogram *h, int binmode);
void cow_histogram_setspacing(cow_histogram *h, int spacing);
void cow_histogram_setnbins(cow_histogram *h, int dim, int nbinsx);
void cow_histogram_setlower(cow_histogram *h, int dim, double v0);
void cow_histogram_setupper(cow_histogram *h, int dim, double v1);
void cow_histogram_setfullname(cow_histogram *h, char *fullname);
void cow_histogram_setnickname(cow_histogram *h, char *nickname);
void cow_histogram_setdomaincomm(cow_histogram *h, cow_domain *d);
void cow_histogram_addsample1(cow_histogram *h, double x, double w);
void cow_histogram_addsample2(cow_histogram *h, double x, double y, double w);
void cow_histogram_dumpascii(cow_histogram *h, char *fn);
void cow_histogram_dumphdf5(cow_histogram *h, char *fn, char *dn);
void cow_histogram_seal(cow_histogram *h);
int cow_histogram_getsealed(cow_histogram *h);
long cow_histogram_gettotalcounts(cow_histogram *h);
void cow_histogram_populate(cow_histogram *h, cow_dfield *f, cow_transform op);
void cow_histogram_getbinlocx(cow_histogram *h, double *x);
void cow_histogram_getbinlocy(cow_histogram *h, double *x);
void cow_histogram_getbinvalv(cow_histogram *h, double *x);
double cow_histogram_getbinval(cow_histogram *h, int i, int j);
char *cow_histogram_getname(cow_histogram *h);

void cow_fft_forward(cow_dfield *f, cow_dfield *fkre, cow_dfield *fkim);
void cow_fft_reverse(cow_dfield *f, cow_dfield *fkre, cow_dfield *fkim);
void cow_fft_pspecscafield(cow_dfield *f, cow_histogram *h);
void cow_fft_pspecvecfield(cow_dfield *f, cow_histogram *h);
void cow_fft_helmholtzdecomp(cow_dfield *f, int mode);
void cow_fft_solvepoisson(cow_dfield *rho, cow_dfield *phi);

void cow_trans_rot5(double *result, double **args, int **s, void *u);
void cow_trans_div5(double *result, double **args, int **s, void *u);
void cow_trans_divcorner(double *result, double **args, int **s, void *u);
void cow_trans_laplacian(double *result, double **args, int **s, void *u);
void cow_trans_component(double *result, double **args, int **s, void *u);
void cow_trans_magnitude(double *result, double **args, int **s, void *u);
void cow_trans_cross(double *result, double **args, int **s, void *u);
void cow_trans_dot3(double *result, double **args, int **s, void *u);


#ifdef COW_PRIVATE_DEFS

void _io_domain_commit(cow_domain *d);
void _io_domain_del(cow_domain *d);

struct cow_domain
{
  double glb_lower[3]; // lower coordinates of global physical domain
  double glb_upper[3]; // upper " "
  double loc_lower[3]; // lower coordinates of local physical domain
  double loc_upper[3]; // upper " "
  double dx[3]; // grid spacing along each direction
  int L_nint[3]; // interior zones on the local subgrid
  int L_ntot[3]; // total " ", including guard zones
  int L_strt[3]; // starting index of interior zones on local subgrid
  int G_ntot[3]; // global domain size
  int G_strt[3]; // starting index into global domain
  int n_dims; // number of dimensions: 1, 2, 3
  int n_ghst; // number of guard zones: >= 0
  int balanced; // true when all subgrids have the same size
  int committed; // true after cow_domain_commit called, locks out size changes
#if (COW_MPI)
  int comm_rank; // rank with respect to MPI_COMM_WORLD communicator
  int comm_size; // size " "
  int cart_rank; // rank with respect to the cartesian communicator
  int cart_size; // size " "
  int proc_sizes[3]; // number of subgrids along each dimension
  int proc_index[3]; // coordinates of local subgrid in cartesian communicator
  int num_neighbors; // 3, 9, or 27 depending on the domain dimensionality
  int *neighbors; // cartesian ranks of the neighboring processors
  int *send_tags; // tag used to on send calls with respective neighbor
  int *recv_tags; // " "            recv " "
  MPI_Comm mpi_cart; // the cartesian communicator
#endif
#if (COW_HDF5)
  hsize_t L_nint_h5[3]; // HDF5 versions of the variables with the same name
  hsize_t L_ntot_h5[3];
  hsize_t L_strt_h5[3];
  hsize_t G_ntot_h5[3];
  hsize_t G_strt_h5[3];
  hid_t fapl; // file access property list
  hid_t dcpl; // data set creation property list
  hid_t dxpl; // data set transfer property list
#endif
} ;

struct cow_dfield
{
  char *name; // name of the data field
  char **members; // list of labels for the data members
  int member_iter; // maintains an index into the last dimension
  int n_members; // size of last dimension
  void *data; // data buffer
  int stride[3]; // strides describing memory layout: C ordering
  int committed; // true after cow_dfield_commit called, locks out most changes
  int ownsdata; // client code can own the data: see setdatabuffer function
  cow_domain *domain; // pointer to an associated domain
  cow_transform transform; // used only by internal code
  cow_dfield **transargs; // list of arguments for transform, used internally
  void *userdata; // shallow pointer to user-supplied data item
  int transargslen;
  int iparam; // extra parameters that might be used by transform functions
  double dparam;
  double *samplecoords;
  double *sampleresult;
  int samplecoordslen;
  int samplemode;
#if (COW_MPI)
  MPI_Datatype *send_type; // chunk of data to be sent to respective neighbor
  MPI_Datatype *recv_type; // " "                 received from " "
#endif
} ;

struct cow_histogram
{
  int nbinsx;
  int nbinsy;
  double x0;
  double x1;
  double y0;
  double y1;
  double *bedgesx;
  double *bedgesy;
  double *weight;
  long totcounts;
  long *counts;
  char *nickname;
  char *fullname;
  int binmode;
  int spacing;
  int n_dims;
  int committed;
  int sealed; // once sealed, is sync'ed and does not accept more samples
  cow_transform transform;
  double *binlocx; // Pointers to these arrays are returned by the getbinlocx
  double *binlocy; // getbinlocy, and getbinval functions.
  double *binvalv;
#if (COW_MPI)
  MPI_Comm comm;
#endif
} ;

#endif // COW_PRIVATE_DEFS
#endif // COW_HEADER_INCLUDED
