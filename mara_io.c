

#include "config.h"
#define __MARA_IO_INCL_PRIVATE_DEFS
#if (__MARA_USE_HDF5)
#define __MARA_IO_HDF5_PRIVATE_DEFS
#include <hdf5.h>
#endif
#if (__MARA_USE_MPI)
#include <mpi.h>
#endif
#include "mara_io.h"

#if (__MARA_USE_HDF5)
struct StatusAttribute_t *Attributes = NULL;
size_t NumAttributes = 0;
hid_t MeasurementType;
hsize_t MeasurementSize;
hsize_t *ChunkSize = NULL;
hsize_t *A_nint = NULL;
hsize_t *L_ntot = NULL;
hsize_t *L_strt = NULL;
hsize_t *G_ntot = NULL;
hsize_t *G_strt = NULL;
#endif

int DiskBlockSize = 1;
int AlignThreshold = 0;

enum MaraIoFunction OutputFunction;
enum MaraIoFunction  InputFunction;

size_t n_dims;
size_t n_prim;
int mpi_rank;
int mpi_size;
FILE *iolog = NULL;
int TotalLocalZones;
int BinaryFileHeader[9];
int EnableChunking;
int EnableAlignment;


void Mara_io_init(const size_t measure_size,
                  const int n_dims_,
                  const int n_prim_,
                  const int *A_nint_,
                  const int *L_ntot_,
                  const int *L_strt_,
                  const int *G_ntot_,
                  const int *G_strt_)
{
  int i;
  _io_set_mpi_rank_size();
  TotalLocalZones = 1;

  for (i=0; i<3; ++i) {
    if (i<n_dims_) {
      BinaryFileHeader[i + 0] = G_strt_[i];
      BinaryFileHeader[i + 3] = G_ntot_[i];
      BinaryFileHeader[i + 6] = A_nint_[i];

      TotalLocalZones *= L_ntot_[i];
    }
    else {
      BinaryFileHeader[i + 0] = 0;
      BinaryFileHeader[i + 3] = 1;
      BinaryFileHeader[i + 6] = 1;
    }
  }

  n_dims = n_dims_;
  n_prim = n_prim_;

#if (__MARA_USE_HDF5)

  EnableChunking = 0;
  EnableAlignment = 0;

  MeasurementType = measure_size > 0 ? H5Tcreate(H5T_COMPOUND, measure_size) : 0;
  MeasurementSize = measure_size;
  Attributes =   (struct StatusAttribute_t*)
    malloc(sizeof(struct StatusAttribute_t));

  A_nint = (hsize_t*) malloc(n_dims*sizeof(hsize_t));
  L_ntot = (hsize_t*) malloc(n_dims*sizeof(hsize_t));
  L_strt = (hsize_t*) malloc(n_dims*sizeof(hsize_t));
  G_ntot = (hsize_t*) malloc(n_dims*sizeof(hsize_t));
  G_strt = (hsize_t*) malloc(n_dims*sizeof(hsize_t));
  ChunkSize = (hsize_t*) malloc(n_dims*sizeof(hsize_t));

  for (i=0; i<n_dims; ++i) {
    A_nint[i] = A_nint_[i]; // Selection size, target and destination
    L_ntot[i] = L_ntot_[i]; // Memory space total size
    L_strt[i] = L_strt_[i]; // Memory space selection start
    G_ntot[i] = G_ntot_[i]; // Global space total size
    G_strt[i] = G_strt_[i]; // Global space selection start
    ChunkSize[i] = G_ntot_[i];
  }
#endif
}
void Mara_io_free()
{
#if (__MARA_USE_HDF5)
  if (MeasurementType) {
    H5Tclose(MeasurementType);
  }
  NumAttributes = 0;
  free(Attributes);
  free(ChunkSize);
  free(A_nint);
  free(L_ntot);
  free(L_strt);
  free(G_ntot);
  free(G_strt);
#endif
}
void Mara_io_set_logfile(FILE *iolog_)
{
  iolog = iolog_;
}
void Mara_io_set_chunk_size(const int *s)
{
#if (__MARA_USE_HDF5)
  int i;
  for (i=0; i<n_dims; ++i) {
    ChunkSize[i] = s[i];
  }
#endif
}
void Mara_io_set_enable_alignment(int s)
{
  EnableAlignment = s;
}
void Mara_io_set_enable_chunking(int s)
{
  EnableChunking = s;
}
void Mara_io_set_disk_block_size(int s)
{
  DiskBlockSize = s;
}
void Mara_io_set_disk_align_threshold(int s)
{
  AlignThreshold = s;
}
void Mara_io_set_output_function(enum MaraIoFunction s)
{
  OutputFunction = s;
}
void Mara_io_set_input_function(enum MaraIoFunction s)
{
  InputFunction = s;
}
void Mara_io_write_prim(const char *fname, const char **pnames, const double *data)
{
  switch (OutputFunction) {

  case MARA_IO_FUNC_H5SER:
    _io_write_prim_h5ser(fname, pnames, data);
    break;
  case MARA_IO_FUNC_H5MPI:
    _io_write_prim_h5mpi(fname, pnames, data);
    break;
  case MARA_IO_FUNC_BINARY:
    _io_write_prim_binary(fname, pnames, data);
    break;
  default:
    _io_write_prim_h5ser(fname, pnames, data);
    break;
  }
}
void Mara_io_read_prim(const char *fname, const char **pnames, double *data)
{
  switch (InputFunction) {

  case MARA_IO_FUNC_H5SER:
    _io_read_prim_h5ser(fname, pnames, data);
    break;
  case MARA_IO_FUNC_H5MPI:
    _io_read_prim_h5mpi(fname, pnames, data);
    break;
  case MARA_IO_FUNC_BINARY:
    _io_read_prim_binary(fname, pnames, data);
    break;
  default:
    _io_read_prim_h5ser(fname, pnames, data);
    break;
  }
}
void Mara_io_make_name_meta(char *fname, const char *dir, const char *base, int num)
{
  sprintf(fname, "%s/%s.%04d.h5", dir, base, num);
}
void Mara_io_make_name_out(char *fname, const char *dir, const char *base, int num)
{
  switch (OutputFunction) {

  case MARA_IO_FUNC_H5SER:
    sprintf(fname, "%s/%s.%04d.h5", dir, base, num);
    break;
  case MARA_IO_FUNC_H5MPI:
    sprintf(fname, "%s/%s.%04d.h5", dir, base, num);
    break;
  case MARA_IO_FUNC_BINARY:
    sprintf(fname, "%s/%s.%04d.%05d.bin", dir, base, num, mpi_rank);
    break;
  default:
    sprintf(fname, "chkpt");
    break;
  }
}
void Mara_io_make_name_in(char *fname, const char *dir, const char *base, int num)
{
  switch (InputFunction) {

  case MARA_IO_FUNC_H5SER:
    sprintf(fname, "%s/%s.%04d.h5", dir, base, num);
    break;
  case MARA_IO_FUNC_H5MPI:
    sprintf(fname, "%s/%s.%04d.h5", dir, base, num);
    break;
  case MARA_IO_FUNC_BINARY:
    sprintf(fname, "%s/%s.%04d.%05d.bin", dir, base, mpi_rank, num);
    break;
  default:
    sprintf(fname, "chkpt");
    break;
  }
}

// The only explicit MPI calls needed by IO routines are wrapped below. This is
// done so that Mara can be compiled with HDF5 but not MPI routines.
// -----------------------------------------------------------------------------
void _io_barrier()
{
#if (__MARA_USE_MPI)
  int run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  if (run_uses_mpi) {
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif
}
void _io_set_mpi_rank_size()
{
#if (__MARA_USE_MPI)
  int run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  if (run_uses_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  }
  else {
    mpi_rank = 0;
    mpi_size = 1;
  }
#else
  mpi_rank = 0;
  mpi_size = 1;
#endif
}
// -----------------------------------------------------------------------------
