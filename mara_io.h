
#ifdef __cplusplus
extern "C" {
#endif

#ifndef __MaraIoModule_HEADER__
#define __MaraIoModule_HEADER__

#include <stdio.h>
#include <stdlib.h>

// Begin MaraIoModule public interface
// -----------------------------------------------------------------------------
enum MaraIoFunction { MARA_IO_FUNC_H5SER,
		      MARA_IO_FUNC_H5MPI,
		      MARA_IO_FUNC_BINARY };

void Mara_io_free();
void Mara_io_init(const size_t measure_size,
                  const int n_dims_,
                  const int n_prim_,
                  const int *A_nint_,
                  const int *L_ntot_,
                  const int *L_strt_,
                  const int *G_ntot_,
                  const int *G_strt_);

void Mara_io_set_logfile(FILE *iolog_);
void Mara_io_set_chunk_size(const int *s);
void Mara_io_set_disk_block_size(int s);
void Mara_io_set_disk_align_threshold(int s);
void Mara_io_set_output_function(enum MaraIoFunction s);
void Mara_io_set_input_function(enum MaraIoFunction s);
void Mara_io_set_enable_chunking(int s);
void Mara_io_set_enable_alignment(int s);

size_t Mara_io_get_config_size(const char *fname);
size_t Mara_io_get_measlog_size(const char *fname);
size_t Mara_io_get_bits_size(const char *fname, const char *dname);

void Mara_io_add_int_measurement(const char *name, const int *offset, int num);
void Mara_io_add_int_status(const char *name, int *data);
void Mara_io_add_dbl_measurement(const char *name, const double *offset, int num);
void Mara_io_add_dbl_status(const char *name, double *data);

void Mara_io_create_file(const char *fname);
void Mara_io_write_config(const char *fname, const char *config);
void Mara_io_write_prim(const char *fname, const char **pnames, const double *data);
void Mara_io_write_measlog(const char *fname, size_t size, const void *buffer);
void Mara_io_write_status(const char *fname);
void Mara_io_write_bits(const char *fname, const char *dname, size_t size, const void *buffer);

void Mara_io_read_config(const char *fname, char *config);
void Mara_io_read_prim(const char *fname, const char **pnames, double *data);
void Mara_io_read_measlog(const char *fname, void *buffer);
void Mara_io_read_status(const char *fname);
void Mara_io_read_bits(const char *fname, const char *dname, void *buffer);

void Mara_io_make_name_meta(char *fname, const char *dir, const char *base, int num);
void Mara_io_make_name_out(char *fname, const char *dir, const char *base, int num);
void Mara_io_make_name_in(char *fname, const char *dir, const char *base, int num);



// Begin MaraIoModule private interface
// -----------------------------------------------------------------------------
#ifdef __MARA_IO_INCL_PRIVATE_DEFS
extern int DiskBlockSize;
extern int AlignThreshold;

extern enum MaraIoFunction OutputFunction;
extern enum MaraIoFunction  InputFunction;

extern size_t n_dims;
extern size_t n_prim;
extern int mpi_rank;
extern int mpi_size;
extern FILE *iolog;
extern int TotalLocalZones;
extern int BinaryFileHeader[9];
extern int EnableChunking;
extern int EnableAlignment;

void _io_barrier();
void _io_set_mpi_rank_size();

void _io_write_prim_h5mpi(const char *fname, const char **pnames, const double *data);
void _io_write_prim_h5ser(const char *fname, const char **pnames, const double *data);
void _io_write_prim_binary(const char *fname, const char **pnames, const double *data);
void _io_read_prim_h5mpi(const char *fname, const char **pnames, double *data);
void _io_read_prim_h5ser(const char *fname, const char **pnames, double *data);
void _io_read_prim_binary(const char *fname, const char **pnames, double *data);

#endif // __MARA_IO_INCL_PRIVATE_DEFS


#ifdef __MARA_IO_HDF5_PRIVATE_DEFS

struct StatusAttribute_t
{
  const char *name;
  void *data;
  hid_t type;
} ;

extern struct StatusAttribute_t *Attributes;
extern size_t NumAttributes;
extern hid_t MeasurementType;
extern hsize_t MeasurementSize;
extern hsize_t *ChunkSize;
extern hsize_t *A_nint;
extern hsize_t *L_ntot;
extern hsize_t *L_strt;
extern hsize_t *G_ntot;
extern hsize_t *G_strt;


#endif // __MARA_IO_HDF5_PRIVATE_DEFS
#endif // __MaraIoModule_HEADER__

#ifdef __cplusplus
}
#endif
