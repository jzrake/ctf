

#define __MARA_IO_INCL_PRIVATE_DEFS
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include "mara_io.h"



void _io_write_prim_binary(const char *fname, const char **pnames, const double *data)
{
  // fname format: data/dir/chkpt.0010.%05d.bin
  // ---------------------------------------------------------------------------
  char fullname[512];
  sprintf(fullname, fname, mpi_rank);
  FILE *outf = fopen(fname, "wb");
  fwrite(BinaryFileHeader, 9, sizeof(int), outf);
  fwrite(data, TotalLocalZones*n_prim, sizeof(double), outf);
  fclose(outf);
}
void _io_read_prim_binary(const char *fname, const char **pnames, double *data)
{
  // This function is only valid for serial reads
  // ---------------------------------------------------------------------------
  assert(mpi_size == 1);
  FILE *outf = fopen(fname, "rb");
  int header_in[9];
  fread(header_in, 9, sizeof(int), outf);
  fread(data, TotalLocalZones*n_prim, sizeof(double), outf);
  fclose(outf);
}
