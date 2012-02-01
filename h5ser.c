
#include "config.h"
#if (__MARA_USE_HDF5)
#define __MARA_IO_INCL_PRIVATE_DEFS
#define __MARA_IO_HDF5_PRIVATE_DEFS

#include <time.h>
#include <hdf5.h>
#include "mara_io.h"


void _io_write_prim_h5ser(const char *fname, const char **pnames, const double *data)
// -----------------------------------------------------------------------------
// This function uses a sequential IO procedure to write the contents of 'data'
// to the HDF5 file named 'fname', which is assumed to have been created
// already. The dataset with name 'dname', which is being written to, must not
// exist already. Chunking is enabled as per the module-wide ChunkSize variable,
// and is disabled by default. Recommended chunk size is local subdomain
// size. This will result in optimized read/write on the same decomposition
// layout, but poor performance for different access patterns, for example the
// slabs used by cluster-FFT functions.
//
//                                   WARNING!
//
// All processors must define the same chunk size, the behavior of this function
// is not defined otherwise. This implies that chunking should be disabled when
// running on a strange number of cores, and subdomain sizes are non-uniform.
// -----------------------------------------------------------------------------
{
  hsize_t ndp1 = n_dims + 1;
  hsize_t *a_nint = (hsize_t*) malloc(ndp1*sizeof(hsize_t));
  hsize_t *l_ntot = (hsize_t*) malloc(ndp1*sizeof(hsize_t));
  hsize_t *l_strt = (hsize_t*) malloc(ndp1*sizeof(hsize_t));
  hsize_t *stride = (hsize_t*) malloc(ndp1*sizeof(hsize_t));

  int i, rank;
  for (i=0; i<n_dims; ++i) {
    a_nint[i] = A_nint[i]; // Selection size, target and destination
    l_ntot[i] = L_ntot[i]; // Memory space total size
    l_strt[i] = L_strt[i]; // Memory space selection start
    stride[i] = 1;
  }
  a_nint[ndp1 - 1] = 1;
  l_ntot[ndp1 - 1] = n_prim;
  stride[ndp1 - 1] = n_prim;

  // Here we create the following property lists:
  //
  // file access property list   ........ for the call to H5Fopen
  // dset creation property list ........ for the call to H5Dcreate
  // dset transfer property list ........ for the call to H5Dwrite
  // ---------------------------------------------------------------------------
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  if (EnableChunking) {
    H5Pset_chunk(dcpl, n_dims, ChunkSize);
  }
  if (EnableAlignment) {
    H5Pset_alignment(fapl, AlignThreshold, DiskBlockSize);
  }

  if (mpi_rank == 0) {
    // Have only rank zero create the prim group and data sets.
    // -------------------------------------------------------------------------
    hid_t file = H5Fopen(fname, H5F_ACC_RDWR, fapl);

    if (H5Lexists(file, "prim", H5P_DEFAULT)) {
      // If the prim group already exists, assume the datasets do as well, and
      // move on.
      H5Fclose(file);
    }
    else {
      hid_t prim = H5Gcreate(file, "prim", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      hid_t fspc = H5Screate_simple(n_dims, G_ntot, NULL);

      for (i=0; i<n_prim; ++i) {
        hid_t dset = H5Dcreate(prim, pnames[i], H5T_NATIVE_DOUBLE, fspc,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Dclose(dset);
      }

      H5Sclose(fspc);
      H5Gclose(prim);
      H5Fclose(file);
    }
  }
  const clock_t start_all = clock();

  for (rank=0; rank<mpi_size; ++rank) {
    const clock_t start = clock();

    if (rank == mpi_rank) {
      hid_t file = H5Fopen(fname, H5F_ACC_RDWR, fapl);
      hid_t prim = H5Gopen(file, "prim", H5P_DEFAULT);
      hid_t mspc = H5Screate_simple(ndp1  , l_ntot, NULL);
      hid_t fspc = H5Screate_simple(n_dims, G_ntot, NULL);

      // Call signature to H5Sselect_hyperslab is (start, stride, count, chunk)
      // -----------------------------------------------------------------------
      for (i=0; i<n_prim; ++i) {
        hid_t dset = H5Dopen(prim, pnames[i], H5P_DEFAULT);
        l_strt[ndp1 - 1] = i;
        H5Sselect_hyperslab(mspc, H5S_SELECT_SET, l_strt, stride, a_nint, NULL);
        H5Sselect_hyperslab(fspc, H5S_SELECT_SET, G_strt,   NULL, A_nint, NULL);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspc, fspc, dxpl, data);
        H5Dclose(dset);
      }
      H5Sclose(fspc);
      H5Sclose(mspc);
      H5Gclose(prim);
      H5Fclose(file);
    }
    _io_barrier();
    if (0 && iolog) { // disabled
      const double sec = (double)(clock() - start) / CLOCKS_PER_SEC;
      fprintf(iolog, "[h5ser] rank %d wrote to %s in %f seconds\n",
              rank, fname, sec);
      fflush(iolog);
    }
  }
  if (iolog) {
    const double sec = (double)(clock() - start_all) / CLOCKS_PER_SEC;
    fprintf(iolog, "[h5ser] write to %s took %f minutes\n",
            fname, sec/60.0);
    fflush(iolog);
  }
  free(a_nint);
  free(l_ntot);
  free(l_strt);
  free(stride);

  // Always close the hid_t handles in the reverse order they were opened in.
  // ---------------------------------------------------------------------------
  H5Pclose(dxpl);
  H5Pclose(dcpl);
  H5Pclose(fapl);
}
void _io_read_prim_h5ser(const char *fname, const char **pnames, double *data)
{
  hsize_t ndp1 = n_dims + 1;
  hsize_t *a_nint = (hsize_t*) malloc(ndp1*sizeof(hsize_t));
  hsize_t *l_ntot = (hsize_t*) malloc(ndp1*sizeof(hsize_t));
  hsize_t *l_strt = (hsize_t*) malloc(ndp1*sizeof(hsize_t));
  hsize_t *stride = (hsize_t*) malloc(ndp1*sizeof(hsize_t));

  int i, rank;
  for (i=0; i<n_dims; ++i) {
    a_nint[i] = A_nint[i]; // Selection size, target and destination
    l_ntot[i] = L_ntot[i]; // Memory space total size
    l_strt[i] = L_strt[i]; // Memory space selection start
    stride[i] = 1;
  }
  a_nint[ndp1 - 1] = 1;
  l_ntot[ndp1 - 1] = n_prim;
  stride[ndp1 - 1] = n_prim;

  // Here we create the following property lists:
  //
  // file access property list   ........ for the call to H5Fopen
  // dset transfer property list ........ for the call to H5Dread
  // ---------------------------------------------------------------------------
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);

  H5Pset_alignment(fapl, AlignThreshold, DiskBlockSize);

  hid_t file = H5Fopen(fname, H5F_ACC_RDONLY, fapl);
  hid_t prim = H5Gopen(file, "prim", H5P_DEFAULT);
  hid_t mspc = H5Screate_simple(ndp1  , l_ntot, NULL);
  hid_t fspc = H5Screate_simple(n_dims, G_ntot, NULL);

  // Call signature to H5Sselect_hyperslab is (start, stride, count, chunk)
  // ---------------------------------------------------------------------------
  const clock_t start_all = clock();

  for (rank=0; rank<mpi_size; ++rank) {
    const clock_t start = clock();

    if (rank == mpi_rank) {
      for (i=0; i<n_prim; ++i) {
        hid_t dset = H5Dopen(prim, pnames[i], H5P_DEFAULT);
        l_strt[ndp1 - 1] = i;
        H5Sselect_hyperslab(mspc, H5S_SELECT_SET, l_strt, stride, a_nint, NULL);
        H5Sselect_hyperslab(fspc, H5S_SELECT_SET, G_strt,   NULL, A_nint, NULL);
        H5Dread(dset, H5T_NATIVE_DOUBLE, mspc, fspc, dxpl, data);
        H5Dclose(dset);
      }
    }

    _io_barrier();
    if (0 && iolog) { // disabled
      const double sec = (double)(clock() - start) / CLOCKS_PER_SEC;
      fprintf(iolog, "[h5ser] rank %d read from %s in %f seconds\n",
              rank, fname, sec);
      fflush(iolog);
    }
  }
  if (iolog) {
    const double sec = (double)(clock() - start_all) / CLOCKS_PER_SEC;
    fprintf(iolog, "[h5ser] read from %s took %f minutes\n",
            fname, sec/60.0);
    fflush(iolog);
  }

  free(a_nint);
  free(l_ntot);
  free(l_strt);
  free(stride);

  // Always close the hid_t handles in the reverse order they were opened in.
  // ---------------------------------------------------------------------------
  H5Sclose(fspc);
  H5Sclose(mspc);
  H5Gclose(prim);
  H5Fclose(file);
  H5Pclose(dxpl);
  H5Pclose(fapl);
}


#else // No HDF5 available
void _io_write_prim_h5ser(const char *fname, const char **pnames, const double *data)
{

}
void _io_read_prim_h5ser(const char *fname, const char **pnames, double *data)
{

}
#endif // (__MARA_USE_HDF5)
