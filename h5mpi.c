
/*------------------------------------------------------------------------------
 * FILE: h5mpi.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 * (1) http://www.hdfgroup.org/HDF5/doc/RM/CollectiveCalls.html
 *
 * (2) http://www.hdfgroup.org/HDF5/Tutor/parallel.html
 *
 * (3) http://www.hdfgroup.org/HDF5/Tutor/property.html
 *
 * (4) http://www.nics.tennessee.edu/io-tips#subsettingio
 *
 * (5) http://www.nersc.gov/nusers/help/tutorials/io/6.php (see bottom of file)
 *
 *
 *
 * (1) Describes which HDF5 functions must be called in parallel, in a reading
 *     or writing context.
 *
 *------------------------------------------------------------------------------
 */

#include "config.h"
#if (__MARA_USE_HDF5_PAR)
#define __MARA_IO_INCL_PRIVATE_DEFS
#define __MARA_IO_HDF5_PRIVATE_DEFS


#include <time.h>
#include <mpi.h>
#include <hdf5.h>
#include "mara_io.h"


void _io_write_prim_h5mpi(const char *fname, const char **pnames, const double *data)
// -----------------------------------------------------------------------------
// This function uses a collective MPI-IO procedure to write the contents of
// 'data' to the HDF5 file named 'fname', which is assumed to have been created
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

  int i;
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

  // Here we define collective (MPI) access to the file with alignment
  // properties optimized for the local file system, according to DiskBlockSize.
  // ---------------------------------------------------------------------------
  if (EnableChunking) {
    H5Pset_chunk(dcpl, n_dims, ChunkSize);
  }
  if (EnableAlignment) {
    H5Pset_alignment(fapl, AlignThreshold, DiskBlockSize);
  }

  H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

  hid_t file = H5Fopen(fname, H5F_ACC_RDWR, fapl);
  const int overwrite = H5Lexists(file, "prim", H5P_DEFAULT);
  hid_t prim = overwrite ? H5Gopen(file, "prim", H5P_DEFAULT) :
    H5Gcreate(file, "prim", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t mspc = H5Screate_simple(ndp1  , l_ntot, NULL);
  hid_t fspc = H5Screate_simple(n_dims, G_ntot, NULL);
  // Call signature to H5Sselect_hyperslab is (start, stride, count, chunk)
  // ---------------------------------------------------------------------------
  const clock_t start_all = clock();

  for (i=0; i<n_prim; ++i) {

    hid_t dset = overwrite ? H5Dopen(prim, pnames[i], H5P_DEFAULT) : 
      H5Dcreate(prim, pnames[i], H5T_NATIVE_DOUBLE, fspc,
		H5P_DEFAULT, dcpl, H5P_DEFAULT);

    l_strt[ndp1 - 1] = i;
    H5Sselect_hyperslab(mspc, H5S_SELECT_SET, l_strt, stride, a_nint, NULL);
    H5Sselect_hyperslab(fspc, H5S_SELECT_SET, G_strt,   NULL, A_nint, NULL);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspc, fspc, dxpl, data);
    H5Dclose(dset);
  }
  if (iolog) {
    const double sec = (double)(clock() - start_all) / CLOCKS_PER_SEC;
    fprintf(iolog, "[h5mpi] write to %s took %f minutes\n", fname, sec/60.0);
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
  H5Pclose(dcpl);
  H5Pclose(fapl);
}
void _io_read_prim_h5mpi(const char *fname, const char **pnames, double *data)
{
  hsize_t ndp1 = n_dims + 1;
  hsize_t *a_nint = (hsize_t*) malloc(ndp1*sizeof(hsize_t));
  hsize_t *l_ntot = (hsize_t*) malloc(ndp1*sizeof(hsize_t));
  hsize_t *l_strt = (hsize_t*) malloc(ndp1*sizeof(hsize_t));
  hsize_t *stride = (hsize_t*) malloc(ndp1*sizeof(hsize_t));

  int i;
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

  // Here we define collective (MPI) access to the file with alignment
  // properties optimized for the local file system, according to DiskBlockSize.
  // ---------------------------------------------------------------------------
  H5Pset_alignment(fapl, AlignThreshold, DiskBlockSize);
  H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

  hid_t file = H5Fopen(fname, H5F_ACC_RDONLY, fapl);
  hid_t prim = H5Gopen(file, "prim", H5P_DEFAULT);
  hid_t mspc = H5Screate_simple(ndp1  , l_ntot, NULL);
  hid_t fspc = H5Screate_simple(n_dims, G_ntot, NULL);

  // Call signature to H5Sselect_hyperslab is (start, stride, count, chunk)
  // ---------------------------------------------------------------------------
  const clock_t start_all = clock();

  for (i=0; i<n_prim; ++i) {
    hid_t dset = H5Dopen(prim, pnames[i], H5P_DEFAULT);
    l_strt[ndp1 - 1] = i;
    H5Sselect_hyperslab(mspc, H5S_SELECT_SET, l_strt, stride, a_nint, NULL);
    H5Sselect_hyperslab(fspc, H5S_SELECT_SET, G_strt,   NULL, A_nint, NULL);
    H5Dread(dset, H5T_NATIVE_DOUBLE, mspc, fspc, dxpl, data);
    H5Dclose(dset);
  }
  if (iolog) {
    const double sec = (double)(clock() - start_all) / CLOCKS_PER_SEC;
    fprintf(iolog, "[h5mpi] read from %s took %f minutes\n", fname, sec/60.0);
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


#else // No parallel HDF5 available
void _io_write_prim_h5mpi(const char *fname, const char **pnames, const double *data)
{

}
void _io_read_prim_h5mpi(const char *fname, const char **pnames, double *data)
{

}
#endif // (__MARA_USE_HDF5_PAR)




/*
  The contents of the following URL is so pertinent, that I am including it
  explicitly here.

  http://www.nersc.gov/nusers/help/tutorials/io/6.php

  6. Optimizations for HDF5 on Lustre

  The Lustre file system performs best when writes are aligned to stripe
  boundaries. This minimizes the overhead of the distributed lock manager which
  has to serialize access when a stripe is accessed by more than one client. Two
  methods to achieve stripe alignment on writes are described below.

  Chunking and Alignment for Balanced 3D Grids When every MPI task has the same
  amount of data to write into a 3D grid, a data layout techique called chunking
  can be enabled in HDF5. This layout maps the 3D block on each task to its own
  contiguous 1D section of the file on disk. Figure 6.1 shows an example in 2D
  where 4 MPI tasks each write their own 4x4 block contiguously into a file, and
  an index in the HDF5 metadata keeps track of where each block resides in the
  file. A second HDF5 tuning parameter called the alignment can be set to pad
  out the size of each chunk to a multiple of the stripe size. Although padding
  wastes space in the file on disk, the gain in write bandwidth often outweighs
  this, especially on HPC systems like Franklin and Hopper with large scratch
  file systems. A more problematic restriction of the chunked layout, however,
  is that future parallel accesses to the file (for instance with an analysis or
  visualization tool) are optimized for reads that are multiples of the chunk
  size.


  Figure 6.1. Four tasks write to a contiguous 2D grid on the left versus a
  chunked grid on the right.

  In our 3D code example, chunking and alignment can be enabled using the
  following:


  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_alignment(fapl, 0, stripe_size);
  file = H5Fcreate("myparfile.h5", H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl, 3, chunk_dims);
  H5Dcreate(file, "mydataset", type, filespace, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  HDF5 uses a B-tree to index the location of chunks inside the file. The
  default size for this B-tree occupies only a few kilobytes in the HDF5
  metadata entry. When this space is exceeded, additional B-trees are created,
  leading to many small metadata writes when thousands of chunks are written to
  a file. The small metadata writes, even when they are padded and aligned using
  the alignment parameter still adversely affect write performance. HDF5
  provides the following mechanism for increasing the default size of the B-tree
  so that it is roughly the same size as a stripe:


  btree_ik = (stripe_size - 4096) / 96;
  fcpl = H5Pcreate(H5P_FILE_CREATE);
  H5Pset_istore_k(fcpl, btree_ik);
  file = H5Fcreate("myparfile.h5", H5F_ACC_TRUNC, fcpl, fapl);


  HDF5 also maintains a metadata cache for open files, and if enough metadata
  fills the cache, it can cause an eviction that interrupts other file
  accesses. Usually, it is better to disabled these evictions

  Metadata in HDF5 files is cached by the HDF5 library to improve access times
  for frequently accessed items. When operating in a sequential application,
  individual metadata items are flushed to the file (if dirty) and evicted from
  the metadata cache. However, when operating in a parallel application, these
  operations are deferred and batched together into eviction epochs, to reduce
  communication and synchronization overhead. At the end of an eviction epoch
  (measured by the amount of dirty metadata produced), the processes in the
  application are synchronized and the oldest dirty metadata items are flushed
  to the file.

  To reduce the frequency of performing small I/O operations, it is possible to
  put the eviction of items from the HDF5 library's metadata cache entirely
  under the application's control with the following:


  mdc_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
  H5Pget_mdc_config(file, &mdc_config)
  mdc_config.evictions_enabled = FALSE;
  mdc_config.incr_mode = H5C_incr__off;
  mdc_config.decr_mode = H5C_decr__off;
  H5Pset_mdc_config(file, &mdc_config);


  This sequence of calls disables evictions from the metadata cache, unless
  H5Fflush is called or the file is closed. Suspending automatic eviction of
  cached metadata items also prevents frequently dirtied items from being
  written to the file repeatedly. Suspending metadata evictions may not be
  appropriate for all applications however, because if the application crashes
  before the cached metadata is written to the file, the HDF5 file will be
  unusable.

  Collective Buffering

  In the case where each task writes different amounts of data, chunking cannot
  be used. Instead, a different optimization called Collective buffering (or
  two-phase IO can be used if collective mode has been enabled, as described in
  Section 4. In fact, collective buffering can also be used in the balanced 3D
  grid case discussed above instead of chunking, and the optimization in HDF5
  for disabling cache evictions also applies to collective buffering.

  Collective buffering works by breaking file accesses down into to stages. For
  a collective read, the first stage uses a subset of MPI tasks (called
  aggregators) to communicate with the IO servers (OSTs in Lustre) and read a
  large chunk of data into a temporary buffer. In the second stage, the
  aggregators ship the data from the buffer to its destination among the
  remaining MPI tasks using point-to-point MPI calls. A collective write does
  the reverse, aggregating the data through MPI into buffers on the aggregator
  nodes, then writing from the aggregator nodes to the IO servers. The advantage
  of collective buffering is that fewer nodes are communicating with the IO
  servers, which reduces contention. In fact, Lustre prefers a one-to-one
  mapping of aggregator nodes to OSTs.

  Since the release of mpt/3.3, Cray has included a Lustre-aware implementation
  of the MPI-IO collective buffering algorithm. This implementation is able to
  buffer data on the aggregator nodes into stripe-sized chunks so that all read
  and writes to the Lustre filesystem are automatically stripe aligned without
  requiring any padding or manual alignment from the developer. Because of the
  way Lustre is designed, alignment is a key factor in achieving optimal
  performance.

  Several environment variables can be used to control the behavior of
  collective buffering on Franklin and Hopper. The MPIIO_MPICH_HINTS variable
  specifies hints to the MPI-IO library that can, for instance, override the
  built-in heuristic and force collective buffering on:


  % setenv MPIIO_MPICH_HINTS="*:romio_cb_write=enable:romio_ds_write=disable"


  Placing this command in your batch file before calling aprun will cause your
  program to use these hints. The * indicates that the hint applies to any file
  opened by MPI-IO, while romio_cb_write controls collective buffering for
  writes and romio_ds_write controls data sieving for writes, an older
  collective mode optimization that is no longer used and can interfere with
  collective buffering. The options for these hints are enabled, disabled, or
  automatic (the default value, which uses the built-in heuristic).

  It is also possible to control the number of aggregator nodes using the
  cb_nodes hint, although the MPI-IO library will automatically set this to the
  stripe count of your file.

  When set to 1, the MPICH_MPIIO_HINTS_DISPLAY variable causes your program to
  dump a summary of the current MPI-IO hints to stderr each time a file is
  opened. This is useful for debugging and as a sanity check againt spelling
  errors in your hints.

  The MPICH_MPIIO_XSTATS variable enables profiling of the collective buffering
  algorithm, including binning of read/write calls and timings for the two
  phases. Setting it to 1 provides summary data, while setting it to 2 or 3
  provides more detail.

  More detail on MPICH runtime environment variables, including a full list and
  description of MPI-IO hints, is available from the intro_mpi man page on
  Franklin or Hopper.
*/
