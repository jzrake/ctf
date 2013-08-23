
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define COW_PRIVATE_DEFS
#include "cow.h"
#define KILOBYTES (1<<10)
#define MODULE "hdf5"

#if (COW_HDF5)
static int _io_write(cow_dfield *f, char *fname);
static int _io_read(cow_dfield *f, char *fname);
static int _io_check_file_exists(char *fname);
#endif

void _io_domain_commit(cow_domain *d)
{
#if (COW_HDF5)
  for (int n=0; n<d->n_dims; ++n) {
    d->L_nint_h5[n] = d->L_nint[n]; // Selection size, target and destination
    d->L_ntot_h5[n] = d->L_ntot[n]; // Memory space total size
    d->L_strt_h5[n] = d->L_strt[n]; // Memory space selection start
    d->G_ntot_h5[n] = d->G_ntot[n]; // Global space total size
    d->G_strt_h5[n] = d->G_strt[n]; // Global space selection start
  }
  // Here we create the following property lists:
  //
  // file access property list   ........ for the call to H5Fopen
  // dset creation property list ........ for the call to H5Dcreate
  // dset transfer property list ........ for the call to H5Dwrite
  // ---------------------------------------------------------------------------
  d->fapl = H5Pcreate(H5P_FILE_ACCESS);
  d->dcpl = H5Pcreate(H5P_DATASET_CREATE);
  d->dxpl = H5Pcreate(H5P_DATASET_XFER);
#if (COW_HDF5_MPI && COW_MPI)
  if (cow_mpirunning()) {
    H5Pset_fapl_mpio(d->fapl, d->mpi_cart, MPI_INFO_NULL);
  }
#endif // COW_HDF5_MPI && COW_MPI
#endif // COW_HDF5
}
void _io_domain_del(cow_domain *d)
{
#if (COW_HDF5)
  H5Pclose(d->fapl);
  H5Pclose(d->dcpl);
  H5Pclose(d->dxpl);
#endif
}

void cow_domain_setcollective(cow_domain *d, int mode)
{
#if (COW_HDF5 && COW_HDF5_MPI)
  if (mode && !cow_mpirunning()) {
    printf("[%s] requested collective without MPI running: "
	   "revert to independent\n", MODULE);
    mode = 0;
  }
  if (mode) {
    printf("[%s] setting HDF5 io mode to collective\n", MODULE);
    H5Pset_dxpl_mpio(d->dxpl, H5FD_MPIO_COLLECTIVE);
  }
  else {
    printf("[%s] setting HDF5 io mode to independent\n", MODULE);
    H5Pset_dxpl_mpio(d->dxpl, H5FD_MPIO_INDEPENDENT);
  }
#endif
}
void cow_domain_setchunk(cow_domain *d, int mode)
{
#if (COW_HDF5)
  if (mode) {
    if (d->balanced) {
      printf("[%s] enabled chunking on HDF5 files\n", MODULE);
      H5Pset_chunk(d->dcpl, d->n_dims, d->L_nint_h5);
    }
    else {
      printf("[%s] chunking could not be enabled because the domain is not "
	     "balanced\n", MODULE);
    }
  }
  else {
    printf("[%s] disabled chunking on HDF5 files\n", MODULE);
  }
#endif
}
void cow_domain_setalign(cow_domain *d, int alignthreshold, int diskblocksize)
{
#if (COW_HDF5)
  printf("[%s] align threshold: %d kB, disk block size: %d kB\n",
	 MODULE, alignthreshold/KILOBYTES, diskblocksize/KILOBYTES);
  H5Pset_alignment(d->fapl, alignthreshold, diskblocksize);
#endif
}
void cow_domain_readsize(cow_domain *d, char *fname, char *dname)
{
#if (COW_HDF5)
  if (_io_check_file_exists(fname)) return;
  hid_t file = H5Fopen(fname, H5F_ACC_RDONLY, d->fapl);
  hid_t dset = H5Dopen(file, dname, H5P_DEFAULT);
  hid_t fspc = H5Dget_space(dset);
  hsize_t dims[3] = { 1, 1, 1 };
  int ndims = H5Sget_simple_extent_dims(fspc, dims, NULL);
  H5Dclose(dset);
  H5Fclose(file);
  cow_domain_setndim(d, ndims);
  for (int n=0; n<ndims; ++n) {
    cow_domain_setsize(d, n, dims[n]);
  }
  printf("[%s] inferred global domain size of (%lld %lld %lld) from %s/%s\n",
	 MODULE, dims[0], dims[1], dims[2], fname, dname);
#endif
}
int cow_dfield_write(cow_dfield *f, char *fname)
{
#if (COW_HDF5)
#if (COW_MPI)
  if (f->domain->cart_rank == 0) {
#endif
    // -------------------------------------------------------------------------
    // The write functions assume the file is already created. Have master
    // create the file if it's not there already.
    // -------------------------------------------------------------------------
    FILE *testf = fopen(fname, "r");
    hid_t fid;
    if (testf == NULL) {
      fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    else {
      fclose(testf);
      fid = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
    }
    if (H5Lexists(fid, f->name, H5P_DEFAULT)) {
      H5Gunlink(fid, f->name);
    }
    hid_t memb = H5Gcreate(fid, f->name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(memb);
    H5Fclose(fid);
#if (COW_MPI)
  }
#endif
  clock_t start = clock();
  _io_write(f, fname);
  double sec = (double)(clock() - start) / CLOCKS_PER_SEC;
  printf("[%s] write to %s/%s took %f minutes\n", MODULE, fname, f->name,
	 sec/60.0);
  fflush(stdout);
#endif
  return 0;
}
int cow_dfield_read(cow_dfield *f, char *fname)
{
  int err = 0;
#if (COW_HDF5)
  if (_io_check_file_exists(fname)) return 1;
  clock_t start = clock();
  err = _io_read(f, fname);
  cow_dfield_syncguard(f);
  double sec = (double)(clock() - start) / CLOCKS_PER_SEC;
  printf("[%s] read from %s/%s took %f minutes\n", MODULE, fname, f->name,
	 sec/60.0);
  fflush(stdout);
#endif
  return err;
}

#if (COW_HDF5)
int _io_write(cow_dfield *f, char *fname)
// -----------------------------------------------------------------------------
// This function uses a collective MPI-IO procedure to write the contents of
// 'data' to the HDF5 file named 'fname', which is assumed to have been created
// already. The dataset with name 'dname', which is being written to, must not
// exist already. Chunking is enabled as per the ChunkSize variable, and is
// disabled by default. Recommended chunk size is local subdomain size. This
// will result in optimized read/write on the same decomposition layout, but
// poor performance for different access patterns, for example the slabs used by
// cluster-FFT functions.
//
//                                   WARNING!
//
// All processors must define the same chunk size, the behavior of this function
// is not defined otherwise. This implies that chunking should be disabled when
// running on a strange number of cores, and subdomain sizes are non-uniform.
// -----------------------------------------------------------------------------
{
  cow_domain *d = f->domain;
  char **pnames = f->members;
  void *data = f->data;
  char *gname = f->name;
  int n_memb = f->n_members;
  int n_dims = d->n_dims;
  hsize_t *L_nint = d->L_nint_h5;
  hsize_t *G_strt = d->G_strt_h5;
  hsize_t *G_ntot = d->G_ntot_h5;

  hsize_t ndp1 = n_dims + 1;
  hsize_t l_nint[4];
  hsize_t l_ntot[4];
  hsize_t l_strt[4];
  hsize_t stride[4];

  for (int i=0; i<n_dims; ++i) {
    l_nint[i] = d->L_nint[i]; // Selection size, target and destination
    l_ntot[i] = d->L_ntot[i]; // Memory space total size
    l_strt[i] = d->L_strt[i]; // Memory space selection start
    stride[i] = 1;
  }
  l_nint[ndp1 - 1] = 1;
  l_ntot[ndp1 - 1] = n_memb;
  stride[ndp1 - 1] = n_memb;

  // The loop over processors is needed if COW_MPI support is enabled and
  // COW_HDF5_MPI is not. If either COW_MPI is disabled, or COW_HDF5_MPI is
  // enabled, then the write calls occur without the loop.
  // ---------------------------------------------------------------------------
  int sequential = 0;
  int rank = 0;
#if (!COW_HDF5_MPI && COW_MPI)
  sequential = 1;
  for (rank=0; rank<d->cart_size; ++rank) {
    if (rank == d->cart_rank) {
#endif
      hid_t file = H5Fopen(fname, H5F_ACC_RDWR, d->fapl);
      hid_t memb = H5Gopen(file, gname, H5P_DEFAULT);
      hid_t mspc = H5Screate_simple(ndp1, l_ntot, NULL);
      hid_t fspc = H5Screate_simple(n_dims, G_ntot, NULL);
      hid_t dset;
      for (int n=0; n<n_memb; ++n) {
	if (sequential && rank != 0) {
	  dset = H5Dopen(memb, pnames[n], H5P_DEFAULT);
	}
	else {
	  dset = H5Dcreate(memb, pnames[n], H5T_NATIVE_DOUBLE, fspc,
			   H5P_DEFAULT, d->dcpl, H5P_DEFAULT);
	}
	l_strt[ndp1 - 1] = n;
	H5Sselect_hyperslab(mspc, H5S_SELECT_SET, l_strt, stride, l_nint, NULL);
	H5Sselect_hyperslab(fspc, H5S_SELECT_SET, G_strt, NULL, L_nint, NULL);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspc, fspc, d->dxpl, data);
	H5Dclose(dset);
      }
      H5Sclose(fspc);
      H5Sclose(mspc);
      H5Gclose(memb);
      H5Fclose(file);
#if (!COW_HDF5_MPI && COW_MPI)
    }
    if (cow_mpirunning()) {
      MPI_Barrier(d->mpi_cart);
    }
  }
#endif // !COW_HDF5_MPI && COW_MPI
  return 0;
}
#endif

#if (COW_HDF5)
int _io_read(cow_dfield *f, char *fname)
{
  cow_domain *d = f->domain;
  char **pnames = f->members;
  void *data = f->data;
  char *gname = f->name;
  int n_memb = f->n_members;
  int n_dims = d->n_dims;
  hsize_t *L_nint = d->L_nint_h5;
  hsize_t *G_strt = d->G_strt_h5;
  hsize_t *G_ntot = d->G_ntot_h5;

  hsize_t ndp1 = n_dims + 1;
  hsize_t l_nint[4];
  hsize_t l_ntot[4];
  hsize_t l_strt[4];
  hsize_t stride[4];

  for (int i=0; i<n_dims; ++i) {
    l_nint[i] = d->L_nint[i]; // Selection size, target and destination
    l_ntot[i] = d->L_ntot[i]; // Memory space total size
    l_strt[i] = d->L_strt[i]; // Memory space selection start
    stride[i] = 1;
  }
  l_nint[ndp1 - 1] = 1;
  l_ntot[ndp1 - 1] = n_memb;
  stride[ndp1 - 1] = n_memb;

  // The loop over processors is needed if COW_MPI support is enabled and
  // COW_HDF5_MPI is not. If either COW_MPI is disabled, or COW_HDF5_MPI is
  // enabled, then the write calls occur without the loop.
  // ---------------------------------------------------------------------------
#if (!COW_HDF5_MPI && COW_MPI)
  for (int rank=0; rank<d->cart_size; ++rank) {
    if (rank == d->cart_rank) {
#endif
      hid_t file = H5Fopen(fname, H5F_ACC_RDONLY, d->fapl);
      hid_t memb = H5Gopen(file, gname, H5P_DEFAULT);
      hid_t mspc = H5Screate_simple(ndp1, l_ntot, NULL);
      hid_t fspc = H5Screate_simple(n_dims, G_ntot, NULL);
      for (int n=0; n<n_memb; ++n) {
	hid_t dset = H5Dopen(memb, pnames[n], H5P_DEFAULT);
	l_strt[ndp1 - 1] = n;
	H5Sselect_hyperslab(mspc, H5S_SELECT_SET, l_strt, stride, l_nint, NULL);
	H5Sselect_hyperslab(fspc, H5S_SELECT_SET, G_strt, NULL, L_nint, NULL);
	H5Dread(dset, H5T_NATIVE_DOUBLE, mspc, fspc, d->dxpl, data);
	H5Dclose(dset);
      }
      H5Sclose(fspc);
      H5Sclose(mspc);
      H5Gclose(memb);
      H5Fclose(file);
#if (!COW_HDF5_MPI && COW_MPI)
    }
    if (cow_mpirunning()) {
      MPI_Barrier(d->mpi_cart);
    }
  }
#endif // !COW_HDF5_MPI && COW_MPI
  return 0;
}
#endif

int _io_check_file_exists(char *fname)
{
  FILE *testf = fopen(fname, "r");
  if (testf == NULL) {
    printf("[%s] error: file does not exist %s\n", MODULE, fname);
    return 1;
  }
  else {
    fclose(testf);
    return 0;
  }
}



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
