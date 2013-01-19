


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "histogram.hpp"



Histogram1d::Histogram1d(int nbins, double x0, double x1,
                         enum Histogram::SpacingType spc)
  : nbins(nbins), nickname("histogram")
{
  bedges = new double[nbins+1];
  weight = new double[nbins];
  counts = new long[nbins];
  const double dx = (x1-x0) / nbins;

  for (int n=0; n<nbins+1; ++n) {
    if (spc == Histogram::Logspace) {
      bedges[n] = x0 * pow(x1/x0, double(n) / nbins);
    }
    else {
      bedges[n] = x0 + n * dx;
    }
  }
  for (int n=0; n<nbins; ++n) {
    counts[n] = 0;
    weight[n] = 0.0;
  }
}
Histogram2d::Histogram2d(int nbinsX, int nbinsY, double x0, double x1,
                         double y0, double y1, enum Histogram::SpacingType spc)
  : nbinsX(nbinsX), nbinsY(nbinsY), nickname("histogram")
{
  nbins = nbinsX*nbinsY;
  bedgesX = new double[nbinsX+1];
  bedgesY = new double[nbinsY+1];
  weight = new double[nbins];
  counts = new long[nbins];

  const double dx = (y1-y0) / nbinsX;
  const double dy = (y1-y0) / nbinsY;

  for (int n=0; n<nbinsX+1; ++n) {
    if (spc == Histogram::Logspace) {
      bedgesX[n] = x0 * pow(x1/x0, double(n) / nbinsX);
    }
    else {
      bedgesX[n] = x0 + n * dx;
    }
  }
  for (int n=0; n<nbinsY+1; ++n) {
    if (spc == Histogram::Logspace) {
      bedgesY[n] = y0 * pow(y1/y0, double(n) / nbinsY);
    }
    else {
      bedgesY[n] = y0 + n * dy;
    }
  }
  for (int n=0; n<nbinsX*nbinsY; ++n) {
    counts[n] = 0;
    weight[n] = 0.0;
  }
}



Histogram1d::~Histogram1d()
{
  delete [] bedges;
  delete [] weight;
  delete [] counts;
}
Histogram2d::~Histogram2d()
{
  delete [] bedgesX;
  delete [] bedgesY;
  delete [] weight;
  delete [] counts;
}



int Histogram1d::add_sample(double x, double w)
{
  for (int n=0; n<nbins; ++n) {
    if (bedges[n] < x && x < bedges[n+1]) {
      weight[n] += w;
      counts[n] += 1;
      return 0;
    }
  }
  return 1;
}

int Histogram2d::add_sample(double x, double y, double w)
{
  int nx=-1, ny=-1;
  for (int n=0; n<nbinsX; ++n) {
    if (bedgesX[n] < x && x < bedgesX[n+1]) {
      nx = n;
      break;
    }
  }
  for (int n=0; n<nbinsY; ++n) {
    if (bedgesY[n] < y && y < bedgesY[n+1]) {
      ny = n;
      break;
    }
  }
  if (nx == -1 || ny == -1) {
    return 1;
  }
  else {
    counts[nx * nbinsY + ny] += 1;
    weight[nx * nbinsY + ny] += w;
    return 0;
  }
}


void Histogram1d::dump_ascii(FILE *file)
{
  for (int n=0; n<nbins; ++n) {
    fprintf(file, "%f %f\n", 0.5*(bedges[n] + bedges[n+1]), weight[n]/counts[n]);
  }
}
void Histogram2d::dump_ascii(FILE *file)
{
  for (int nx=0; nx<nbinsX; ++nx) {
    for (int ny=0; ny<nbinsY; ++ny) {
      fprintf(file, "%f %f %f\n",
              0.5*(bedgesX[nx] + bedgesX[nx+1]),
              0.5*(bedgesY[ny] + bedgesY[ny+1]),
	      counts[nx * nbinsY + ny] / weight[nx * nbinsY + ny]);
    }
  }
}


#include "config.h"
// -----------------------------------------------------------------------------
// The code below is implemented when HDF5 or MPI are available on the system.
// -----------------------------------------------------------------------------
#if (__MARA_USE_HDF5)
#include <hdf5.h>


void Histogram1d::dump_hdf5(int base)
{
  if (base == 0) return;

  // Create a group to represent this histogram, and an attribute to name it
  // ---------------------------------------------------------------------------
  hid_t hgrp = H5Gcreate(base, nickname.c_str(),
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (!fullname.empty()) {

    hid_t aspc = H5Screate(H5S_SCALAR);
    hid_t strn = H5Tcopy(H5T_C_S1);
    H5Tset_size(strn, fullname.length());
    hid_t attr = H5Acreate(hgrp, "fullname", strn, aspc,
                           H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr, strn, fullname.c_str()); // write the full name

    H5Aclose(attr);
    H5Tclose(strn);
    H5Sclose(aspc);
  }

  // Create a temporary array for the bin centers
  // ---------------------------------------------------------------------------
  double *binval = new double[nbins];
  double *binloc = new double[nbins];
  for (int i=0; i<nbins; ++i) {
    const long c = counts[i];
    binloc[i] = 0.5*(bedges[i] + bedges[i+1]);

    switch (binning_mode) {
    case Histogram::BinAverage:
      binval[i] = (c == 0) ? 0.0 : weight[i] / c;
      break;
    case Histogram::BinDensity:
      binval[i] = weight[i] / (bedges[i+1] - bedges[i]);
      break;
    default:
      binval[i] = 0.0;
      break;
    }
  }

  // Create the data sets in the group: binloc (bin centers) and binval (values)
  // ---------------------------------------------------------------------------
  hsize_t size = nbins;
  hid_t fspc = H5Screate_simple(1, &size, NULL);
  hid_t dsetval = H5Dcreate(hgrp, "binval", H5T_NATIVE_DOUBLE, fspc,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t dsetbin = H5Dcreate(hgrp, "binloc", H5T_NATIVE_DOUBLE, fspc,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Dwrite(dsetbin, H5T_NATIVE_DOUBLE, fspc, fspc, H5P_DEFAULT, binloc);
  H5Dwrite(dsetval, H5T_NATIVE_DOUBLE, fspc, fspc, H5P_DEFAULT, binval);

  H5Dclose(dsetval);
  H5Dclose(dsetbin);
  H5Sclose(fspc);
  H5Gclose(hgrp);

  delete [] binloc;
  delete [] binval;
}

void Histogram2d::dump_hdf5(int base)
{
  if (base == 0) return;

  // Create a group to represent this histogram, and an attribute to name it
  // ---------------------------------------------------------------------------
  hid_t hgrp = H5Gcreate(base, nickname.c_str(),
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (!fullname.empty()) {

    hid_t aspc = H5Screate(H5S_SCALAR);
    hid_t strn = H5Tcopy(H5T_C_S1);
    H5Tset_size(strn, fullname.length());
    hid_t attr = H5Acreate(hgrp, "fullname", strn, aspc,
                           H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr, strn, fullname.c_str()); // write the full name

    H5Aclose(attr);
    H5Tclose(strn);
    H5Sclose(aspc);
  }

  // Create a temporary array for the bin centers
  // ---------------------------------------------------------------------------
  double *binlocX = new double[nbinsX];
  double *binlocY = new double[nbinsY];
  double *binval  = new double[nbins];
  for (int i=0; i<nbinsX; ++i) {
    binlocX[i] = 0.5*(bedgesX[i] + bedgesX[i+1]);
  }
  for (int j=0; j<nbinsY; ++j) {
    binlocY[j] = 0.5*(bedgesY[j] + bedgesY[j+1]);
  }
  for (int i=0; i<nbinsX; ++i) {
    for (int j=0; j<nbinsY; ++j) {
      const long c = counts[i*nbinsY + j];

      switch (binning_mode) {
      case Histogram::BinAverage:
	binval[i*nbinsY + j] = (c == 0) ? 0.0 : weight[i*nbinsY + j] / c;
	break;
      case Histogram::BinDensity:
	binval[i*nbinsY + j] = weight[i*nbinsY + j] /
	  ((bedgesX[i+1] - bedgesX[i])*(bedgesY[j+1] - bedgesY[j]));
	break;
      default:
	binval[i*nbinsY + j] = 0.0;
	break;
      }
    }
  }

  // Create the data sets in the group: binloc (bin centers) and binval (values)
  // ---------------------------------------------------------------------------
  hsize_t sizeX[2] = { nbinsX };
  hsize_t sizeY[2] = { nbinsY };
  hsize_t sizeZ[2] = { nbinsX, nbinsY };
  hid_t fspcX = H5Screate_simple(1, sizeX, NULL);
  hid_t fspcY = H5Screate_simple(1, sizeY, NULL);
  hid_t fspcZ = H5Screate_simple(2, sizeZ, NULL);

  hid_t dsetbinX = H5Dcreate(hgrp, "binlocX", H5T_NATIVE_DOUBLE, fspcX,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t dsetbinY = H5Dcreate(hgrp, "binlocY", H5T_NATIVE_DOUBLE, fspcY,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t dsetval = H5Dcreate(hgrp, "binval", H5T_NATIVE_DOUBLE, fspcZ,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Dwrite(dsetbinX, H5T_NATIVE_DOUBLE, fspcX, fspcX, H5P_DEFAULT, binlocX);
  H5Dwrite(dsetbinY, H5T_NATIVE_DOUBLE, fspcY, fspcY, H5P_DEFAULT, binlocY);
  H5Dwrite(dsetval , H5T_NATIVE_DOUBLE, fspcZ, fspcZ, H5P_DEFAULT, binval);

  H5Dclose(dsetval);
  H5Dclose(dsetbinX);
  H5Dclose(dsetbinY);
  H5Sclose(fspcX);
  H5Sclose(fspcY);
  H5Sclose(fspcZ);
  H5Gclose(hgrp);

  delete [] binlocX;
  delete [] binlocY;
  delete [] binval;
}


#else
void Histogram1d::dump_hdf5(int hid) { }
void Histogram2d::dump_hdf5(int hid) { }
#endif // (__MARA_USE_HDF5)



#if (__MARA_USE_MPI)
#include <mpi.h>

void Histogram1d::synchronize()
{
  int run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  if (!run_uses_mpi) {
    return;
  }
  MPI_Allreduce(MPI_IN_PLACE, weight, nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, counts, nbins, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
}
void Histogram2d::synchronize()
{
  int run_uses_mpi;
  MPI_Initialized(&run_uses_mpi);
  if (!run_uses_mpi) {
    return;
  }
  MPI_Allreduce(MPI_IN_PLACE, weight, nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, counts, nbins, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
}

#else
void Histogram1d::synchronize() { }
void Histogram2d::synchronize() { }
#endif // (__MARA_USE_MPI)
