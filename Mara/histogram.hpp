


#ifndef __Histogram_HEADER__
#define __Histogram_HEADER__

#include <string>

namespace Histogram {
  enum SpacingType { Linspace, Logspace };
  enum BinningMode { BinAverage, BinDensity };
} ;

class Histogram1d
{
private:
  int nbins;
  double *bedges;
  double *weight;
  long *counts;

public:
  std::string nickname;
  std::string fullname;
  enum Histogram::BinningMode binning_mode;

  Histogram1d(int nbins, double x0, double x1, enum Histogram::SpacingType spc=Histogram::Linspace);
  ~Histogram1d();

  int add_sample(double x, double w);
  void dump_ascii(FILE *file);
  void dump_hdf5(int hid);
  void synchronize();
} ;

class Histogram2d
{
private:
  int nbinsX, nbinsY, nbins;
  double *bedgesX, *bedgesY;
  double *weight;
  long *counts;

public:
  std::string nickname;
  std::string fullname;
  enum Histogram::BinningMode binning_mode;

  Histogram2d(int nbinsX, int nbinsY, double x0, double x1,
	      double y0, double y1, enum Histogram::SpacingType spc=Histogram::Linspace);
  ~Histogram2d();

  int add_sample(double x, double y, double w);
  void dump_ascii(FILE *file);
  void dump_hdf5(int hid);
  void synchronize();
} ;


#endif // __Histogram_HEADER__
