/* parallel FFT functions - 1998, 1999

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs
   (505) 845-7873
   sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level directory of the distribution.
*/

/* User-settable FFT precision */

/* FFT_PRECISION = 1 is single-precision complex (4-byte real, 4-byte imag) */
/* FFT_PRECISION = 2 is double-precision complex (8-byte real, 8-byte imag) */

#ifndef FFT_3D_HEADER

#define FFT_PRECISION 2
#define FFT_FFTW3

#include "fftw3.h"
typedef fftw_complex FFT_DATA;

/* details of how to do a 3d FFT */

struct fft_plan_3d {
  struct remap_plan_3d *pre_plan;       /* remap from input -> 1st FFTs */
  struct remap_plan_3d *mid1_plan;      /* remap from 1st -> 2nd FFTs */
  struct remap_plan_3d *mid2_plan;      /* remap from 2nd -> 3rd FFTs */
  struct remap_plan_3d *post_plan;      /* remap from 3rd FFTs -> output */
  FFT_DATA *copy;                   /* memory for remap results (if needed) */
  FFT_DATA *scratch;                /* scratch space for remaps */
  int total1,total2,total3;         /* # of 1st,2nd,3rd FFTs (times length) */
  int length1,length2,length3;      /* length of 1st,2nd,3rd FFTs */
  int pre_target;                   /* where to put remap results */
  int mid1_target,mid2_target;
  int scaled;                       /* whether to scale FFT results */
  int normnum;                      /* # of values to rescale */
  double norm;                      /* normalization factor for rescaling */
};

/* function prototypes */

void fft_3d(FFT_DATA *, FFT_DATA *, int, struct fft_plan_3d *);
struct fft_plan_3d *fft_3d_create_plan
(MPI_Comm, int, int, int,
 int, int, int, int, int, int, int, int, int, int, int, int,
 int, int, int *);
void fft_3d_destroy_plan(struct fft_plan_3d *);
void factor(int, int *, int *);
void bifactor(int, int *, int *);

#endif // FFT_3D_HEADER
