
/*------------------------------------------------------------------------------
 * FILE: matinv.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP (adapted from C++ code by Mike Dinolfo)
 *
 * DESCRIPTION: Calculates the inverse of a square matrix

 A    : Pointer to input square double matrix
 B    : Pointer to output memory space with size of A
 rank : The number of rows/columns

 * NOTES:
 - the matrix must be invertible
 - there's no pivoting of rows or columns,
 -         accuracy might not be adequate for your needs.

 * REFERENCES:

 http://en.wikipedia.org/wiki/LU_decomposition
 *------------------------------------------------------------------------------
 */

#include "matrix.h"

static const double SMALL = 1e-12;
static double Determinant = 0.0;

double matrix_determinant()
{
  return Determinant;
}

int matrix_inverse(const double *A, double *B, int rank)
{
  int i,j,k;
  double sum,x;

  if (rank < 1) { // matrix size must be larger than zero
    return 1;
  }
  else if (rank == 1) {
    Determinant = A[0];
    B[0] = 1.0 / A[0];
    return 0;
  }
  else if (rank == 2) {
    Determinant = A[0*rank+0]*A[1*rank+1] - A[1*rank+0]*A[0*rank+1];

    B[0*rank+0] =  A[1*rank+1] / Determinant;
    B[1*rank+1] =  A[0*rank+0] / Determinant;
    B[0*rank+1] = -A[0*rank+1] / Determinant;
    B[1*rank+0] = -A[1*rank+0] / Determinant;

    return 0;
  }

  for (i=0; i<rank*rank; ++i) { // copy the input matrix to output matrix
    B[i] = A[i];
  }

  for (i=0; i<rank; ++i) { // add small value to diagonal if diagonal is zero
    j = i*rank + i;
    if (B[j] < SMALL && B[j] > -SMALL) {
      B[j] = SMALL;
    }
  }

  for (i=1; i<rank; ++i) { // normalize row 0
    B[i] /= B[0];
  }

  for (i=1; i<rank; ++i) {
    for (j=i; j<rank; ++j) { // do a column of L
      sum = 0.0;
      for (k=0; k<i; ++k) {
        sum += B[j*rank+k] * B[k*rank+i];
      }
      B[j*rank+i] -= sum;
    }
    if (i == rank-1) continue;
    for (j=i+1; j<rank; ++j) { // do a row of U
      sum = 0.0;
      for (k=0; k<i; ++k) {
        sum += B[i*rank+k]*B[k*rank+j];
      }
      B[i*rank+j] = (B[i*rank+j]-sum) / B[i*rank+i];
    }
  }

  Determinant = 1.0; // compute the determinant, product of diag(U)
  for (i=0; i<rank; ++i) {
    Determinant *= B[i*rank+i];
  }

  for (i=0; i<rank; ++i) { // invert L
    for (j=i; j<rank; ++j) {
      x = 1.0;
      if (i != j) {
        x = 0.0;
        for (k=i; k<j; ++k) {
          x -= B[j*rank+k]*B[k*rank+i];
        }
      }
      B[j*rank+i] = x / B[j*rank+j];
    }
  }

  for (i=0; i<rank; ++i) { // invert U
    for (j=i; j<rank; ++j) {
      if (i == j) continue;
      sum = 0.0;
      for (k=i; k<j; ++k) {
        sum += B[k*rank+j]*((i==k) ? 1.0 : B[i*rank+k]);
      }
      B[i*rank+j] = -sum;
    }
  }

  for (i=0; i<rank; ++i) { // final inversion
    for (j=0; j<rank; ++j) {
      sum = 0.0;
      for (k=((i>j) ? i : j); k<rank; ++k) {
        sum += ((j==k) ? 1.0 : B[j*rank+k])*B[k*rank+i];
      }
      B[j*rank+i] = sum;
    }
  }

  return 0;
}
