
/*------------------------------------------------------------------------------
 * FILE: matrix.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION:
 *
 * - Calculates the matrix-matrix and matrix-vector product
 * - Calculates the inverse of a square matrix (adapted from Mike Dinolfo)
 *
 * A    : Pointer to input square double matrix
 * B    : Pointer to output memory space with size of A
 * rank : The number of rows/columns
 *
 * NOTES:
 * - the matrix must be invertible
 * - there's no pivoting of rows or columns, accuracy might not be adequate for
 *   your needs
 *
 * REFERENCES:
 *
 * http://en.wikipedia.org/wiki/LU_decomposition
 *------------------------------------------------------------------------------
 */
#include "matrix.h"
#define SMALL 1e-12

int matrix_matrix_product(const double *A, const double *B, double *C,
			  int ni, int nj, int nk)
{
  /* ---------------------------------------------------------------------------
   * Multiplies matrices C = A*B according to
   *
   * C_ij = Sum_k{ A_ik * B_kj }
   *
   * A is (ni x nk)
   * B is (nk x nj)
   * C is (ni x nj)
   * ---------------------------------------------------------------------------
   */
  int i,j,k;
  for (i=0; i<ni; ++i) {
    for (j=0; j<nj; ++j) {
      C[i*nj+j] = 0.0;
      for (k=0; k<nk; ++k) {
	C[i*nj+j] += A[i*nk+k] * B[k*nj+j];
      }
    }
  }
  return 0;
}

int matrix_vector_product(const double *A, const double *x, double *b,
			  int ni, int nj)
{
  /* ---------------------------------------------------------------------------
   * Multiplies matrix to vector b = A*x according to
   *
   * b_i = Sum_j{ A_ij * x_j }
   *
   * A is (ni x nj)
   * x is (nj x 1)
   * b is (ni x 1)
   * ---------------------------------------------------------------------------
   */
  int i,j;
  for (i=0; i<ni; ++i) {
    b[i] = 0.0;
    for (j=0; j<nj; ++j) {
      b[i] += A[i*nj+j] * x[j];
    }
  }
  return 0;
}

int vector_matrix_product(const double *A, const double *x, double *b,
			  int ni, int nj)
{
  /* ---------------------------------------------------------------------------
   * Multiplies matrix to vector b = x*A according to
   *
   * b_j = Sum_i{ x_i * A_ij }
   *
   * A is (ni x nj)
   * x is (1  x nj)
   * b is (1  x ni)
   * ---------------------------------------------------------------------------
   */
  int i,j;
  for (j=0; j<nj; ++j) {
    b[j] = 0.0;
    for (i=0; i<ni; ++i) {
      b[j] += x[i] * A[i*nj+j];
    }
  }
  return 0;
}


double matrix_inverse(const double *A, double *B, int rank)
{
  int i, j, k;
  double det, sum, x;

  if (rank < 1) { // matrix size must be larger than zero
    return 1;
  }
  else if (rank == 1) {
    det = A[0];
    B[0] = 1.0 / A[0];
    return det;
  }
  else if (rank == 2) {
    det = A[0*rank+0]*A[1*rank+1] - A[1*rank+0]*A[0*rank+1];
    B[0*rank+0] =  A[1*rank+1] / det;
    B[1*rank+1] =  A[0*rank+0] / det;
    B[0*rank+1] = -A[0*rank+1] / det;
    B[1*rank+0] = -A[1*rank+0] / det;
    return det;
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

  det = 1.0; // compute the determinant, product of diag(U)
  for (i=0; i<rank; ++i) {
    det *= B[i*rank+i];
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
  return det;
}
