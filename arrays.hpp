
#ifndef __ArrayTools_HEADER__
#define __ArrayTools_HEADER__

// i0: first index, which can be negative
// ni: number of elements in i-direction

template <class T> T *array_1d(int i0, int ni)
{
  T *p;

  p = new T[ni];
  p -= i0;
  return p;
}

template <class T> T **array_2d(int i0, int ni, int j0, int nj)
{
  T **pp;

  pp = new T*[ni];
  pp[0] = new T[ni*nj];

  pp -= i0;
  pp[i0] -= j0;

  for (int i=i0+1; i<ni+i0; ++i) {
    pp[i] = pp[i0] + (i-i0)*nj;
  }
  return pp;
}

template <class T> T ***array_3d(int i0, int ni, int j0, int nj, int k0, int nk)
{
  T ***ppp;

  ppp = new T**[ni];
  ppp[0] = new T*[ni*nj];
  ppp[0][0] = new T[ni*nj*nk];

  ppp -= i0;
  ppp[i0] -= j0;
  ppp[i0][j0] -= k0;

  for (int i=i0+1; i<ni+i0; ++i) {
    ppp[i] = ppp[i0] + (i-i0)*nj;
  }

  for (int i=i0; i<ni+i0; ++i) {
    for (int j=j0; j<nj+j0; ++j) {
      ppp[i][j] = ppp[i0][j0] + (i-i0)*nj*nk + (j-j0)*nk ;
    }
  }
  return ppp;
}


template <class T> T ****array_4d(int i0, int ni, int j0, int nj,
                                  int k0, int nk, int l0, int nl)
{
  T ****p4;
  p4 = new T***[ni];
  p4[0] = new T**[ni*nj];
  p4[0][0] = new T*[ni*nj*nk];
  p4[0][0][0] = new T[ni*nj*nk*nl];

  p4 -= i0;
  p4[i0] -= j0;
  p4[i0][j0] -= k0;
  p4[i0][j0][k0] -= l0;

  for (int i=i0+1; i<ni+i0; ++i) {
    p4[i] = p4[i0] + (i-i0)*nj;
  }

  for (int i=i0; i<ni+i0; ++i) {
    for (int j=j0; j<nj+j0; ++j) {
      p4[i][j] = p4[i0][j0] + (i-i0)*nj*nk + (j-j0)*nk;
    }
  }

  for (int i=i0; i<ni+i0; ++i) {
    for (int j=j0; j<nj+j0; ++j) {
      for (int k=k0; k<nk+k0; ++k) {
        p4[i][j][k] = p4[i0][j0][k0] + (i-i0)*nj*nk*nl + (j-j0)*nk*nl + (k-k0)*nl;
      }
    }
  }

  return p4;
}

template <class T> void delete_array_1d(T *p, int i0)
{
  p += i0;
  delete [] p;
}

template <class T> void delete_array_2d(T **pp, int i0, int j0)
{
  pp[i0] += j0;
  pp += i0;

  delete [] pp[0];
  delete [] pp;
}

template <class T> void delete_array_3d(T ***ppp, int i0, int j0, int k0)
{
  ppp[i0][j0] += k0;
  ppp[i0] += j0;
  ppp += i0;

  delete [] ppp[0][0];
  delete [] ppp[0];
  delete [] ppp;
}

template <class T> void delete_array_4d(T ****p4, int i0, int j0, int k0, int l0)
{
  p4[i0][j0][k0] += l0;
  p4[i0][j0] += k0;
  p4[i0] += j0;
  p4 += i0;

  delete [] p4[0][0][0];
  delete [] p4[0][0];
  delete [] p4[0];
  delete [] p4;
}

#endif // __ArrayTools_HEADER__
