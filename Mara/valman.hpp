
#ifndef __ValarrayManager_HEADER__
#define __ValarrayManager_HEADER__

#include <functional>
#include <numeric>
#include <valarray>
#include <vector>

class ValarrayManager
{
public:
  ValarrayManager() { }
  ValarrayManager(std::vector<int> N, int Q) : N(N), Q(Q)
  {
    R = N.size();
    S = std::accumulate(N.begin(), N.end(), Q, std::multiplies<int>());
  }
  const std::vector<int> &get_N() const { return N; }

  int get_Q() const { return Q; }
  int get_R() const { return R; }
  int get_S() const { return S; }

  std::slice operator()(int i)
  {
    return std::slice(Q*(i), Q, 1);
  }
  std::slice operator()(int i, int j)
  {
    return std::slice(Q*(i*N[1] + j), Q, 1);
  }
  std::slice operator()(int i, int j, int k)
  {
    return std::slice(Q*(i*N[1]*N[2] + j*N[2] + k), Q, 1);
  }
  std::slice operator[](int q)
  {
    return std::slice(q, S/Q, Q);
  }

private:
  std::vector<int> N;
  int Q, R, S;
} ;

class ValarrayIndexer
{
public:
  ValarrayIndexer() { }
  ValarrayIndexer(std::vector<int> N_)
  {
    N = N_;
    R = N.size();
    S = std::accumulate(N.begin(), N.end(), 1, std::multiplies<int>());
  }
  const std::vector<int> &get_N() { return N; }

  int get_R() { return R; }
  int get_S() { return S; }

  int operator()(int i)
  {
    return i;
  }
  int operator()(int i, int j)
  {
    return i*N[1] + j;
  }
  int operator()(int i, int j, int k)
  {
    return i*N[1]*N[2] + j*N[2] + k;
  }

private:
  std::vector<int> N;
  int R, S;
} ;

#endif // __ValarrayManager_HEADER__
