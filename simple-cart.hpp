
#ifndef __SimpleCartesianDomain_HEADER__
#define __SimpleCartesianDomain_HEADER__


#include "hydro.hpp"

class SimpleCartesianDomain : public PhysicalDomain
{
private:
  int mpi_size, old_rank, crt_rank;
  int mpi_index[3], mpi_sizes[3];

  double glb_x0[3], glb_x1[3]; // glb upper & lower domain limits
  double loc_x0[3], loc_x1[3]; // loc upper & lower domain limits
  double dx[3], min_dx;

  int Ng, Nq, num_dims; // ghost zones, fluid components, # of dimensions
  int glb_start[3];     // starting index in global array
  int glb_shape[3];     // glb number of zones (not including guard)
  int loc_shape[3];     // loc number of zones (not including guard)
  int ttl_zones;

public:
  SimpleCartesianDomain(const double *x0, const double *x1, const int *N,
			int Nd, int Nq, int Ng);
  ~SimpleCartesianDomain();
  void Synchronize(std::valarray<double> &A) const;
  int SubgridRank() const;
  int SubgridSize() const;
  int GetNumberOfZones() const { return ttl_zones; }
  int GetGlobalNumberOfZones() const { return glb_shape[0]*glb_shape[1]*glb_shape[2]; }
  int GetSubgridIndex(int i) const { return mpi_index[i]; }
  int GetSubgridSizes(int i) const { return mpi_sizes[i]; }

  double get_dx(int d) const { return dx[d-1]         ; }
  double get_min_dx()  const { return min_dx          ; }
  int get_N(int d)     const { return loc_shape[d-1]  ; } // not including guard
  int get_Ng()         const { return Ng              ; }
  int get_Nq()         const { return Nq              ; }
  int get_Nd()         const { return num_dims        ; } // number of dimensions
  double x_at(int i) const { return (num_dims > 0) ? loc_x0[0] + dx[0]*(i-Ng+0.5) : 0.0; }
  double y_at(int j) const { return (num_dims > 1) ? loc_x0[1] + dx[1]*(j-Ng+0.5) : 0.0; }
  double z_at(int k) const { return (num_dims > 2) ? loc_x0[2] + dx[2]*(k-Ng+0.5) : 0.0; }

  const int *GetGlobalShape() const { return glb_shape; }
  const int *GetGlobalStart() const { return glb_start; }
  const int *GetLocalShape() const { return loc_shape; } // not including guard
  const double *GetGlobalX0() const { return glb_x0; }
  const double *GetGlobalX1() const { return glb_x1; }
  std::vector<int> aug_shape() const; // including guard
  std::vector<double> get_x0() const { return std::vector<double>(loc_x0, loc_x0+num_dims); }
  std::vector<double> get_x1() const { return std::vector<double>(loc_x1, loc_x1+num_dims); }


  int SubgridAtPosition(const double *r) const;
  int IndexAtPosition(const double *r, int d) const; // d: 0,1,2 for x,y,z
} ;

#endif // __SimpleCartesianDomain_HEADER__
