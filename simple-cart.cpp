
#include <cstdio>
#include <iostream>
#include "valman.hpp"
#include "simple-cart.hpp"

typedef SimpleCartesianDomain Domain;

Domain::SimpleCartesianDomain(const double *x0, const double *x1, const int *N,
			      int Nd, int Nq, int Ng)
  : Ng(Ng), Nq(Nq), num_dims(Nd)
{

  for (int d=0; d<3; ++d) {
    glb_shape[d] = (Nd >= d+1) ? N [d] : 1;
    glb_x0   [d] = (Nd >= d+1) ? x0[d] : 0.0;
    glb_x1   [d] = (Nd >= d+1) ? x1[d] : 1.0;
  }

  min_dx = glb_x1[0] - glb_x0[0];
  for (int i=0; i < num_dims; ++i) {
    dx[i] = (glb_x1[i] - glb_x0[i]) / glb_shape[i];
    if (dx[i] < min_dx) min_dx = dx[i];
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  old_rank = crt_rank = 0;
  mpi_index[0] = mpi_index[1] = mpi_index[2] = 0;
  mpi_size = mpi_sizes[0] = mpi_sizes[1] = mpi_sizes[2] = 1;
  // ---------------------------------------------------------------------------

  for (int i=0; i<num_dims; ++i) {
    const int R = glb_shape[i] % mpi_sizes[i];
    const int normal_size = glb_shape[i] / mpi_sizes[i];
    const int augmnt_size = normal_size + 1;
    const int thisdm_size = (mpi_index[i]<R) ? augmnt_size : normal_size;

    loc_shape[i] = thisdm_size;
    glb_start[i] = 0;
    for (int j=0; j<mpi_index[i]; ++j) {
      glb_start[i] += (j<R) ? augmnt_size : normal_size;
    }

    loc_x0[i] = glb_x0[i] + dx[i] *  glb_start[i];
    loc_x1[i] = glb_x0[i] + dx[i] * (glb_start[i] + thisdm_size);
  }

  ttl_zones = 1;
  for (int i=0; i < num_dims; ++i) {
    ttl_zones *= loc_shape[i] + 2*Ng;
  }

  // ---------------------------------------------------------------------------
  Specs.n_dims = num_dims;
  Specs.n_prim = Nq;
  Specs.A_nint = (int*) malloc(num_dims*sizeof(int));
  Specs.L_ntot = (int*) malloc(num_dims*sizeof(int));
  Specs.L_strt = (int*) malloc(num_dims*sizeof(int));
  Specs.G_ntot = (int*) malloc(num_dims*sizeof(int));
  Specs.G_strt = (int*) malloc(num_dims*sizeof(int));

  for (int i=0; i<num_dims; ++i) {
    Specs.A_nint[i] = loc_shape[i];
    Specs.L_ntot[i] = loc_shape[i] + 2*Ng;
    Specs.L_strt[i] = Ng;
    Specs.G_ntot[i] = glb_shape[i];
    Specs.G_strt[i] = glb_start[i];
  }
  // ---------------------------------------------------------------------------
}
Domain::~SimpleCartesianDomain()
{
  free(Specs.A_nint);
  free(Specs.L_ntot);
  free(Specs.L_strt);
  free(Specs.G_ntot);
  free(Specs.G_strt);
}

void Domain::Synchronize(std::valarray<double> &A) const
{
  const int &Nx = loc_shape[0];
  const int &Ny = loc_shape[1];
  const int &Nz = loc_shape[2];
  // ---------------------------------------------------------------------------
  // Periodic BC's implemented here
  ValarrayManager M(this->aug_shape(), Nq);
  switch (num_dims) {

  case 1: // ---------------------------------------------------------------- 1d
    for (int i=0; i<Ng; ++i) {
      A[ M(i) ] = A[ M(i+Nx) ]; // set_bc_x0_wall
    }
    for (int i=Nx+Ng; i<Nx+2*Ng; ++i) {
      A[ M(i) ] = A[ M(i-Nx) ]; // set_bc_x1_wall
    }
    break;
    // -------------------------------------------------------------------------

  case 2: // ---------------------------------------------------------------- 2d
    for (int i=0; i<Ng; ++i) {
      for (int j=0; j<Ny+2*Ng; ++j) {
        A[ M(i,j) ] = A[ M(i+Nx,j) ]; // set_bc_x0_wall
      }
    }
    for (int i=Nx+Ng; i<Nx+2*Ng; ++i) {
      for (int j=0; j<Ny+2*Ng; ++j) {
        A[ M(i,j) ] = A[ M(i-Nx,j) ]; // set_bc_x1_wall
      }
    }

    for (int i=0; i<Nx+2*Ng; ++i) {
      for (int j=0; j<Ng; ++j) {
        A[ M(i,j) ] = A[ M(i,j+Ny) ]; // set_bc_y0_wall
      }
    }
    for (int i=0; i<Nx+2*Ng; ++i) {
      for (int j=Ny+Ng; j<Ny+2*Ng; ++j) {
        A[ M(i,j) ] = A[ M(i,j-Ny) ]; // set_bc_y1_wall
      }
    }
    break;
    // -------------------------------------------------------------------------

  case 3: // ---------------------------------------------------------------- 3d
    for (int i=0; i<Ng; ++i) {
      for (int j=0; j<Ny+2*Ng; ++j) {
        for (int k=0; k<Nz+2*Ng; ++k) {
          A[ M(i,j,k) ] = A[ M(i+Nx,j,k) ]; // set_bc_x0_wall
        }
      }
    }
    for (int i=Nx+Ng; i<Nx+2*Ng; ++i) {
      for (int j=0; j<Ny+2*Ng; ++j) {
        for (int k=0; k<Nz+2*Ng; ++k) {
          A[ M(i,j,k) ] = A[ M(i-Nx,j,k) ]; // set_bc_x1_wall
        }
      }
    }

    for (int i=0; i<Nx+2*Ng; ++i) {
      for (int j=0; j<Ng; ++j) {
        for (int k=0; k<Nz+2*Ng; ++k) {
          A[ M(i,j,k) ] = A[ M(i,j+Ny,k) ]; // set_bc_y0_wall
        }
      }
    }
    for (int i=0; i<Nx+2*Ng; ++i) {
      for (int j=Ny+Ng; j<Ny+2*Ng; ++j) {
        for (int k=0; k<Nz+2*Ng; ++k) {
          A[ M(i,j,k) ] = A[ M(i,j-Ny,k) ]; // set_bc_y1_wall
        }
      }
    }

    for (int i=0; i<Nx+2*Ng; ++i) {
      for (int j=0; j<Ny+2*Ng; ++j) {
        for (int k=0; k<Ng; ++k) {
          A[ M(i,j,k) ] = A[ M(i,j,k+Nz) ]; // set_bc_z0_wall
        }
      }
    }
    for (int i=0; i<Nx+2*Ng; ++i) {
      for (int j=0; j<Ny+2*Ng; ++j) {
        for (int k=Nz+Ng; k<Nz+2*Ng; ++k) {
          A[ M(i,j,k) ] = A[ M(i,j,k-Nz) ]; // set_bc_z1_wall
        }
      }
    }
    break;
    // -------------------------------------------------------------------------
  }
}
int Domain::SubgridAtPosition(const double *r) const { return 0; }
int Domain::IndexAtPosition(const double *r, int d) const
{
  // d: 0,1,2 for x,y,z

  // r is a global position, i.e. not relative to this subgrid
  // The return value is the integer index, which is relative to this subgrid

  // Within the zone i+Ng, the value (x-x0)/dx ranges from i to i+1.
  // Then we correct for ghosts by adding Ng.

  return Ng + int((r[d] - loc_x0[d]) / dx[d]);
}

int Domain::SubgridRank() const { return 0; }
int Domain::SubgridSize() const { return 1; }

std::vector<int> Domain::aug_shape() const // including guard
{
  std::vector<int> shape(loc_shape, loc_shape+num_dims);
  for (int i=0; i<num_dims; ++i) {
    shape[i] += 2*Ng;
  }
  return shape;
}
