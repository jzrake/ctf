
#include "config.h"
#if (__MARA_USE_MPI)

#include "subgrids.hpp"
#include "valman.hpp"

typedef DecomposedCartesianDomain Domain;

Domain::DecomposedCartesianDomain(const double *x0, const double *x1, const int *N,
				  int Nd, int Nq, int Ng)
  : Ng(Ng), Nq(Nq), num_dims(Nd)
{
  const int *dims_request = NULL; // not used presently

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
  MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Type_contiguous(Nq, MPI_DOUBLE, &mpi_type);

  int dims[3] = { 0,0,0 };
  int wrap[3] = { 1,1,1 };
  int reorder = 0; // Change this to 1 soon!

  if (dims_request == NULL) {
    MPI_Dims_create(mpi_size, num_dims, dims);
  }
  else {
    dims[0] = dims_request[0];
    dims[1] = dims_request[1];
    dims[2] = dims_request[2];
  }
  if (Nd < 2) dims[1] = 1;
  if (Nd < 3) dims[2] = 1;

  std::copy(dims, dims + num_dims, mpi_sizes);

  if (dims[0]*dims[1]*dims[2] != mpi_size) {
    printf("[subgrids] warning: requested subgrid dimensions not equal to mpi size.\n");
  }

  MPI_Cart_create(MPI_COMM_WORLD, num_dims, dims, wrap, reorder, &mpi_cart);
  MPI_Comm_rank(mpi_cart, &crt_rank);
  MPI_Cart_coords(mpi_cart, crt_rank, num_dims, mpi_index);

  printf("[subgrids] layout is (%d %d %d)\n", dims[0], dims[1], dims[2]);

  // ---------------------------------------------------------------------------

  for (int i=0; i < num_dims; ++i) {
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

  if (Ng == 0) return; // No need for boundary synchronization if no ghost zones

  switch(num_dims) {
  case(1):
    create_cart_1d();
    break;
  case(2):
    create_cart_2d();
    break;
  case(3):
    create_cart_3d();
    break;
  }
}
Domain::~DecomposedCartesianDomain()
{
  for (size_t i=0; i<send_type.size(); ++i) MPI_Type_free(&send_type[i]);
  for (size_t i=0; i<recv_type.size(); ++i) MPI_Type_free(&recv_type[i]);
  MPI_Type_free(&mpi_type);
  MPI_Comm_free(&mpi_cart);

  free(Specs.A_nint);
  free(Specs.L_ntot);
  free(Specs.L_strt);
  free(Specs.G_ntot);
  free(Specs.G_strt);
}

void Domain::Synchronize(std::valarray<double> &A) const
{
  std::vector<MPI_Request> requests;
  std::vector<MPI_Status> statuses;
  for (size_t i=0; i < neighbors.size(); ++i) {
    MPI_Request req1, req2;

    MPI_Isend(&A[0], 1, send_type[i], neighbors[i], send_tags[i], mpi_cart, &req1);
    MPI_Irecv(&A[0], 1, recv_type[i], neighbors[i], recv_tags[i], mpi_cart, &req2);
    requests.push_back(req1);
    requests.push_back(req2);
  }
  statuses.resize(requests.size());
  MPI_Waitall(requests.size(), &requests[0], &statuses[0]);
}
int Domain::SubgridRank() const
{
  return crt_rank;
}
int Domain::SubgridSize() const
{
  return mpi_size;
}
std::vector<int> Domain::aug_shape() const // including guard
{
  std::vector<int> shape(loc_shape, loc_shape+num_dims);
  for (int i=0; i<num_dims; ++i) {
    shape[i] += 2*Ng;
  }
  return shape;
}

void Domain::create_cart_1d()
{
  for (int i=-1; i<=1; ++i) {
    if (!i) continue; // if neighbor is me, don't do anything

    int nint[1] = { loc_shape[0]      };
    int size[1] = { loc_shape[0]+2*Ng };

    int Plx[3] = { Ng,Ng,nint[0] };

    int Qlx[3] = { 0,Ng,nint[0]+Ng };

    int rel_index [1] = { i };
    int start_send[1] = { Plx[i+1] };
    int start_recv[1] = { Qlx[i+1] };

    int subsize[1] = { (1-abs(i))*nint[0] + abs(i)*Ng };
    int index[1];
    for (int d=0; d<1; ++d) {
      index[d] = mpi_index[d] + rel_index[d];
    }
    int their_rank;
    MPI_Cart_rank(mpi_cart, index, &their_rank);
    neighbors.push_back(their_rank);
    send_tags.push_back(1*(+i+5));
    recv_tags.push_back(1*(-i+5));

    MPI_Datatype send, recv;
    MPI_Type_create_subarray(1, size, subsize, start_send, MPI_ORDER_C, mpi_type, &send);
    MPI_Type_create_subarray(1, size, subsize, start_recv, MPI_ORDER_C, mpi_type, &recv);
    MPI_Type_commit(&send);
    MPI_Type_commit(&recv);

    send_type.push_back(send);
    recv_type.push_back(recv);
  }
}
void Domain::create_cart_2d()
{
  for (int i=-1; i<=1; ++i) {
    for (int j=-1; j<=1; ++j) {
      if (!i && !j) continue; // if neighbor is me, don't do anything

      int nint[2] = { loc_shape[0]     , loc_shape[1]      };
      int size[2] = { loc_shape[0]+2*Ng, loc_shape[1]+2*Ng };

      int Plx[3] = { Ng,Ng,nint[0] };
      int Ply[3] = { Ng,Ng,nint[1] };

      int Qlx[3] = { 0,Ng,nint[0]+Ng };
      int Qly[3] = { 0,Ng,nint[1]+Ng };

      int rel_index [2] = { i, j };
      int start_send[2] = { Plx[i+1], Ply[j+1] };
      int start_recv[2] = { Qlx[i+1], Qly[j+1] };

      int subsize[2] = { (1-abs(i))*nint[0] + abs(i)*Ng,
                         (1-abs(j))*nint[1] + abs(j)*Ng };
      int index[2];
      for (int d=0; d<2; ++d) {
        index[d] = mpi_index[d] + rel_index[d];
      }
      int their_rank;
      MPI_Cart_rank(mpi_cart, index, &their_rank);
      neighbors.push_back(their_rank);
      send_tags.push_back(10*(+i+5) + 1*(+j+5));
      recv_tags.push_back(10*(-i+5) + 1*(-j+5));

      MPI_Datatype send, recv;
      MPI_Type_create_subarray(2, size, subsize, start_send, MPI_ORDER_C, mpi_type, &send);
      MPI_Type_create_subarray(2, size, subsize, start_recv, MPI_ORDER_C, mpi_type, &recv);
      MPI_Type_commit(&send);
      MPI_Type_commit(&recv);

      send_type.push_back(send);
      recv_type.push_back(recv);
    }
  }
}
void Domain::create_cart_3d()
{
  for (int i=-1; i<=1; ++i) {
    for (int j=-1; j<=1; ++j) {
      for (int k=-1; k<=1; ++k) {
        if (!i && !j && !k) continue; // if neighbor is me, don't do anything

        int nint[3] = { loc_shape[0]     , loc_shape[1]     , loc_shape[2]      };
        int size[3] = { loc_shape[0]+2*Ng, loc_shape[1]+2*Ng, loc_shape[2]+2*Ng };

        int Plx[3] = { Ng,Ng,nint[0] };
        int Ply[3] = { Ng,Ng,nint[1] };
        int Plz[3] = { Ng,Ng,nint[2] };

        int Qlx[3] = { 0,Ng,nint[0]+Ng };
        int Qly[3] = { 0,Ng,nint[1]+Ng };
        int Qlz[3] = { 0,Ng,nint[2]+Ng };

        int rel_index [3] = { i, j, k };
        int start_send[3] = { Plx[i+1], Ply[j+1], Plz[k+1] };
        int start_recv[3] = { Qlx[i+1], Qly[j+1], Qlz[k+1] };

        int subsize[3] = { (1-abs(i))*nint[0] + abs(i)*Ng,
                           (1-abs(j))*nint[1] + abs(j)*Ng,
                           (1-abs(k))*nint[2] + abs(k)*Ng };
        int index[3];
        for (int d=0; d<3; ++d) {
          index[d] = mpi_index[d] + rel_index[d];
        }

        int their_rank;
        MPI_Cart_rank(mpi_cart, index, &their_rank);
        neighbors.push_back(their_rank);
        send_tags.push_back(100*(+i+5) + 10*(+j+5) + 1*(+k+5));
        recv_tags.push_back(100*(-i+5) + 10*(-j+5) + 1*(-k+5));


        MPI_Datatype send, recv;
        MPI_Type_create_subarray(3, size, subsize, start_send, MPI_ORDER_C, mpi_type, &send);
        MPI_Type_create_subarray(3, size, subsize, start_recv, MPI_ORDER_C, mpi_type, &recv);
        MPI_Type_commit(&send);
        MPI_Type_commit(&recv);

        send_type.push_back(send);
        recv_type.push_back(recv);
      }
    }
  }
}
int Domain::SubgridAtPosition(const double *r) const
{
  int index[3];

  for (int i=0; i<num_dims; ++i) {
    index[i] = mpi_sizes[i] * (r[i] - glb_x0[i]) / (glb_x1[i] -  glb_x0[i]);
  }
  int their_rank;
  MPI_Cart_rank(mpi_cart, index, &their_rank);

  return their_rank;
}
int Domain::IndexAtPosition(const double *r, int d) const
{
  // d: 0,1,2 for x,y,z

  // r is a global position, i.e. not relative to this subgrid
  // The return value is the integer index, which is relative to this subgrid

  // Within the zone i+Ng, the value (x-x0)/dx ranges from i to i+1.
  // Then we correct for ghosts by adding Ng.

  return Ng + int((r[d] - loc_x0[d]) / dx[d]);
}

#else
// To stop ar from complaining that no symbols exist
void __subgrids_stub_object() { }
#endif
