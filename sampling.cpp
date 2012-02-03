
#include "config.h"
#if (__MARA_USE_MPI)

#include <vector>
#include <mpi.h>

#include "valman.hpp"
#include "hydro.hpp"
#include "mara_mpi.h"



void Mara_prim_at_point(const double *r0, double *P1)
{
  const PhysicalDomain &domain = *HydroModule::Mara->domain;

  const int rank = domain.SubgridRank();
  const int size = domain.SubgridSize();
  const double Lx = domain.GetGlobalX1()[0] - domain.GetGlobalX0()[0];
  const double Ly = domain.GetGlobalX1()[1] - domain.GetGlobalX0()[1];
  const double Lz = domain.GetGlobalX1()[2] - domain.GetGlobalX0()[2];

  const int Nq = domain.get_Nq();
  const int Ng = domain.get_Ng();
  const int Ny = domain.GetLocalShape()[1];
  const int Nz = domain.GetLocalShape()[2];

  double r1[3];
  memcpy(r1, r0, 3*sizeof(double));

  // Ensure that the target point is in the domain.
  // -------------------------------------------------------------------------
  if (r1[0] > domain.GetGlobalX1()[0]) r1[0] -= Lx;
  if (r1[1] > domain.GetGlobalX1()[1]) r1[1] -= Ly;
  if (r1[2] > domain.GetGlobalX1()[2]) r1[2] -= Lz;

  if (r1[0] < domain.GetGlobalX0()[0]) r1[0] += Lx;
  if (r1[1] < domain.GetGlobalX0()[1]) r1[1] += Ly;
  if (r1[2] < domain.GetGlobalX0()[2]) r1[2] += Lz;


  // Determine on which process that point resides, set up some MPI variables.
  // -------------------------------------------------------------------------
  const int tag = 123;
  const int dest_rank = domain.SubgridAtPosition(r1);
  const MPI_Comm comm = MPI_COMM_WORLD;


  // Figure out how many of other people's points lie in my domain.
  // -------------------------------------------------------------------------
  int *number_of_queries = new int[size];
  for (int i=0; i<size; ++i) { // ensures that all entries are initialized
    number_of_queries[i] = (i == dest_rank);
  }
  MPI_Allreduce(MPI_IN_PLACE, number_of_queries, size, MPI_INT, MPI_SUM, comm);
  const int queries_to_answer = number_of_queries[rank];
  delete [] number_of_queries;


  // Send out a query for the target position, non-blocking.
  // -------------------------------------------------------------------------
  MPI_Request request, *Prequest = new MPI_Request[queries_to_answer];
  MPI_Status status, *Pstatus = new MPI_Status[queries_to_answer];
  MPI_Isend(r1, 3, MPI_DOUBLE, dest_rank, tag, comm, &request);

  for (int n=0; n<queries_to_answer; ++n) {

    double r_query[3];
    MPI_Recv(r_query, 3, MPI_DOUBLE, MPI_ANY_SOURCE, tag, comm, &status);

    const int i = domain.IndexAtPosition(r_query, 0);
    const int j = domain.IndexAtPosition(r_query, 1);
    const int k = domain.IndexAtPosition(r_query, 2);

    const int mm = (i*(Ny+2*Ng)*(Nz+2*Ng) + j*(Nz+2*Ng) + k)*Nq;
    double *P_answer = &HydroModule::Mara->PrimitiveArray[mm];

    MPI_Isend(P_answer, Nq, MPI_DOUBLE, status.MPI_SOURCE, tag+1, comm, &Prequest[n]);
  }

  MPI_Recv(P1, Nq, MPI_DOUBLE, dest_rank, tag+1, comm, &status);
  MPI_Waitall(queries_to_answer, Prequest, Pstatus);
  MPI_Waitall(1, &request, &status);
  MPI_Barrier(comm);

  delete [] Prequest;
  delete [] Pstatus;
}




#else
void Mara_prim_at_point(const double *r0, double *P1) { }
#endif // __MARA_USE_MPI
