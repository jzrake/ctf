
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
  const int Nd = domain.get_Nd();
  const int Nq = domain.get_Nq();

  const double *gx0 = domain.GetGlobalX0();
  const double *gx1 = domain.GetGlobalX1();

  const double Lx = (Nd>=1) ? gx1[0] - gx0[0] : 0.0;
  const double Ly = (Nd>=2) ? gx1[1] - gx0[1] : 0.0;
  const double Lz = (Nd>=3) ? gx1[2] - gx0[2] : 0.0;

  double r1[3];
  memcpy(r1, r0, Nd*sizeof(double));

  // Ensure that the target point is in the domain.
  // -------------------------------------------------------------------------
  if (Nd>=1) if (r1[0] > gx1[0]) r1[0] -= Lx;
  if (Nd>=2) if (r1[1] > gx1[1]) r1[1] -= Ly;
  if (Nd>=3) if (r1[2] > gx1[2]) r1[2] -= Lz;

  if (Nd>=1) if (r1[0] < gx0[0]) r1[0] += Lx;
  if (Nd>=2) if (r1[1] < gx0[1]) r1[1] += Ly;
  if (Nd>=3) if (r1[2] < gx0[2]) r1[2] += Lz;


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
  MPI_Isend(r1, Nd, MPI_DOUBLE, dest_rank, tag, comm, &request);

  ValarrayManager M(domain.aug_shape(), Nq);
  std::valarray<double> Panswer(Nq);

  for (int n=0; n<queries_to_answer; ++n) {

    double r_query[3];
    MPI_Recv(r_query, Nd, MPI_DOUBLE, MPI_ANY_SOURCE, tag, comm, &status);

    const int i = (Nd>=1) ? domain.IndexAtPosition(r_query, 0) : 0;
    const int j = (Nd>=2) ? domain.IndexAtPosition(r_query, 1) : 0;
    const int k = (Nd>=3) ? domain.IndexAtPosition(r_query, 2) : 0;

    if (Nd == 1) Panswer = HydroModule::Mara->PrimitiveArray[M(i)];
    if (Nd == 2) Panswer = HydroModule::Mara->PrimitiveArray[M(i,j)];
    if (Nd == 3) Panswer = HydroModule::Mara->PrimitiveArray[M(i,j,k)];

    MPI_Isend(&Panswer[0], Nq, MPI_DOUBLE, status.MPI_SOURCE, tag+1, comm,
	      &Prequest[n]);
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
