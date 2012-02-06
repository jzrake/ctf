
#include "config.h"
#if (__MARA_USE_MPI)

#include <vector>
#include <cstring>
#include <mpi.h>

#include "valman.hpp"
#include "hydro.hpp"
#include "sampling.hpp"
#include "mara_mpi.h"


static std::valarray<double> unilinear_interp(const double *r);
static std::valarray<double>  bilinear_interp(const double *r);
static std::valarray<double> trilinear_interp(const double *r);


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


  // Since non-blocking sends are being used below, it is nesessary to manually
  // buffer the send buffer data. Each 'Panswer' will have its own place in this
  // vector which is unchanged throughout this function, and is thus guarenteed
  // to be in the correct place when the matching receive is posted.
  // ---------------------------------------------------------------------------
  std::valarray<double> *Pvector = new std::valarray<double>[queries_to_answer];


  for (int n=0; n<queries_to_answer; ++n) {

    double r_query[3];
    MPI_Recv(r_query, Nd, MPI_DOUBLE, MPI_ANY_SOURCE, tag, comm, &status);

    std::valarray<double> &Panswer = Pvector[n];
    Panswer.resize(Nq);


    if (Nd == 1) {
      Panswer = unilinear_interp(r_query);
    }
    else if (Nd == 2) {
      Panswer = bilinear_interp(r_query);
    }
    else if (Nd == 3) {
      Panswer = trilinear_interp(r_query);
    }

    MPI_Isend(&Panswer[0], Nq, MPI_DOUBLE, status.MPI_SOURCE, tag+1, comm,
              &Prequest[n]);
  }

  MPI_Recv(P1, Nq, MPI_DOUBLE, dest_rank, tag+1, comm, &status);
  MPI_Waitall(queries_to_answer, Prequest, Pstatus);
  MPI_Waitall(1, &request, &status);
  MPI_Barrier(comm);

  delete [] Prequest;
  delete [] Pstatus;
  delete [] Pvector;
}




void Mara_prim_at_point_many(const double *Rin, double *Rlist, double *Plist, int Nsamp)
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

  std::vector<double> *remote_r1 = new std::vector<double>[size];
  std::vector<double> *remote_P1 = new std::vector<double>[size];

  for (int m=0; m<Nsamp; ++m) {

    double r1[3];
    r1[0] = Rin[3*m + 0];
    r1[1] = Rin[3*m + 1];
    r1[2] = Rin[3*m + 2];

    // Ensure that the target point is in the domain.
    // -------------------------------------------------------------------------
    if (Nd>=1) if (r1[0] > gx1[0]) r1[0] -= Lx;
    if (Nd>=2) if (r1[1] > gx1[1]) r1[1] -= Ly;
    if (Nd>=3) if (r1[2] > gx1[2]) r1[2] -= Lz;

    if (Nd>=1) if (r1[0] < gx0[0]) r1[0] += Lx;
    if (Nd>=2) if (r1[1] < gx0[1]) r1[1] += Ly;
    if (Nd>=3) if (r1[2] < gx0[2]) r1[2] += Lz;

    const int remote = domain.SubgridAtPosition(r1);

    for (int d=0; d<3; ++d) {
      remote_r1[remote].push_back(r1[d]);
    }

    for (int q=0; q<Nq; ++q) {
      remote_P1[remote].push_back(0.0); // not known yet; will be looked up on remote
    }
  }


  // We will be sampling pairs between ourselves and the remote process, rank +
  // dn. That process is referred to as 'lawyer' because they will work for us,
  // obtaining the records we have determined live on their domain. Similarly,
  // we will be the lawyer for process rank-dn, so that process is called
  // 'client'.
  // ---------------------------------------------------------------------------
  MPI_Status status;
  MPI_Comm comm = MPI_COMM_WORLD;
  int queries_satisfied = 0;

  for (int dn=0; dn<size; ++dn) {
    // This loop contains three pairs of matching Send/Recv's. For the first
    // two, the send is place to the lawyer process, and the receive comes from
    // the client. We call our lawyer and let him know to expect a message of
    // length 'num_lawyer' double[3]'s. Those are the positions of the remote
    // points we have chosen, and then determined live on his domain. At the
    // same time, we receive a message from our client, which contains the
    // length 'num_client' of double[3]'s he will be asking us to fetch. The
    // next pair of Send/Recv's is the list of coordinates themselves. The last
    // one is the list of primitive quantities at those locations.
    // -------------------------------------------------------------------------

    const int lawyer = (rank + size + dn) % size;
    const int client = (rank + size - dn) % size;

    int num_lawyer = remote_r1[lawyer].size() / 3;
    int num_client;

    MPI_Sendrecv(&num_lawyer, 1, MPI_INT, lawyer, 123,
                 &num_client, 1, MPI_INT, client, 123, comm, &status);

    double *you_get_for_me_r = &remote_r1[lawyer][0];
    double *you_get_for_me_P = &remote_P1[lawyer][0];
    double *I_find_for_you_r = new double[num_client*3];
    double *I_find_for_you_P = new double[num_client*Nq];

    MPI_Sendrecv(you_get_for_me_r, num_lawyer*3, MPI_DOUBLE, lawyer, 123,
                 I_find_for_you_r, num_client*3, MPI_DOUBLE, client, 123,
                 comm, &status);

    for (int s=0; s<num_client; ++s) {

      const double *r_query = &I_find_for_you_r[3*s];
      std::valarray<double> Panswer(Nq);

      if (Nd == 1) {
	Panswer = unilinear_interp(r_query);
      }
      else if (Nd == 2) {
	Panswer = bilinear_interp(r_query);
      }
      else if (Nd == 3) {
	Panswer = trilinear_interp(r_query);
      }
      std::memcpy(&I_find_for_you_P[s*Nq], &Panswer[0], Nq*sizeof(double));
    }

    MPI_Sendrecv(I_find_for_you_P, num_client*Nq, MPI_DOUBLE, client, 123,
                 you_get_for_me_P, num_lawyer*Nq, MPI_DOUBLE, lawyer, 123,
                 comm, &status);

    std::memcpy(&Rlist[queries_satisfied* 3], &remote_r1[lawyer][0],
		remote_r1[lawyer].size()*sizeof(double));
    std::memcpy(&Plist[queries_satisfied*Nq], &remote_P1[lawyer][0],
		remote_P1[lawyer].size()*sizeof(double));

    queries_satisfied += num_lawyer;

    delete [] I_find_for_you_r;
    delete [] I_find_for_you_P;
  }

  delete [] remote_r1;
  delete [] remote_P1;
}







std::valarray<double> unilinear_interp(const double *r)
{
  const PhysicalDomain &domain = *HydroModule::Mara->domain;
  const int Nq = domain.get_Nq();
  ValarrayManager M(domain.aug_shape(), Nq);
  std::valarray<double> Panswer(Nq);

  const int i = domain.IndexAtPosition(r, 0);
  const std::valarray<double> &Prim = HydroModule::Mara->PrimitiveArray;

  std::valarray<double> P0 = Prim[M(i-1)];
  std::valarray<double> P1 = Prim[M(i+1)];

  double delx[1];
  delx[0] = 0.5 * (r[0] - domain.x_at(i-1)) / domain.get_dx(1);

  for (int q=0; q<Nq; ++q) {
    double b1 = P0[q];
    double b2 = P1[q] - P0[q];

    Panswer[q] = b1 + b2*delx[0];
  }
  return Panswer;
}

std::valarray<double> bilinear_interp(const double *r)
{
  const PhysicalDomain &domain = *HydroModule::Mara->domain;
  const int Nq = domain.get_Nq();
  ValarrayManager M(domain.aug_shape(), Nq);
  std::valarray<double> Panswer(Nq);

  const int i = domain.IndexAtPosition(r, 0);
  const int j = domain.IndexAtPosition(r, 1);
  const std::valarray<double> &Prim = HydroModule::Mara->PrimitiveArray;

  std::valarray<double> P00 = Prim[M(i-1,j-1)];
  std::valarray<double> P01 = Prim[M(i-1,j+1)];
  std::valarray<double> P10 = Prim[M(i+1,j-1)];
  std::valarray<double> P11 = Prim[M(i+1,j+1)];

  double delx[2];
  delx[0] = 0.5 * (r[0] - domain.x_at(i-1)) / domain.get_dx(1);
  delx[1] = 0.5 * (r[1] - domain.y_at(j-1)) / domain.get_dx(2);

  for (int q=0; q<Nq; ++q) {
    double b1 = P00[q];
    double b2 = P10[q] - P00[q];
    double b3 = P01[q] - P00[q];
    double b4 = P00[q] - P10[q] - P01[q] + P11[q];

    Panswer[q] = b1 + b2*delx[0] + b3*delx[1] + b4*delx[0]*delx[1];
  }
  return Panswer;
}

std::valarray<double> trilinear_interp(const double *r)
{
  const PhysicalDomain &domain = *HydroModule::Mara->domain;
  const int Nq = domain.get_Nq();
  ValarrayManager M(domain.aug_shape(), Nq);
  std::valarray<double> Panswer(Nq);

  const int i = domain.IndexAtPosition(r, 0);
  const int j = domain.IndexAtPosition(r, 1);
  const int k = domain.IndexAtPosition(r, 2);
  const std::valarray<double> &Prim = HydroModule::Mara->PrimitiveArray;

  std::valarray<double> P000 = Prim[M(i-1,j-1,k-1)];
  std::valarray<double> P001 = Prim[M(i-1,j-1,k+1)];
  std::valarray<double> P010 = Prim[M(i-1,j+1,k-1)];
  std::valarray<double> P011 = Prim[M(i-1,j+1,k+1)];
  std::valarray<double> P100 = Prim[M(i+1,j-1,k-1)];
  std::valarray<double> P101 = Prim[M(i+1,j-1,k+1)];
  std::valarray<double> P110 = Prim[M(i+1,j+1,k-1)];
  std::valarray<double> P111 = Prim[M(i+1,j+1,k+1)];

  double delx[3];
  delx[0] = 0.5 * (r[0] - domain.x_at(i-1)) / domain.get_dx(1);
  delx[1] = 0.5 * (r[1] - domain.y_at(j-1)) / domain.get_dx(2);
  delx[2] = 0.5 * (r[2] - domain.z_at(k-1)) / domain.get_dx(3);

  for (int q=0; q<Nq; ++q) {
    // http://en.wikipedia.org/wiki/Trilinear_interpolation

    double i1 = P000[q] * (1.0 - delx[2]) + P001[q] * delx[2];
    double i2 = P010[q] * (1.0 - delx[2]) + P011[q] * delx[2];
    double j1 = P100[q] * (1.0 - delx[2]) + P101[q] * delx[2];
    double j2 = P110[q] * (1.0 - delx[2]) + P111[q] * delx[2];

    double w1 = i1 * (1.0 - delx[1]) + i2 * delx[1];
    double w2 = j1 * (1.0 - delx[1]) + j2 * delx[1];

    Panswer[q] = w1 * (1.0 - delx[0]) + w2 * delx[0];
  }
  return Panswer;
}


#else
void Mara_prim_at_point(const double *r0, double *P1) { }
void Mara_prim_at_point_many(const double *Rin, double *Rlist, double *Plist, int Nsamp) { }
#endif // __MARA_USE_MPI
