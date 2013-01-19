
#include <Python.h>
#if (COW_MPI)
#include <mpi.h>
#endif

void init_cow(void);

int main(int argc, char **argv)
{
#if (COW_MPI)
  {
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) freopen("/dev/null", "w", stdout);
    printf("was compiled with MPI support\n");
  }
#endif
  Py_Initialize();
  PySys_SetArgv(argc, argv);
  Py_SetProgramName("/Users/jzrake/Work/cow/cowpy");
  init_cow();

  if (argc > 1) {
    FILE *fp = fopen(argv[1], "r");
    if (fp) {
      PyRun_SimpleFileExFlags(fp, argv[1], 1, NULL);
      fclose(fp);
    }
    else {
      printf("No such file %s\n", argv[1]);
    }
  }

  Py_Finalize();
  printf("finished...\n");
#if (COW_MPI)
  MPI_Finalize();
#endif
  return 0;
}
