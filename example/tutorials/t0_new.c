#include <p4est.h>

int main (int argc, char ** argv) {
  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  P4EST_GLOBAL_PRODUCTION ("Hello World!\n");
  P4EST_PRODUCTION ("Hello World from the parallel process!\n");

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
