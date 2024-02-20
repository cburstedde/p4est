/* please add the usual GPL license header */

/* what about renaming?
   t0_init
   t1_new
   t1_file
   t2_refine
 */

#include <p4est.h>   /* Include the p4est library header for parallel adaptive mesh refinement. */

int main (int argc, char ** argv) {
  /* Initialize the MPI communicator to the default world communicator, which includes all MPI processes. */
  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;

  /* Initialize the MPI environment with arguments from the command line. */
  int mpiret = sc_MPI_Init (&argc, &argv);

  /* Check the return status of MPI initialization and abort if it failed. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the SC library. */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  /* Initialize the p4est library with default logging priority. */
  p4est_init (NULL, SC_LP_DEFAULT);

  /* Print a global production-level message saying "Hello World!". */
  P4EST_GLOBAL_PRODUCTION ("Hello World!\n");

  /* Print a production-level message from the current MPI process. */
  P4EST_PRODUCTION ("Hello World from the parallel process!\n");

  /* Finalize the MPI environment and clean up all MPI resources. */
  mpiret = sc_MPI_Finalize ();

  /* Check the return status of MPI finalization and abort if an error occurred. */
  SC_CHECK_MPI (mpiret);

  /* Return 0 from main to indicate that the program has finished successfully. */
  return 0;
}
