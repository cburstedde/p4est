#include <p4est_extended.h>
#include <p4est_vtk.h>

/* also please run p4estindent on all main programs */

int main (int argc, char **argv) {
  /* document me */
  const int minlevel = 3;

  /* Declare the MPI communicator and initialize the MPI environment */
  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;
  int mpiret = sc_MPI_Init(&argc, &argv);
  /* Check for successful MPI initialization */
  SC_CHECK_MPI(mpiret);

  /* Initialize the SC library with 1-byte alignment and default log priority */
  sc_init(mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  /* Initialize the p4est library with default parameters */
  p4est_init(NULL, SC_LP_DEFAULT);

  /* Create a new connectivity for a unit square domain */
  p4est_connectivity_t *conn;
  conn = p4est_connectivity_new_unitsquare();

  /* Create a new p4est structure (forest of quadtrees) with the created connectivity */
  p4est_t *p4est;
  /* please comment on minlevel, fill_uniform always 1, min_quadrants always 0
     (to make the mesh independend from the number of processes */
  p4est = p4est_new_ext (mpicomm, conn, 0, minlevel, 1, 0, NULL, NULL);

  /* Write the forest structure to a VTK file for visualization purposes */
  /* The filename will be prefixed with "p4est_a_unitsquare" */
  p4est_vtk_write_file(p4est, NULL, P4EST_STRING "_a_unitsquare");

  /* Destroy the p4est structure to free memory */
  p4est_destroy(p4est);
  /* Destroy the connectivity structure to free memory */
  p4est_connectivity_destroy(conn);

  /* Finalize the MPI environment and check for errors */
  mpiret = sc_MPI_Finalize();
  SC_CHECK_MPI(mpiret);

  return 0;
}
