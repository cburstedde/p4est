#include <p4est.h>
#include <p4est_vtk.h>

int main(int argc, char **argv) {
  p4est_connectivity_t *conn_out, *conn_in; /* Pointers to p4est connectivity structures */
  p4est_t *p4est;                          /* Pointer to the p4est structure */
  int save;                                /* Variable to store the status of save operation. 0: a successfull save*/
  size_t *bytes = NULL;                    /* Pointer to store the size of the connectivity in bytes */
  const char *filename = P4EST_STRING "_a_unitsquare"; /* File name for saving/loading the connectivity */

  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;
  int mpiret = sc_MPI_Init(&argc, &argv);
  SC_CHECK_MPI(mpiret);
  sc_init(mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init(NULL, SC_LP_DEFAULT);

  /* Create a new unit square connectivity which describes how the quadrants are connected */
  conn_out = p4est_connectivity_new_unitsquare();

  /* Create a new p4est (a forest of quadtrees) structure using the unit square connectivity */
  p4est = p4est_new(mpicomm, conn_out, 0, NULL, NULL);

  /* Print a message indicating the start of the connectivity file writing process */
  P4EST_PRODUCTION("Writing connectivity file\n");
  /* Save the connectivity to a file */
  save = p4est_connectivity_save(filename, conn_out);
  if (save) {
    /* Print an error message if saving the connectivity failed */
    P4EST_PRODUCTION("Error writing connectivity file\n");
  }
  /* Destroy the p4est structure to free memory */
  p4est_destroy(p4est);
  /* Destroy the output connectivity structure to free memory */
  p4est_connectivity_destroy(conn_out);

  /* Print a message indicating the start of the connectivity file loading process */
  P4EST_PRODUCTION("Loading connectivity file\n");
  /* Load the connectivity from a file */
  conn_in = p4est_connectivity_load(filename, bytes);
  if (conn_in == NULL) {
    /* Print an error message if loading the connectivity failed */
    P4EST_LERRORF("Could not read file %s\n", filename);
    return 1;
  }
  /* Create a new p4est structure using the loaded connectivity */
  p4est = p4est_new(mpicomm, conn_in, 0, NULL, NULL);
  /* Destroy the p4est structure to free memory */
  p4est_destroy(p4est);
  /* Destroy the input connectivity structure to free memory */
  p4est_connectivity_destroy(conn_in);

  /* Finalize the MPI environment */
  mpiret = sc_MPI_Finalize();
  SC_CHECK_MPI(mpiret);
  return 0; /* Return 0 to indicate successful execution */
}
