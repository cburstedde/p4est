#include <p4est.h>
#include <p4est_vtk.h>
#include <p4est_extended.h>

/* refinement level initialization */
static int          refine_level = 6;

/* refinement function */
static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  double x_center, y_center, radius;

  /* Calculate the center coordinates of the quadrant */
  x_center = (double)(quadrant->x + (P4EST_QUADRANT_LEN(quadrant->level) / 2)) / P4EST_ROOT_LEN;
  y_center = (double)(quadrant->y + (P4EST_QUADRANT_LEN(quadrant->level) / 2)) / P4EST_ROOT_LEN;

  /* Define the radius of the circle */
  radius = 0.4;
  /* If the refinement level  */
  if ((int) quadrant->level <= (refine_level)) {

    /* Check if the center of the quadrant lies within the circle of radius 0.1 */
    if ((x_center - 0.5) * (x_center - 0.5) + (y_center - 0.5) * (y_center - 0.5) < radius * radius) {
      /* The center is within the circle, so refine this quadrant */
      return 1;
    }
  }
  return 0;
}


int main (int argc, char **argv) {
  p4est_connectivity_t *conn;
  p4est_t *p4est;
  /* Usage parameters for exercise 3*/


  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;
  int mpiret = sc_MPI_Init(&argc, &argv);
  SC_CHECK_MPI(mpiret);

  sc_init(mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init(NULL, SC_LP_DEFAULT);

  /* Exercise 1*/
  conn = p4est_connectivity_new_unitsquare();
  p4est = p4est_new(mpicomm, conn, 0, NULL, NULL);
  p4est_vtk_write_file(p4est, NULL, P4EST_STRING "_unitsquare_new");

  /* Exercise 2*/
  p4est_refine (p4est, 1, refine_fn, NULL);
  p4est_vtk_write_file(p4est, NULL, P4EST_STRING "_unitsquare_refine");

  // p4est_partition (p4est, 0, NULL);
  // p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_unitsquare_partition");

  /* Exercise 3
  Here there are two options */
  // p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  // p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_unitsquare_balance");
  
  p4est_partition_ext(p4est, 1, NULL);
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_unitsquare_partition_ext");

  p4est_destroy(p4est);
  p4est_connectivity_destroy(conn);


  mpiret = sc_MPI_Finalize();
  SC_CHECK_MPI(mpiret);

  return 0;
}
