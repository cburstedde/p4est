
#include <p4est_algorithms.h>
#include <p4est_base.h>

static int
refine_fn (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * q)
{
  return (q->level < 5);
}

int
main (void)
{
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est1;
  p4est_t            *p4est2;

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_unitsquare ();
  p4est1 = p4est_new (MPI_COMM_NULL, stdout, connectivity, 0, NULL);
  p4est2 = p4est_new (MPI_COMM_NULL, stdout, connectivity, 0, NULL);

  /* refine the second tree to a uniform level */
  p4est_refine (p4est2, refine_fn, NULL);

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est1);
  p4est_destroy (p4est2);
  p4est_connectivity_destroy (connectivity);

  p4est_memory_check ();

  return 0;
}

/* EOF test_quadrants.c */
