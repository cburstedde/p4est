
#include <p4est_algorithms.h>
#include <p4est_base.h>

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

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est1);
  p4est_destroy (p4est2);
  p4est_connectivity_destroy (connectivity);

  p4est_memory_check ();

  return 0;
}

/* EOF test_quadrants.c */
