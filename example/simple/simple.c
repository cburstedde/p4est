
#include <p4est.h>
#include <p4est_base.h>

typedef struct
{
  int32_t             a;
}
user_data_t;

static void
init_fn (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * quadrant)
{
  user_data_t        *data = quadrant->user_data;

  data->a = which_tree;
}

int
main (int argc, char **argv)
{
#ifdef HAVE_MPI
  int                 use_mpi = 1;
  int                 mpiret;
#endif
  MPI_Comm            mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;

  mpicomm = MPI_COMM_NULL;
#ifdef HAVE_MPI
  if (use_mpi) {
    mpiret = MPI_Init (&argc, &argv);
    P4EST_CHECK_MPI (mpiret);
    mpicomm = MPI_COMM_WORLD;
  }
#endif

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_unitsquare ();
  p4est = p4est_new (mpicomm, stdout, connectivity,
                     sizeof (user_data_t), init_fn);

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  p4est_memory_check ();

#ifdef HAVE_MPI
  if (use_mpi) {
    mpiret = MPI_Finalize ();
    P4EST_CHECK_MPI (mpiret);
  }
#endif

  return 0;
}

/* EOF simple.c */
