
#include <p4est.h>
#include <p4est_base.h>

static const int32_t initial_quadrants_per_processor = 15;

p4est_connectivity_t *
p4est_connectivity_new (int32_t num_trees, int32_t num_vertices)
{
  p4est_connectivity_t *connectivity;

  connectivity = P4EST_ALLOC_ZERO (p4est_connectivity_t, 1);
  P4EST_CHECK_ALLOC (connectivity);

  connectivity->num_trees = num_trees;
  connectivity->num_vertices = num_vertices;

  connectivity->tree_to_vertex = P4EST_ALLOC (int32_t, 4 * num_trees);
  P4EST_CHECK_ALLOC (connectivity->tree_to_vertex);

  connectivity->tree_to_tree = P4EST_ALLOC (int32_t, 4 * num_trees);
  P4EST_CHECK_ALLOC (connectivity->tree_to_tree);

  connectivity->tree_to_face = P4EST_ALLOC (int8_t, 4 * num_trees);
  P4EST_CHECK_ALLOC (connectivity->tree_to_face);

  return connectivity;
}

void
p4est_connectivity_destroy (p4est_connectivity_t * connectivity)
{
  P4EST_FREE (connectivity->tree_to_face);
  P4EST_FREE (connectivity->tree_to_tree);
  P4EST_FREE (connectivity->tree_to_vertex);

  P4EST_FREE (connectivity);
}

p4est_t            *
p4est_new (MPI_Comm mpicomm, FILE * nout,
           p4est_connectivity_t * connectivity, int data_size)
{
#ifdef HAVE_MPI
  int                 mpiret;
#endif
  int8_t              level;
  int32_t             j, num_trees;
  int32_t             num_quadrants;
  int32_t             local_num_quadrants;
  int64_t             global_num_quadrants;
  p4est_t            *p4est;
  p4est_tree_t       *tree;

  p4est = P4EST_ALLOC_ZERO (p4est_t, 1);
  P4EST_CHECK_ALLOC (p4est);

  /* assign some data members */
  p4est->nout = nout;
  p4est->data_size = data_size;
  p4est->connectivity = connectivity;
  num_trees = connectivity->num_trees;

  p4est->mpicomm = mpicomm;
  p4est->mpisize = 1;
  p4est->mpirank = 0;
#ifdef HAVE_MPI
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Comm_size (p4est->mpicomm, &p4est->mpisize);
    P4EST_CHECK_MPI (mpiret);

    mpiret = MPI_Comm_rank (p4est->mpicomm, &p4est->mpirank);
    P4EST_CHECK_MPI (mpiret);
  }
#endif

  /* determine uniform level of initial tree */
  num_quadrants = 1;
  for (level = 0; level < 16; ++level) {
    if (num_quadrants >=
        (p4est->mpisize * initial_quadrants_per_processor) / num_trees) {
      break;
    }
    num_quadrants *= 4;
  }
  P4EST_ASSERT (level < 16);

  /* compute index of first tree for this processor */
  global_num_quadrants = (int64_t) num_quadrants *(int64_t) num_trees;
  local_num_quadrants = (int32_t) (global_num_quadrants / p4est->mpisize);

  if (p4est->mpirank == 0 && p4est->nout != NULL) {
    fprintf (p4est->nout, "New forest: %d trees %d processors\n",
             num_trees, p4est->mpisize);
    fprintf (p4est->nout,
             "   initial level %d global quadrants %lld local %d\n", level,
             (long long) global_num_quadrants, local_num_quadrants);
  }

  /* allocate memory pools */
  p4est->user_data_pool = p4est_mempool_new (p4est->data_size);

  /* allocate trees and quadrants */
  p4est->trees = p4est_array_new (sizeof (p4est_tree_t));
  p4est_array_resize (p4est->trees, num_trees);
  for (j = 0; j < num_trees; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    tree->quadrants = p4est_array_new (sizeof (p4est_quadrant_t));
  }

  return p4est;
}

void
p4est_destroy (p4est_t * p4est)
{
  int32_t             j;
  p4est_tree_t       *tree;

  for (j = 0; j < p4est->connectivity->num_trees; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    p4est_array_destroy (tree->quadrants);
  }
  p4est_array_destroy (p4est->trees);

  p4est_mempool_destroy (p4est->user_data_pool);

  p4est_connectivity_destroy (p4est->connectivity);

  P4EST_FREE (p4est);
}

/* EOF p4est.c */
