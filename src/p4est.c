
#include <p4est.h>
#include <p4est_base.h>

static const int64_t initial_quadrants_per_processor = 15;

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
  int64_t             tree_num_quadrants, global_num_quadrants;
  int64_t             first_tree, first_quadrant, first_tree_quadrant;
  int64_t             last_tree, last_quadrant, last_tree_quadrant;
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

  /* allocate memory pools */
  p4est->user_data_pool = p4est_mempool_new (p4est->data_size);

  /* determine uniform level of initial tree */
  tree_num_quadrants = 1;
  for (level = 0; level < 16; ++level) {
    if (tree_num_quadrants >=
        (p4est->mpisize * initial_quadrants_per_processor) / num_trees) {
      break;
    }
    tree_num_quadrants *= 4;
  }
  P4EST_ASSERT (level < 16 && tree_num_quadrants <= INT32_MAX);

  /* compute index of first tree for this processor */
  global_num_quadrants = tree_num_quadrants * num_trees;
  first_quadrant = (global_num_quadrants * p4est->mpirank) / p4est->mpisize;
  first_tree = first_quadrant / tree_num_quadrants;
  first_tree_quadrant = first_quadrant - first_tree * tree_num_quadrants;
  last_quadrant =
    (global_num_quadrants * (p4est->mpirank + 1)) / p4est->mpisize - 1;
  last_tree = last_quadrant / tree_num_quadrants;
  last_tree_quadrant = last_quadrant - last_tree * tree_num_quadrants;

  /* check ranges of various integers */
  P4EST_ASSERT (first_tree <= last_tree && last_tree < num_trees);
  P4EST_ASSERT (0 <= first_tree_quadrant && 0 <= last_tree_quadrant);
  P4EST_ASSERT (first_tree_quadrant < tree_num_quadrants);
  P4EST_ASSERT (last_tree_quadrant < tree_num_quadrants);

  /* print some diagnostics */
  if (p4est->nout != NULL) {
    if (p4est->mpirank == 0) {
      fprintf (p4est->nout, "New forest: %d trees %d processors\n",
               num_trees, p4est->mpisize);
      fprintf (p4est->nout,
               "   initial level %d global quadrants %lld per tree %lld\n",
               level, (long long) global_num_quadrants,
               (long long) tree_num_quadrants);
    }
    fprintf (p4est->nout, "[%d] first tree %lld first quadrant %lld"
             " global quadrant %lld\n", p4est->mpirank,
             (long long) first_tree, (long long) first_tree_quadrant,
             (long long) first_quadrant);
    fprintf (p4est->nout, "[%d] last tree %lld last quadrant %lld"
             " global quadrant %lld\n", p4est->mpirank,
             (long long) last_tree, (long long) last_tree_quadrant,
             (long long) last_quadrant);
  }

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
