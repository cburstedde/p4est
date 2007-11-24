
#include <p4est.h>
#include <p4est_base.h>

p4est_connectivity_t *
p4est_connectivity_new (int32_t num_trees, int32_t num_vertices)
{
  p4est_connectivity_t *connectivity;

  connectivity = P4EST_ALLOC (p4est_connectivity_t, 1);

  connectivity->num_trees = num_trees;
  connectivity->num_vertices = num_vertices;

  connectivity->tree_to_vertex = P4EST_ALLOC (int32_t, 4 * num_trees);
  connectivity->tree_to_tree = P4EST_ALLOC (int32_t, 4 * num_trees);
  connectivity->tree_to_face = P4EST_ALLOC (int8_t, 4 * num_trees);

  return connectivity;
}

void
p4est_connectivity_destroy (p4est_connectivity_t * connectivity)
{
  P4EST_FREE (connectivity->tree_to_face);
  P4EST_FREE (connectivity->tree_to_tree);
  P4EST_FREE (connectivity->tree_to_vertex);

  P4EST_FREE (connectivity);

  return;
}

/* EOF p4est.c */
