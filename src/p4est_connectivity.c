
#include <p4est_connectivity.h>
#include <p4est_base.h>

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

p4est_connectivity_t *
p4est_connectivity_new_unitsquare (void)
{
  p4est_connectivity_t *connectivity;

  connectivity = p4est_connectivity_new (1, 4);

  /* assign vertex numbers */
  connectivity->tree_to_vertex[0] = 0;
  connectivity->tree_to_vertex[1] = 1;
  connectivity->tree_to_vertex[2] = 2;
  connectivity->tree_to_vertex[3] = 3;

  /* we do not have neighbors, put in ourself */
  connectivity->tree_to_tree[0] = 0;
  connectivity->tree_to_tree[1] = 0;
  connectivity->tree_to_tree[2] = 0;
  connectivity->tree_to_tree[3] = 0;

  /* we do not share faces, put in our own face */
  connectivity->tree_to_face[0] = 0;
  connectivity->tree_to_face[1] = 1;
  connectivity->tree_to_face[2] = 2;
  connectivity->tree_to_face[3] = 3;

  return connectivity;
}

/* EOF p4est_connectivity.c */
