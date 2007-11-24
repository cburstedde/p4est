
#include <p4est.h>

int
main (void)
{
  int32_t             i;
  p4est_t            *p4est;
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

  /* ownership of the connectivity structure transfers to p4est */
  p4est = p4est_new (MPI_COMM_NULL, connectivity);

  /* destroy the 4est and its connectivity structure */
  p4est_destroy (p4est);

  return 0;
}

/* EOF simple.c */
