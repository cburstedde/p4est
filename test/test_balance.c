/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007 Carsten Burstedde, Lucas Wilcox.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <p4est_base.h>
#include <p4est_algorithms.h>

int
main (int argc, char **argv)
{
  size_t              k;
  int8_t              l;
  p4est_t            *p4est;
  p4est_tree_t        stree, *tree = &stree;
  p4est_quadrant_t   *q;
  p4est_connectivity_t *connectivity;

  p4est_init (stdout, 0, NULL, NULL);

  connectivity = p4est_connectivity_new_unitsquare ();
  p4est = p4est_new (MPI_COMM_NULL, connectivity, 4, NULL);

  /* build empty tree */
  sc_array_init (&tree->quadrants, sizeof (p4est_quadrant_t));
  for (l = 0; l <= P4EST_MAXLEVEL; ++l) {
    tree->quadrants_per_level[l] = 0;
  }
  tree->maxlevel = 0;

  /* insert two quadrants */
  sc_array_resize (&tree->quadrants, 4);
  q = sc_array_index (&tree->quadrants, 0);
  p4est_quadrant_set_morton (q, 3, 13);
  q = sc_array_index (&tree->quadrants, 1);
  p4est_quadrant_set_morton (q, 1, 1);
  q = sc_array_index (&tree->quadrants, 2);
  p4est_quadrant_set_morton (q, 1, 2);
  q = sc_array_index (&tree->quadrants, 3);
  p4est_quadrant_set_morton (q, 1, 3);
  for (k = 0; k < tree->quadrants.elem_count; ++k) {
    q = sc_array_index (&tree->quadrants, k);
    q->p.user_data = sc_mempool_alloc (p4est->user_data_pool);
    ++tree->quadrants_per_level[q->level];
    tree->maxlevel = (int8_t) SC_MAX (tree->maxlevel, q->level);
  }

  /* balance the tree, print and destroy */
  p4est_balance_subtree (p4est, tree, 0, NULL);
  p4est_tree_print (SC_LP_INFO, tree);
  for (k = 0; k < tree->quadrants.elem_count; ++k) {
    q = sc_array_index (&tree->quadrants, k);
    sc_mempool_free (p4est->user_data_pool, q->p.user_data);
  }
  sc_array_reset (&tree->quadrants);

  /* balance the forest */
  p4est_balance (p4est, NULL);
  tree = sc_array_index (p4est->trees, 0);
  p4est_tree_print (SC_LP_INFO, tree);

  /* clean up memory */
  P4EST_ASSERT (p4est->user_data_pool->elem_count ==
                p4est->local_num_quadrants);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  sc_memory_check ();

  return 0;
}

/* EOF test_balance.c */
