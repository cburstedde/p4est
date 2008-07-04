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

#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_mesh.h>

static const int    refine_level = 5;

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && p4est_quadrant_child_id (quadrant) == 3) {
    return 1;
  }
  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->y >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

int
main (int argc, char **argv)
{
  MPI_Comm            mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  unsigned            crc;
  size_t              k;
  int8_t              l;
  p4est_t            *p4est;
  p4est_tree_t        stree, *tree = &stree;
  p4est_quadrant_t   *q;
  p4est_connectivity_t *connectivity;

  /* initialize MPI */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpirank, NULL, NULL, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  connectivity = p4est_connectivity_new_unitsquare ();
  p4est = p4est_new (mpicomm, connectivity, 4, NULL);

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
#if 0
  p4est_balance_subtree (p4est, 0, NULL);
  p4est_tree_print (SC_LP_INFO, tree);
#endif
  for (k = 0; k < tree->quadrants.elem_count; ++k) {
    q = sc_array_index (&tree->quadrants, k);
    sc_mempool_free (p4est->user_data_pool, q->p.user_data);
  }
  sc_array_reset (&tree->quadrants);

  /* refine and balance the forest */
  p4est_refine (p4est, refine_fn, NULL);
  SC_CHECK_ABORT (!p4est_is_balanced (p4est), "!Balance");
  p4est_balance (p4est, NULL);
  SC_CHECK_ABORT (p4est_is_balanced (p4est), "Balance");

  /* checksum and rebalance */
  crc = p4est_checksum (p4est);
  p4est_balance (p4est, NULL);
  SC_CHECK_ABORT (p4est_checksum (p4est) == crc, "Rebalance");

  /* clean up and exit */
  P4EST_ASSERT (p4est->user_data_pool->elem_count ==
                (size_t) p4est->local_num_quadrants);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_balance.c */
