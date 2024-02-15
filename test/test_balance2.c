/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef P4_TO_P8
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_nodes.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_nodes.h>
#endif

#ifndef P4_TO_P8
static const int    refine_level = 5;
#else
static const int    refine_level = 3;
#endif

static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
}

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 cid;

  if (which_tree == 2 || which_tree == 3) {
    return 0;
  }

  cid = p4est_quadrant_child_id (quadrant);

  if (cid == P4EST_CHILDREN - 1 ||
      (quadrant->x >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2) &&
       quadrant->y >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#ifdef P4_TO_P8
       && quadrant->z >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#endif
      )) {
    return 1;
  }
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && cid == 2) {
    return 1;
  }
  if (quadrant->x == P4EST_QUADRANT_LEN (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->y >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

static unsigned
test_checksum (p4est_t * p4est, int have_zlib)
{
  return have_zlib ? p4est_checksum (p4est) : 0;
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  int                 have_zlib;
  unsigned            crc;
#ifndef P4_TO_P8
  size_t              kz;
  int8_t              l;
  p4est_quadrant_t   *q;
  p4est_tree_t        stree, *tree = &stree;
#endif
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* establish parallel logging */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* check for ZLIB usability */
  if (!(have_zlib = p4est_have_zlib ())) {
    P4EST_GLOBAL_LERROR
      ("Not found a working ZLIB installation: ignoring CRCs\n");
  }

#ifndef P4_TO_P8
  connectivity = p4est_connectivity_new_star ();
#else
  connectivity = p8est_connectivity_new_rotcubes ();
#endif
  p4est = p4est_new_ext (mpicomm, connectivity, 0, 0, 0, 4, NULL, NULL);

#ifndef P4_TO_P8
  /* build empty tree */
  sc_array_init (&tree->quadrants, sizeof (p4est_quadrant_t));
  for (l = 0; l <= P4EST_MAXLEVEL; ++l) {
    tree->quadrants_per_level[l] = 0;
  }
  tree->maxlevel = 0;

  /* insert two quadrants */
  sc_array_resize (&tree->quadrants, 4);
  q = p4est_quadrant_array_index (&tree->quadrants, 0);
  p4est_quadrant_set_morton (q, 3, 13);
  q = p4est_quadrant_array_index (&tree->quadrants, 1);
  p4est_quadrant_set_morton (q, 1, 1);
  q = p4est_quadrant_array_index (&tree->quadrants, 2);
  p4est_quadrant_set_morton (q, 1, 2);
  q = p4est_quadrant_array_index (&tree->quadrants, 3);
  p4est_quadrant_set_morton (q, 1, 3);
  for (kz = 0; kz < tree->quadrants.elem_count; ++kz) {
    q = p4est_quadrant_array_index (&tree->quadrants, kz);
    q->p.user_data = sc_mempool_alloc (p4est->user_data_pool);
    ++tree->quadrants_per_level[q->level];
    tree->maxlevel = (int8_t) SC_MAX (tree->maxlevel, q->level);
  }

  /* balance the tree, print and destroy */
#if 0
  p4est_balance_subtree (p4est, P4EST_CONNECT_FULL, 0, NULL);
  p4est_tree_print (SC_LP_INFO, tree);
#endif
  for (kz = 0; kz < tree->quadrants.elem_count; ++kz) {
    q = p4est_quadrant_array_index (&tree->quadrants, kz);
    sc_mempool_free (p4est->user_data_pool, q->p.user_data);
  }
  sc_array_reset (&tree->quadrants);
#endif /* !P4_TO_P8 */

  /* check reset data function */
  p4est_reset_data (p4est, 0, init_fn, NULL);
  p4est_reset_data (p4est, 0, NULL, NULL);

  /* refine and balance the forest */
  SC_CHECK_ABORT (p4est_is_balanced (p4est, P4EST_CONNECT_FULL), "Balance 1");
  p4est_refine (p4est, 1, refine_fn, NULL);
  SC_CHECK_ABORT (!p4est_is_balanced (p4est, P4EST_CONNECT_FULL),
                  "Balance 2");
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  SC_CHECK_ABORT (p4est_is_balanced (p4est, P4EST_CONNECT_FULL), "Balance 3");

  /* check reset data function */
  p4est_reset_data (p4est, 17, NULL, NULL);
  p4est_reset_data (p4est, 8, init_fn, NULL);

  /* checksum and partition */
  crc = test_checksum (p4est, have_zlib);
  p4est_partition (p4est, 0, NULL);
  SC_CHECK_ABORT (test_checksum (p4est, have_zlib) == crc, "Partition");
  SC_CHECK_ABORT (p4est_is_balanced (p4est, P4EST_CONNECT_FULL), "Balance 4");

  /* check reset data function */
  p4est_reset_data (p4est, 3, NULL, NULL);
  p4est_reset_data (p4est, 3, NULL, NULL);

  /* checksum and rebalance */
  crc = test_checksum (p4est, have_zlib);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  SC_CHECK_ABORT (test_checksum (p4est, have_zlib) == crc, "Rebalance");

  /* clean up and exit */
  P4EST_ASSERT (p4est->user_data_pool->elem_count ==
                (size_t) p4est->local_num_quadrants);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
