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
#include <p4est_wrap.h>
#else
#include <p8est_wrap.h>
#endif

static int
wrap_adapt_partition (p4est_wrap_t * wrap, int weight_exponent)
{
  p4est_locidx_t      uf, ul;

  if (p4est_wrap_adapt (wrap)) {
    if (p4est_wrap_partition (wrap, weight_exponent, &uf, &ul, NULL)) {

      SC_CHECK_ABORT (uf >= 0 && ul >= 0, "Invalid post window");
      SC_CHECK_ABORT (uf + ul <= wrap->p4est->local_num_quadrants,
                      "Invalid post count");

      p4est_wrap_complete (wrap);
    }
    return 1;
  }

  return 0;
}

static void
test_coarsen_delay (p4est_wrap_t * wrap)
{
  p4est_locidx_t      jl;
  p4est_wrap_leaf_t  *leaf;
  p4est_wrap_t       *copy1, *copy2;

  p4est_wrap_set_coarsen_delay (wrap, 2, 0);

  for (jl = 0, leaf = p4est_wrap_leaf_first (wrap, 1); leaf != NULL;
       jl++, leaf = p4est_wrap_leaf_next (leaf)) {
    if (leaf->which_quad % 4 == 0) {
      p4est_wrap_mark_refine (wrap, leaf->which_tree, leaf->which_quad);
    }
  }

  copy1 = p4est_wrap_new_copy (wrap, 17, NULL, NULL);
  copy2 = p4est_wrap_new_copy (copy1, 0, NULL, NULL);

  wrap_adapt_partition (wrap, 1);

  /* copies must be destroyed before the original is destroyed, in any order */
  p4est_wrap_destroy (copy1);

  for (jl = 0, leaf = p4est_wrap_leaf_first (wrap, 1); leaf != NULL;
       jl++, leaf = p4est_wrap_leaf_next (leaf)) {
    p4est_wrap_mark_coarsen (wrap, leaf->which_tree, leaf->which_quad);
  }
  wrap_adapt_partition (wrap, 1);

  /* copies must be destroyed before the original is destroyed, anywhere */
  p4est_wrap_destroy (copy2);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 changed;
  int                 loop;
#ifdef P4EST_ENABLE_DEBUG
  int                 lp = SC_LP_DEFAULT;
#else
  int                 lp = SC_LP_PRODUCTION;
#endif
  p4est_topidx_t      treecount;
  p4est_locidx_t      jl;
  p4est_wrap_leaf_t  *leaf;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
  sc_MPI_Comm         mpicomm;
  p4est_wrap_t       *wrap;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  sc_init (mpicomm, 0, 0, NULL, lp);
  p4est_init (NULL, lp);

#ifndef P4_TO_P8
  wrap = p4est_wrap_new_rotwrap (mpicomm, 0);
#else
  wrap = p8est_wrap_new_rotwrap (mpicomm, 0);
#endif
  ghost = p4est_wrap_get_ghost (wrap);
  SC_CHECK_ABORT (ghost != NULL, "Get ghost");
  ghost = NULL;
  mesh = p4est_wrap_get_mesh (wrap);
  SC_CHECK_ABORT (mesh != NULL, "Get mesh");
  mesh = NULL;

  for (loop = 0; loop < 3; ++loop) {
    /* mark for refinement */
    for (jl = 0, leaf = p4est_wrap_leaf_first (wrap, 1); leaf != NULL;
         jl++, leaf = p4est_wrap_leaf_next (leaf)) {
      if (leaf->which_quad % 3 == 0) {
        p4est_wrap_mark_refine (wrap, leaf->which_tree, leaf->which_quad);
      }
    }
    SC_CHECK_ABORT (jl == wrap->p4est->local_num_quadrants, "Iterator");

    changed = wrap_adapt_partition (wrap, 1);
    SC_CHECK_ABORT (changed, "Wrap refine");
  }

  for (loop = 0; loop < 2; ++loop) {
    /* mark some elements for coarsening that does not effect anything */
    treecount = 0;
    for (jl = 0, leaf = p4est_wrap_leaf_first (wrap, 0); leaf != NULL;
         jl++, leaf = p4est_wrap_leaf_next (leaf)) {
      if (P4EST_LEAF_IS_FIRST_IN_TREE (leaf)) {
        ++treecount;
      }
      if (leaf->which_quad % 5 == 0) {
        p4est_wrap_mark_refine (wrap, leaf->which_tree, leaf->which_quad);
        p4est_wrap_mark_coarsen (wrap, leaf->which_tree, leaf->which_quad);
      }
    }
    SC_CHECK_ABORT (jl == wrap->p4est->local_num_quadrants, "Iterator");
    /* this test should also be fine with empty processors */
    SC_CHECK_ABORT (treecount == wrap->p4est->last_local_tree -
                    wrap->p4est->first_local_tree + 1, "Iterator");

    changed = wrap_adapt_partition (wrap, 0);
    SC_CHECK_ABORT (!changed, "Wrap noop");
  }

  for (loop = 0; loop < 2; ++loop) {
    /* mark for coarsening */
    for (jl = 0, leaf = p4est_wrap_leaf_first (wrap, 1); leaf != NULL;
         jl++, leaf = p4est_wrap_leaf_next (leaf)) {
      if ((leaf->which_quad / 13) % 17 != 3) {
        p4est_wrap_mark_coarsen (wrap, leaf->which_tree, leaf->which_quad);
      }
    }
    SC_CHECK_ABORT (jl == wrap->p4est->local_num_quadrants, "Iterator");

    (void) wrap_adapt_partition (wrap, 0);
  }

  test_coarsen_delay (wrap);

  p4est_wrap_destroy (wrap);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
