/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2012 Carsten Burstedde

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
  if (p4est_wrap_adapt (wrap)) {
    if (p4est_wrap_partition (wrap, weight_exponent)) {
      p4est_wrap_complete (wrap);
    }
    return 1;
  }

  return 0;
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
    for (jl = 0, leaf = p4est_wrap_leaf_first (wrap); leaf != NULL;
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
    for (jl = 0, leaf = p4est_wrap_leaf_first (wrap); leaf != NULL;
         jl++, leaf = p4est_wrap_leaf_next (leaf)) {
      if (leaf->which_quad % 5 == 0) {
        p4est_wrap_mark_refine (wrap, leaf->which_tree, leaf->which_quad);
        p4est_wrap_mark_coarsen (wrap, leaf->which_tree, leaf->which_quad);
      }
    }
    SC_CHECK_ABORT (jl == wrap->p4est->local_num_quadrants, "Iterator");

    changed = wrap_adapt_partition (wrap, 0);
    SC_CHECK_ABORT (!changed, "Wrap noop");
  }
  
  for (loop = 0; loop < 2; ++loop) {
    /* mark for coarsening */
    for (jl = 0, leaf = p4est_wrap_leaf_first (wrap); leaf != NULL;
         jl++, leaf = p4est_wrap_leaf_next (leaf)) {
      if ((leaf->which_quad / 13) % 17 != 3) {
        p4est_wrap_mark_coarsen (wrap, leaf->which_tree, leaf->which_quad);
      }
    }
    SC_CHECK_ABORT (jl == wrap->p4est->local_num_quadrants, "Iterator");

    (void) wrap_adapt_partition (wrap, 0);
  }

  p4est_wrap_destroy (wrap);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
