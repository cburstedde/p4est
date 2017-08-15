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
#include <p4est_communication.h>
#include <p4est_extended.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#endif

typedef struct
{
  int                 maxlevel;
  int                 counter;
  int                 wrapper;
}
test_build_t;

static int
test_build_refine (p4est_t * p4est, p4est_topidx_t which_tree,
                   p4est_quadrant_t * quadrant)
{
  test_build_t       *tb;

  tb = (test_build_t *) p4est->user_pointer;

  if (quadrant->level >= tb->maxlevel) {
    return 0;
  }
  return !(tb->counter = (tb->counter + 1) % tb->wrapper);
}

/* specify in which tree we add which quadrant at all */
static p4est_topidx_t local_id[2][5] = {
  {2, 48, 94, -1, -1},
  {-1, -1, 2, 25, 71}
};

static void
test_build_local (sc_MPI_Comm mpicomm)
{
  p4est_topidx_t      treeid;
  p4est_locidx_t      lid;
  p4est_locidx_t      correct;
  p4est_connectivity_t *conn;
  p4est_t            *p4est, *built;
  p4est_tree_t       *ptree;
  p4est_tree_t        stree, *subtree = &stree;
  p4est_quadrant_t   *quadrant;
  test_build_t        stb, *tb = &stb;

  /* 0. prepare data that we will reuse */
  tb->maxlevel = 7 - P4EST_DIM;
  tb->counter = -1;
  tb->wrapper = 3;
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_moebius ();
#else
  SC_ABORT_NOT_REACHED ();
  conn = p8est_connectivity_new_rotcubes ();
#endif /* P4_TO_P8 */
  p4est = p4est_new_ext (mpicomm, conn, 0, 0, 0, 0, NULL, tb);
  p4est_refine (p4est, 1, test_build_refine, NULL);
  p4est_partition (p4est, 0, NULL);

  /* Create a minimal p4est to enable the use of complete_subtree */
  built = p4est_copy (p4est, 0);
  correct = 0;
  for (treeid = 0; treeid < 5; ++treeid) {
    /* grab subtree to keep current as much as possible */
    subtree = p4est_tree_array_index (built->trees, treeid);
    subtree->quadrants_offset += correct;

    /* now construct one quadrant to challenge complete_subtree */
    lid = -1;
    if (p4est->mpisize <= 2 && conn->num_trees <= 5) {
      lid = local_id[p4est->mpirank][treeid];
    }
    if (lid >= 0) {
      ptree = p4est_tree_array_index (p4est->trees, treeid);
      correct -= (p4est_locidx_t) ptree->quadrants.elem_count;

      /* we leak the quadrant's user data but there is none allocated */
      sc_array_resize (&subtree->quadrants, 0);
      memset (subtree->quadrants_per_level, 0,
              sizeof (p4est_locidx_t) * P4EST_MAXLEVEL);
      subtree->quadrants_per_level[P4EST_MAXLEVEL] = -1;
      subtree->maxlevel = 0;

      /* use lid to construct the prescribed quadrant */
      quadrant = p4est_quadrant_array_index
        (&ptree->quadrants, lid - ptree->quadrants_offset);
      P4EST_ASSERT (p4est_quadrant_is_valid (quadrant));

      /* add to emptied subtree */
      *(p4est_quadrant_t *) sc_array_push (&subtree->quadrants) = *quadrant;
      subtree->quadrants_per_level[quadrant->level] = 1;
      subtree->maxlevel = quadrant->level;

      /* and call the offending subtree routine */
      p4est_complete_subtree (built, treeid, NULL);
      correct += (p4est_locidx_t) subtree->quadrants.elem_count;
    }
  }
  built->local_num_quadrants += correct;
  p4est_comm_count_quadrants (built);
  P4EST_ASSERT (p4est_is_valid (built));

  /* clean up */
  p4est_destroy (built);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;

  /* Initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  /* Initialize packages */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* Test complete_subtree */
  test_build_local (mpicomm);

  /* Finalize */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
