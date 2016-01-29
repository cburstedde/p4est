/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
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
#include <p4est_communication.h>
#include <p4est_search_build.h>
#else
#include <p8est_algorithms.h>
#include <p8est_communication.h>
#include <p8est_search_build.h>
#endif

p4est_search_build_t *
p4est_search_build_new (p4est_t * from)
{
  p4est_topidx_t      jt, num_trees;
  p4est_t            *p4est;
  p4est_tree_t       *ftree, *ptree;
  p4est_search_build_t *build;

  P4EST_ASSERT (p4est_is_valid (from));

  /* create a new p4est structure to be populated */
  build = P4EST_ALLOC (p4est_search_build_t, 1);
  build->from = from;
  build->p4est = p4est = P4EST_ALLOC (p4est_t, 1);
  memcpy (p4est, from, sizeof (p4est_t));
  num_trees = from->connectivity->num_trees;

  /* remove anything that we will not use from the template forest */
  p4est->mpicomm_owned = 0;
  p4est->data_size = 0;
  p4est->user_pointer = NULL;
  p4est->local_num_quadrants = 0;
  p4est->global_num_quadrants = 0;
  p4est->global_first_quadrant = NULL;
  p4est->global_first_position = NULL;
  p4est->trees = NULL;
  p4est->user_data_pool = NULL;
  p4est->quadrant_pool = NULL;
  p4est->inspect = NULL;

  /* start populating missing members */
  p4est->global_first_quadrant =
    P4EST_ALLOC (p4est_gloidx_t, p4est->mpisize + 1);
  p4est->global_first_position =
    P4EST_ALLOC (p4est_quadrant_t, p4est->mpisize + 1);
  memcpy (p4est->global_first_position, from->global_first_position,
          (p4est->mpisize + 1) * sizeof (p4est_quadrant_t));
  p4est->trees = sc_array_new_size (sizeof (p4est_tree_t), num_trees);
  for (jt = 0; jt < num_trees; ++jt) {
    ftree = p4est_tree_array_index (from->trees, jt);
    ptree = p4est_tree_array_index (p4est->trees, jt);
    sc_array_init (&ptree->quadrants, sizeof (p4est_quadrant_t));
    ptree->first_desc = ftree->first_desc;
    ptree->last_desc = ftree->last_desc;
    ptree->quadrants_offset = 0;
    memset (ptree->quadrants_per_level, 0,
            (P4EST_MAXLEVEL + 1) * sizeof (p4est_locidx_t));
    ptree->maxlevel = 0;
  }
  p4est->quadrant_pool = sc_mempool_new (sizeof (p4est_quadrant_t));

  /*
   * what remains to be filled:
   * local_num_quadrants
   * trees:
   *   quadrants
   *   quadrants_offset
   *   quadrants_per_level
   *   maxlevel
   * global_first_quadrant, global_num_quadrants:
   *   filled after local_num_quadrants has been set correctly
   */

  return build;
}

#ifndef P4_TO_P8

int
p4est_search_build_local (p4est_search_build_t * build,
                          p4est_topidx_t which_tree,
                          p4est_quadrant_t * quadrant,
                          p4est_locidx_t local_num)
{
  p4est_t            *p4est;

  P4EST_ASSERT (build != NULL);
  P4EST_ASSERT (build->from != NULL);
  P4EST_ASSERT (build->p4est != NULL);

  p4est = build->p4est;

  /*           *** Notes ***
   * 0. The search goes through the nonempty local trees.
   * 1. The search may begin below the root at the top of a branch.
   * 2. The search may skip intermediate levels in the tree.
   */

  return 0;
}

#endif /* !P4_TO_P8 */

p4est_t            *
p4est_search_build_complete (p4est_search_build_t * build)
{
  p4est_topidx_t      jt, last_local_tree, num_trees;
  p4est_t            *p4est;
  p4est_tree_t       *ptree;

  P4EST_ASSERT (build != NULL);
  P4EST_ASSERT (build->from != NULL);
  P4EST_ASSERT (build->p4est != NULL);

  p4est = build->p4est;
  P4EST_FREE (build);

  if (p4est->first_local_tree <= p4est->last_local_tree) {
    /* fix quadrants_offset in empty trees > last_local_tree */
    num_trees = p4est->connectivity->num_trees;
    last_local_tree = p4est->last_local_tree;
    P4EST_ASSERT (0 <= last_local_tree && last_local_tree < num_trees);
    P4EST_ASSERT (p4est->local_num_quadrants > 0);
    for (jt = last_local_tree + 1; jt < num_trees; ++jt) {
      ptree = p4est_tree_array_index (p4est->trees, jt);
      ptree->quadrants_offset = p4est->local_num_quadrants;
    }
  }
  else {
    /* this is an empty processor */
    P4EST_ASSERT (p4est->first_local_tree == -1);
    P4EST_ASSERT (p4est->last_local_tree == -2);
    P4EST_ASSERT (p4est->local_num_quadrants == 0);
  }

  /* fix global cumulative quadrant count per processor */
  p4est_comm_count_quadrants (p4est);
  P4EST_ASSERT (p4est_is_valid (p4est));

  return p4est;
}
