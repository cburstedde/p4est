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
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_search_build.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_search_build.h>
#endif

/** Context object for building a new p4est from search callbacks.
 */
struct p4est_search_build
{
  p4est_t            *p4est;    /**< New forest being built. */

  /* Context for the tree walk */
  p4est_topidx_t      cur_tree;         /**< Current tree under examination. */
  p4est_tree_t       *tree;             /**< Pointer to current tree. */
  sc_array_t         *tquadrants;       /**< Points to the tree's quadrants. */

  /* TODO: Does it make sense to add an init_fn callback? */
  p4est_init_t        init_local;       /**< By p4est_search_build_local. */
  p4est_init_t        init_complete;    /**< Used on tree completion. */
};

static void
p4est_search_build_begin_tree (p4est_search_build_t * build,
                               p4est_topidx_t which_tree,
                               p4est_locidx_t quadrants_offset)
{
  p4est_t            *p4est;

  /* check sanity of call */
  P4EST_ASSERT (build != NULL);
  P4EST_ASSERT (build->p4est != NULL);

  /* check sanity of build object */
  p4est = build->p4est;
  P4EST_ASSERT (p4est->first_local_tree <= which_tree &&
                which_tree <= p4est->last_local_tree);
  P4EST_ASSERT (0 <= quadrants_offset);

  /* initialize context for a new tree */
  build->cur_tree = which_tree;
  build->tree = p4est_tree_array_index (p4est->trees, build->cur_tree);
  build->tree->quadrants_offset = quadrants_offset;
  build->tquadrants = &build->tree->quadrants;
  P4EST_ASSERT (build->tquadrants->elem_size == sizeof (p4est_quadrant_t));
  P4EST_ASSERT (build->tquadrants->elem_count == 0);
}

p4est_search_build_t *
p4est_search_build_new (p4est_t * from, size_t data_size)
{
  p4est_topidx_t      jt, num_trees;
  p4est_t            *p4est;
  p4est_tree_t       *ftree, *ptree;
  p4est_search_build_t *build;

  P4EST_ASSERT (p4est_is_valid (from));

  /* create a new p4est structure to be populated */
  build = P4EST_ALLOC (p4est_search_build_t, 1);
  build->p4est = p4est = P4EST_ALLOC (p4est_t, 1);
  memcpy (p4est, from, sizeof (p4est_t));
  num_trees = from->connectivity->num_trees;

  /* remove anything that we will not use from the template forest */
  p4est->mpicomm_owned = 0;
  p4est->data_size = data_size;
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
  if (p4est->data_size > 0) {
    p4est->user_data_pool = sc_mempool_new (p4est->data_size);
  }
  p4est->quadrant_pool = sc_mempool_new (sizeof (p4est_quadrant_t));

  /* initialize context structure */
  build->init_local = NULL;
  build->init_complete = NULL;
  p4est_search_build_begin_tree (build, p4est->first_local_tree, 0);

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

static              p4est_locidx_t
p4est_search_build_end_tree (p4est_search_build_t * build)
{
  int8_t              ell;
  p4est_t            *p4est;

  /* check sanity of call */
  P4EST_ASSERT (build != NULL);
  P4EST_ASSERT (build->p4est != NULL);
  P4EST_ASSERT (build->tree != NULL);

  p4est = build->p4est;
  P4EST_ASSERT (build->tree->quadrants_per_level[P4EST_MAXLEVEL] == 0);
  P4EST_ASSERT (build->tree->maxlevel == 0);
  for (ell = P4EST_QMAXLEVEL; ell > 0; --ell) {
    if (build->tree->quadrants_per_level[ell] > 0) {
      build->tree->maxlevel = ell;
      break;
    }
  }
  P4EST_ASSERT (build->tquadrants->elem_size == sizeof (p4est_quadrant_t));

  /* do the heavy lifting: complete this tree as coarsely as possible */
  p4est_complete_subtree (p4est, build->cur_tree, build->init_complete);
  P4EST_ASSERT (&build->tree->quadrants == build->tquadrants);

  return
    build->tree->quadrants_offset +
    (p4est_locidx_t) build->tquadrants->elem_count;
}

#ifndef P4_TO_P8

int
p4est_search_build_local (p4est_search_build_t * build,
                          p4est_topidx_t which_tree,
                          p4est_quadrant_t * quadrant,
                          p4est_locidx_t local_num)
{
  p4est_t            *p4est;
  p4est_quadrant_t   *q;
  p4est_locidx_t      quadrants_offset;

  /* check sanity of call */
  P4EST_ASSERT (build != NULL);
  P4EST_ASSERT (build->p4est != NULL);

  /* check sanity of build object */
  p4est = build->p4est;
  P4EST_ASSERT (p4est->first_local_tree <= which_tree &&
                which_tree <= p4est->last_local_tree);
  P4EST_ASSERT (p4est->first_local_tree <= build->cur_tree &&
                build->cur_tree <= p4est->last_local_tree);
  P4EST_ASSERT (which_tree == build->cur_tree ||
                which_tree == build->cur_tree + 1);

  /* finish up previous tree if we are entering a new one */
  if (which_tree > build->cur_tree) {
    quadrants_offset = p4est_search_build_end_tree (build);
    p4est_search_build_begin_tree (build, which_tree, quadrants_offset);
  }

  /* we do nothing if we are not at a leaf of the tree */
  if (local_num < 0) {
    return 0;
  }

  /*           *** Notes ***
   * 0. The search goes through the nonempty local trees.
   * 1. The search may begin below the root at the top of a branch.
   * 2. The search may skip intermediate levels in the tree.
   */

  /* insert only relevant leaves */
  P4EST_ASSERT (p4est_quadrant_is_valid (quadrant));
  P4EST_ASSERT (build->tquadrants->elem_size == sizeof (p4est_quadrant_t));
  q = (p4est_quadrant_t *) sc_array_push (build->tquadrants);
  *q = *quadrant;
  p4est_quadrant_init_data (p4est, which_tree, q, build->init_local);
  ++build->tree->quadrants_per_level[q->level];

  /* TODO: figure out if we need a return value */
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
  P4EST_ASSERT (build->p4est != NULL);

  p4est = build->p4est;
  P4EST_FREE (build);

  if (p4est->first_local_tree <= p4est->last_local_tree) {
    /* finish last tree of the iteration */
    P4EST_ASSERT (build->cur_tree == p4est->last_local_tree);
    p4est->local_num_quadrants = p4est_search_build_end_tree (build);

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
