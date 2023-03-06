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

#include <p4est_iterate.h>
#include <p4est_trimesh.h>

typedef struct trimesh_meta
{
  p4est_locidx_t      lenum;
  p4est_locidx_t      num_owned;
  p4est_locidx_t      num_shared;
  p4est_locidx_t      szero[25];
  p4est_lnodes_code_t *face_code;
  p4est_trimesh_t    *tm;
}
trimesh_meta_t;

static const int vertex_seq[9] = { 4, 3, 5, 1, 7, 0, 2, 6, 8 };

static void
iter_volume1 (p4est_iter_volume_info_t * vi, void *user_data)
{
  trimesh_meta_t     *me = (trimesh_meta_t *) user_data;
  p4est_lnodes_t     *ln = me->tm->lnodes;
  p4est_locidx_t      le;
#ifdef P4EST_ENABLE_DEBUG
  p4est_tree_t       *tree;

  tree = p4est_tree_array_index (vi->p4est->trees, vi->treeid);
  P4EST_ASSERT (tree->quadrants_offset + vi->quadid == me->lenum);

#endif
  le = me->lenum++;
  P4EST_ASSERT (ln->face_code[le] == 0);
  P4EST_ASSERT (!memcmp (ln->element_nodes + ln->vnodes * le,
                         me->szero, sizeof (p4est_locidx_t) * ln->vnodes));

  /* place owned node in center of quadrant */
  ln->element_nodes[ln->vnodes * le + vertex_seq[0]] = me->num_owned++;
}

static void
iter_face1 (p4est_iter_face_info_t * fi, void *user_data)
{
}

static void
iter_corner1 (p4est_iter_corner_info_t * ci, void *user_data)
{
}

p4est_trimesh_t    *
p4est_trimesh_new (p4est_t * p4est, p4est_ghost_t * ghost, int with_edge)
{
  int                 mpiret;
  int                 p, q, s;
  int                 vn;
  p4est_locidx_t      le;
  p4est_gloidx_t      gc;
  p4est_trimesh_t    *tm;
  p4est_lnodes_t     *ln;
  trimesh_meta_t      tmeta, *me = &tmeta;

  P4EST_ASSERT (p4est_is_balanced (p4est, P4EST_CONNECT_FACE));

  memset (me, 0, sizeof (trimesh_meta_t));
  tm = me->tm = P4EST_ALLOC_ZERO (p4est_trimesh_t, 1);
  ln = tm->lnodes = P4EST_ALLOC_ZERO (p4est_lnodes_t, 1);

  p = p4est->mpirank;
  s = p4est->mpisize;
  ln->mpicomm = p4est->mpicomm;
  ln->sharers = sc_array_new (sizeof (p4est_lnodes_rank_t));
  ln->degree = 0;
  vn = ln->vnodes = 9 + (with_edge ? 16 : 0);
  le = ln->num_local_elements = p4est->local_num_quadrants;
  ln->face_code = P4EST_ALLOC_ZERO (p4est_lnodes_code_t, le);
  ln->element_nodes = P4EST_ALLOC_ZERO (p4est_locidx_t, le * vn);

  /* determine the face_code for each element */
  me->lenum = 0;
  me->face_code = ln->face_code;
  p4est_iterate (p4est, ghost, me, iter_volume1, iter_face1, iter_corner1);
  P4EST_ASSERT (me->lenum == le);
  P4EST_INFOF ("p4est_trimesh_new: owned %ld shared %ld\n",
               (long) me->num_owned, (long) me->num_shared);

  /* send messages */

  /* share owned count */
  ln->global_owned_count = P4EST_ALLOC (p4est_locidx_t, s);
  mpiret = sc_MPI_Allgather (&me->num_owned, 1, P4EST_MPI_LOCIDX,
                             ln->global_owned_count, 1, P4EST_MPI_LOCIDX,
                             p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  ln->global_offset = 0;
  for (q = 0; q < p; ++q) {
    ln->global_offset += ln->global_owned_count[q];
  }
  for (gc = ln->global_offset; q < s; ++q) {
    gc += ln->global_owned_count[q];
  }
  P4EST_GLOBAL_PRODUCTIONF ("p4est_trimesh_new: global owned %lld\n",
                            (long long) gc);

  /* receive messages */
  /* populate nflags */
  /* finalize lnodes */

  return tm;
}

void
p4est_trimesh_destroy (p4est_trimesh_t * tm)
{
  P4EST_ASSERT (tm != NULL);
  P4EST_ASSERT (tm->lnodes != NULL);

  p4est_lnodes_destroy (tm->lnodes);
  P4EST_FREE (tm->nflags);
  P4EST_FREE (tm);
}
