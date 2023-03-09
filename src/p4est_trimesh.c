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
  int                 with_faces;
  int                 mpisize, mpirank;
  int                *ghost_rank;
  p4est_locidx_t      lenum;
  p4est_locidx_t      num_owned;
  p4est_locidx_t      num_shared;
  p4est_locidx_t      szero[25];
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;
  p4est_trimesh_t    *tm;
}
trimesh_meta_t;

static void
set_lnodes_corner_center (p4est_lnodes_t * ln, p4est_locidx_t le,
                          p4est_locidx_t lni)
{
  p4est_locidx_t      lpos;

  P4EST_ASSERT (ln != NULL);
  P4EST_ASSERT (ln->vnodes == 9 || ln->vnodes == 25);
  P4EST_ASSERT (0 <= le && le < ln->num_local_elements);

  lpos = le * ln->vnodes + 0;
  P4EST_ASSERT (ln->element_nodes[lpos] == 0);
  ln->element_nodes[lpos] = lni;
}

static void
set_lnodes_face_full (p4est_lnodes_t * ln, p4est_locidx_t le,
                      int face,  p4est_locidx_t lni)
{
  p4est_locidx_t      lpos;

  P4EST_ASSERT (ln != NULL);
  P4EST_ASSERT (ln->vnodes == 25);
  P4EST_ASSERT (0 <= le && le < ln->num_local_elements);
  P4EST_ASSERT (0 <= face && face < P4EST_FACES);

  lpos = le * ln->vnodes + (9 + 8) + 2 * face;
  P4EST_ASSERT (ln->element_nodes[lpos] == 0);
  P4EST_ASSERT (ln->element_nodes[lpos + 1] == 0);
  ln->element_nodes[lpos] = lni;
  ln->element_nodes[lpos + 1] = -1;
}

static void
iter_volume1 (p4est_iter_volume_info_t * vi, void *user_data)
{
  trimesh_meta_t     *me = (trimesh_meta_t *) user_data;
  p4est_lnodes_t     *ln = me->tm->lnodes;
  p4est_locidx_t      le;
#ifdef P4EST_ENABLE_DEBUG
  p4est_tree_t       *tree;

  /* initial checks  */
  P4EST_ASSERT (vi->p4est == me->p4est);
  tree = p4est_tree_array_index (vi->p4est->trees, vi->treeid);
  P4EST_ASSERT (tree->quadrants_offset + vi->quadid == me->lenum);

#endif
  /* create owned node */
  le = me->lenum++;
  P4EST_ASSERT (ln->face_code[le] == 0);
  P4EST_ASSERT (!memcmp (ln->element_nodes + ln->vnodes * le,
                         me->szero, sizeof (p4est_locidx_t) * ln->vnodes));

  /* place owned node at quadrant midpoint */
  set_lnodes_corner_center (ln, le, me->num_owned++);
}

static void
iter_face1 (p4est_iter_face_info_t * fi, void *user_data)
{
  trimesh_meta_t     *me = (trimesh_meta_t *) user_data;
  p4est_lnodes_t     *ln = me->tm->lnodes;
  int                 i, j;
  int                 q;
  /* each face connection produces at most 3 nodes: 1 corner, 2 face */
  int                 nunodes;          /**< nodes on interface */
  int                 codim[3];         /**< codimension of a node */
  int                 is_owned[3];      /**< is that node locally owned */
  int                 is_shared[3];     /**< does the node have sharers */
  int                 owners[3][3];     /**< owner processes for each node */
  p4est_locidx_t      le;               /**< local element number */
  p4est_locidx_t      lo;               /**< owned node number */
  p4est_tree_t       *tree;             /**< tree within forest */
  p4est_iter_face_side_t *fs, *fss[2];
  p4est_iter_face_side_full_t *fu;

  /* initial checks  */
  P4EST_ASSERT (fi->p4est == me->p4est);

  /* a boundary face is the easiest case */
  if (fi->sides.elem_count == 1) {
    P4EST_ASSERT (fi->orientation == 0);
    P4EST_ASSERT (fi->tree_boundary == P4EST_CONNECT_FACE);
    fs = (p4est_iter_face_side_t *) sc_array_index_int (&fi->sides, 0);
    P4EST_ASSERT (!fs->is_hanging);
    P4EST_ASSERT (!fs->is.full.is_ghost);
    if (me->with_faces) {
      /* place owned node at boundary face midpoint */
      tree = p4est_tree_array_index (fi->p4est->trees, fs->treeid);
      le = tree->quadrants_offset + fs->is.full.quadid;
      set_lnodes_face_full (ln, le, fs->face, me->num_owned++);
    }
    return;
  }

  /* find ownership of all nodes on this face connection */
  nunodes = 0;
  for (i = 0; i < 3; ++i) {
    codim[i] = -1;
    is_owned[i] = is_shared[i] = 0;
    for (j = 0; j < 3; ++j) {
      owners[i][j] = -1;
    }
  }

  /* we have two sides to the face connection */
  P4EST_ASSERT (fi->sides.elem_count == 2);
  fss[0] = (p4est_iter_face_side_t *) sc_array_index_int (&fi->sides, 0);
  fss[1] = (p4est_iter_face_side_t *) sc_array_index_int (&fi->sides, 1);
  P4EST_ASSERT (!fss[0]->is_hanging || !fss[1]->is_hanging);
  if (!fss[0]->is_hanging && !fss[1]->is_hanging) {
    if (me->with_faces) {
      /* one face node on same-size connection */
      nunodes = 1;
      codim[0] = 1;
      is_owned[0] = 1;
      for (i = 0; i < 2; ++i) {
        fu = &fss[i]->is.full;

        /* examine ownership situation */
        q = -1;
        if (!fu->is_ghost) {
          q = owners[0][i] = me->mpirank;
        }
        else if (fu->quadid >= 0) {
          P4EST_ASSERT (me->ghost != NULL);
          q = owners[0][i] = me->ghost_rank[fu->quadid];
          is_shared[0] = 1;
        }
        if (q >= 0) {
          /* this side face is local or found in ghost layer */
          if (q < me->mpirank) {
            is_owned[0] = 0;
          }
        }
      }
      if (!is_shared[0]) {
        P4EST_ASSERT (is_owned[0]);
        lo = me->num_owned++;
        for (i = 0; i < 2; ++i) {
          if (owners[0][i] >= 0) {
            tree = p4est_tree_array_index (fi->p4est->trees, fss[i]->treeid);
            le = tree->quadrants_offset + fss[i]->is.full.quadid;
            set_lnodes_face_full (ln, le, fss[i]->face, lo);
          }
        }
      }
      else {
        /* the node is shared with one other process */

      }

    }
    return;
  }

  /* this is a hanging face connection */
  nunodes = 1 + (me->with_faces ? 2 : 0);
  codim[0] = P4EST_DIM;
  if (me->with_faces) {
    codim[1] = codim[2] = 1;
  }
  for (j = 0, i = 0; i < 2; ++i) {
    fs = (p4est_iter_face_side_t *) sc_array_index_int (&fi->sides, i);
    if (!fs->is_hanging) {
      /* add midface corner and possibly two half face nodes */
      fu = &fs->is.full;
      q = -1;
      if (!fu->is_ghost) {
        q = owners[0][j] = me->mpirank;
      }
      else if (fu->quadid >= 0) {
        P4EST_ASSERT (me->ghost != NULL);
        q = owners[0][j] = me->ghost_rank[fu->quadid];
        is_shared[0] = 1;
      }
      if (q >= 0) {
        if (me->with_faces) {
          owners[1][j] = owners[2][j] = q;
        }
        if (me->with_faces) {

        }
        for (j = 0; j < nunodes; ++j) {
          if (q < me->mpirank) {
            is_owned[j] = 0;
          }
        }

      }

    }

  }
}

static void
iter_corner1 (p4est_iter_corner_info_t * ci, void *user_data)
{
}

p4est_trimesh_t    *
p4est_trimesh_new (p4est_t * p4est, p4est_ghost_t * ghost, int with_faces)
{
  int                 mpiret;
  int                 p, q, s;
  int                 vn;
  p4est_locidx_t      le, lg, ng;
  p4est_gloidx_t      gc;
  p4est_trimesh_t    *tm;
  p4est_lnodes_t     *ln;
  trimesh_meta_t      tmeta, *me = &tmeta;

  P4EST_ASSERT (p4est_is_balanced (p4est, P4EST_CONNECT_FACE));

  /* basic assignment of members */
  memset (me, 0, sizeof (trimesh_meta_t));
  me->p4est = p4est;
  me->with_faces = with_faces;
  s = me->mpisize = p4est->mpisize;
  p = me->mpirank = p4est->mpirank;
  tm = me->tm = P4EST_ALLOC_ZERO (p4est_trimesh_t, 1);
  ln = tm->lnodes = P4EST_ALLOC_ZERO (p4est_lnodes_t, 1);

  /* lookup structure for ghost owner rank */
  if ((me->ghost = ghost) != NULL) {
    P4EST_ASSERT (ghost->proc_offsets[0] == 0);
    P4EST_ASSERT (ghost->proc_offsets[s] ==
                  (p4est_locidx_t) ghost->ghosts.elem_count);
    me->ghost_rank = P4EST_ALLOC (int, ghost->ghosts.elem_count);
    lg = 0;
    for (q = 0; q < s; ++q) {
      ng = ghost->proc_offsets[q + 1];
      for (; lg < ng; ++lg) {
        me->ghost_rank[lg] = q;
      }
    }
    P4EST_ASSERT (lg == (p4est_locidx_t) ghost->ghosts.elem_count);
  }

  /* prepare node information */
  ln->mpicomm = p4est->mpicomm;
  ln->sharers = sc_array_new (sizeof (p4est_lnodes_rank_t));
  ln->degree = 0;
  vn = ln->vnodes = 9 + (with_faces ? 16 : 0);
  le = ln->num_local_elements = p4est->local_num_quadrants;
  ln->face_code = P4EST_ALLOC_ZERO (p4est_lnodes_code_t, le);
  ln->element_nodes = P4EST_ALLOC_ZERO (p4est_locidx_t, le * vn);

  /* determine node count and ownership */
  me->lenum = 0;
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

  /* free memory */
  P4EST_FREE (me->ghost_rank);

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
