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
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#else
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#endif

#ifndef P4_TO_P8
#define P4EST_LN_C_OFFSET 4
#else
#define P8EST_LN_E_OFFSET 6
#define P4EST_LN_C_OFFSET 18
#endif

static int
p4est_locidx_offset_compare (const void *key, const void *elem)
{
  const p4est_locidx_t *start = (p4est_locidx_t *) elem;
  const p4est_locidx_t *stop = (p4est_locidx_t *) start + 1;
  int                 comp = p4est_locidx_compare (key, start);
  if (comp < 0) {
    return -1;
  }
  comp = p4est_locidx_compare (key, stop);
  if (comp >= 0) {
    return 1;
  }
  return 0;
}

/** dep: dependent quads and processes.
 * Suppose quadrants q0, q1, q2, and q3 share the same face neighbor p.
 *
 *       _______ o       +____________
 *      /       /        /            /
 *     /______ /| q3    /            /|
 *    /       / |      /            / |
 *   /______ /| /     /___________ /  |
 *   |      | |/| q1  |           |   |
 *   |  q2  | / |     |           |   |
 *   |______|/| /     |     p     |   /
 *   |      | |/      |           |  /
 *   |  q0  | /       |           | /
 *   |______|/        |___________|/
 *
 * Even though quads q0, q1, and q2 do not touch the corner marked by "+",
 * they share the node that is created there.  When the callback that creates
 * that node is run, it needs to know about all quads that share the node.
 * Quadrant q3 has one dep: during the face callback between face p and the q
 * faces, the dep associated with the marked corner "o" will have one of its
 * face values record the presence of q0, on edge value record the presence of
 * q1, and one edge value record the presence of q2.  If the quads are local,
 * their indices are stored; if they are ghosts, their owner is encoded as
 * -(proc + 3).  This allows -1 to stand for "uninitialized", and -2 to stand
 *  for "no dependencies".
 */

/* new idea: quad index */
typedef struct p4est_lnodes_dep
{
  p4est_locidx_t      face[P4EST_DIM];
#ifdef P4_TO_P8
  p4est_locidx_t      edge[P4EST_DIM];
#endif
}
p4est_lnodes_dep_t;

/** buf_info: encodes/decodes the transmission of node information.
 * share_offset and share_count index into the inode_sharers array
 * for a list of all processes that share the nodes (local process included).
 */
typedef struct p4est_lnodes_buf_info
{
  int8_t              type;     /* which nodes it shares */
  int8_t              send_sharers;     /* whether the sharers are included in
                                           the message */
  p4est_locidx_t      first_index;      /* inodes array, first node to/from */
  p4est_locidx_t      share_offset;
  int8_t              share_count;
}
p4est_lnodes_buf_info_t;

typedef struct p4est_lnodes_data
{
  p4est_lnodes_dep_t *local_dep;        /* num local quads */
  p4est_lnodes_dep_t *ghost_dep;        /* num ghost quads */
  p4est_locidx_t     *local_elem_nodes; /* num local quads * nodes per q */
  p4est_locidx_t     *poff;     /* mpisize + 1 */
  sc_array_t         *inodes;   /* 2 * p4est_locidx_t */
  sc_array_t         *inode_sharers;    /* int */
  sc_array_t         *send_buf_info;    /* one for each proc: type buf_info_t */
  sc_array_t         *recv_buf_info;    /* one for each proc: type buf_info_t */
  p4est_lnodes_code_t *face_codes;
  int                 nodes_per_elem;
  int                 nodes_per_volume;
  int                *volume_nodes;
  int                 nodes_per_face;
  int                *face_nodes[P4EST_FACES];
#ifdef P4_TO_P8
  int                 nodes_per_edge;
  int                *edge_nodes[P8EST_EDGES];
#endif
  int                 nodes_per_corner;
  int                *corner_nodes[P4EST_CHILDREN];
  sc_array_t          send_requests;
  sc_array_t         *send_buf;
  sc_array_t         *touching_procs;
  sc_array_t         *all_procs;
}
p4est_lnodes_data_t;

static inline int
fside_get_fields (p4est_iter_face_side_t * fside, int *is_hanging,
                  p4est_topidx_t * tid, int *f, int8_t ** is_ghost,
                  p4est_locidx_t ** quadid, p4est_quadrant_t *** quad)
{
  int                 limit;

  *is_hanging = fside->is_hanging;
  *tid = fside->treeid;
  *f = (int) fside->face;
  if (fside->is_hanging) {
    limit = P4EST_HALF;
    *is_ghost = fside->is.hanging.is_ghost;
    *quadid = fside->is.hanging.quadid;
    *quad = fside->is.hanging.quad;
  }
  else {
    limit = 1;
    *is_ghost = &fside->is.full.is_ghost;
    *quadid = &fside->is.full.quadid;
    *quad = &fside->is.full.quad;
  }

  return limit;
}

/** lnodes_face_simple_callback: runs even if there are no face nodes.  If a
 * side of the face is not hanging, then there are no other quadrants that are
 * facewise dependent on its corner or edge nodes.  If a side of the face is
 * hanging, we set up the facewise dep values.  We also update the face_codes
 * for local quads.  Store a list of all touching processors */
static void
p4est_lnodes_face_simple_callback (p4est_iter_face_info_t * info, void *Data)
{
  int                 i, f, fdir, limit, cid, xind, *ip;
  sc_array_t         *sides = &(info->sides);
  size_t              zz, count = sides->elem_count;
  p4est_lnodes_data_t *data = (p4est_lnodes_data_t *) Data;
  p4est_iter_face_side_t *fside;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_topidx_t      tid;
  sc_array_t         *touching_procs = data->touching_procs;
  sc_array_t          proc_offsets;
  p4est_lnodes_dep_t *local_dep = data->local_dep;
  p4est_lnodes_dep_t *ghost_dep = data->ghost_dep;
  p4est_lnodes_dep_t *dep;
  p4est_locidx_t      quadrants_offset;
  int                 rank = info->p4est->mpirank;
  p4est_lnodes_code_t *face_codes = data->face_codes;
  int8_t             *is_ghost;
  p4est_locidx_t     *quadid;
  p4est_quadrant_t  **quad;
  int                 is_hanging;
  int                 procs[P4EST_HALF];
  p4est_locidx_t      qid[P4EST_HALF];
#ifdef P4_TO_P8
  int                 j;
#endif
  int                 k, c;
  p4est_quadrant_t    tempq;

  P4EST_ASSERT (touching_procs->elem_size == sizeof (int));
  sc_array_truncate (touching_procs);
  /* even though the original is size mpisize+1, proc_offsets uses
   * p4est_locidx_offset_compare, and we don't want to read past the end of the
   * array */
  sc_array_init_data (&proc_offsets, info->ghost_layer->proc_offsets,
                      sizeof (p4est_locidx_t), (size_t) info->p4est->mpisize);

  for (zz = 0; zz < count; zz++) {
    fside = p4est_iter_fside_array_index (sides, zz);
    limit = fside_get_fields (fside, &is_hanging, &tid, &f, &is_ghost,
                              &quadid, &quad);
    fdir = f / 2;
    tree = p4est_tree_array_index (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    k = -1;
    for (i = 0; i < limit; i++) {
      qid[i] = quadid[i];
      if (qid[i] < 0) {
        P4EST_ASSERT (limit == P4EST_HALF);
        if (k < 0) {
          for (k = 0; k < P4EST_HALF; k++) {
            if (quadid[k] >= 0) {
              P4EST_ASSERT (quad[k]);
              break;
            }
          }
        }
        P4EST_ASSERT (k >= 0 && k < P4EST_HALF);
        c = p4est_face_corners[f][i];
        p4est_quadrant_sibling (quad[k], &tempq, c);
        procs[i] = p4est_comm_find_owner (info->p4est, tid,
                                          &tempq, info->p4est->mpirank);
        P4EST_ASSERT (procs[i] >= 0 && procs[i] != rank);
        ip = (int *) sc_array_push (touching_procs);
        *ip = procs[i];
      }
      else if (is_ghost[i]) {
        procs[i] = (int) sc_array_bsearch (&proc_offsets, &(qid[i]),
                                           p4est_locidx_offset_compare);
        P4EST_ASSERT (procs[i] >= 0 && procs[i] != rank);
        ip = (int *) sc_array_push (touching_procs);
        *ip = procs[i];
      }
      else {
        qid[i] += quadrants_offset;
        procs[i] = rank;
        /* update face code */
        if (is_hanging) {
          face_codes[qid[i]] |= ((p4est_lnodes_code_t)
                                 p4est_face_corners[f][i]);
          face_codes[qid[i]] |= ((p4est_lnodes_code_t)
                                 1 << (P4EST_DIM + f / 2));
        }
      }
    }
    if (!data->nodes_per_corner &&
#ifdef P4_TO_P8
        !data->nodes_per_edge &&
#endif
        1) {
      continue;
    }
    for (i = 0; i < limit; i++) {
      dep = !is_ghost[i] ? &(local_dep[qid[i]]) : &(ghost_dep[qid[i]]);
      if (is_hanging) {
#ifdef P4_TO_P8
        int                 ndir[2];

        ndir[0] = SC_MIN (((fdir + 1) % 3), ((fdir + 2) % 3));
        ndir[1] = SC_MAX (((fdir + 1) % 3), ((fdir + 2) % 3));

        for (j = 0; j < 2; j++) {
          xind = i ^ (j + 1);
          if (is_ghost[xind]) {
            dep->edge[ndir[j]] = -((p4est_locidx_t) procs[xind] + 3);
          }
          else {
            dep->edge[ndir[j]] = qid[xind];
          }
        }
#endif
        xind = i ^ (P4EST_HALF - 1);
        if (is_ghost[xind]) {
          dep->face[fdir] = -((p4est_locidx_t) procs[xind] + 3);
        }
        else {
          dep->face[fdir] = qid[xind];
        }
      }
      else {
        cid = p4est_quadrant_child_id (quad[i]);
        if (p4est_corner_face_corners[cid][f] >= 0) {
          dep->face[fdir] = -2;
        }
      }
    }
  }
}

#ifdef P4_TO_P8

static inline int
eside_get_fields (p8est_iter_edge_side_t * eside, int *is_hanging,
                  p4est_topidx_t * tid, int *e, int *o, int8_t ** is_ghost,
                  p4est_locidx_t ** quadid, p4est_quadrant_t *** quad)
{
  int                 limit;

  *is_hanging = eside->is_hanging;
  *tid = eside->treeid;
  *e = (int) eside->edge;
  *o = (int) eside->orientation;
  if (eside->is_hanging) {
    limit = 2;
    *is_ghost = eside->is.hanging.is_ghost;
    *quadid = eside->is.hanging.quadid;
    *quad = eside->is.hanging.quad;
  }
  else {
    limit = 1;
    *is_ghost = &eside->is.full.is_ghost;
    *quadid = &eside->is.full.quadid;
    *quad = &eside->is.full.quad;
  }

  return limit;
}

/** lnodes_edge_simple_callback: runs even if there are no edge nodes.  If a
 * side of the face is not hanging, then there are no other quadrants that are
 * facewise dependent on its corner or edge nodes.  If a side of the face is
 * hanging, we set up the facewise dep values.  We also update the face_codes
 * for local quads.  Store a list of all touching procs.  return whether there
 * is a local touching quadrant.
 */
static int
p8est_lnodes_edge_simple_callback (p8est_iter_edge_info_t * info, void *Data)
{
  int                 i, limit, e, edir, cid, *ip;
  sc_array_t         *sides = &(info->sides);
  size_t              zz, count = sides->elem_count;
  p4est_lnodes_data_t *data = (p4est_lnodes_data_t *) Data;
  p8est_iter_edge_side_t *eside;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_topidx_t      tid;
  sc_array_t          proc_offsets;
  sc_array_t         *touching_procs = data->touching_procs;
  p4est_lnodes_dep_t *local_dep = data->local_dep;
  p4est_lnodes_dep_t *ghost_dep = data->ghost_dep;
  p4est_lnodes_dep_t *dep;
  int8_t             *is_ghost;
  int                 procs[2];
  int                 rank = info->p4est->mpirank;
  p4est_locidx_t     *quadid;
  p4est_quadrant_t  **quad;
  p4est_locidx_t      qid[2];
  p4est_locidx_t      quadrants_offset;
  p4est_lnodes_code_t *face_codes = data->face_codes;
  int                 is_hanging, o, has_local = 0, c;
  p4est_quadrant_t    tempq;

  P4EST_ASSERT (touching_procs->elem_size == sizeof (int));
  sc_array_truncate (touching_procs);
  /* even though the original is size mpisize+1, proc_offsets uses
   * p4est_locidx_offset_compare, and we don't want to read past the end of the
   * array */
  sc_array_init_data (&proc_offsets, info->ghost_layer->proc_offsets,
                      sizeof (p4est_locidx_t), (size_t) info->p4est->mpisize);

  for (zz = 0; zz < count; zz++) {
    eside = p8est_iter_eside_array_index (sides, zz);
    limit = eside_get_fields (eside, &is_hanging, &tid, &e, &o, &is_ghost,
                              &quadid, &quad);
    edir = e / 4;
    tree = p4est_tree_array_index (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    for (i = 0; i < limit; i++) {
      qid[i] = quadid[i];
      if (qid[i] < 0) {
        if (limit == 2 && quadid[i ^ 1] >= 0) {
          P4EST_ASSERT (quad[i ^ 1]);
          c = p8est_edge_corners[e][i];
          p4est_quadrant_sibling (quad[i ^ 1], &tempq, c);
          procs[i] = p4est_comm_find_owner (info->p4est, tid,
                                            &tempq, info->p4est->mpirank);
          P4EST_ASSERT (procs[i] >= 0 && procs[i] != rank);
          ip = (int *) sc_array_push (touching_procs);
          *ip = procs[i];
        }
      }
      else if (is_ghost[i]) {
        procs[i] = (int) sc_array_bsearch (&proc_offsets, &(qid[i]),
                                           p4est_locidx_offset_compare);
        P4EST_ASSERT (procs[i] >= 0 && procs[i] != rank);
        ip = (int *) sc_array_push (touching_procs);
        *ip = procs[i];
      }
      else {
        has_local = 1;
        qid[i] += quadrants_offset;
        procs[i] = rank;
        if (is_hanging) {
          /* update face code */
          face_codes[qid[i]] |=
            ((p4est_lnodes_code_t) p8est_edge_corners[e][i]);
          face_codes[qid[i]] |= ((p4est_lnodes_code_t) 1 << (6 + e / 4));
        }
      }
    }
    for (i = 0; i < limit; i++) {
      if (qid[i] < 0) {
        continue;
      }
      dep = !is_ghost[i] ? &(local_dep[qid[i]]) : &(ghost_dep[qid[i]]);
      if (is_hanging) {
        if (!has_local && qid[i ^ 1] < 0) {
          dep->edge[edir] = -1;
        }
        else if (is_ghost[i ^ 1]) {
          dep->edge[edir] = -((p4est_locidx_t) procs[i ^ 1] + 3);
        }
        else {
          dep->edge[edir] = qid[i ^ 1];
        }
      }
      else {
        cid = p4est_quadrant_child_id (quad[i]);
        if (p8est_edge_corners[e][0] == cid ||
            p8est_edge_corners[e][1] == cid) {
          dep->edge[edir] = -2;
        }
      }
    }
  }

  return has_local;
}

static void
p8est_lnodes_edge_simple_callback_void (p8est_iter_edge_info_t * info,
                                        void *Data)
{
  (void) p8est_lnodes_edge_simple_callback (info, Data);
}
#endif

static inline void
cside_get_fields (p4est_iter_corner_side_t * cside,
                  p4est_topidx_t * tid, int *c, int8_t * is_ghost,
                  p4est_locidx_t * quadid, p4est_quadrant_t ** quad)
{
  *tid = cside->treeid;
  *c = cside->corner;
  *is_ghost = cside->is_ghost;
  *quadid = cside->quadid;
  *quad = cside->quad;
}

/** Once we have found the quadrant (\a q, \a tid, \a type) that owns a set of
 * nodes, push the info describing the owner quadrant on the appropriate
 * send/recv lists.
 */
static inline void
p4est_lnodes_push_binfo (sc_array_t * touch, sc_array_t * all,
                         sc_array_t * send, sc_array_t * recv,
                         sc_array_t * share, int owner, int rank,
                         int mpisize, int is_remote,
                         int8_t type, p4est_locidx_t nin)
{
  size_t              zz, count = all->elem_count;
  int                *ip, proc;
  p4est_lnodes_buf_info_t *binfo;
  int8_t              scount = -1;
  p4est_locidx_t      offset = (p4est_locidx_t) share->elem_count;

  if (!is_remote) {
    ip = (int *) sc_array_push (share);
    *ip = rank;
    scount = (int8_t) (count + 1);
  }
  for (zz = 0; zz < count; zz++) {
    proc = *((int *) sc_array_index (all, zz));
    if (!is_remote) {
      ip = (int *) sc_array_push (share);
      *ip = proc;
    }
    if (owner == rank) {
      P4EST_ASSERT (proc != rank);
      P4EST_ASSERT (!is_remote);
      P4EST_ASSERT (0 <= proc && proc < mpisize);
      binfo = (p4est_lnodes_buf_info_t *) sc_array_push (&(send[proc]));
      binfo->send_sharers = 1;
      if (touch == NULL ||
          sc_array_bsearch (touch, &proc, sc_int_compare) >= 0) {
        binfo->send_sharers = 0;
      }
    }
    else if (proc == owner) {
      P4EST_ASSERT (0 <= proc && proc < mpisize);
      binfo = (p4est_lnodes_buf_info_t *) sc_array_push (&(recv[proc]));
      if (!is_remote) {
        binfo->send_sharers = 0;
      }
      else {
        binfo->send_sharers = 1;
      }
    }
    else {
      continue;
    }
    binfo->type = type;
    binfo->first_index = nin;
    if (!is_remote) {
      binfo->share_offset = offset;
      binfo->share_count = scount;
    }
    else {
      binfo->share_offset = -1;
      binfo->share_count = -1;
    }
  }
}

/** p4est_lnodes_missing_proc_corner: figure out processors that may share a
 * corner node remotely
 */
static int
p4est_lnodes_missing_proc_corner (p4est_iter_corner_info_t * info, int side,
                                  int b)
{
  sc_array_t         *sides = &(info->sides);
  int                 i, nsides = (int) sides->elem_count;
  p4est_iter_corner_side_t *thisside = p4est_iter_cside_array_index_int
    (sides, side);
  p4est_iter_corner_side_t *cside;
  p4est_quadrant_t   *q, tempq;
  int                 l, j, f, fc, c2, nproc, key, test;
  int                 c = thisside->corner;
#ifdef P4_TO_P8
  int                 e;
#endif

  q = thisside->quad;
  P4EST_ASSERT (q != NULL);
  l = q->level;
  if (l == 0) {
    return -1;
  }

  if (b < P4EST_DIM) {
    key = thisside->faces[b];
    f = p4est_corner_faces[c][b];
    fc = p4est_corner_face_corners[c][f];
    c2 = p4est_face_corners[f][fc ^ (P4EST_HALF - 1)];
  }
#ifndef P4_TO_P8
  else {
    key = -1;
    c2 = -1;
    SC_ABORT_NOT_REACHED ();
  }
#else
  else {
    key = thisside->edges[b - 3];
    e = p8est_corner_edges[c][b - 3];
    if (p8est_edge_corners[e][0] == c) {
      c2 = p8est_edge_corners[e][1];
    }
    else {
      c2 = p8est_edge_corners[e][0];
    }
  }
#endif
  p4est_quadrant_sibling (q, &tempq, c2);

  for (i = 0; i < nsides; i++) {
    if (i == side) {
      continue;
    }
    cside = p4est_iter_cside_array_index_int (sides, i);
    for (j = 0; j < P4EST_DIM; j++) {
      test = cside->faces[j];
#ifdef P4_TO_P8
      if (b >= P4EST_DIM) {
        test = cside->edges[j];
      }
#endif
      if (test == key) {
        P4EST_ASSERT (cside->quad != NULL);
        if (cside->quad->level < l) {
          P4EST_ASSERT (cside->quad->level == l - 1);
          nproc =
            p4est_comm_find_owner (info->p4est, thisside->treeid, &tempq,
                                   info->p4est->mpirank);
          P4EST_ASSERT (nproc >= 0);
          return nproc;
        }
      }
    }
  }
  return -1;
}

/* p4est_lnodes_corner_callback:
 *
 * Create a new independent node at a corner.  The quad in the first side is
 * the owner: determine the owning proc.  If the node isn't remote, compute
 * all processes that share the node.  For every local quad that shares the
 * node, point local_elem_nodes at the new node.  If the node is locally
 * owned, add info describing the node to the send buffer of all processes
 * that need the node.  If the node is not locally owned, add info describing
 * the node to the receive buffer of the owner.
 */
static void
p4est_lnodes_corner_callback (p4est_iter_corner_info_t * info, void *Data)
{
  int                 i, j, limit;
  sc_array_t         *sides = &(info->sides);
  size_t              zz, count = sides->elem_count;
  p4est_lnodes_data_t *data = (p4est_lnodes_data_t *) Data;
  p4est_iter_corner_side_t *cside, *owner_cside;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *inode, *lp;
  sc_array_t         *inode_sharers = data->inode_sharers;
  p4est_lnodes_dep_t *local_dep = data->local_dep;
  p4est_lnodes_dep_t *ghost_dep = data->ghost_dep;
  p4est_lnodes_dep_t *dep;
  p4est_locidx_t     *local_elem_nodes = data->local_elem_nodes;
  sc_array_t         *send_buf_info = data->send_buf_info;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  sc_array_t         *touching_procs = data->touching_procs;
  sc_array_t         *all_procs = data->all_procs;
  int                *ip;
  p4est_topidx_t      tid, owner_tid;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_ghost_t      *ghost_layer = info->ghost_layer;
  sc_array_t          proc_offsets;
  p4est_locidx_t      qid, owner_qid, nqid;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) inodes->elem_count;
  int                 npc = data->nodes_per_corner;
  int8_t              is_ghost, owner_is_ghost;
  p4est_locidx_t      nid;
  int                 proc, owner_proc, nproc;
  int                 rank = info->p4est->mpirank;
  int                 c, owner_c;
  p4est_quadrant_t    ownq, tempq, tempr;
  p4est_quadrant_t   *q, *owner_q;
  int               **corner_nodes = data->corner_nodes;
  int                 nodes_per_elem = data->nodes_per_elem;
  p4est_locidx_t      quadrants_offset;
  int8_t              type;
  p4est_connectivity_t *conn = info->p4est->connectivity;
  int                 is_remote;
  int                 has_local;

  /* even though the original is size mpisize+1, proc_offsets uses
   * p4est_locidx_offset_compare, and we don't want to read past the end of the
   * array */
  sc_array_init_data (&proc_offsets, ghost_layer->proc_offsets,
                      sizeof (p4est_locidx_t), (size_t) info->p4est->mpisize);

  sc_array_truncate (touching_procs);
  sc_array_truncate (all_procs);

  /* figure out which proc owns the node */
  owner_cside = p4est_iter_cside_array_index (sides, 0);
  cside_get_fields (owner_cside, &owner_tid, &owner_c, &owner_is_ghost,
                    &owner_qid, &owner_q);
  if (owner_q == NULL) {
    p4est_qcoord_t      x, y, h;
#ifdef P4_TO_P8
    p4est_qcoord_t      l, z;
#endif
    /* if this is a remote node, we don't have a quad available for
     * determining ownership, so we have to create it */
    P4EST_ASSERT (count > 1);
    cside = NULL;
    for (zz = 1; zz < count; zz++) {
      /* find a nonempty side */
      cside = p4est_iter_cside_array_index (sides, zz);
      if (cside->quad) {
        break;
      }
    }
    P4EST_ASSERT (zz < count);
    cside_get_fields (cside, &tid, &c, &is_ghost, &qid, &q);
    p4est_quadrant_corner_descendant (q, &tempr, c, P4EST_QMAXLEVEL);
    q = &tempr;
    P4EST_ASSERT (p4est_quadrant_child_id (q) == c);
    /* we want the coordinates of the common corner as the owning process sees
     * them.  transform the quad across the corner: the status of the transformed
     * quad tells us how to proceed */
    p4est_quadrant_corner_neighbor (q, c, &tempq);
    if (p4est_quadrant_is_inside_root (&tempq)) {
      /* inside the root, only one set of coordinates */
      h = P4EST_QUADRANT_LEN (q->level);
      x = q->x + h * (c & 1);
      y = q->y + h * ((c & 2) >> 1);
#ifdef P4_TO_P8
      z = q->z + h * ((c & 4) >> 2);
#endif
    }
    else if (p4est_quadrant_is_outside_corner (&tempq)) {
      /* outside a corner, trivially set the coordinates to the appropriate
       * corner */
      h = P4EST_QUADRANT_LEN (0);
      x = h * (owner_c & 1);
      y = h * ((owner_c & 2) >> 1);
#ifdef P4_TO_P8
      z = h * ((owner_c & 4) >> 2);
#endif
    }
#ifdef P4_TO_P8
    else if (p8est_quadrant_is_outside_edge (&tempq)) {
      /* outside an edge: use some knowledge about how p4est_iterate orders
       * the sides around a corner that is in the middle of an edge */
      int                 owner_e, owner_c2, e, c2, this_c, this_c2;

      P4EST_ASSERT (count % 2 == 0);
      /* side[count/2] should be on the same edge as side[0] */
      cside = p4est_iter_cside_array_index (sides, count / 2);
      P4EST_ASSERT (cside->treeid == owner_tid);
      P4EST_ASSERT (p8est_child_corner_edges[owner_c][cside->corner] >= 0);
      /* we now have two corners, which determines the edge from the owner's
       * point of view */
      owner_c2 = cside->corner;
      owner_e = p8est_child_corner_edges[owner_c][owner_c2];
      /* side[zz] is on the same edge as side[zz +- count / 2] */
      if (zz < count / 2) {
        cside = p4est_iter_cside_array_index (sides, zz + count / 2);
        /* the order of this_c and this_c2 reflects the orientation of the
         * edge that the corners touch */
        this_c = c;
        this_c2 = cside->corner;
      }
      else {
        cside = p4est_iter_cside_array_index (sides, zz - count / 2);
        /* the order of this_c and this_c2 reflects the orientation of the
         * edge that the corners touch */
        this_c = cside->corner;
        this_c2 = c;
      }
      /* we now have two corners, which determines the edge from zz's point of
       * view */
      P4EST_ASSERT (cside->treeid == tid);
      c2 = cside->corner;
      e = p8est_child_corner_edges[c][c2];
      P4EST_ASSERT (e >= 0);
      h = P4EST_QUADRANT_LEN (q->level);
      /* get the coordinate of the corner along the edge from zz's point of
       * view */
      if (e / 4 == 0) {
        l = q->x + h * (c & 1);
      }
      else if (e / 4 == 1) {
        l = q->y + h * ((c & 2) >> 1);
      }
      else {
        l = q->z + h * ((c & 4) >> 2);
      }
      /* if the two edges are oppositely oriented, get the complement
       * coordinate */
      if ((owner_c > owner_c2) != (this_c > this_c2)) {
        l = P4EST_ROOT_LEN - l;
      }
      /* we combine the knowledge about which edge we are looking for with the
       * coordinate distance along the edge into the coordinates of the corner
       */
      h = P4EST_QUADRANT_LEN (0);
      if (owner_e / 4 == 0) {
        x = l;
        y = h * (owner_e & 1);
        z = h * ((owner_e & 2) >> 1);
      }
      else if (owner_e / 4 == 1) {
        x = h * (owner_e & 1);
        y = l;
        z = h * ((owner_e & 2) >> 1);
      }
      else {
        x = h * (owner_e & 1);
        y = h * ((owner_e & 2) >> 1);
        z = l;
      }
    }
#endif
    else {
      int                 owner_c2, owner_f, c2;
      p4est_topidx_t      nt;
      int                 nf;

      /* this uses some knowledge about how iterate orders the sides of a
       * corner that is in the middle of a face */
      P4EST_ASSERT (count == P4EST_HALF || count == P4EST_CHILDREN);
      /* side[count - (count / P4EST_HALF)] is the same side of the face,
       * opposite corner */
      cside =
        p4est_iter_cside_array_index (sides, count - (count / P4EST_HALF));
      P4EST_ASSERT (cside->treeid == owner_tid);

      owner_c2 = cside->corner;
      /* the two coordinates determine the face */
      owner_f = p4est_child_corner_faces[owner_c][owner_c2];
      P4EST_ASSERT (owner_f >= 0);

      /* figure out which tree is on the other side of the face */
      nt = conn->tree_to_tree[P4EST_FACES * owner_tid + owner_f];
      nf = conn->tree_to_face[P4EST_FACES * owner_tid + owner_f];

      nf %= P4EST_FACES;

      if ((nt == owner_tid && nf == owner_f) || (zz % 2) == 0) {
        /* one-sided face: q must be on the same face: the corner is in the
         * same coordinates */
        h = P4EST_QUADRANT_LEN (q->level);
        x = q->x + h * (c & 1);
        y = q->y + h * ((c & 2) >> 1);
#ifdef P4_TO_P8
        z = q->z + h * ((c & 4) >> 2);
#endif
      }
      else {
        /* transform q across the face */
        p4est_topidx_t      nnt;

        P4EST_ASSERT (nt == tid);
        nnt = p4est_quadrant_face_neighbor_extra (q, tid, nf, &tempq, NULL,
                                                  conn);
        P4EST_ASSERT (nnt == owner_tid);
        c2 = p4est_quadrant_child_id (&tempq);
        h = P4EST_QUADRANT_LEN (tempq.level);
        x = tempq.x + h * (c2 & 1);
        y = tempq.y + h * ((c2 & 2) >> 1);
#ifdef P4_TO_P8
        z = tempq.z + h * ((c2 & 4) >> 2);
#endif
      }
    }
    /* turn the coordinates of the corner into the coordinates of a smallest
     * quad that touches the corner from the correct side */
    h = P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
    ownq.x = x - h * (owner_c & 1);
    ownq.y = y - h * ((owner_c & 2) >> 1);
#ifdef P4_TO_P8
    ownq.z = z - h * ((owner_c & 4) >> 2);
#endif
    ownq.level = P4EST_QMAXLEVEL;
    /* find the owner */
    owner_proc = p4est_comm_find_owner (info->p4est, owner_tid, &ownq, rank);
    P4EST_ASSERT (owner_proc >= 0 && owner_proc != rank);
  }
  else {
    if (owner_is_ghost) {
      owner_proc = (int) sc_array_bsearch (&proc_offsets, &(owner_qid),
                                           p4est_locidx_offset_compare);
      P4EST_ASSERT (owner_proc >= 0 && owner_proc != rank);
    }
    else {
      tree = p4est_tree_array_index (trees, owner_tid);
      quadrants_offset = tree->quadrants_offset;
      owner_qid += quadrants_offset;
      owner_proc = rank;
    }
  }
  /* create the new node */
  for (j = 0; j < npc; j++) {
    inode = (p4est_locidx_t *) sc_array_push (inodes);
    inode[0] = owner_proc;
    inode[1] = owner_qid;
  }

  /* figure out if this is a remote corner or one for which we can determing
   * all touching and sharing procs */
  has_local = 0;
  for (zz = 0; zz < count; zz++) {
    cside = p4est_iter_cside_array_index (sides, zz);
    if (!cside->is_ghost) {
      has_local = 1;
    }
  }
  is_remote = !has_local;
  if (is_remote) {
    ip = (int *) sc_array_push (all_procs);
    *ip = owner_proc;
  }

  for (zz = 0; zz < count; zz++) {
    cside = p4est_iter_cside_array_index (sides, zz);
    cside_get_fields (cside, &tid, &c, &is_ghost, &qid, &q);
    if (q == NULL) {
      P4EST_ASSERT (is_ghost);
      continue;
    }

    tree = p4est_tree_array_index (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    if (!is_ghost) {
      proc = rank;
      qid += quadrants_offset;
      for (j = 0; j < npc; j++) {
        nid = qid * nodes_per_elem + corner_nodes[c][j];
        P4EST_ASSERT (local_elem_nodes[nid] == -1);
        local_elem_nodes[nid] = num_inodes + j;
      }
    }
    else if (!is_remote) {
      P4EST_ASSERT (qid >= 0);
      proc = (int) sc_array_bsearch (&proc_offsets, &qid,
                                     p4est_locidx_offset_compare);
      P4EST_ASSERT (proc >= 0 && proc != rank);
      ip = (int *) sc_array_push (touching_procs);
      *ip = proc;
      ip = (int *) sc_array_push (all_procs);
      *ip = proc;
    }
    else {
      proc = -1;
    }
    if (p4est_quadrant_child_id (q) != c) {
      /* there can be no remote quads / processes */
      continue;
    }
    P4EST_ASSERT (qid >= 0);
    if (is_ghost) {
      P4EST_ASSERT ((size_t) qid < info->ghost_layer->ghosts.elem_count);
    }
    else {
      P4EST_ASSERT (qid < info->p4est->local_num_quadrants);
    }
    dep = !is_ghost ? &(local_dep[qid]) : &(ghost_dep[qid]);
#ifndef P4_TO_P8
    limit = P4EST_DIM;
#else
    limit = 2 * P4EST_DIM;
#endif
    for (i = 0; i < limit; i++) {
#ifndef P4_TO_P8
      lp = &dep->face[i];
#else
      lp = (i < P4EST_DIM) ? &dep->face[i] : &dep->edge[i - 3];
#endif
      nqid = *lp;
      if (nqid >= 0) {
        has_local = 1;
        /* remote local quad */
        for (j = 0; j < npc; j++) {
          nid = nqid * nodes_per_elem + corner_nodes[c][j];
          P4EST_ASSERT (local_elem_nodes[nid] == -1);
          local_elem_nodes[nid] = num_inodes + j;
        }
      }
      else if (!is_remote) {
        nproc = nqid;
        if (nproc == -1) {
          P4EST_ASSERT (is_ghost);
          nproc = p4est_lnodes_missing_proc_corner (info, zz, i);
          P4EST_ASSERT (nproc != rank);
          P4EST_ASSERT (nproc >= -1);
          P4EST_ASSERT (nproc < info->p4est->mpisize);
          *lp = -((p4est_locidx_t) nproc + 3);
        }
        else {
          nproc = -(nproc + 3);
        }
        P4EST_ASSERT (nproc >= -1 && nproc != rank);
        if (nproc >= 0) {
          P4EST_ASSERT (nproc != rank);
          ip = (int *) sc_array_push (all_procs);
          *ip = nproc;
        }
      }
    }
  }
  P4EST_ASSERT (has_local);
  sc_array_sort (touching_procs, sc_int_compare);
  sc_array_uniq (touching_procs, sc_int_compare);
  sc_array_sort (all_procs, sc_int_compare);
  sc_array_uniq (all_procs, sc_int_compare);
  count = all_procs->elem_count;
  if (count) {
    type = (int8_t) (P4EST_LN_C_OFFSET + owner_c);
    p4est_lnodes_push_binfo (touching_procs, all_procs, send_buf_info,
                             recv_buf_info, inode_sharers, owner_proc, rank,
                             info->p4est->mpisize, is_remote, type,
                             num_inodes);
  }
  else {
    P4EST_ASSERT (owner_proc == rank);
  }
}

#ifdef P4_TO_P8
/** p8est_lnodes_missing_proc_edge: figure out processors that may share an
 * edge node remotely
 */
static void
p8est_lnodes_missing_proc_edge (p8est_iter_edge_info_t * info, int side,
                                int b, int *mproc)
{
  sc_array_t         *sides = &(info->sides);
  int                 i, nsides = (int) sides->elem_count;
  p8est_iter_edge_side_t *thisside = p8est_iter_eside_array_index_int
    (sides, side);
  p8est_iter_edge_side_t *eside;
  p4est_quadrant_t   *q, tempq, tempr;
  int                 key, test;
  int                 j;
  int                 e = thisside->edge, f;
  int                 edir = e / 4;
  int                 missdir = 3 - edir - b;
  int                 c, c2;

  P4EST_ASSERT (edir != b);
  P4EST_ASSERT (thisside->is_hanging);
  P4EST_ASSERT (p8est_edge_faces[e][b < missdir ? 0 : 1] / 2 == b);
  q = thisside->is.hanging.quad[0];
  if (!q) {
    q = thisside->is.hanging.quad[1];
    P4EST_ASSERT (q);
  }
  key = thisside->faces[b < missdir ? 0 : 1];
  c = p8est_edge_corners[e][0];
  f = p8est_corner_faces[c][b];
  c = p8est_corner_face_corners[c][f];
  c = p8est_face_corners[f][c ^ 3];
  c2 = p8est_edge_corners[e][1];
  c2 = p8est_corner_face_corners[c2][f];
  c2 = p8est_face_corners[f][c2 ^ 3];
  p4est_quadrant_sibling (q, &tempq, c);
  p4est_quadrant_sibling (q, &tempr, c2);
  for (i = 0; i < nsides; i++) {
    if (i == side) {
      continue;
    }
    eside = p8est_iter_eside_array_index_int (sides, i);
    for (j = 0; j < 2; j++) {
      test = eside->faces[j];
      if (test == key) {
        if (!eside->is_hanging && eside->is.full.quad != NULL) {
          mproc[0] =
            p4est_comm_find_owner (info->p4est, thisside->treeid, &tempq,
                                   info->p4est->mpirank);
          P4EST_ASSERT (mproc[0] >= 0);
          mproc[1] =
            p4est_comm_find_owner (info->p4est, thisside->treeid, &tempr,
                                   mproc[0]);
          P4EST_ASSERT (mproc[1] >= 0);
          return;
        }
        else {
          mproc[0] = -1;
          mproc[1] = -1;
          return;
        }
      }
    }
  }
  mproc[0] = -1;
  mproc[1] = -1;
}

/* p8est_lnodes_edge_callback:
 *
 * Create new independent nodes on an edge.
 * Set all touching element nodes to point to the newly created independent
 * nodes.
 * Compute all processes that share the nodes.
 * If the nodes are locally owned, add info describing the nodes to the send
 * buffer of all processes that share the nodes.
 * If the nodes are not locally owned, add info describing the nodes to the
 * receive buffer of the owner.
 */
static void
p8est_lnodes_edge_callback (p8est_iter_edge_info_t * info, void *Data)
{
  int                 i, j, k, xdir[2];
  sc_array_t         *sides = &(info->sides);
  size_t              zz, count = sides->elem_count;
  p4est_lnodes_data_t *data = (p4est_lnodes_data_t *) Data;
  p8est_iter_edge_side_t *eside, *owner_eside;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *inode;
  sc_array_t         *inode_sharers = data->inode_sharers;
  p4est_lnodes_dep_t *local_dep = data->local_dep;
  p4est_lnodes_dep_t *ghost_dep = data->ghost_dep;
  p4est_lnodes_dep_t *dep;
  p4est_locidx_t     *local_elem_nodes = data->local_elem_nodes;
  sc_array_t         *send_buf_info = data->send_buf_info;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  sc_array_t         *touching_procs = data->touching_procs;
  sc_array_t         *all_procs = data->all_procs;
  int                *ip;
  p4est_topidx_t      tid, owner_tid;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_locidx_t      quadrants_offset;
  p4est_locidx_t     *qids;
  p4est_locidx_t      qid, owner_qid, nqid;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) inodes->elem_count;
  int8_t             *is_ghost, owner_is_ghost;
  int                 e, edir, owner_e, owner_c, o;
  p4est_locidx_t      nid;
  int                 owner_proc, nproc;
  int                 rank = info->p4est->mpirank;
  p4est_quadrant_t   *owner_q = NULL;
  p4est_quadrant_t   *q;
  p4est_quadrant_t  **quad;
  p4est_quadrant_t    tempq, tempr, ownq;
  int                 nodes_per_edge = data->nodes_per_edge;
  int                 nodes_per_elem = data->nodes_per_elem;
  int               **edge_nodes = data->edge_nodes;
  int                 is_hanging;
  int                 limit;
  int                 stride;
  p4est_locidx_t      start_node;
  int8_t              type;
  int                 is_remote, has_local;
  p4est_connectivity_t *conn = info->p4est->connectivity;
  int                 mproc[2][2] = { {-1, -1}, {-1, -1} };

  sc_array_truncate (touching_procs);
  sc_array_truncate (all_procs);
  has_local = p8est_lnodes_edge_simple_callback (info, data);

  owner_eside = p8est_iter_eside_array_index (sides, 0);
  owner_e = owner_eside->edge;
  owner_tid = owner_eside->treeid;
  if (owner_eside->is_hanging) {
    owner_qid = owner_eside->is.hanging.quadid[0];
    owner_is_ghost = owner_eside->is.hanging.is_ghost[0];
    owner_q = owner_eside->is.hanging.quad[0];
  }
  else {
    owner_qid = owner_eside->is.full.quadid;
    owner_is_ghost = owner_eside->is.full.is_ghost;
    owner_q = owner_eside->is.full.quad;
  }
  P4EST_ASSERT (!owner_eside->orientation);
  owner_c = p8est_edge_corners[owner_e][0];
  if (owner_q == NULL) {
    int                 c;
    p4est_qcoord_t      x, y, z, h, l;

    P4EST_ASSERT (count > 1);
    eside = NULL;
    for (zz = 1; zz < count; zz++) {
      eside = p8est_iter_eside_array_index (sides, zz);
      if ((!eside->is_hanging) && eside->is.full.quad) {
        break;
      }
    }
    P4EST_ASSERT (zz < count);
    e = eside->edge;
    tid = eside->treeid;
    o = eside->orientation;
    c = p8est_edge_corners[e][o];
    q = eside->is.full.quad;
    P4EST_ASSERT (q->level < P4EST_QMAXLEVEL);
    p4est_quadrant_corner_descendant (q, &tempr, c, q->level + 1);
    q = &tempr;
    P4EST_ASSERT (c == p4est_quadrant_child_id (q));
    p8est_quadrant_edge_neighbor (q, e, &tempq);
    /* get the coordinates of the edge */
    if (p4est_quadrant_is_inside_root (&tempq)) {
      h = P4EST_QUADRANT_LEN (q->level);
      x = q->x + h * (c & 1);
      y = q->y + h * ((c & 2) >> 1);
      z = q->z + h * ((c & 4) >> 2);
    }
    else if (p8est_quadrant_is_outside_edge (&tempq)) {
      h = P4EST_QUADRANT_LEN (q->level);
      if (e / 4 == 0) {
        l = q->x + h * (c & 1);
      }
      else if (e / 4 == 1) {
        l = q->y + h * ((c & 2) >> 1);
      }
      else {
        l = q->z + h * ((c & 4) >> 2);
      }
      if (o) {
        l = P4EST_ROOT_LEN - l;
      }
      h = P4EST_QUADRANT_LEN (0);
      if (owner_e / 4 == 0) {
        x = l;
        y = h * (owner_e & 1);
        z = h * ((owner_e & 2) >> 1);
      }
      else if (owner_e / 4 == 1) {
        x = h * (owner_e & 1);
        y = l;
        z = h * ((owner_e & 2) >> 1);
      }
      else {
        x = h * (owner_e & 1);
        y = h * ((owner_e & 2) >> 1);
        z = l;
      }
    }
    else {
      /* outside face */
      int                 owner_f, c1, c2, nf;
      p4est_topidx_t      nt;
      /* this uses some knowledge about how iterate orders the sides of a
       * corner that is in the middle of a face */
      P4EST_ASSERT (count == 2 || count == 4);
      eside = p8est_iter_eside_array_index (sides, count / 2);
      P4EST_ASSERT (eside->treeid == owner_tid);
      P4EST_ASSERT (eside->edge != owner_e);

      c1 = p8est_edge_corners[owner_e][0];
      c2 = p8est_edge_corners[eside->edge][1];
      owner_f = p4est_child_corner_faces[c1][c2];
      P4EST_ASSERT (owner_f >= 0);

      nt = conn->tree_to_tree[P4EST_FACES * owner_tid + owner_f];
      nf = conn->tree_to_face[P4EST_FACES * owner_tid + owner_f];

      /* o2 = nf / P4EST_FACES; */
      nf %= P4EST_FACES;

      if ((nt == owner_tid && nf == owner_f) || (zz % 2) == 0) {
        /* q must be on the same side: the corner is in the same coordinates
         */
        P4EST_ASSERT (!o);
        h = P4EST_QUADRANT_LEN (q->level);
        x = q->x + h * (c & 1);
        y = q->y + h * ((c & 2) >> 1);
        z = q->z + h * ((c & 4) >> 2);
      }
      else {
        int                 c2;
        p4est_topidx_t      nnt;

        P4EST_ASSERT (nt == tid);
        P4EST_ASSERT (count == 4);
        nnt =
          p4est_quadrant_face_neighbor_extra (q, tid, nf, &tempq, NULL, conn);
        P4EST_ASSERT (nnt == owner_tid);
        c2 = p4est_quadrant_child_id (&tempq);
        P4EST_ASSERT (p4est_corner_face_corners[c2][owner_f] >= 0);
        h = P4EST_QUADRANT_LEN (tempq.level);
        x = tempq.x + h * (c2 & 1);
        y = tempq.y + h * ((c2 & 2) >> 1);
        z = tempq.z + h * ((c2 & 4) >> 2);
      }
    }
    h = P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
    ownq.x = x - h * (owner_c & 1);
    ownq.y = y - h * ((owner_c & 2) >> 1);
#ifdef P4_TO_P8
    ownq.z = z - h * ((owner_c & 4) >> 2);
#endif
    ownq.level = P4EST_QMAXLEVEL;
    owner_proc = p4est_comm_find_owner (info->p4est, owner_tid, &ownq, rank);
  }
  else {
    if (owner_is_ghost) {
      owner_proc = *((int *) sc_array_index (touching_procs, 0));
      P4EST_ASSERT (owner_proc >= 0 && owner_proc != rank);
    }
    else {
      owner_proc = rank;
      tree = p4est_tree_array_index (trees, owner_tid);
      quadrants_offset = tree->quadrants_offset;
      owner_qid += quadrants_offset;
    }
  }
  if (has_local) {
    sc_array_sort (touching_procs, sc_int_compare);
    sc_array_uniq (touching_procs, sc_int_compare);
  }
  /* create nodes */
  for (i = 0; i < nodes_per_edge; i++) {
    inode = (p4est_locidx_t *) sc_array_push (inodes);
    P4EST_ASSERT (inodes->elem_count <= (size_t)
                  (nodes_per_elem * info->p4est->local_num_quadrants));
    inode[0] = owner_proc;
    inode[1] = owner_qid;
  }
  /* point element nodes at created nodes; find all sharing procs */
  is_remote = !has_local;
  if (!is_remote) {
    sc_array_copy (all_procs, touching_procs);
  }
  else {
    ip = (int *) sc_array_push (all_procs);
    *ip = owner_proc;
  }
  for (zz = 0; zz < count; zz++) {
    eside = p8est_iter_eside_array_index (sides, zz);
    limit = eside_get_fields (eside, &is_hanging, &tid, &e, &o, &is_ghost,
                              &qids, &quad);
    tree = p4est_tree_array_index (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    if (!is_hanging && quad[0] == NULL) {
      continue;
    }
    mproc[0][0] = -2;
    mproc[0][1] = -2;
    mproc[1][0] = -2;
    mproc[1][1] = -2;
    for (i = 0; i < limit; i++) {
      qid = qids[i];
      if (qid < 0) {
        continue;
      }
      stride = (o ? -1 : 1);
      if (!is_ghost[i]) {
        qid += quadrants_offset;
        P4EST_ASSERT (qid < info->p4est->local_num_quadrants);
        start_node = num_inodes + (o ? nodes_per_edge - 1 : 0);
        for (k = 0; k < nodes_per_edge; k++, start_node += stride) {
          nid = qid * nodes_per_elem + edge_nodes[e][k];
          P4EST_ASSERT (local_elem_nodes[nid] == -1);
          local_elem_nodes[nid] = start_node;
        }
      }
      if (!is_hanging) {
        continue;
      }
      if (is_remote && qids[i ^ 1] < 0) {
        continue;
      }
      /* get quads that may be dependent because of hanging faces */
      dep = !is_ghost[i] ? &local_dep[qid] : &ghost_dep[qid];
      edir = e / 4;
      for (j = 0; j < 2; j++) {
        xdir[0] = (edir + j + 1) % 3;
        xdir[1] = (edir + 2 - j) % 3;
        if (dep->face[xdir[1]] == -1) {
          P4EST_ASSERT (is_ghost[i]);
          if (!is_ghost[i ^ 1]) {
            P4EST_ASSERT (local_dep[qids[i ^ 1] + quadrants_offset].face
                          [xdir[1]] == -2);
            dep->face[xdir[1]] = -2;
          }
          else if (mproc[j][i] == -2) {
            p8est_lnodes_missing_proc_edge (info, zz, xdir[1],
                                            &(mproc[j][0]));
            P4EST_ASSERT (mproc[j][0] != -2 && mproc[j][1] != -2);
            P4EST_ASSERT (mproc[j][0] != rank && mproc[j][1] != rank);
            for (k = 0; k < 2; k++) {
              if (mproc[j][k] >= 0) {
                ip = (int *) sc_array_push (all_procs);
                *ip = mproc[j][k];
              }
            }
            dep->face[xdir[1]] = -((p4est_locidx_t) mproc[j][i] + 3);
          }
          else {
            dep->face[xdir[1]] = -((p4est_locidx_t) mproc[j][i] + 3);
          }
        }
        if (dep->face[xdir[1]] == -2) {
          continue;
        }
        nqid = dep->edge[xdir[0]];
        if (nqid >= 0) {
          has_local = 1;
          start_node = num_inodes + (o ? nodes_per_edge - 1 : 0);
          for (k = 0; k < nodes_per_edge; k++, start_node += stride) {
            nid = nqid * nodes_per_elem + edge_nodes[e][k];
            P4EST_ASSERT (local_elem_nodes[nid] == -1);
            local_elem_nodes[nid] = start_node;
          }
        }
        else if (!is_remote) {
          nproc = nqid;
          if (nproc == -1) {
            nproc = mproc[j][i ^ 1];
            dep->edge[xdir[0]] = -((p4est_locidx_t) nproc + 3);
          }
          else {
            nproc = -(nproc + 3);
          }
          P4EST_ASSERT (nproc >= -1);
          if (nproc >= 0 && nproc != rank) {
            ip = (int *) sc_array_push (all_procs);
            *ip = nproc;
          }
        }
      }
    }
  }
  P4EST_ASSERT (has_local);
  sc_array_sort (all_procs, sc_int_compare);
  sc_array_uniq (all_procs, sc_int_compare);

  count = all_procs->elem_count;
  if (count) {
    type = (int8_t) (P8EST_LN_E_OFFSET + owner_e);
    p4est_lnodes_push_binfo (touching_procs, all_procs, send_buf_info,
                             recv_buf_info, inode_sharers, owner_proc, rank,
                             info->p4est->mpisize, is_remote,
                             type, num_inodes);
  }
  else {
    P4EST_ASSERT (owner_proc == rank);
  }
}

/** p8est_lnodes_face_node_transform:
 *
 * Compute the transformation from an independent node's position on a face
 * to the position on the face of the touching element.
 */
static inline void
p8est_lnodes_face_node_transform (int orig_f, int f, int8_t orientation,
                                  int8_t * flipj, int8_t * flipk,
                                  int8_t * swapjk)
{
  int                 ref = p8est_face_permutation_refs[f][orig_f];
  int                 set = p8est_face_permutation_sets[ref][orientation];
  int                 c0 = p8est_face_permutations[set][0];
  int                 c1 = p8est_face_permutations[set][1];
  int                 c2 = p8est_face_permutations[set][2];
  *flipj = (c1 < c0);
  *flipk = (c2 < c0);
  *swapjk = ((c0 ^ c2) == 1);
}
#endif

/* p4est_lnodes_face_callback:
 *
 * Create new independent nodes on a face.
 * Set all touching element nodes to point to the newly created independent
 * nodes.
 * Compute all processes that share the nodes.
 * If the nodes are locally owned, add info describing the nodes to the send
 * buffer of all processes that share the nodes.
 * If the nodes are not locally owned, add info describing the nodes to the
 * receive buffer of the owner.
 */
static void
p4est_lnodes_face_callback (p4est_iter_face_info_t * info, void *Data)
{
  sc_array_t         *sides = &(info->sides);
  size_t              zz, count = sides->elem_count;
  p4est_lnodes_data_t *data = (p4est_lnodes_data_t *) Data;
  p4est_iter_face_side_t *fside;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *inode;
  sc_array_t         *inode_sharers = data->inode_sharers;
  p4est_locidx_t     *local_elem_nodes = data->local_elem_nodes;
  sc_array_t         *send_buf_info = data->send_buf_info;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  sc_array_t         *touching_procs = data->touching_procs;
  p4est_topidx_t      tid;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_locidx_t      quadrants_offset;
  p4est_locidx_t     *qids;
  p4est_locidx_t      qid, owner_qid;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) inodes->elem_count;
  int8_t             *is_ghost, owner_is_ghost;
  int                 f, owner_f;
  p4est_locidx_t      nid;
  int                 owner_proc;
  int                 rank = info->p4est->mpirank;
  p4est_quadrant_t  **q;
  /* p4est_quadrant_t   *owner_q; */
  int                 nodes_per_face = data->nodes_per_face;
  int                 nodes_per_elem = data->nodes_per_elem;
  int               **face_nodes = data->face_nodes;
  int                 is_hanging;
  int                 i, j, limit;
#ifndef P4_TO_P8
  int                 stride;
#else
  int                 nodes_per_edge = SC_MAX (1, data->nodes_per_edge);
  int8_t              flipj, flipk, swapjk;
  int                 k, l, jind, kind, lind;
#endif
  p4est_locidx_t      start_node;
  int8_t              type;

  sc_array_truncate (touching_procs);
  p4est_lnodes_face_simple_callback (info, data);

  /* the first touching quad is the owner */
  fside = p4est_iter_fside_array_index (sides, 0);
  if (fside->is_hanging) {
    /* owner_q = fside->is.hanging.quad[0]; */
    owner_is_ghost = fside->is.hanging.is_ghost[0];
    owner_qid = fside->is.hanging.quadid[0];
    owner_f = fside->face;
  }
  else {
    /* owner_q = fside->is.full.quad; */
    owner_is_ghost = fside->is.full.is_ghost;
    owner_qid = fside->is.full.quadid;
    owner_f = fside->face;
  }
  if (!(owner_is_ghost)) {
    owner_proc = rank;
    tid = fside->treeid;
    tree = p4est_tree_array_index (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    owner_qid += quadrants_offset;
  }
  else {
    owner_proc = *((int *) sc_array_index (touching_procs, 0));
    P4EST_ASSERT (owner_proc >= 0 && owner_proc != rank);
  }
  sc_array_sort (touching_procs, sc_int_compare);
  sc_array_uniq (touching_procs, sc_int_compare);
  /* create the nodes */
  for (i = 0; i < nodes_per_face; i++) {
    inode = (p4est_locidx_t *) sc_array_push (inodes);
    P4EST_ASSERT (inodes->elem_count <= (size_t)
                  (nodes_per_elem * info->p4est->local_num_quadrants));
    inode[0] = owner_proc;
    inode[1] = owner_qid;
  }

  /* point element nodes to created nodes */
  for (zz = 0; zz < count; zz++) {
    fside = p4est_iter_fside_array_index (sides, zz);
    limit = fside_get_fields (fside, &is_hanging, &tid, &f, &is_ghost, &qids,
                              &q);
    tree = p4est_tree_array_index (trees, tid);
    quadrants_offset = tree->quadrants_offset;

    for (i = 0; i < limit; i++) {
      qid = qids[i];
      if (!is_ghost[i]) {
        qid += quadrants_offset;
#ifndef P4_TO_P8
        start_node = num_inodes + (!zz ? 0 :
                                   !info->orientation ? 0 :
                                   nodes_per_face - 1);
        stride = !zz ? 1 : !info->orientation ? 1 : -1;
        for (j = 0; j < nodes_per_face; j++, start_node += stride) {
          nid = qid * nodes_per_elem + face_nodes[f][j];
          P4EST_ASSERT (local_elem_nodes[nid] == -1);
          local_elem_nodes[nid] = start_node;
        }
#else
        if (!zz) {
          flipj = 0;
          flipk = 0;
          swapjk = 0;
        }
        else {
          p8est_lnodes_face_node_transform (owner_f, f, info->orientation,
                                            &flipj, &flipk, &swapjk);
        }
        start_node = num_inodes;
        for (l = 0, k = 0; k < nodes_per_edge; k++) {
          for (j = 0; j < nodes_per_edge; j++, l++) {
            nid = qid * nodes_per_elem + face_nodes[f][l];
            jind = flipj ? (nodes_per_edge - 1 - j) : j;
            kind = flipk ? (nodes_per_edge - 1 - k) : k;
            lind = swapjk ? (nodes_per_edge * jind + kind) :
              (nodes_per_edge * kind + jind);
            P4EST_ASSERT (local_elem_nodes[nid] == -1);
            local_elem_nodes[nid] = start_node + lind;
          }
        }
#endif
      }
    }
  }

  count = touching_procs->elem_count;
  if (count) {
    type = (int8_t) owner_f;
    p4est_lnodes_push_binfo (NULL, touching_procs, send_buf_info,
                             recv_buf_info, inode_sharers, owner_proc, rank,
                             info->p4est->mpisize, 0, type, num_inodes);
  }
}

/* p4est_lnodes_volume_callback:
 *
 * Create independent nodes and set volume nodes to point to them.
 */
static void
p4est_lnodes_volume_callback (p4est_iter_volume_info_t * info, void *Data)
{
  p4est_lnodes_data_t *data = (p4est_lnodes_data_t *) Data;
  p4est_tree_t       *tree = p4est_tree_array_index (info->p4est->trees,
                                                     info->treeid);
  p4est_locidx_t      qid = info->quadid + tree->quadrants_offset;
  p4est_locidx_t     *elem_nodes = data->local_elem_nodes;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) inodes->elem_count;
  p4est_locidx_t      nid;
  p4est_locidx_t     *inode;
  int                 nodes_per_volume = data->nodes_per_volume;
  int                *volume_nodes = data->volume_nodes;
  int                 nodes_per_elem = data->nodes_per_elem;
  int                 i;
  int                 rank = info->p4est->mpirank;

  for (i = 0; i < nodes_per_volume; i++) {
    nid = qid * nodes_per_elem + volume_nodes[i];
    P4EST_ASSERT (elem_nodes[nid] == -1);
    elem_nodes[nid] = num_inodes + (p4est_locidx_t) i;
    inode = (p4est_locidx_t *) sc_array_push (inodes);
    P4EST_ASSERT (inodes->elem_count <= (size_t)
                  (nodes_per_elem * info->p4est->local_num_quadrants));
    inode[0] = rank;
    inode[1] = qid;
  }
}

static void
p4est_lnodes_init_data (p4est_lnodes_data_t * data, int p, p4est_t * p4est,
                        p4est_ghost_t * ghost_layer, p4est_lnodes_t * lnodes)
{
  int                 i, j, n;
  int                 npv;
  int                 vcount;
  int                 npf, npc;
  int                 fcount[P4EST_FACES];
  int                 ccount[P4EST_CHILDREN];
  int                 f;
  int                 bcount;
  int                 c;
#ifdef P4_TO_P8
  int                 e;
  int                 eshift;
  int                 k;
  int                 npe;
  int                 ecount[P8EST_EDGES];
#endif
  p4est_locidx_t      nlq = p4est->local_num_quadrants;
  p4est_locidx_t      ngq = (p4est_locidx_t) ghost_layer->ghosts.elem_count;
  p4est_locidx_t      nldep = nlq;
  p4est_locidx_t      ngdep = ngq;
  int                 mpisize = p4est->mpisize;

  if (p == -1) {
    data->nodes_per_elem = P4EST_FACES;
    npv = data->nodes_per_volume = 0;
    npf = data->nodes_per_face = 1;
#ifdef P4_TO_P8
    npe = data->nodes_per_edge = 0;
#endif
    npc = data->nodes_per_corner = 0;
  }
#ifdef P4_TO_P8
  else if (p == -2) {
    data->nodes_per_elem = P4EST_FACES + P8EST_EDGES;
    npv = data->nodes_per_volume = 0;
    npf = data->nodes_per_face = 1;
    npe = data->nodes_per_edge = 1;
    npc = data->nodes_per_corner = 0;
  }
#endif
  else if (p == -P4EST_DIM) {
    data->nodes_per_elem = P4EST_FACES +
#ifdef P4_TO_P8
      P8EST_EDGES +
#endif
      P4EST_CHILDREN;

    npv = data->nodes_per_volume = 0;
    npf = data->nodes_per_face = 1;
#ifdef P4_TO_P8
    npe = data->nodes_per_edge = 1;
#endif
    npc = data->nodes_per_corner = 1;
  }
  else {
#ifndef P4_TO_P8
    data->nodes_per_elem = (p + 1) * (p + 1);
    npv = data->nodes_per_volume = (p - 1) * (p - 1);
    npf = data->nodes_per_face = p - 1;
#else
    data->nodes_per_elem = (p + 1) * (p + 1) * (p + 1);
    npv = data->nodes_per_volume = (p - 1) * (p - 1) * (p - 1);
    npf = data->nodes_per_face = (p - 1) * (p - 1);
    npe = data->nodes_per_edge = (p - 1);
#endif
    npc = data->nodes_per_corner = 1;
  }
#ifndef P4_TO_P8
  fcount[0] = fcount[1] = fcount[2] = fcount[3] = 0;
  ccount[0] = ccount[1] = ccount[2] = ccount[3] = 0;
#else
  fcount[0] = fcount[1] = fcount[2] = fcount[3] = fcount[4] = fcount[5] = 0;
  ccount[0] = ccount[1] = ccount[2] = ccount[3] = 0;
  ccount[4] = ccount[5] = ccount[6] = ccount[7] = 0;
  ecount[0] = ecount[1] = ecount[2] = ecount[3] = ecount[4] = ecount[5] = 0;
  ecount[6] = ecount[7] = ecount[8] = ecount[9] = ecount[10] = ecount[11] = 0;
#endif
  vcount = 0;

  data->volume_nodes = P4EST_ALLOC (int, npv);
  for (i = 0; i < P4EST_FACES; i++) {
    data->face_nodes[i] = P4EST_ALLOC (int, npf);
  }
#ifdef P4_TO_P8
  for (i = 0; i < P8EST_EDGES; i++) {
    data->edge_nodes[i] = P4EST_ALLOC (int, npe);
  }
#endif
  for (i = 0; i < P4EST_CHILDREN; i++) {
    data->corner_nodes[i] = P4EST_ALLOC (int, npc);
  }

  if (p > 0) {
    /* figure out which nodes live on which parts of the quadrants */
    n = 0;
#ifdef P4_TO_P8
    for (k = 0; k < p + 1; k++) {
#endif
      for (j = 0; j < p + 1; j++) {
        for (i = 0; i < p + 1; i++, n++) {
          bcount = f = c = 0;
#ifdef P4_TO_P8
          e = 0;
          eshift = -1;
          switch (k == 0 ? 0 : k == p ? 1 : 2) {
          case 0:
            f = 4;
            bcount++;
            break;
          case 1:
            f = 5;
            c |= 4;
            e++;
            bcount++;
            break;
          default:
            eshift = 8;
            break;
          }
#endif
          switch (j == 0 ? 0 : j == p ? 1 : 2) {
          case 0:
            f = 2;
#ifdef P4_TO_P8
            e <<= 1;
#endif
            bcount++;
            break;
          case 1:
            f = 3;
            c |= 2;
#ifdef P4_TO_P8
            e <<= 1;
            e++;
#endif
            bcount++;
            break;
          default:
#ifdef P4_TO_P8
            eshift = 4;
#endif
            break;
          }
          switch (i == 0 ? 0 : i == p ? 1 : 2) {
          case 0:
            bcount++;
#ifdef P4_TO_P8
            e <<= 1;
#endif
            break;
          case 1:
            f = 1;
            c |= 1;
#ifdef P4_TO_P8
            e <<= 1;
            e++;
#endif
            bcount++;
            break;
          default:
#ifdef P4_TO_P8
            eshift = 0;
#endif
            break;
          }
          switch (bcount) {
          case 0:
            data->volume_nodes[vcount++] = n;
            break;
          case 1:
            data->face_nodes[f][fcount[f]++] = n;
            break;
#ifdef P4_TO_P8
          case 2:
            P4EST_ASSERT (eshift >= 0);
            e += eshift;
            data->edge_nodes[e][ecount[e]++] = n;
            break;
#endif
          default:
            data->corner_nodes[c][ccount[c]++] = n;
            break;
          }
        }
      }
#ifdef P4_TO_P8
    }
#endif
  }
  else {
    int                 offset = 0;

    for (i = 0; i < npv; i++) {
      data->volume_nodes[vcount++] = offset++;
    }
    for (f = 0; f < P4EST_FACES; f++) {
      for (i = 0; i < npf; i++) {
        data->face_nodes[f][fcount[f]++] = offset++;
      }
    }
#ifdef P4_TO_P8
    for (e = 0; e < P8EST_EDGES; e++) {
      for (i = 0; i < npe; i++) {
        data->edge_nodes[e][ecount[e]++] = offset++;
      }
    }
#endif
    for (c = 0; c < P4EST_CHILDREN; c++) {
      for (i = 0; i < npc; i++) {
        data->corner_nodes[c][ccount[c]++] = offset++;
      }
    }
  }

  data->local_dep = P4EST_ALLOC (p4est_lnodes_dep_t, nldep);
  memset (data->local_dep, -1, nldep * sizeof (p4est_lnodes_dep_t));
  data->ghost_dep = P4EST_ALLOC (p4est_lnodes_dep_t, ngdep);
  memset (data->ghost_dep, -1, ngdep * sizeof (p4est_lnodes_dep_t));

  data->local_elem_nodes = lnodes->element_nodes;

  data->inodes = sc_array_new (2 * sizeof (p4est_locidx_t));
  data->inode_sharers = sc_array_new (sizeof (int));
  data->send_buf_info = P4EST_ALLOC (sc_array_t, mpisize);
  data->recv_buf_info = P4EST_ALLOC (sc_array_t, mpisize);
  for (i = 0; i < mpisize; i++) {
    sc_array_init (&(data->send_buf_info[i]),
                   sizeof (p4est_lnodes_buf_info_t));
    sc_array_init (&(data->recv_buf_info[i]),
                   sizeof (p4est_lnodes_buf_info_t));
  }
  data->face_codes = lnodes->face_code;
  data->poff = P4EST_ALLOC_ZERO (p4est_locidx_t, mpisize + 1);
  data->touching_procs = sc_array_new (sizeof (int));
  data->all_procs = sc_array_new (sizeof (int));
}

static void
p4est_lnodes_reset_data (p4est_lnodes_data_t * data, p4est_t * p4est)
{
  int                 mpisize = p4est->mpisize;
  int                 i;

  sc_array_destroy (data->touching_procs);
  sc_array_destroy (data->all_procs);
  P4EST_FREE (data->poff);
  P4EST_FREE (data->volume_nodes);
  for (i = 0; i < P4EST_FACES; i++) {
    P4EST_FREE (data->face_nodes[i]);
  }
#ifdef P4_TO_P8
  for (i = 0; i < P8EST_EDGES; i++) {
    P4EST_FREE (data->edge_nodes[i]);
  }
#endif
  for (i = 0; i < P4EST_CHILDREN; i++) {
    P4EST_FREE (data->corner_nodes[i]);
  }

  sc_array_destroy (data->inodes);
  sc_array_destroy (data->inode_sharers);
  for (i = 0; i < mpisize; i++) {
    sc_array_reset (&(data->send_buf_info[i]));
    sc_array_reset (&(data->recv_buf_info[i]));
  }
  P4EST_FREE (data->send_buf_info);
  P4EST_FREE (data->recv_buf_info);
  P4EST_FREE (data->local_dep);
  P4EST_FREE (data->ghost_dep);
  /* do not free face_codes: controlled by lnodes_t */
}

/* p4est_lnodes_count_send:
 *
 * Coming out of the main iteration that finds the independent nodes, but before
 * hanging faces and hanging edges have been fixed, we assign numbers to
 * independent nodes based on greedy numbering: looping through the quadrants,
 * a quadrant numbers of node if it touches it (touches the
 * volume/face/edge/corner containing it) and if it hasn't received a number
 * yet.  This is consistent with the node ownership scheme.
 *
 * At this point each process knows where and what it needs to send, so the
 * sends are initiated.
 */
static void
p4est_lnodes_count_send (p4est_lnodes_data_t * data, p4est_t * p4est,
                         p4est_lnodes_t * lnodes)
{
  p4est_locidx_t      nlq = p4est->local_num_quadrants;
  p4est_locidx_t      nlen, nln;
  p4est_locidx_t      li, *lp;
  p4est_locidx_t      inidx;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *inode;
  p4est_topidx_t     *local_en = data->local_elem_nodes;
  int                 i, j;
  int                 rank = p4est->mpirank;
  int                 mpisize = p4est->mpisize;
  int                 npe = data->nodes_per_elem;
  p4est_locidx_t      count = 0;
  sc_array_t         *send_buf_info = data->send_buf_info;
  sc_array_t         *send_info;
  sc_array_t         *send;
  size_t              zz;
  size_t              zindex;
  p4est_lnodes_buf_info_t *binfo;
  int8_t              type;
  int                 limit;
  int                 nodes_per_face = data->nodes_per_face;
#ifdef P4_TO_P8
  int                 nodes_per_edge = data->nodes_per_edge;
#endif
  int                 nodes_per_corner = data->nodes_per_corner;
  int                 share_count;
  int                 share_proc;
  sc_array_t         *inode_sharers = data->inode_sharers;
  size_t              send_count;
  sc_MPI_Request     *send_request;
  int                 num_send_procs;
  size_t              total_sent;
  int                 mpiret;
  size_t              countz;
  p4est_locidx_t     *poff = data->poff;
  p4est_locidx_t      pcount;

  nlen = ((p4est_locidx_t) npe) * nlq;
  for (li = 0; li < nlen; li++) {
    inidx = local_en[li];
    P4EST_ASSERT (inidx >= 0);
    inode = (p4est_locidx_t *) sc_array_index (inodes, (size_t) inidx);
    /* if this quadrant owns the node */
    if (inode[0] == rank && inode[1] == li / npe) {
      inode[0] = -1;
      inode[1] = count++;
    }
  }
  for (zz = 0; zz < inodes->elem_count; zz++) {
    inode = (p4est_locidx_t *) sc_array_index (inodes, zz);
    if (inode[0] >= 0) {
      P4EST_ASSERT (inode[0] != rank);
      poff[inode[0]]++;
    }
  }

  pcount = 0;
  for (i = 0; i < mpisize; i++) {
    p4est_topidx_t      temp = pcount;

    pcount += poff[i];
    poff[i] = temp;
  }
  poff[mpisize] = pcount;

  lnodes->owned_count = count;
  lnodes->num_local_nodes = nln = (p4est_locidx_t) inodes->elem_count;
  lnodes->nonlocal_nodes = P4EST_ALLOC (p4est_gloidx_t, nln - count);
  memset (lnodes->nonlocal_nodes, -1,
          (nln - count) * sizeof (p4est_gloidx_t));

  num_send_procs = 0;
  total_sent = 0;
  sc_array_init (&(data->send_requests), sizeof (sc_MPI_Request));
  data->send_buf = P4EST_ALLOC (sc_array_t, mpisize);
  for (i = 0; i < mpisize; i++) {
    sc_array_init (&(data->send_buf[i]), sizeof (p4est_locidx_t));
  }
  for (i = 0; i < mpisize; i++) {
    send_info = &(send_buf_info[i]);
    countz = send_info->elem_count;
    if (countz > 0) {
      P4EST_ASSERT (i != p4est->mpirank);
      send = &(data->send_buf[i]);
      for (zz = 0; zz < countz; zz++) {
        binfo = (p4est_lnodes_buf_info_t *) sc_array_index (send_info, zz);
        zindex = (size_t) binfo->first_index;
        type = binfo->type;
        if (type >= P4EST_LN_C_OFFSET) {
          limit = nodes_per_corner;
        }
#ifdef P4_TO_P8
        else if (type >= P8EST_LN_E_OFFSET) {
          limit = nodes_per_edge;
        }
#endif
        else {
          P4EST_ASSERT (0 <= type && type < P4EST_FACES);
          limit = nodes_per_face;
        }
        for (j = 0; j < limit; j++) {
          lp = (p4est_locidx_t *) sc_array_push (send);
          inode = (p4est_locidx_t *) sc_array_index (inodes, zindex++);
          P4EST_ASSERT (inode[0] == -1 && inode[1] >= 0 && inode[1] < count);
          *lp = inode[1];
        }
        if (binfo->send_sharers) {
          lp = (p4est_locidx_t *) sc_array_push (send);
          *lp = (p4est_locidx_t) binfo->share_count;
          P4EST_ASSERT (binfo->share_count > 0);
          zindex = (size_t) binfo->share_offset;
          share_count = (int) binfo->share_count;
          for (j = 0; j < share_count; j++) {
            lp = (p4est_locidx_t *) sc_array_push (send);
            share_proc = *((int *) sc_array_index (inode_sharers, zindex++));
            *lp = (p4est_locidx_t) share_proc;
            P4EST_ASSERT (0 <= share_proc && share_proc < mpisize);
          }
        }
      }
      send_count = send->elem_count;
      send_request = (sc_MPI_Request *) sc_array_push (&data->send_requests);
      mpiret = sc_MPI_Isend (send->array,
                             (int) (send_count * sizeof (p4est_locidx_t)),
                             sc_MPI_BYTE, i, P4EST_COMM_LNODES_PASS,
                             p4est->mpicomm, send_request);
      SC_CHECK_MPI (mpiret);
      num_send_procs++;
      total_sent += (send_count * sizeof (p4est_locidx_t));
    }
  }
  P4EST_VERBOSEF ("Total of %lld bytes sent to %d processes\n",
                  (unsigned long long) total_sent, num_send_procs);
}

#ifdef P4EST_ENABLE_DEBUG
/* p4est_lnodes_test_comm:
 *
 * If the buf_info_t array is the same on both ends of a communication, then the
 * information sent will be properly decoded.
 *
 */
static              int8_t
p4est_lnodes_test_comm (p4est_t * p4est, p4est_lnodes_data_t * data)
{
  int                 mpisize = p4est->mpisize;
  int                 i, j;
  sc_array_t         *send;
  sc_array_t         *send_buf_info = data->send_buf_info;
  sc_array_t         *recv, *recv2;
  sc_array_t         *recv_buf;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  size_t              zz;
  int                 mpiret;
  sc_array_t          send_requests;
  sc_MPI_Request     *send_request;
  sc_MPI_Status       probe_status, recv_status;
  int                 num_recv_procs = 0;
  int                *num_recv_expect = P4EST_ALLOC_ZERO (int, mpisize);
  int                 byte_count;
  size_t              count, elem_count;
  p4est_lnodes_buf_info_t *binfo, *binfo2;

  sc_array_init (&send_requests, sizeof (sc_MPI_Request));
  for (i = 0; i < mpisize; i++) {
    send = &(send_buf_info[i]);
    count = send->elem_count;
    if (count > 0) {
      P4EST_ASSERT (i != p4est->mpirank);
      send_request = (sc_MPI_Request *) sc_array_push (&send_requests);
      mpiret = sc_MPI_Isend (send->array,
                             (int) (count * sizeof (p4est_lnodes_buf_info_t)),
                             sc_MPI_BYTE, i, P4EST_COMM_LNODES_TEST,
                             p4est->mpicomm, send_request);
      SC_CHECK_MPI (mpiret);
    }
    recv = &(recv_buf_info[i]);
    count = recv->elem_count;
    if (count) {
      P4EST_ASSERT (i != p4est->mpirank);
      num_recv_procs++;
      num_recv_expect[i]++;
    }
  }

  recv_buf = P4EST_ALLOC (sc_array_t, mpisize);
  for (i = 0; i < mpisize; i++) {
    sc_array_init (&(recv_buf[i]), sizeof (p4est_lnodes_buf_info_t));
  }
  for (i = 0; i < num_recv_procs; i++) {
    mpiret = sc_MPI_Probe (sc_MPI_ANY_SOURCE, P4EST_COMM_LNODES_TEST,
                           p4est->mpicomm, &probe_status);
    SC_CHECK_MPI (mpiret);
    j = probe_status.MPI_SOURCE;
    P4EST_ASSERT (j != p4est->mpirank && num_recv_expect[j] == 1);
    recv = &(recv_buf[j]);
    mpiret = sc_MPI_Get_count (&probe_status, sc_MPI_BYTE, &byte_count);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (byte_count % ((int) sizeof (p4est_lnodes_buf_info_t)) == 0);
    elem_count = ((size_t) byte_count) / sizeof (p4est_lnodes_buf_info_t);
    sc_array_resize (recv, elem_count);
    mpiret = sc_MPI_Recv (recv->array, byte_count, sc_MPI_BYTE, j,
                          P4EST_COMM_LNODES_TEST, p4est->mpicomm,
                          &recv_status);
    SC_CHECK_MPI (mpiret);
    num_recv_expect[j]--;

    recv2 = &(recv_buf_info[j]);
    P4EST_ASSERT (recv2->elem_count == recv->elem_count);
    for (zz = 0; zz < elem_count; zz++) {
      binfo2 = (p4est_lnodes_buf_info_t *) sc_array_index (recv2, zz);
      binfo = (p4est_lnodes_buf_info_t *) sc_array_index (recv, zz);
      P4EST_ASSERT (binfo->type == binfo2->type);
      P4EST_ASSERT (binfo->send_sharers == binfo2->send_sharers);
      if (!binfo->send_sharers) {
        P4EST_ASSERT (binfo->share_count == binfo2->share_count);
      }
    }
  }

  if (send_requests.elem_count > 0) {
    mpiret = sc_MPI_Waitall ((int) send_requests.elem_count,
                             (sc_MPI_Request *) send_requests.array,
                             sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  sc_array_reset (&send_requests);
  for (i = 0; i < mpisize; i++) {
    sc_array_reset (&(recv_buf[i]));
  }

  P4EST_FREE (recv_buf);
  P4EST_FREE (num_recv_expect);
  return 1;
}
#endif

/* p4est_lnodes_recv:
 *
 * Each process has its sorted receive lists.
 * Each process knows who to expect nodes from.
 * When the nodes are received, they are put into a sorter list dedicated just
 * to the process that the nodes come from, this is then sorted in ascending
 * node order.
 */
static void
p4est_lnodes_recv (p4est_t * p4est, p4est_lnodes_data_t * data,
                   p4est_lnodes_t * lnodes)
{
  int                 mpisize = p4est->mpisize;
  int                 i, j, k;
  int                 limit;
  sc_array_t         *recv, *recv_info;
  sc_array_t         *recv_buf;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  size_t              count, info_count, zz;
  int                 mpiret;
  sc_MPI_Status       probe_status, recv_status;
  int                 num_recv_procs = 0;
  size_t              total_recv = 0;
  int                *num_recv_expect = P4EST_ALLOC_ZERO (int, mpisize);
  int                 byte_count;
  size_t              elem_count;
  p4est_lnodes_buf_info_t *binfo;
  size_t              zindex;
  int                 nodes_per_face = data->nodes_per_face;
#ifdef P4_TO_P8
  int                 nodes_per_edge = data->nodes_per_edge;
#endif
  int                 nodes_per_corner = data->nodes_per_corner;
  p4est_locidx_t     *lp;
  int                *ip;
  p4est_locidx_t     *inode;
  sc_array_t         *inode_sharers = data->inode_sharers;
  sc_array_t         *inodes = data->inodes;
  int                 share_count;
  sc_array_t         *sorter;
  p4est_gloidx_t     *nonlocal_nodes = lnodes->nonlocal_nodes;
  p4est_locidx_t     *poff = data->poff;

  for (i = 0; i < mpisize; i++) {
    recv_info = &(recv_buf_info[i]);
    count = recv_info->elem_count;
    if (count) {
      P4EST_ASSERT (i != p4est->mpirank);
      P4EST_ASSERT (poff[i + 1] - poff[i] > 0);
      num_recv_procs++;
      num_recv_expect[i]++;
    }
  }

  sorter = sc_array_new (2 * sizeof (p4est_locidx_t));

  recv_buf = P4EST_ALLOC (sc_array_t, mpisize);
  for (i = 0; i < mpisize; i++) {
    sc_array_init (&(recv_buf[i]), sizeof (p4est_locidx_t));
  }
  for (i = 0; i < num_recv_procs; i++) {
    mpiret = sc_MPI_Probe (sc_MPI_ANY_SOURCE, P4EST_COMM_LNODES_PASS,
                           p4est->mpicomm, &probe_status);
    SC_CHECK_MPI (mpiret);
    j = probe_status.MPI_SOURCE;
    P4EST_ASSERT (j != p4est->mpirank && num_recv_expect[j] == 1);
    recv = &(recv_buf[j]);
    recv_info = &(recv_buf_info[j]);
    mpiret = sc_MPI_Get_count (&probe_status, sc_MPI_BYTE, &byte_count);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (byte_count % ((int) sizeof (p4est_locidx_t)) == 0);
    elem_count = ((size_t) byte_count) / sizeof (p4est_locidx_t);
    sc_array_resize (recv, elem_count);
    mpiret = sc_MPI_Recv (recv->array, byte_count, sc_MPI_BYTE, j,
                          P4EST_COMM_LNODES_PASS, p4est->mpicomm,
                          &recv_status);
    SC_CHECK_MPI (mpiret);
    num_recv_expect[j]--;

    info_count = recv_info->elem_count;
    count = 0;
    for (zz = 0; zz < info_count; zz++) {
      binfo = (p4est_lnodes_buf_info_t *) sc_array_index (recv_info, zz);
      if (binfo->type >= P4EST_LN_C_OFFSET) {
        limit = nodes_per_corner;
      }
#ifdef P4_TO_P8
      else if (binfo->type >= P8EST_LN_E_OFFSET) {
        limit = nodes_per_edge;
      }
#endif
      else {
        limit = nodes_per_face;
      }
      zindex = (size_t) binfo->first_index;
      for (k = 0; k < limit; k++) {
        inode = (p4est_locidx_t *) sc_array_index (inodes, zindex);
        lp = (p4est_locidx_t *) sc_array_index (recv, count++);
        P4EST_ASSERT (inode[0] == j);
        P4EST_ASSERT (*lp >= 0);
        inode[1] = *lp;
        lp = (p4est_locidx_t *) sc_array_push (sorter);
        lp[0] = (p4est_locidx_t) inode[1];
        lp[1] = (p4est_locidx_t) zindex++;
      }
      if (binfo->send_sharers) {
        lp = (p4est_locidx_t *) sc_array_index (recv, count++);
        share_count = (int) (*lp);
        P4EST_ASSERT (share_count > 0);
        P4EST_ASSERT (binfo->share_count == -1);
        P4EST_ASSERT (binfo->share_offset == -1);
        binfo->share_count = (int8_t) share_count;
        binfo->share_offset = (p4est_locidx_t) inode_sharers->elem_count;
        ip = (int *) sc_array_push_count (inode_sharers, share_count);
        for (k = 0; k < share_count; k++) {
          lp = (p4est_locidx_t *) sc_array_index (recv, count++);
          ip[k] = (int) (*lp);
          P4EST_ASSERT (0 <= *ip && *ip < mpisize);
        }
      }
    }
    P4EST_ASSERT (count == elem_count);
    total_recv += byte_count;
    P4EST_ASSERT ((p4est_locidx_t) sorter->elem_count ==
                  poff[j + 1] - poff[j]);
    sc_array_sort (sorter, p4est_locidx_compare);
    for (zz = 0; zz < sorter->elem_count; zz++) {
      lp = (p4est_locidx_t *) sc_array_index (sorter, zz);
      nonlocal_nodes[poff[j] + zz] = lp[1];
    }
    sc_array_reset (sorter);
  }

  if (data->send_requests.elem_count > 0) {
    mpiret = sc_MPI_Waitall ((int) data->send_requests.elem_count,
                             (sc_MPI_Request *) data->send_requests.array,
                             sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  sc_array_reset (&data->send_requests);
  for (i = 0; i < mpisize; i++) {
    sc_array_reset (&(data->send_buf[i]));
    sc_array_reset (&(recv_buf[i]));
  }

  P4EST_VERBOSEF ("Total of %lld bytes received from %d processes\n",
                  (unsigned long long) total_recv, num_recv_procs);
  P4EST_FREE (data->send_buf);
  P4EST_FREE (recv_buf);
  P4EST_FREE (num_recv_expect);
  sc_array_destroy (sorter);
}

/* p4est_lnodes_global_and_sharers:
 *
 * Each process that owns a node shared by the local process has a list of nodes
 * that is sorted in increasing local_index, while the local element nodes refer
 * to indices in the inodes array.  A map between the two is created, which
 * allows local element nodes to point to global nodes.  After this is done, the
 * sharers can be created.
 */
static              p4est_gloidx_t
p4est_lnodes_global_and_sharers (p4est_lnodes_data_t * data,
                                 p4est_lnodes_t * lnodes, p4est_t * p4est)
{
  int                 i, j, k, l;
  int                 mpisize = p4est->mpisize;
  p4est_gloidx_t     *gnodes = lnodes->nonlocal_nodes, gtotal;
  size_t              count, zz;
  p4est_locidx_t     *lp, li, *inode;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *elnodes = lnodes->element_nodes;
  p4est_locidx_t      nlen = lnodes->num_local_elements * lnodes->vnodes;
#ifdef P4EST_ENABLE_DEBUG
  p4est_locidx_t      num_inodes = (p4est_locidx_t) data->inodes->elem_count;
#endif
  p4est_locidx_t      inidx;
  int                *comm_proc;
  int                 comm_proc_count;
  sc_array_t         *inode_sharers = data->inode_sharers;
  sc_array_t         *sharers;
  p4est_lnodes_rank_t *lrank;
  sc_array_t         *binfo_array;
  p4est_lnodes_buf_info_t *binfo;
  p4est_locidx_t      share_offset, owned_count = lnodes->owned_count;
  int                 share_count;
  int                 limit;
  int                 nodes_per_face = data->nodes_per_face;
  size_t              zindex;
#ifdef P4_TO_P8
  int                 nodes_per_edge = data->nodes_per_edge;
#endif
  int                 proc;
  int                 shareidx;
  p4est_locidx_t      gidx;
  sc_array_t         *shared_nodes;
  p4est_locidx_t     *global_num_indep;
  p4est_gloidx_t     *global_offsets = P4EST_ALLOC (p4est_gloidx_t,
                                                    mpisize + 1);
  p4est_locidx_t     *poff = data->poff;

  global_num_indep = lnodes->global_owned_count = P4EST_ALLOC (p4est_locidx_t,
                                                               mpisize);
  sc_MPI_Allgather (&owned_count, 1, P4EST_MPI_LOCIDX, global_num_indep, 1,
                    P4EST_MPI_LOCIDX, p4est->mpicomm);

  global_offsets[0] = 0;
  for (i = 0; i < mpisize; i++) {
    global_offsets[i + 1] = global_offsets[i] +
      (p4est_gloidx_t) global_num_indep[i];
  }
  lnodes->global_offset = global_offsets[p4est->mpirank];
  gtotal = global_offsets[p4est->mpisize];

  i = p4est->mpirank;
  for (i = 0; i < mpisize; i++) {
    if (i == p4est->mpirank) {
      continue;
    }
    for (j = poff[i]; j < poff[i + 1]; j++) {
      li = gnodes[j];
      inode = (p4est_locidx_t *) sc_array_index (inodes, li);
      P4EST_ASSERT (inode[0] == i);
      gnodes[j] = inode[1] + global_offsets[i];
      inode[1] = j + owned_count;
    }
  }

  for (li = 0; li < nlen; li++) {
    inidx = elnodes[li];
    P4EST_ASSERT (0 <= inidx && inidx < num_inodes);
    inode = (p4est_locidx_t *) sc_array_index (inodes, (size_t) inidx);
    if (inode[0] == -1) {
      P4EST_ASSERT (0 <= inode[1] && inode[1] < lnodes->owned_count);
      elnodes[li] = inode[1];
    }
    else {
      P4EST_ASSERT (inode[0] >= 0 && inode[0] != p4est->mpirank &&
                    inode[0] < mpisize);
      P4EST_ASSERT (inode[1] >= poff[inode[0]] + owned_count &&
                    inode[1] < poff[inode[0] + 1] + owned_count);
      elnodes[li] = inode[1];
    }
  }

  /* figure out all nodes that also share nodes shared by the local process */
  comm_proc = P4EST_ALLOC_ZERO (int, mpisize);
  count = inode_sharers->elem_count;
  for (zz = 0; zz < count; zz++) {
    i = *((int *) sc_array_index (inode_sharers, zz));
    comm_proc[i] = 1;
  }
  /* create an entry in sharers for each such process, providing a map from
   * process id to sharer index */
  comm_proc_count = 0;
  lnodes->sharers = sharers = sc_array_new (sizeof (p4est_lnodes_rank_t));
  for (i = 0; i < mpisize; i++) {
    if (comm_proc[i]) {
      lrank = (p4est_lnodes_rank_t *) sc_array_push (sharers);
      lrank->rank = i;
      sc_array_init (&(lrank->shared_nodes), sizeof (p4est_locidx_t));
      comm_proc[i] = comm_proc_count++;
    }
    else {
      comm_proc[i] = -1;
    }
  }

  /* for every node in a send or receive list, figure out which global node it
   * is, and which processes share it, and add the index in global nodes to that
   * sharer's element_nodes array.
   */
  for (i = 0; i < mpisize; i++) {
    for (j = 0; j < 2; j++) {
      if (j == 0) {
        binfo_array = &(data->send_buf_info[i]);
      }
      else {
        binfo_array = &(data->recv_buf_info[i]);
      }
      count = binfo_array->elem_count;
      for (zz = 0; zz < count; zz++) {
        binfo = (p4est_lnodes_buf_info_t *) sc_array_index (binfo_array, zz);
        if (binfo->type >= P4EST_LN_C_OFFSET) {
          limit = 1;
        }
#ifdef P4_TO_P8
        else if (binfo->type >= P8EST_LN_E_OFFSET) {
          limit = nodes_per_edge;
        }
#endif
        else {
          limit = nodes_per_face;
        }
        zindex = (size_t) binfo->first_index;
        share_offset = binfo->share_offset;
        share_count = (int) binfo->share_count;
        for (k = 0; k < limit; k++) {
          inode = (p4est_locidx_t *) sc_array_index (inodes, zindex++);
          gidx = inode[1];
          if (j == 0) {
            P4EST_ASSERT (inode[0] == -1);
            P4EST_ASSERT (gidx < owned_count);
            shareidx = comm_proc[i];
            P4EST_ASSERT (shareidx >= 0);
            lrank = p4est_lnodes_rank_array_index_int (sharers, shareidx);
            P4EST_ASSERT (lrank->rank == i);
            lp = (p4est_locidx_t *) sc_array_push (&(lrank->shared_nodes));
            *lp = gidx;

            P4EST_ASSERT (share_count >= 2);
            proc = *((int *) sc_array_index (inode_sharers,
                                             (size_t) share_offset + 1));
            P4EST_ASSERT (proc != p4est->mpirank);
            if (proc == i) {
              shareidx = comm_proc[p4est->mpirank];
              P4EST_ASSERT (shareidx >= 0);
              lrank = p4est_lnodes_rank_array_index_int (sharers, shareidx);
              P4EST_ASSERT (lrank->rank == p4est->mpirank);
              lp = (p4est_locidx_t *) sc_array_push (&(lrank->shared_nodes));
              *lp = gidx;
            }
          }
          else {
            P4EST_ASSERT (inode[0] == i);
            P4EST_ASSERT (poff[i] <= inode[1] - owned_count &&
                          inode[1] - owned_count < poff[i + 1]);
            for (l = 0; l < share_count; l++) {
              proc = *((int *) sc_array_index (inode_sharers,
                                               (size_t) share_offset +
                                               (size_t) l));
              shareidx = comm_proc[proc];
              P4EST_ASSERT (shareidx >= 0);
              lrank = p4est_lnodes_rank_array_index_int (sharers, shareidx);
              P4EST_ASSERT (lrank->rank == proc);
              lp = (p4est_locidx_t *) sc_array_push (&(lrank->shared_nodes));
              *lp = gidx;
            }
          }
        }
      }
    }
  }

  /* for each sharer, figure out which entries in element_nodes are owned by
   * the current process, and which are owned by the sharer's rank */
  for (i = 0; i < comm_proc_count; i++) {
    lrank = p4est_lnodes_rank_array_index_int (sharers, i);
    shared_nodes = &(lrank->shared_nodes);
    count = shared_nodes->elem_count;
    if (count) {
      sc_array_t         *sortshared =
        sc_array_new_size (2 * sizeof (p4est_gloidx_t),
                           count);
      for (zz = 0; zz < count; zz++) {
        p4est_gloidx_t     *gp;

        gidx = *((p4est_locidx_t *) sc_array_index (shared_nodes, zz));
        gp = (p4est_gloidx_t *) sc_array_index (sortshared, zz);
        gp[0] = p4est_lnodes_global_index (lnodes, gidx);
        gp[1] = gidx;
      }
      sc_array_sort (sortshared, p4est_gloidx_compare);
      for (zz = 0; zz < count; zz++) {
        p4est_gloidx_t     *gp;

        gp = (p4est_gloidx_t *) sc_array_index (sortshared, zz);
        *((p4est_locidx_t *) sc_array_index (shared_nodes, zz)) = gp[1];
      }
      sc_array_destroy (sortshared);
    }
    proc = lrank->rank;
    lrank->shared_mine_offset = -1;
    lrank->shared_mine_count = 0;
    for (zz = 0; zz < count; zz++) {
      gidx = *((p4est_locidx_t *) sc_array_index (shared_nodes, zz));
      if (gidx < lnodes->owned_count) {
        if (lrank->shared_mine_count == 0) {
          lrank->shared_mine_offset = (p4est_locidx_t) zz;
        }
        lrank->shared_mine_count++;
      }
    }
    if (proc == p4est->mpirank) {
      lrank->owned_count = lnodes->owned_count;
      lrank->owned_offset = 0;
    }
    else {
      lrank->owned_offset = poff[proc] + owned_count;
      lrank->owned_count = poff[proc + 1] - poff[proc];
      P4EST_VERBOSEF ("Processor %d shares %llu nodes with processor %d\n",
                      p4est->mpirank, (unsigned long long) count,
                      lrank->rank);
      P4EST_VERBOSEF ("Processor %d owns %d nodes used by processor %d\n",
                      p4est->mpirank, lrank->shared_mine_count, lrank->rank);
      P4EST_VERBOSEF ("Processor %d borrows %d nodes from processor %d\n",
                      p4est->mpirank, lrank->owned_count, lrank->rank);
    }
  }
  P4EST_FREE (comm_proc);
  P4EST_FREE (global_offsets);

  return gtotal;
}

p4est_lnodes_t     *
p4est_lnodes_new (p4est_t * p4est, p4est_ghost_t * ghost_layer, int degree)
{
  p4est_iter_face_t   fiter;
  p4est_iter_volume_t viter;
  p4est_iter_corner_t citer;
#ifdef P4_TO_P8
  p8est_iter_edge_t   eiter;
#endif
  p4est_lnodes_data_t data;
  p4est_locidx_t      nel;
  p4est_locidx_t      nlen;
#ifdef P4EST_ENABLE_DEBUG
  p4est_locidx_t      lj;
#endif
  p4est_lnodes_t     *lnodes = P4EST_ALLOC (p4est_lnodes_t, 1);
  p4est_gloidx_t      gtotal;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING "_lnodes_new, degree %d\n",
                            degree);
  p4est_log_indent_push ();

#ifndef P4_TO_P8
  P4EST_ASSERT (degree >= 1 || degree == -1 || degree == -P4EST_DIM);
#else
  P4EST_ASSERT (degree >= 1 || degree == -1 ||
                degree == -2 || degree == -P4EST_DIM);
#endif

  lnodes->mpicomm = p4est->mpicomm;
  lnodes->degree = degree;
  lnodes->num_local_elements = nel = p4est->local_num_quadrants;
  if (degree > 0) {
#ifndef P4_TO_P8
    lnodes->vnodes = (degree + 1) * (degree + 1);
#else
    lnodes->vnodes = (degree + 1) * (degree + 1) * (degree + 1);
#endif
  }
  else if (degree == -1) {
    lnodes->vnodes = P4EST_FACES;
  }
#ifdef P4_TO_P8
  else if (degree == -2) {
    lnodes->vnodes = P4EST_FACES + P8EST_EDGES;
  }
#endif
  else if (degree == -P4EST_DIM) {
    lnodes->vnodes = P4EST_FACES +
#ifdef P4_TO_P8
      P8EST_EDGES +
#endif
      P4EST_CHILDREN;
  }
  lnodes->face_code = P4EST_ALLOC_ZERO (p4est_lnodes_code_t, nel);
  nlen = nel * lnodes->vnodes;
  lnodes->element_nodes = P4EST_ALLOC (p4est_locidx_t, nlen);
  memset (lnodes->element_nodes, -1, nlen * sizeof (p4est_locidx_t));

  p4est_lnodes_init_data (&data, degree, p4est, ghost_layer, lnodes);
  viter = data.nodes_per_volume ? p4est_lnodes_volume_callback : NULL;
  fiter = data.nodes_per_face ? p4est_lnodes_face_callback :
    ((data.nodes_per_corner ||
#ifdef P4_TO_P8
      data.nodes_per_edge ||
#endif
      0) ? p4est_lnodes_face_simple_callback : NULL);
#ifdef P4_TO_P8
  eiter = data.nodes_per_edge ? p8est_lnodes_edge_callback :
    (data.nodes_per_corner ? p8est_lnodes_edge_simple_callback_void : NULL);
#endif
  citer = data.nodes_per_corner ? p4est_lnodes_corner_callback : NULL;

  p4est_iterate_ext (p4est, ghost_layer, &data, viter, fiter,
#ifdef P4_TO_P8
                     eiter,
#endif
                     citer, 1);

#ifdef P4EST_ENABLE_DEBUG
  for (lj = 0; lj < nlen; lj++) {
    P4EST_ASSERT (lnodes->element_nodes[lj] >= 0);
  }
#endif

  P4EST_ASSERT (p4est_lnodes_test_comm (p4est, &data));

  p4est_lnodes_count_send (&data, p4est, lnodes);

  p4est_lnodes_recv (p4est, &data, lnodes);

  gtotal = p4est_lnodes_global_and_sharers (&data, lnodes, p4est);

  p4est_lnodes_reset_data (&data, p4est);

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING "_lnodes_new with"
                            " %lld global nodes\n",
                            (unsigned long long) gtotal);
  return lnodes;
}

void
p4est_lnodes_destroy (p4est_lnodes_t * lnodes)
{
  size_t              zz, count;
  p4est_lnodes_rank_t *lrank;

  P4EST_FREE (lnodes->element_nodes);
  P4EST_FREE (lnodes->nonlocal_nodes);
  P4EST_FREE (lnodes->global_owned_count);
  P4EST_FREE (lnodes->face_code);

  count = lnodes->sharers->elem_count;
  for (zz = 0; zz < count; zz++) {
    lrank = p4est_lnodes_rank_array_index (lnodes->sharers, zz);
    sc_array_reset (&(lrank->shared_nodes));
  }
  sc_array_destroy (lnodes->sharers);

  P4EST_FREE (lnodes);
}

#ifdef P4EST_ENABLE_MPI

static              size_t
ghost_tree_type (sc_array_t * array, size_t zindex, void *data)
{
  p4est_quadrant_t   *q;

  P4EST_ASSERT (array->elem_size == sizeof (p4est_quadrant_t));

  q = (p4est_quadrant_t *) sc_array_index (array, zindex);
  return (size_t) q->p.which_tree;
}

#endif /* P4EST_ENABLE_MPI */

void
p4est_ghost_support_lnodes (p4est_t * p4est, p4est_lnodes_t * lnodes,
                            p4est_ghost_t * ghost)
{
#ifdef P4EST_ENABLE_MPI
  sc_array_t         *ghosts = &ghost->ghosts;
  sc_array_t         *mirrors = &ghost->mirrors;
  p4est_locidx_t     *proc_offsets = ghost->proc_offsets;
  p4est_locidx_t     *tree_offsets = ghost->tree_offsets;
  p4est_locidx_t     *mirror_proc_offsets = ghost->mirror_proc_offsets;
  p4est_locidx_t     *mirror_proc_mirrors = ghost->mirror_proc_mirrors;
  p4est_locidx_t     *mirror_tree_offsets = ghost->mirror_tree_offsets;
  sc_array_t         *new_ghosts, *new_mirrors;
  int                 mpisize = p4est->mpisize;
  int                 self = p4est->mpirank;
  int                 mpiret;
  int                 n_comm;
  p4est_locidx_t      old_num_mirrors, *newmpcount;
  sc_array_t         *recv_counts, *send_counts;
  sc_array_t         *recv_all, *send_all;
  sc_array_t         *recv_requests;
  sc_array_t         *send_requests;
  p4est_connectivity_t *conn = p4est->connectivity;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING "_ghost_support_lnodes %s\n",
                            p4est_connect_type_string (ghost->btype));
  p4est_log_indent_push ();

  /* this should only be done with an unexpanded ghost layer */
  P4EST_ASSERT (ghost->mirror_proc_fronts == ghost->mirror_proc_mirrors &&
                ghost->mirror_proc_front_offsets ==
                ghost->mirror_proc_front_offsets);

  /* get the topological nodes */
  n_comm = lnodes->sharers ? (int) lnodes->sharers->elem_count : 0;
  recv_counts = sc_array_new_size (sizeof (p4est_locidx_t), n_comm);
  send_counts = sc_array_new_size (sizeof (p4est_locidx_t), n_comm);
  recv_all = sc_array_new_size (sizeof (sc_array_t), n_comm);
  send_all = sc_array_new_size (sizeof (sc_array_t), n_comm);
  recv_requests = sc_array_new (sizeof (sc_MPI_Request));
  send_requests = sc_array_new (sizeof (sc_MPI_Request));

  new_mirrors =
    sc_array_new_size (sizeof (p4est_quadrant_t), mirrors->elem_count);
  old_num_mirrors = (p4est_locidx_t) mirrors->elem_count;
  sc_array_copy (new_mirrors, mirrors);

  /* figure out which quads to send and send them */
  {
    int                *quads_per_node;
    int                *tmp_count;
    p4est_quadrant_t   *node_to_quad;
    p4est_locidx_t     *qpn_offsets;
    int                 V, vid;
    p4est_locidx_t      nid, elid, N, K;
    p4est_lnodes_rank_t *lrank;
    int                 i, p;
    p4est_topidx_t      flt, llt, t;

    N = lnodes->num_local_nodes;
    K = lnodes->num_local_elements;
    V = lnodes->vnodes;
    /* count the quadrants per node */
    quads_per_node = P4EST_ALLOC_ZERO (int, N);
    for (elid = 0; elid < K; elid++) {
      for (vid = 0; vid < V; vid++) {
        nid = lnodes->element_nodes[V * elid + vid];
        quads_per_node[nid]++;
      }
    }

    /* convert counts to offsets */
    qpn_offsets = P4EST_ALLOC (int, N + 1);
    qpn_offsets[0] = 0;
    for (nid = 0; nid < N; nid++) {
      qpn_offsets[nid + 1] = qpn_offsets[nid] + quads_per_node[nid];
    }

    /* create and fill the node to quadrant array */
    node_to_quad = P4EST_ALLOC (p4est_quadrant_t, qpn_offsets[N]);
    tmp_count = P4EST_ALLOC_ZERO (int, N);

    flt = p4est->first_local_tree;
    llt = p4est->last_local_tree;

    for (elid = 0, t = flt; t <= llt; t++) {
      p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, t);
      sc_array_t         *quadrants = &tree->quadrants;
      p4est_locidx_t      il, nquads = (p4est_locidx_t) quadrants->elem_count;

      for (il = 0; il < nquads; il++, elid++) {
        for (vid = 0; vid < V; vid++) {
          p4est_quadrant_t   *q;

          nid = lnodes->element_nodes[V * elid + vid];
          q = &(node_to_quad[qpn_offsets[nid] + tmp_count[nid]++]);
          *q = *p4est_quadrant_array_index (quadrants, (size_t) il);
          q->p.piggy3.which_tree = t;
          q->p.piggy3.local_num = elid;
        }
      }
    }

    newmpcount = P4EST_ALLOC (p4est_locidx_t, p4est->mpisize);
    for (i = 0; i < p4est->mpisize; i++) {
      newmpcount[i] = mirror_proc_offsets[i + 1] - mirror_proc_offsets[i];
    }

    /* for each communication partner */
    for (i = 0; i < n_comm; i++) {
      p4est_locidx_t     *count_recv, *count_send;
      sc_MPI_Request     *req;
      p4est_locidx_t      nshared, il, jl, nquads;
      sc_array_t         *shared, *send_quads;
      sc_array_t          mirror_view;
      p4est_locidx_t      startmirror, endmirror;

      lrank = p4est_lnodes_rank_array_index_int (lnodes->sharers, i);

      p = lrank->rank;
      if (p == self) {
        continue;
      }
      /* post the count receive */
      count_recv = (p4est_locidx_t *) sc_array_index_int (recv_counts, i);
      req = (sc_MPI_Request *) sc_array_push (recv_requests);
      mpiret = sc_MPI_Irecv (count_recv, 1, P4EST_MPI_LOCIDX, p,
                             P4EST_COMM_GHOST_SUPPORT_COUNT, p4est->mpicomm,
                             req);
      SC_CHECK_MPI (mpiret);

      /* create the send buffer */
      shared = &(lrank->shared_nodes);
      nshared = (p4est_locidx_t) shared->elem_count;
      send_quads = (sc_array_t *) sc_array_index_int (send_all, i);
      sc_array_init (send_quads, sizeof (p4est_quadrant_t));
      for (il = 0; il < nshared; il++) {
        p4est_locidx_t      startquad, endquad, qid;

        nid = *((p4est_locidx_t *) sc_array_index (shared, il));
        startquad = qpn_offsets[nid];
        endquad = qpn_offsets[nid + 1];
        for (qid = startquad; qid < endquad; qid++) {
          p4est_quadrant_t   *q;

          q = p4est_quadrant_array_push (send_quads);
          *q = node_to_quad[qid];
        }
      }
      sc_array_sort (send_quads, p4est_quadrant_compare_piggy);
      sc_array_uniq (send_quads, p4est_quadrant_compare_piggy);

      nquads = (p4est_locidx_t) send_quads->elem_count;

      startmirror = mirror_proc_offsets[p];
      endmirror = mirror_proc_offsets[p + 1];
      sc_array_init_data (&mirror_view,
                          &(mirror_proc_mirrors[startmirror]),
                          sizeof (p4est_locidx_t),
                          (size_t) (endmirror - startmirror));

      for (il = 0, jl = 0; (size_t) il < send_quads->elem_count; il++) {
        p4est_quadrant_t   *q;
        ssize_t             idx;
        int                 already_sent = 0;

        q = p4est_quadrant_array_index (send_quads, (size_t) il);
        idx = sc_array_bsearch (mirrors, q, p4est_quadrant_compare_piggy);
        if (idx >= 0) {
          p4est_locidx_t      key = (p4est_locidx_t) idx;
          idx = sc_array_bsearch (&mirror_view, &key, p4est_locidx_compare);
          if (idx >= 0) {
            already_sent = 1;
          }
        }
        else {
          p4est_quadrant_t   *q2;

          q2 = p4est_quadrant_array_push (new_mirrors);
          *q2 = *q;
        }
        if (already_sent) {
          nquads--;
        }
        else {
          size_t              jz;

          jz = (size_t) jl++;
          if (jz != (size_t) il) {
            p4est_quadrant_t   *q2;

            q2 = p4est_quadrant_array_index (send_quads, jz);
            *q2 = *q;
          }
        }
      }
      sc_array_resize (send_quads, nquads);
      newmpcount[p] += nquads;

      P4EST_LDEBUGF ("ghost layer support nodes sending %lld new to %d\n",
                     (long long) nquads, p);

      count_send = (p4est_locidx_t *) sc_array_index_int (send_counts, i);
      *count_send = nquads;
      req = (sc_MPI_Request *) sc_array_push (send_requests);
      mpiret = sc_MPI_Isend (count_send, 1, P4EST_MPI_LOCIDX, p,
                             P4EST_COMM_GHOST_SUPPORT_COUNT, p4est->mpicomm,
                             req);
      SC_CHECK_MPI (mpiret);

      if (nquads) {
        req = (sc_MPI_Request *) sc_array_push (send_requests);
        mpiret = sc_MPI_Isend (send_quads->array, nquads * sizeof
                               (p4est_quadrant_t), sc_MPI_BYTE, p,
                               P4EST_COMM_GHOST_SUPPORT_LOAD, p4est->mpicomm,
                               req);
        SC_CHECK_MPI (mpiret);
      }
    }

    P4EST_FREE (quads_per_node);
    P4EST_FREE (tmp_count);
    P4EST_FREE (node_to_quad);
    P4EST_FREE (qpn_offsets);
  }

  /* update mirror structures */
  {
    p4est_locidx_t      new_num_mirrors;
    p4est_locidx_t     *newmpoffset, *new_mirror_proc_mirrors;
    int                 i, p;
    p4est_lnodes_rank_t *lrank;

    newmpoffset = P4EST_ALLOC (p4est_locidx_t, mpisize + 1);
    newmpoffset[0] = 0;
    for (i = 0; i < mpisize; i++) {
      newmpoffset[i + 1] = newmpoffset[i] + newmpcount[i];
    }
    new_mirror_proc_mirrors =
      P4EST_ALLOC (p4est_locidx_t, newmpoffset[mpisize]);

    sc_array_sort (new_mirrors, p4est_quadrant_compare_piggy);
    sc_array_uniq (new_mirrors, p4est_quadrant_compare_piggy);
    new_num_mirrors = (p4est_locidx_t) new_mirrors->elem_count;
    P4EST_ASSERT (new_num_mirrors >= old_num_mirrors);

    if (new_num_mirrors > old_num_mirrors) {
      sc_array_t          split;
      p4est_topidx_t      t;

      /* update mirror_tree_offsets */
      sc_array_init (&split, sizeof (size_t));
      sc_array_split (new_mirrors, &split,
                      (size_t) conn->num_trees, ghost_tree_type, NULL);
      P4EST_ASSERT (split.elem_count == (size_t) conn->num_trees + 1);
      for (t = 0; t <= conn->num_trees; ++t) {
        size_t             *ppz;

        ppz = (size_t *) sc_array_index (&split, (size_t) t);
        mirror_tree_offsets[t] = (p4est_locidx_t) (*ppz);
      }
      sc_array_reset (&split);
    }

    /* update mirror_proc_mirrors */
    for (p = 0; p < mpisize; p++) {
      p4est_locidx_t      old_offset = mirror_proc_offsets[p];
      p4est_locidx_t      old_count = mirror_proc_offsets[p + 1] - old_offset;
      p4est_locidx_t      il;

      for (il = 0; il < old_count; il++) {
        ssize_t             idx;
        p4est_quadrant_t   *q;

        q = p4est_quadrant_array_index (mirrors, (size_t)
                                        mirror_proc_mirrors[old_offset + il]);

        idx = sc_array_bsearch (new_mirrors, q, p4est_quadrant_compare_piggy);
        P4EST_ASSERT (idx >= 0);
        new_mirror_proc_mirrors[newmpoffset[p] + il] = (p4est_locidx_t) idx;
      }
    }

    for (i = 0; i < n_comm; i++) {
      sc_array_t         *send_quads;
      int                 p;

      lrank = p4est_lnodes_rank_array_index_int (lnodes->sharers, i);

      p = lrank->rank;
      if (p == self) {
        continue;
      }
      send_quads = (sc_array_t *) sc_array_index_int (send_all, i);
      if (send_quads->elem_count) {
        p4est_locidx_t      il, startidx, endidx;
        p4est_locidx_t      old_offset = mirror_proc_offsets[p];
        p4est_locidx_t      old_count =
          mirror_proc_offsets[p + 1] - old_offset;
        sc_array_t          pview;

        startidx = newmpoffset[p];
        endidx = newmpoffset[p + 1];

        for (il = 0; (size_t) il < send_quads->elem_count; il++) {
          p4est_quadrant_t   *q;
          ssize_t             idx;
          p4est_locidx_t      idxl;

          q = p4est_quadrant_array_index (send_quads, (size_t) il);
          idx =
            sc_array_bsearch (new_mirrors, q, p4est_quadrant_compare_piggy);
          P4EST_ASSERT (idx >= 0);
          idxl = (p4est_locidx_t) idx;
          new_mirror_proc_mirrors[startidx + old_count + il] = idxl;
        }

        sc_array_init_data (&pview, &new_mirror_proc_mirrors[startidx],
                            sizeof (p4est_locidx_t),
                            (size_t) (endidx - startidx));
        sc_array_sort (&pview, p4est_locidx_compare);
        sc_array_reset (&pview);
      }
    }
    P4EST_FREE (mirror_proc_mirrors);
    P4EST_FREE (mirror_proc_offsets);
    ghost->mirror_proc_mirrors = mirror_proc_mirrors =
      new_mirror_proc_mirrors;
    ghost->mirror_proc_offsets = mirror_proc_offsets = newmpoffset;
    ghost->mirror_proc_fronts = mirror_proc_mirrors;
    ghost->mirror_proc_front_offsets = mirror_proc_offsets;

    sc_array_resize (mirrors, new_mirrors->elem_count);
    sc_array_copy (mirrors, new_mirrors);
    sc_array_destroy (new_mirrors);
    P4EST_FREE (newmpcount);
  }

  /* update the ghost layer */
  {
    int                 i, p;
    p4est_locidx_t     *new_proc_counts, *new_proc_offsets;
    p4est_lnodes_rank_t *lrank;

    mpiret =
      sc_MPI_Waitall ((int) recv_requests->elem_count,
                      (sc_MPI_Request *) recv_requests->array,
                      sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    sc_array_reset (recv_requests);
    sc_array_init (recv_requests, sizeof (sc_MPI_Request));
    new_proc_counts = P4EST_ALLOC (p4est_locidx_t, mpisize);
    for (p = 0; p < mpisize; p++) {
      new_proc_counts[p] = proc_offsets[p + 1] - proc_offsets[p];
    }

    for (i = 0; i < n_comm; i++) {
      sc_array_t         *recv_quads;
      p4est_locidx_t      recv_count;

      lrank = p4est_lnodes_rank_array_index_int (lnodes->sharers, i);

      p = lrank->rank;
      if (p == self) {
        continue;
      }

      recv_count = *((p4est_locidx_t *) sc_array_index_int (recv_counts, i));

      P4EST_LDEBUGF
        ("ghost layer support nodes receiving  %lld new from %d\n",
         (long long) recv_count, p);

      new_proc_counts[p] += recv_count;

      recv_quads = (sc_array_t *) sc_array_index (recv_all, i);
      sc_array_init_size (recv_quads, sizeof (p4est_quadrant_t),
                          (size_t) recv_count);
      if (recv_count) {
        sc_MPI_Request     *req;

        req = (sc_MPI_Request *) sc_array_push (recv_requests);

        mpiret =
          sc_MPI_Irecv (recv_quads->array,
                        recv_count * sizeof (p4est_quadrant_t), sc_MPI_BYTE,
                        p, P4EST_COMM_GHOST_SUPPORT_LOAD, p4est->mpicomm,
                        req);
        SC_CHECK_MPI (mpiret);
      }
    }

    new_proc_offsets = P4EST_ALLOC (p4est_locidx_t, mpisize + 1);
    new_proc_offsets[0] = 0;
    for (p = 0; p < mpisize; p++) {
      new_proc_offsets[p + 1] = new_proc_offsets[p] + new_proc_counts[p];
    }

    new_ghosts = sc_array_new_size (sizeof (p4est_quadrant_t),
                                    (size_t) new_proc_offsets[mpisize]);
    for (p = 0; p < mpisize; p++) {
      p4est_quadrant_t   *dest, *src;
      p4est_locidx_t      count;

      count = proc_offsets[p + 1] - proc_offsets[p];
      if (count) {
        dest = p4est_quadrant_array_index (new_ghosts, new_proc_offsets[p]);
        src = p4est_quadrant_array_index (ghosts, proc_offsets[p]);
        memcpy (dest, src, count * sizeof (p4est_quadrant_t));
      }
    }

    mpiret =
      sc_MPI_Waitall ((int) recv_requests->elem_count,
                      (sc_MPI_Request *) recv_requests->array,
                      sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    for (i = 0; i < n_comm; i++) {
      sc_array_t         *recv_quads;

      lrank = p4est_lnodes_rank_array_index_int (lnodes->sharers, i);

      p = lrank->rank;
      if (p == self) {
        continue;
      }

      recv_quads = (sc_array_t *) sc_array_index (recv_all, i);

      if (recv_quads->elem_count) {
        p4est_quadrant_t   *dest, *src;
        p4est_locidx_t      startidx, endidx, oldcount;
        sc_array_t          pview;

        startidx = new_proc_offsets[p];
        endidx = new_proc_offsets[p + 1];

        oldcount = proc_offsets[p + 1] - proc_offsets[p];

        dest = p4est_quadrant_array_index (new_ghosts, (size_t)
                                           (new_proc_offsets[p] + oldcount));

        src = p4est_quadrant_array_index (recv_quads, 0);

        memcpy (dest, src,
                recv_quads->elem_count * sizeof (p4est_quadrant_t));

        sc_array_init_view (&pview, new_ghosts, startidx,
                            (size_t) (endidx - startidx));
        sc_array_sort (&pview, p4est_quadrant_compare_piggy);
        sc_array_reset (&pview);
      }
    }

    if (new_ghosts->elem_count > ghosts->elem_count) {
      p4est_topidx_t      t;
      sc_array_t          split;

      /* update tree_offsets */
      sc_array_init (&split, sizeof (size_t));
      sc_array_split (new_ghosts, &split,
                      (size_t) conn->num_trees, ghost_tree_type, NULL);
      P4EST_ASSERT (split.elem_count == (size_t) conn->num_trees + 1);
      for (t = 0; t <= conn->num_trees; ++t) {
        size_t             *ppz;

        ppz = (size_t *) sc_array_index (&split, (size_t) t);
        tree_offsets[t] = *ppz;
      }
      sc_array_reset (&split);

      P4EST_ASSERT (sc_array_is_sorted
                    (new_ghosts, p4est_quadrant_compare_piggy));
#ifdef P4EST_ENABLE_DEBUG
      {
        size_t              zz;

        for (zz = 0; zz < ghosts->elem_count - 1; zz++) {
          p4est_quadrant_t   *q1 = p4est_quadrant_array_index (ghosts, zz);
          p4est_quadrant_t   *q2 =
            p4est_quadrant_array_index (ghosts, zz + 1);

          if (!p4est_quadrant_compare_piggy (q1, q2)) {
            P4EST_LERROR ("already duplicate in ghost:\n");
            p4est_quadrant_print (SC_LP_ERROR, q1);
          }
        }
        for (zz = 0; zz < new_ghosts->elem_count - 1; zz++) {
          p4est_quadrant_t   *q1 =
            p4est_quadrant_array_index (new_ghosts, zz);
          p4est_quadrant_t   *q2 =
            p4est_quadrant_array_index (new_ghosts, zz + 1);

          if (!p4est_quadrant_compare_piggy (q1, q2)) {
            P4EST_LERROR ("duplicate in new ghost:\n");
            p4est_quadrant_print (SC_LP_ERROR, q1);
          }
        }
      }
#endif
      sc_array_resize (ghosts, new_ghosts->elem_count);
      sc_array_copy (ghosts, new_ghosts);
    }
    P4EST_FREE (proc_offsets);
    ghost->proc_offsets = new_proc_offsets;
    sc_array_destroy (new_ghosts);
    P4EST_FREE (new_proc_counts);

    /* end sends */
    mpiret =
      sc_MPI_Waitall ((int) send_requests->elem_count,
                      (sc_MPI_Request *) send_requests->array,
                      sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    /* clean up */
    sc_array_destroy (recv_counts);
    sc_array_destroy (send_counts);
    for (i = 0; i < n_comm; i++) {
      sc_array_t         *quads;

      lrank = p4est_lnodes_rank_array_index_int (lnodes->sharers, i);

      p = lrank->rank;
      if (p == self) {
        continue;
      }

      quads = (sc_array_t *) sc_array_index (recv_all, i);
      sc_array_reset (quads);
      quads = (sc_array_t *) sc_array_index (send_all, i);
      sc_array_reset (quads);
    }
    sc_array_destroy (recv_all);
    sc_array_destroy (send_all);
    sc_array_destroy (recv_requests);
    sc_array_destroy (send_requests);
  }

  P4EST_ASSERT (p4est_ghost_is_valid (p4est, ghost));

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_ghost_support_lnodes\n");
#endif
}

p4est_lnodes_buffer_t *
p4est_lnodes_share_owned_begin (sc_array_t * node_data,
                                p4est_lnodes_t * lnodes)
{
  int                 mpiret;
  int                 p, proc;
  sc_array_t         *sharers = lnodes->sharers;
  int                 npeers = (int) sharers->elem_count;
  p4est_lnodes_rank_t *lrank;
  sc_array_t         *requests;
  sc_MPI_Request     *request;
  p4est_locidx_t      li, lz;
  void               *dest;
  sc_array_t         *send_bufs;
  sc_array_t         *send_buf;
  p4est_locidx_t      mine_offset, mine_count;
  size_t              elem_size = node_data->elem_size;
  p4est_lnodes_buffer_t *buffer;
  sc_MPI_Comm         comm = lnodes->mpicomm;
  int                 mpirank;

  P4EST_ASSERT (node_data->elem_count == (size_t) lnodes->num_local_nodes);

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);

  buffer = P4EST_ALLOC (p4est_lnodes_buffer_t, 1);

  buffer->requests = requests = sc_array_new (sizeof (sc_MPI_Request));
  buffer->send_buffers = send_bufs = sc_array_new (sizeof (sc_array_t));
  /* in this routine, the values from other processes are written directly
   * into node_data */
  buffer->recv_buffers = NULL;

  for (p = 0; p < npeers; p++) {
    lrank = p4est_lnodes_rank_array_index_int (sharers, p);
    proc = lrank->rank;
    if (proc == mpirank) {
      continue;
    }
    if (lrank->owned_count) {
      request = (sc_MPI_Request *) sc_array_push (requests);
      mpiret =
        sc_MPI_Irecv (node_data->array + elem_size * lrank->owned_offset,
                      (int) (lrank->owned_count * elem_size), sc_MPI_BYTE,
                      proc, P4EST_COMM_LNODES_OWNED, comm, request);
      SC_CHECK_MPI (mpiret);
    }
    mine_count = lrank->shared_mine_count;
    if (mine_count) {
      mine_offset = lrank->shared_mine_offset;
      send_buf = (sc_array_t *) sc_array_push (send_bufs);
      sc_array_init (send_buf, elem_size);
      sc_array_resize (send_buf, mine_count);
      for (li = 0; li < mine_count; li++) {
        lz = *((p4est_locidx_t *) sc_array_index (&lrank->shared_nodes,
                                                  (size_t) (li +
                                                            mine_offset)));
        dest = sc_array_index (send_buf, (size_t) li);
        memcpy (dest, node_data->array + elem_size * lz, elem_size);
      }
      request = (sc_MPI_Request *) sc_array_push (requests);
      mpiret = sc_MPI_Isend (send_buf->array, (int) (mine_count * elem_size),
                             sc_MPI_BYTE, proc, P4EST_COMM_LNODES_OWNED,
                             comm, request);
      SC_CHECK_MPI (mpiret);
    }
  }

  return buffer;
}

void
p4est_lnodes_share_owned_end (p4est_lnodes_buffer_t * buffer)
{
  int                 mpiret;
  size_t              zz;
  sc_array_t         *requests = buffer->requests;
  sc_array_t         *send_bufs = buffer->send_buffers;
  sc_array_t         *send_buf;

  P4EST_ASSERT (buffer->recv_buffers == NULL);

  if (requests->elem_count) {
    mpiret = sc_MPI_Waitall ((int) requests->elem_count,
                             (sc_MPI_Request *) requests->array,
                             sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  sc_array_destroy (requests);
  for (zz = 0; zz < send_bufs->elem_count; zz++) {
    send_buf = (sc_array_t *) sc_array_index (send_bufs, zz);
    sc_array_reset (send_buf);
  }
  sc_array_destroy (send_bufs);
  buffer->requests = NULL;
  buffer->send_buffers = NULL;
}

void
p4est_lnodes_share_owned (sc_array_t * array, p4est_lnodes_t * lnodes)
{
  p4est_lnodes_buffer_t *buffer;

  buffer = p4est_lnodes_share_owned_begin (array, lnodes);
  p4est_lnodes_share_owned_end (buffer);
  p4est_lnodes_buffer_destroy (buffer);
}

p4est_lnodes_buffer_t *
p4est_lnodes_share_all_begin (sc_array_t * node_data, p4est_lnodes_t * lnodes)
{
  int                 mpiret;
  int                 p, proc;
  sc_array_t         *sharers = lnodes->sharers;
  int                 npeers = (int) sharers->elem_count;
  p4est_lnodes_rank_t *lrank;
  p4est_lnodes_buffer_t *buffer;
  sc_array_t         *requests;
  sc_MPI_Request     *request;
  sc_array_t         *send_bufs;
  sc_array_t         *send_buf;
  sc_array_t         *recv_bufs;
  sc_array_t         *recv_buf;
  p4est_locidx_t      lz;
  void               *dest;
  size_t              zz;
  size_t              count;
  size_t              elem_size = node_data->elem_size;
  sc_MPI_Comm         comm = lnodes->mpicomm;
  int                 mpirank;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);

  P4EST_ASSERT (node_data->elem_count == (size_t) lnodes->num_local_nodes);

  buffer = P4EST_ALLOC (p4est_lnodes_buffer_t, 1);
  buffer->requests = requests = sc_array_new (sizeof (sc_MPI_Request));
  buffer->send_buffers = send_bufs = sc_array_new (sizeof (sc_array_t));
  buffer->recv_buffers = recv_bufs = sc_array_new (sizeof (sc_array_t));
  sc_array_resize (recv_bufs, (size_t) npeers);
  sc_array_resize (send_bufs, (size_t) npeers);

  for (p = 0; p < npeers; p++) {
    lrank = p4est_lnodes_rank_array_index_int (sharers, p);
    proc = lrank->rank;
    if (proc == mpirank) {
      /* there is no buffer for the current process: look in node_data
       * for values */
      recv_buf = (sc_array_t *) sc_array_index_int (recv_bufs, p);
      sc_array_init (recv_buf, elem_size);
      send_buf = (sc_array_t *) sc_array_index_int (send_bufs, p);
      sc_array_init (send_buf, elem_size);
      continue;
    }
    count = lrank->shared_nodes.elem_count;
    if (lrank->shared_nodes.elem_count) {
      recv_buf = (sc_array_t *) sc_array_index_int (recv_bufs, p);
      sc_array_init (recv_buf, elem_size);
      sc_array_resize (recv_buf, count);
      request = (sc_MPI_Request *) sc_array_push (requests);
      mpiret = sc_MPI_Irecv (recv_buf->array, (int) (count * elem_size),
                             sc_MPI_BYTE, proc, P4EST_COMM_LNODES_ALL,
                             comm, request);
      SC_CHECK_MPI (mpiret);

      send_buf = (sc_array_t *) sc_array_index_int (send_bufs, p);
      sc_array_init (send_buf, elem_size);
      sc_array_resize (send_buf, count);
      for (zz = 0; zz < count; zz++) {
        lz = *((p4est_locidx_t *) sc_array_index (&lrank->shared_nodes, zz));
        dest = sc_array_index (send_buf, zz);
        memcpy (dest, node_data->array + elem_size * lz, elem_size);
      }
      request = (sc_MPI_Request *) sc_array_push (requests);
      mpiret = sc_MPI_Isend (send_buf->array, (int) (count * elem_size),
                             sc_MPI_BYTE, proc, P4EST_COMM_LNODES_ALL,
                             comm, request);
      SC_CHECK_MPI (mpiret);
    }
  }

  return buffer;
}

void
p4est_lnodes_share_all_end (p4est_lnodes_buffer_t * buffer)
{
  int                 mpiret;
  size_t              zz;
  sc_array_t         *requests = buffer->requests;
  sc_array_t         *send_bufs = buffer->send_buffers;
  sc_array_t         *send_buf;

  if (requests->elem_count) {
    mpiret = sc_MPI_Waitall ((int) requests->elem_count,
                             (sc_MPI_Request *) requests->array,
                             sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  sc_array_destroy (requests);
  for (zz = 0; zz < send_bufs->elem_count; zz++) {
    send_buf = (sc_array_t *) sc_array_index (send_bufs, zz);
    sc_array_reset (send_buf);
  }
  sc_array_destroy (send_bufs);
  buffer->requests = NULL;
  buffer->send_buffers = NULL;
}

p4est_lnodes_buffer_t *
p4est_lnodes_share_all (sc_array_t * node_data, p4est_lnodes_t * lnodes)
{
  p4est_lnodes_buffer_t *buffer;

  buffer = p4est_lnodes_share_all_begin (node_data, lnodes);
  p4est_lnodes_share_all_end (buffer);

  return buffer;
}

void
p4est_lnodes_buffer_destroy (p4est_lnodes_buffer_t * buffer)
{
  int                 i;
  size_t              zz;
  sc_array_t         *requests = buffer->requests;
  sc_array_t         *send_bufs = buffer->send_buffers;
  sc_array_t         *recv_bufs = buffer->recv_buffers;
  sc_array_t         *bufs, *buf;

  if (requests != NULL) {
    P4EST_ASSERT (requests->elem_size == sizeof (sc_MPI_Request));
    sc_array_destroy (requests);
  }
  for (i = 0; i < 2; i++) {
    bufs = (i == 0) ? send_bufs : recv_bufs;
    if (bufs == NULL) {
      continue;
    }
    P4EST_ASSERT (bufs->elem_size == sizeof (sc_array_t));
    for (zz = 0; zz < bufs->elem_count; zz++) {
      buf = (sc_array_t *) sc_array_index (bufs, zz);
      sc_array_reset (buf);
    }
    sc_array_destroy (bufs);
  }
  P4EST_FREE (buffer);
}
