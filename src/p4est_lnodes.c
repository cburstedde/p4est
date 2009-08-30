/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2009 Carsten Burstedde, Lucas Wilcox.

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

#ifndef P4_TO_P8
#include <p4est_lnodes.h>
#include <p4est_iterate.h>
#include <p4est_communication.h>
#include <p4est_bits.h>
#else
#include <p8est_lnodes.h>
#include <p8est_iterate.h>
#include <p8est_communication.h>
#include <p8est_bits.h>
#endif

#ifdef SC_ALLGATHER
#include <sc_allgather.h>
#define MPI_Allgather sc_allgather
#endif

#ifndef P4_TO_P8
#define P4EST_LN_C_OFFSET 4
#else
#define P8EST_LN_E_OFFSET 6
#define P4EST_LN_C_OFFSET 18
#endif

/** cdp: corner dependent processes.
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
 * Suppose process q0, q1, and q2 are owned by process 0, 1, and 2 respectively,
 * while q3 and p are owned by process 3. Even though processes 0, 1, and 2 do
 * not touch the corner marked by "+", they share the node that is created
 * there.  When process 3 runs the callback that creates that node, it needs to
 * know that processes 0, 1, and 2 share that node, even though no quadrant
 * owned by those process is present in the callback.  Quadrant q3 has one cdp
 * for each corner: during the face callback between face p and the q faces,
 * the cdp associated with the marked corner "o" will have one of its face
 * values set to 0, one of its edge values set to 1, and one of its edge values
 * set to 2.  Then, when the corner callback is executed where "o" ans "+" meet,
 * process 3 can check the cdp associated with "o" to see that the other
 * processes share.
 */
typedef struct p4est_lnodes_cdp
{
  int                 face[P4EST_DIM];
#ifdef P4_TO_P8
  int                 edge[P4EST_DIM];
#endif
}
p4est_lnodes_cdp_t;

#ifdef P4_TO_P8
/* edp: edge dependent processes.
 * Using the same setup as above, when the far edge callback is executed,
 * quadrants q1, q3, and p are present, but q0 and q2, and thus processes i and
 * j, also share the nodes created on that edge.  The edp associated with q1's
 * edge will have one value set to 0, and the edp associated with q3's edge
 * will have one value set to 2.
 */
typedef struct p8est_lnodes_edp
{
  int                 face[2];
}
p8est_lnodes_edp_t;
#endif

/** buf_info: encodes/decodes the transmission of node information.
 * share_offset and share_count index into the inode_sharers array
 * for a list of all processes that share the nodes (local process included).
 */
typedef struct p4est_lnodes_buf_info
{
  p4est_quadrant_t    q;        /* p.which_tree filled */
  int8_t              type;     /* which nodes it shares */
  bool                send_sharers;     /* whether the sharers are included in
                                           the message */
  p4est_locidx_t      first_index;      /* inodes array, first node to/from */
  p4est_locidx_t      share_offset;
  int8_t              share_count;
}
p4est_lnodes_buf_info_t;

typedef struct p4est_lnodes_sorter
{
  p4est_locidx_t      local_index;
  p4est_locidx_t      inode_index;
}
p4est_lnodes_sorter_t;

static int
p4est_lnodes_sorter_compare (const void *a, const void *b)
{
  const p4est_lnodes_sorter_t *A = a;
  const p4est_lnodes_sorter_t *B = b;
  return (p4est_locidx_compare (&(A->local_index), &(B->local_index)));
}

typedef struct p4est_lnodes_data
{
  p4est_lnodes_cdp_t *local_cdp;        /* num local quads * corners per quad */
  p4est_lnodes_cdp_t *ghost_cdp;        /* num ghost quads * corners per quad */
#ifdef P4_TO_P8
  p8est_lnodes_edp_t *local_edp;        /* num local quads * edges per quad */
  p8est_lnodes_edp_t *ghost_edp;        /* num ghost quads * edges per quad */
#endif
  p4est_locidx_t     *local_elem_nodes; /* number of local q's * nodes per q */
  p4est_locidx_t     *ghost_elem_nodes; /* number of ghost q's * nodes per q */
  sc_array_t         *hfaces;   /* p4est_iter_face_side_t: hanging faces */
#ifdef P4_TO_P8
  sc_array_t         *hedges;   /* p8est_iter_edge_side_t: hanging edges */
#endif
  sc_array_t         *inodes;   /* p4est_locidx_t */
  sc_array_t         *inode_sharers;    /* int */
  sc_array_t         *send_buf_info;    /* one for each proc: type buf_info_t */
  sc_array_t         *recv_buf_info;    /* one for each proc: type buf_info_t */
  sc_array_t         *sorters;  /* one for each proc: type sorter_t */
#ifndef P4_TO_P8
  int8_t             *face_codes;
#else
  int16_t            *face_codes;
#endif
  int                 nodes_per_elem;
  int                 nodes_per_volume;
  int                *volume_nodes;
  int                 nodes_per_face;
  int                *face_nodes[2 * P4EST_DIM];
#ifdef P4_TO_P8
  int                 nodes_per_edge;
  int                *edge_nodes[12];
#endif
  int                 corner_nodes[P4EST_CHILDREN];
  p4est_gloidx_t     *global_offsets;
}
p4est_lnodes_data_t;

static void
p4est_quadrant_smallest_corner_descendent (p4est_quadrant_t * q,
                                           p4est_quadrant_t * r, int c)
{
  p4est_qcoord_t      shift = P4EST_QUADRANT_LEN (q->level) -
    P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
  r->x = q->x + (c % 2 ? shift : 0);
  r->y = q->y + ((c % 4) / 2 ? shift : 0);
#ifdef P4_TO_P8
  r->z = q->z + (c / 4 ? shift : 0);
#endif
  r->level = P4EST_QMAXLEVEL;
}

/** lnodes_face_simple_callback: runs even if there are no face nodes.
 * If a side of the face is not hanging, then there are no other quadrants that
 * are facewise dependent on its corner or edge nodes, so we set those cdp/edp
 * values to -2.
 * If a side of the face is hanging, we store the hanging face in hfaces, and we
 * set up the facewise cdp/edp values.
 */
static void
p4est_lnodes_face_simple_callback (p4est_iter_face_info_t * info, void *Data)
{
  int                 i;
  int                 c;
  int                 f;
  size_t              zz;
  p4est_lnodes_data_t *data = Data;
  sc_array_t         *sides = &(info->sides);
  size_t              count = sides->elem_count;
  p4est_iter_face_side_t *fside;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_topidx_t      tid;
  p4est_locidx_t      qid;
  p4est_lnodes_cdp_t *local_cdp = data->local_cdp;
  p4est_lnodes_cdp_t *ghost_cdp = data->ghost_cdp;
  p4est_locidx_t      cdpid;
#ifdef P4_TO_P8
  p8est_lnodes_edp_t *local_edp = data->local_edp;
  p8est_lnodes_edp_t *ghost_edp = data->ghost_edp;
  p4est_locidx_t      edpid;
  int                 j;
  int                 e;
  int                 he[2];
  int                 c2;
  int                 k;
  int                 eside;
#endif
  bool                is_local, *h_is_local;
  int                 procs[P4EST_CHILDREN / 2];
  p4est_iter_face_side_t *hface;
  sc_array_t         *hfaces = data->hfaces;
  int                 fdir;
  int                 rank = info->p4est->mpirank;
  p4est_locidx_t      hqid[P4EST_CHILDREN / 2];
  p4est_locidx_t      quadrants_offset;
#ifndef P4_TO_P8
  int8_t             *face_codes = data->face_codes;
#else
  int16_t            *face_codes = data->face_codes;
#endif
  bool                has_local;

  for (zz = 0; zz < count; zz++) {
    fside = sc_array_index (sides, zz);
    tid = fside->treeid;
    f = fside->face;
    fdir = f / 2;
    tree = p4est_array_index_topidx (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    if (fside->is_hanging) {
      has_local = false;
      h_is_local = fside->is.hanging.is_local;
      for (i = 0; i < P4EST_CHILDREN / 2; i++) {
        hqid[i] = fside->is.hanging.quadid[i];
        if (h_is_local[i]) {
          has_local = true;
          procs[i] = rank;
          hqid[i] += quadrants_offset;
#ifndef P4_TO_P8
          face_codes[hqid[i]] |=
            ((int8_t) p4est_face_corners[p4est_zface_to_rface[f]][i]);
          face_codes[hqid[i]] |= (0x01 << (2 + f / 2));
#else
          face_codes[hqid[i]] |= ((int16_t) p8est_face_corners[f][i]);
          face_codes[hqid[i]] |= (0x0001 << (6 + f / 2));
#endif
        }
        /* TODO: replace when the new ghost_layer API is implemented */
        else {
          procs[i] = p4est_comm_find_owner (info->p4est, tid,
                                            fside->is.hanging.quad[i], rank);
        }
      }
      if (has_local) {
        hface = sc_array_push (hfaces);
        *hface = *fside;
      }
      for (i = 0; i < P4EST_CHILDREN / 2; i++) {
#ifndef P4_TO_P8
        c = p4est_face_corners[p4est_zface_to_rface[f]][i];
#else
        c = p4est_face_corners[f][i];
#endif
        cdpid = hqid[i] * P4EST_CHILDREN + c;
        if (h_is_local[i]) {
          local_cdp[cdpid].face[fdir] = procs[P4EST_CHILDREN / 2 - 1 - i];
        }
        else {
          ghost_cdp[cdpid].face[fdir] = procs[P4EST_CHILDREN / 2 - 1 - i];
        }
#ifdef P4_TO_P8
        he[0] = p8est_corner_edges[c][(fdir + 1) % 3];
        he[1] = p8est_corner_edges[c][(fdir + 2) % 3];
        for (j = 0; j < 2; j++) {
          edpid = hqid[i] * 12 + he[j];
          c2 = p8est_edge_corners[he[1 - j]][0];
          if (c2 == c) {
            c2 = p8est_edge_corners[he[1 - j]][1];
          }
          k = p8est_corner_face_corners[c2][f];
          eside = (p8est_edge_faces[he[j]][0] == f) ? 0 : 1;
          if (h_is_local[i]) {
            local_edp[edpid].face[eside] = procs[k];
          }
          else {
            ghost_edp[edpid].face[eside] = procs[k];
          }
        }
#endif
      }
    }
    else {
      is_local = fside->is.full.is_local;
      qid = fside->is.full.quadid;
      if (is_local) {
        qid += quadrants_offset;
      }
      for (i = 0; i < P4EST_CHILDREN / 2; i++) {
#ifndef P4_TO_P8
        c = p4est_face_corners[p4est_zface_to_rface[f]][i];
#else
        c = p4est_face_corners[f][i];
#endif
        cdpid = qid * P4EST_CHILDREN + c;
#ifdef P4_TO_P8
        e = p8est_face_edges[f][i];
        edpid = qid * 12 + e;
        eside = (p8est_edge_faces[e][0] == f) ? 0 : 1;
#endif
        if (is_local) {
          local_cdp[cdpid].face[fdir] = -2;
#ifdef P4_TO_P8
          local_edp[edpid].face[eside] = -2;
#endif
        }
        else {
          ghost_cdp[cdpid].face[fdir] = -2;
#ifdef P4_TO_P8
          ghost_edp[edpid].face[eside] = -2;
#endif
        }
      }
    }
  }
}

#ifdef P4_TO_P8
/** lnodes_edge_simple_callback: runs even if there are no edge nodes.
 * If a side of the face is not hanging, then there are no other quadrants that
 * are facewise dependent on its corner or edge nodes, so we set those cdp/edp
 * values to -2.
 * If a side of the face is hanging, we store the hanging face in hfaces, and
 * we set up the facewise cdp/edp values.
 */
static void
p8est_lnodes_edge_simple_callback (p8est_iter_edge_info_t * info, void *Data)
{
  int                 i;
  int                 c;
  size_t              zz;
  p4est_lnodes_data_t *data = Data;
  sc_array_t         *sides = &(info->sides);
  size_t              count = sides->elem_count;
  p8est_iter_edge_side_t *eside;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_topidx_t      tid;
  p4est_locidx_t      qid;
  p4est_lnodes_cdp_t *local_cdp = data->local_cdp;
  p4est_lnodes_cdp_t *ghost_cdp = data->ghost_cdp;
  p4est_locidx_t      cdpid;
  int                 e;
  bool                is_local, *h_is_local;
  int                 procs[2];
  p8est_iter_edge_side_t *hedge;
  sc_array_t         *hedges = data->hedges;
  int                 edir;
  int                 rank = info->p4est->mpirank;
  p4est_locidx_t      hqid[2];
  p4est_locidx_t      quadrants_offset;
  int16_t            *face_codes = data->face_codes;
  bool                has_local;

  for (zz = 0; zz < count; zz++) {
    eside = sc_array_index (sides, zz);
    tid = eside->treeid;
    e = eside->edge;
    edir = e / 4;
    tree = p4est_array_index_topidx (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    if (eside->is_hanging) {
      has_local = false;
      h_is_local = eside->is.hanging.is_local;
      for (i = 0; i < 2; i++) {
        hqid[i] = eside->is.hanging.quadid[i];
        if (h_is_local[i]) {
          has_local = true;
          procs[i] = rank;
          hqid[i] += quadrants_offset;
          face_codes[hqid[i]] |= ((int16_t) p8est_edge_corners[e][i]);
          face_codes[hqid[i]] |= (0x0001 << (3 + e / 4));
        }
        /* TODO: replace when the new ghost_layer API is implemented */
        else {
          procs[i] = p4est_comm_find_owner (info->p4est, tid,
                                            eside->is.hanging.quad[i], rank);
        }
      }
      if (has_local) {
        hedge = sc_array_push (hedges);
        *hedge = *eside;
      }
      for (i = 0; i < 2; i++) {
        c = p8est_edge_corners[e][i];
        cdpid = hqid[i] * P4EST_CHILDREN + c;
        if (h_is_local[i]) {
          local_cdp[cdpid].edge[edir] = procs[1 - i];
        }
        else {
          ghost_cdp[cdpid].edge[edir] = procs[1 - i];
        }
      }
    }
    else {
      is_local = eside->is.full.is_local;
      qid = eside->is.full.quadid;
      if (is_local) {
        qid += quadrants_offset;
      }
      for (i = 0; i < 2; i++) {
        c = p8est_edge_corners[e][i];
        cdpid = qid * P4EST_CHILDREN + c;
        if (is_local) {
          local_cdp[cdpid].edge[edir] = -2;
        }
        else {
          ghost_cdp[cdpid].edge[edir] = -2;
        }
      }
    }
  }
}
#endif

/** p4est_lnodes_missing_proc_cdp_face:
 * Quadrant q touches corner c. Face f = p4est_corner_faces[c][dir].
 *
 *                         Case 1                   Case 2
 *       _______   |      ____________   |    |      _______   |
 *      /       /  |     /            /  |    |     /       /  |
 *     /______ /|p |    /            /|  |    |    /______ /|? |
 *    /       /f|  |   /            / |  |    |   /       / |  |
 *   /______ /| /  |  /___________ /  |  |    |  /______ /| /  |
 *   |      |f|/|  |  |           |   |  | or |  |      | |/|? |
 *   |      | /f|  |  |           |   |  |    |  |  ?   | / |  |
 *   |______|/| /  |  |     ?     |   /  |    |  |______|/| /  |
 *   |      |f|/   |  |           |  /   |    |  |      | |/   |
 *   |  q   | /    |  |           | /    |    |  |  ?   | /    |
 *   |______|/     |  |___________|/     |    |  |______|/     |
 *          c         *                          *
 *
 * We do not know which process owns p.  If opposite f is case 1, then quadrant
 * p shares in the node at c, and so does the process that owns it.  In case 2,
 * quadrants p does not share in the node at c.
 * \return the process that owns p in case 1, -2 otherwise.
 */
static int
p4est_lnodes_missing_proc_cdp_face (p4est_quadrant_t * q, p4est_topidx_t tid,
                                    int c, int dir,
                                    p4est_iter_corner_info_t * info)
{
  int                 rank = info->p4est->mpirank;
  p4est_quadrant_t    tempq, tempr;
#ifndef P4_TO_P8
  int                 f =
    (p4est_rface_to_zface[p4est_corner_faces[c][dir]] / 2 == dir) ?
    p4est_corner_faces[c][dir] : p4est_corner_faces[c][1 - dir];
  int                 j = p4est_face_corners[f][0] == c ? 0 : 1;
#else
  int                 f = p8est_corner_faces[c][dir];
  int                 j = p8est_corner_face_corners[c][f];
#endif
  int                 c2 = p4est_face_corners[f][P4EST_CHILDREN / 2 - 1 - j];
  int                 nproc;
  int                 nc;
  p4est_locidx_t      ntid;
  sc_array_t         *sides = &(info->sides);
  size_t              zz;
  size_t              count = sides->elem_count;
  p4est_t            *p4est = info->p4est;
  p4est_iter_corner_side_t *cside;
  int8_t              level, nlevel;

  if (q->level == 0) {
    return -2;
  }

  /* tempq == p in the diagram above */
  p4est_quadrant_sibling (q, &tempq, c2);
  nproc = p4est_comm_find_owner (p4est, tid, &tempq, rank);

  ntid = p4est_quadrant_face_neighbor_extra (q, tid, f, &tempr,
                                             p4est->connectivity);
  if (ntid == -1) {
    return -2;
  }
  nc = p4est_quadrant_child_id (&tempr);
  p4est_quadrant_parent (&tempr, &tempq);
  nlevel = tempq.level;

  for (zz = 0; zz < count; zz++) {
    cside = sc_array_index (sides, zz);
    if (cside->treeid != ntid || cside->corner != nc) {
      continue;
    }
    level = cside->quad->level;
    P4EST_ASSERT (nlevel <= level && level <= nlevel + 2);
    return (level == nlevel) ? nproc : -2;
  }

  /* not reached because the face neighbor of q, if it exists, should be in the
   * iterate corner info */
  SC_ABORT_NOT_REACHED ();

  return -2;
}

#ifdef P4_TO_P8
/** p8est_lnodes_missing_proc_cdp_edge:
 * Quadrant q touches corner c. Edge e = p8est_corner_edges[c][dir].
 *
 *                         Case 1                   Case 2
 *       _______   |      ____________   |    |      _______   |
 *      /       /  |     /            /  |    |     /       /  |
 *     /______ /|  |    /            /|  |    |    /______ /|? |
 *    /       / |  |   /            / |  |    |   /       / |  |
 *   /______ /| /  |  /___________ /  |  |    |  /______ /| /  |
 *   |      e |/|  |  +           |   |  | or |  +      | |/|? |
 *   |  p   e / |  |  +           |   |  |    |  +  ?   | / |  |
 *   |______e/| /  |  +     ?     |   /  |    |  +______|/| /  |
 *   |      e |/   |  +           |  /   |    |  +      | |/   |
 *   |  q   e /    |  +           | /    |    |  +  ?   | /    |
 *   |______e/     |  +___________|/     |    |  +______|/     |
 *          c         *                          *
 *
 * We do not know which process owns p.  If also touching e is case 1, then
 * quadrant p shares in the node at c, and so does the process that owns it.
 * In case 2, quadrants p does not share in the node at c.
 * \return the process that owns p in case 1, -2 otherwise.
 */
static int
p8est_lnodes_missing_proc_cdp_edge (p4est_quadrant_t * q, p4est_topidx_t tid,
                                    int c, int dir,
                                    p4est_iter_corner_info_t * info)
{
  int                 i;
  int                 rank = info->p4est->mpirank;
  p4est_quadrant_t    tempq, tempr;
  p4est_quadrant_t   *ptemp;
  int                 f;
  int                 e = p8est_corner_edges[c][dir];
  int                 j = p8est_edge_corners[e][0] == c ? 0 : 1;
  int                 c2 = p8est_edge_corners[e][1 - j];
  int                 nproc;
  int                 nc;
  p4est_locidx_t      ntid;
  sc_array_t         *sides = &(info->sides);
  size_t              zz, zy;
  size_t              count = sides->elem_count;
  size_t              count2;
  p4est_t            *p4est = info->p4est;
  p4est_iter_corner_side_t *cside;
  int8_t              level, nlevel;
  sc_array_t          quads, treeids;

  if (q->level == 0) {
    return -2;
  }

  /* tempq is now p in the diagram */
  p4est_quadrant_sibling (q, &tempq, c2);
  nproc = p4est_comm_find_owner (p4est, tid, &tempq, rank);

  for (i = 0; i < 2; i++) {
    f = p8est_edge_faces[e][i];
    ntid = p4est_quadrant_face_neighbor_extra (q, tid, f, &tempr,
                                               p4est->connectivity);
    if (ntid == -1) {
      continue;
    }
    nc = p4est_quadrant_child_id (&tempr);
    p4est_quadrant_parent (&tempr, &tempq);
    nlevel = tempq.level;

    for (zz = 0; zz < count; zz++) {
      cside = sc_array_index (sides, zz);
      if (cside->treeid != ntid || cside->corner != nc) {
        continue;
      }
      level = cside->quad->level;
      P4EST_ASSERT (nlevel <= level && level <= nlevel + 2);
      if (level == nlevel) {
        return nproc;
      }
    }
  }

  sc_array_init (&quads, sizeof (p4est_quadrant_t));
  sc_array_init (&treeids, sizeof (p4est_topidx_t));
  p8est_quadrant_edge_neighbor_extra (q, tid, e, &quads, &treeids,
                                      p4est->connectivity);
  count2 = quads.elem_count;
  for (zy = 0; zy < count2; zy++) {
    ptemp = sc_array_index (&quads, zy);
    ntid = *((p4est_topidx_t *) sc_array_index (&treeids, zy));
    nc = p4est_quadrant_child_id (ptemp);
    p4est_quadrant_parent (ptemp, &tempq);
    nlevel = tempq.level;
    for (zz = 0; zz < count; zz++) {
      cside = sc_array_index (sides, zz);
      if (cside->treeid != ntid || cside->corner != nc) {
        continue;
      }
      level = cside->quad->level;
      P4EST_ASSERT (nlevel <= level && level <= nlevel + 2);
      if (level == nlevel) {
        sc_array_reset (&quads);
        sc_array_reset (&treeids);
        return nproc;
      }
    }
  }
  sc_array_reset (&quads);
  sc_array_reset (&treeids);

  return -2;
}

/** p8est_lnodes_missing_proc_edp:
 * Quadrant q touches edge e. f = p8est_edge_faces[e][dir], and
 * p8est_edge_corners[e][pos] = child_id(q).
 *
 *                         Case 1                   Case 2
 *       _______   |      ____________   |    |      _______   |
 *      /       /  |     /            /  |    |     /       /  |
 *     /______ /|  |    /            /|  |    |    /______ /|? |
 *    /       /f|  |   /            / |  |    |   /       / |  |
 *   /______ /| /  |  /___________ /  |  |    |  /______ /| /  |
 *   |      ef|/|p |  +           |   |  | or |  +      | |/|? |
 *   |      e /f|  |  +           |   |  |    |  +  ?   | / |  |
 *   |______e/| /  |  +     ?     |   /  |    |  +______|/| /  |
 *   |      ef|/   |  +           |  /   |    |  +      | |/   |
 *   |  q   e /    |  +           | /    |    |  +  ?   | /    |
 *   |______e/     |  +___________|/     |    |  +______|/     |
 *
 * We do not know which process owns p.  If opposite f is case 1, then
 * quadrant p shares in the nodes at e, and so does the process that owns it.
 * In case 2, quadrants p does not share in the nodes at e.
 * \return the process that owns p in case 1, -2 otherwise.
 */
static int
p8est_lnodes_missing_proc_edp (p4est_quadrant_t * q, p4est_topidx_t tid,
                               int e, int dir, int pos,
                               p8est_iter_edge_info_t * info)
{
  int                 rank = info->p4est->mpirank;
  p4est_quadrant_t    tempq, tempr, temps;
  int                 f = p8est_edge_faces[e][dir];
  int                 c = p8est_edge_corners[e][pos];
  int                 j = p8est_corner_face_corners[c][f];
  int                 k =
    p8est_face_corners[f][j ^ 1] ==
    p8est_edge_corners[e][1 - pos] ? (j ^ 2) : (j ^ 1);
  int                 c2 = p8est_face_corners[f][k];
  int                 nproc;
  int                 nc, nc2;
  int                 ne;
  p4est_locidx_t      ntid;
  sc_array_t         *sides = &(info->sides);
  size_t              zz;
  size_t              count = sides->elem_count;
  p4est_t            *p4est = info->p4est;
  p8est_iter_edge_side_t *eside;

  if (q->level == 0) {
    return -2;
  }

  /* tempq is now p in the diagram */
  p4est_quadrant_sibling (q, &tempq, c2);
  nproc = p4est_comm_find_owner (p4est, tid, &tempq, rank);

  ntid = p4est_quadrant_face_neighbor_extra (q, tid, f, &tempr,
                                             p4est->connectivity);
  if (ntid == -1) {
    return -2;
  }
  nc = p4est_quadrant_child_id (&tempr);
  p4est_quadrant_sibling (q, &tempq, p8est_edge_corners[e][1 - pos]);
  ntid = p4est_quadrant_face_neighbor_extra (&tempq, tid, f, &temps,
                                             p4est->connectivity);
  P4EST_ASSERT (ntid != -1);
  nc2 = p4est_quadrant_child_id (&temps);
  ne = p8est_child_corner_edges[nc][nc2];
  for (zz = 0; zz < count; zz++) {
    eside = sc_array_index (sides, zz);
    if (eside->treeid != ntid || eside->edge != ne) {
      continue;
    }
    return (!eside->is_hanging) ? nproc : -2;
  }

  /* not reached because the face neighbor should be in the edge iterate info
   */
  SC_ABORT_NOT_REACHED ();

  return -2;
}
#endif

static inline void
p4est_lnodes_push_binfo (sc_array_t * touch, sc_array_t * all,
                         sc_array_t * send, sc_array_t * recv,
                         sc_array_t * share, int owner, int rank,
                         p4est_quadrant_t * q, p4est_locidx_t tid,
                         int8_t type, p4est_locidx_t nin)
{
  size_t              zz, count = all->elem_count;
  int                *ip, proc;
  p4est_lnodes_buf_info_t *binfo;
  int8_t              scount;
  p4est_locidx_t      offset = (p4est_locidx_t) share->elem_count;

  ip = sc_array_push (share);
  *ip = rank;
  scount = (int8_t) (count + 1);
  for (zz = 0; zz < count; zz++) {
    ip = sc_array_push (share);
    proc = *((int *) sc_array_index (all, zz));
    *ip = proc;
    if (owner == rank) {
      binfo = sc_array_push (&(send[proc]));
      binfo->send_sharers = true;
      if (touch == NULL ||
          sc_array_bsearch (touch, &proc, sc_int_compare) >= 0) {
        binfo->send_sharers = false;
      }
    }
    else if (proc == owner) {
      binfo = sc_array_push (&(recv[proc]));
      binfo->send_sharers = false;
    }
    else {
      continue;
    }
    binfo->q = *q;
    binfo->q.p.which_tree = tid;
    binfo->type = type;
    binfo->first_index = nin;
    binfo->share_offset = offset;
    binfo->share_count = scount;
  }
}

/* p4est_lnodes_corner_callback:
 *
 * Create a new independent node at a corner.  Set all touching element nodes
 * to point to the newly created independent node.
 * Compute all processes that share the node.
 * If the node is locally owned, add info describing the node to the send
 * buffer of all processes that need the node.
 * If the node is not locally owned, add info describing the node to the 
 * receive buffer of the owner.
 */
static void
p4est_lnodes_corner_callback (p4est_iter_corner_info_t * info, void *Data)
{
  int                 i;
  size_t              zz;
  p4est_lnodes_data_t *data = Data;
  sc_array_t         *sides = &(info->sides);
  size_t              count = sides->elem_count;
  p4est_iter_corner_side_t *cside;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *inode;
  sc_array_t         *inode_sharers = data->inode_sharers;
  p4est_lnodes_cdp_t *local_cdp = data->local_cdp;
  p4est_lnodes_cdp_t *ghost_cdp = data->ghost_cdp;
  p4est_lnodes_cdp_t *cdp;
  p4est_locidx_t     *local_elem_nodes = data->local_elem_nodes;
  p4est_locidx_t     *ghost_elem_nodes = data->ghost_elem_nodes;
  p4est_locidx_t     *elnode;
  sc_array_t         *send_buf_info = data->send_buf_info;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  sc_array_t          touching_procs;
  sc_array_t          all_procs;
  int                *ip;
  p4est_topidx_t      tid;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_locidx_t      qid;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) inodes->elem_count;
  bool                is_local;
  int                 c;
  p4est_locidx_t      cdpid;
  p4est_locidx_t      nid;
  int                 proc;
  int                 rank = info->p4est->mpirank;
  p4est_topidx_t      owner_tid = P4EST_TOPIDX_MAX;
  int                 owner_corner = -1;
  p4est_quadrant_t   *owner_quad = NULL;
  p4est_quadrant_t   *q;
  p4est_quadrant_t    tempq;
  int                 owner_proc = INT_MAX;
  int                *corner_nodes = data->corner_nodes;
  int                 nodes_per_elem = data->nodes_per_elem;
  p4est_locidx_t      quadrants_offset;
  bool                intra_tree = info->intra_tree;
  int8_t              type;

  if (intra_tree) {
    cside = sc_array_index (sides, 0);
    owner_corner = cside->corner;
    owner_tid = cside->treeid;
    owner_quad = cside->quad;
    is_local = cside->is_local;
    if (is_local) {
      owner_proc = rank;
    }
    else {
      /* TODO: replace when the new ghost_layer API is implemented */
      owner_proc = p4est_comm_find_owner (info->p4est, owner_tid, owner_quad,
                                          rank);
      P4EST_ASSERT (owner_proc != rank);
    }
  }

  sc_array_init (&touching_procs, sizeof (int));
  sc_array_init (&all_procs, sizeof (int));
  for (zz = 0; zz < count; zz++) {
    cside = sc_array_index (sides, zz);
    c = cside->corner;
    tid = cside->treeid;
    qid = cside->quadid;
    q = cside->quad;
    is_local = cside->is_local;
    tree = p4est_array_index_topidx (trees, tid);
    quadrants_offset = tree->quadrants_offset;

    if (is_local) {
      proc = rank;
      qid += quadrants_offset;
    }
    else {
      /* TODO: replace when the new ghost_layer API is implemented */
      proc = p4est_comm_find_owner (info->p4est, tid, q, rank);
      P4EST_ASSERT (proc != rank);
      ip = sc_array_push (&touching_procs);
      *ip = proc;
      ip = sc_array_push (&all_procs);
      *ip = proc;
    }

    if (!intra_tree && tid <= owner_tid) {
      if (tid < owner_tid || p4est_quadrant_compare (q, owner_quad) < 0) {
        owner_proc = proc;
        owner_tid = tid;
        owner_quad = q;
        owner_corner = c;
      }
    }
#ifdef P4EST_DEBUG
    else {
      P4EST_ASSERT (zz == 0 || owner_tid < tid ||
                    (owner_tid == tid &&
                     p4est_quadrant_compare (owner_quad, q) < 0));
    }
#endif

    cdpid = qid * P4EST_CHILDREN + c;
    cdp = is_local ? &(local_cdp[cdpid]) : &(ghost_cdp[cdpid]);
    for (i = 0; i < P4EST_DIM; i++) {
      proc = cdp->face[i];
      if (proc == -1) {
        proc = p4est_lnodes_missing_proc_cdp_face (q, tid, c, i, info);
      }
      if (proc >= 0 && proc != rank) {
        ip = sc_array_push (&all_procs);
        *ip = proc;
      }
#ifdef P4_TO_P8
      proc = cdp->edge[i];
      if (proc == -1) {
        proc = p8est_lnodes_missing_proc_cdp_edge (q, tid, c, i, info);
      }
      if (proc >= 0 && proc != rank) {
        ip = sc_array_push (&all_procs);
        *ip = proc;
      }
#endif
    }

    nid = qid * nodes_per_elem + corner_nodes[c];
    elnode = is_local ? &(local_elem_nodes[nid]) : &(ghost_elem_nodes[nid]);
    *elnode = num_inodes;
  }
  sc_array_sort (&touching_procs, sc_int_compare);
  sc_array_uniq (&touching_procs, sc_int_compare);
  sc_array_sort (&all_procs, sc_int_compare);
  sc_array_uniq (&all_procs, sc_int_compare);

  inode = sc_array_push (inodes);
  *inode = -((p4est_locidx_t) owner_proc + 1);

  count = all_procs.elem_count;
  if (count) {
    P4EST_QUADRANT_INIT (&tempq);
    p4est_quadrant_smallest_corner_descendent (owner_quad, &tempq,
                                               owner_corner);
    type = (int8_t) (P4EST_LN_C_OFFSET + owner_corner);
    p4est_lnodes_push_binfo (&touching_procs, &all_procs, send_buf_info,
                             recv_buf_info, inode_sharers, owner_proc, rank,
                             &tempq, owner_tid, type, num_inodes);
  }

  sc_array_reset (&touching_procs);
  sc_array_reset (&all_procs);
}

#ifdef P4_TO_P8
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
  int                 i;
  int                 j;
  size_t              zz;
  p4est_lnodes_data_t *data = Data;
  sc_array_t         *sides = &(info->sides);
  size_t              count = sides->elem_count;
  p8est_iter_edge_side_t *eside;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *inode;
  sc_array_t         *inode_sharers = data->inode_sharers;
  p8est_lnodes_edp_t *local_edp = data->local_edp;
  p8est_lnodes_edp_t *ghost_edp = data->ghost_edp;
  p8est_lnodes_edp_t *edp;
  p4est_locidx_t     *local_elem_nodes = data->local_elem_nodes;
  p4est_locidx_t     *ghost_elem_nodes = data->ghost_elem_nodes;
  p4est_locidx_t     *elnode;
  sc_array_t         *send_buf_info = data->send_buf_info;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  sc_array_t          touching_procs;
  sc_array_t          all_procs;
  int                *ip;
  p4est_topidx_t      tid;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_locidx_t      quadrants_offset;
  p4est_locidx_t     *qids;
  p4est_locidx_t      qid;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) inodes->elem_count;
  bool               *is_local;
  int                 e;
  p4est_locidx_t      edpid;
  p4est_locidx_t      nid;
  int                 proc;
  int                 rank = info->p4est->mpirank;
  p4est_topidx_t      owner_tid = P4EST_TOPIDX_MAX;
  int                 owner_orientation = -1;
  int                 owner_edge = -1;
  p4est_quadrant_t   *owner_quad = NULL;
  p4est_quadrant_t  **q;
  p4est_quadrant_t    tempq, tempr;
  int                 owner_proc = INT_MAX;
  int                 nodes_per_edge = data->nodes_per_edge;
  int                 nodes_per_elem = data->nodes_per_elem;
  int               **edge_nodes = data->edge_nodes;
  bool                is_hanging;
  int                 limit;
  int                 stride;
  int                 orientation;
  int8_t              min_level = P4EST_QMAXLEVEL;
  int8_t              max_level = 0;
  p4est_locidx_t      start_node;
  bool                intra_tree = info->intra_tree;
  int8_t              type;

  p8est_lnodes_edge_simple_callback (info, data);

  if (intra_tree) {
    eside = sc_array_index (sides, 0);
    owner_edge = eside->edge;
    owner_tid = eside->treeid;
    owner_orientation = eside->orientation;
    if (eside->is_hanging) {
      owner_quad = eside->is.hanging.quad[0];
      is_local = eside->is.hanging.is_local;
    }
    else {
      owner_quad = eside->is.full.quad;
      is_local = &(eside->is.full.is_local);
    }
    if (*is_local) {
      owner_proc = rank;
    }
    else {
      /* TODO: replace when the new ghost_layer API is implemented */
      owner_proc = p4est_comm_find_owner (info->p4est, owner_tid, owner_quad,
                                          rank);
      P4EST_ASSERT (owner_proc != rank);
    }
  }

  sc_array_init (&touching_procs, sizeof (int));
  sc_array_init (&all_procs, sizeof (int));
  for (zz = 0; zz < count; zz++) {
    eside = sc_array_index (sides, zz);
    e = eside->edge;
    tid = eside->treeid;
    tree = p4est_array_index_topidx (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    is_hanging = eside->is_hanging;
    orientation = eside->orientation;
    if (!is_hanging) {
      limit = 1;
      is_local = &(eside->is.full.is_local);
      qids = &(eside->is.full.quadid);
      q = &(eside->is.full.quad);
      min_level = eside->is.full.quad->level;
      max_level = min_level + 1;
    }
    else {
      limit = 2;
      is_local = eside->is.hanging.is_local;
      qids = eside->is.hanging.quadid;
      q = eside->is.hanging.quad;
    }

    for (i = 0; i < limit; i++) {
      qid = qids[i];
      if (is_local[i]) {
        proc = rank;
        qid += quadrants_offset;
      }
      else {
        /* TODO: replace when the new ghost_layer API is implemented */
        proc = p4est_comm_find_owner (info->p4est, tid, q[i], rank);
        P4EST_ASSERT (proc != rank);
        ip = sc_array_push (&touching_procs);
        *ip = proc;
        ip = sc_array_push (&all_procs);
        *ip = proc;
      }

      if (!intra_tree && i == 0 && tid <= owner_tid) {
        if (tid < owner_tid || p4est_quadrant_compare (q[i], owner_quad) < 0) {
          owner_proc = proc;
          owner_tid = tid;
          owner_quad = q[i];
          owner_edge = e;
          owner_orientation = orientation;
        }
      }
#ifdef P4EST_DEBUG
      else {
        P4EST_ASSERT ((zz == 0 && i == 0) || owner_tid < tid ||
                      (owner_tid == tid &&
                       p4est_quadrant_compare (owner_quad, q[i]) < 0));
      }
#endif

      if (is_hanging) {
        edpid = qid * 12 + e;
        edp = is_local[i] ? &(local_edp[edpid]) : &(ghost_edp[edpid]);
        for (j = 0; j < 2; j++) {
          proc = edp->face[j];
          if (proc == -1) {
            proc = p8est_lnodes_missing_proc_edp (q[i], tid, e, j, i, info);
          }
          if (proc >= 0 && proc != rank) {
            ip = sc_array_push (&all_procs);
            *ip = proc;
          }
        }
      }
    }
  }

  for (zz = 0; zz < count; zz++) {
    eside = sc_array_index (sides, zz);
    e = eside->edge;
    tid = eside->treeid;
    tree = p4est_array_index_topidx (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    is_hanging = eside->is_hanging;
    orientation = eside->orientation;
    if (!is_hanging) {
      limit = 1;
      is_local = &(eside->is.full.is_local);
      qids = &(eside->is.full.quadid);
    }
    else {
      limit = 2;
      is_local = eside->is.hanging.is_local;
      qids = eside->is.hanging.quadid;
    }

    for (i = 0; i < limit; i++) {
      qid = qids[i];
      if (is_local[i]) {
        qid += quadrants_offset;
      }
      start_node = num_inodes + (orientation == owner_orientation ?
                                 0 : nodes_per_edge - 1);
      stride = (orientation == owner_orientation ? 1 : -1);
      for (j = 0; j < nodes_per_edge; j++, start_node += stride) {
        nid = qid * nodes_per_elem + edge_nodes[e][j];
        elnode = is_local[i] ? &(local_elem_nodes[nid]) :
          &(ghost_elem_nodes[nid]);
        *elnode = start_node;
      }
    }
  }
  P4EST_ASSERT (max_level == min_level + 1);
  sc_array_sort (&touching_procs, sc_int_compare);
  sc_array_uniq (&touching_procs, sc_int_compare);
  sc_array_sort (&all_procs, sc_int_compare);
  sc_array_uniq (&all_procs, sc_int_compare);

  count = all_procs.elem_count;
  if (count) {
    P4EST_QUADRANT_INIT (&tempq);
    if (owner_quad->level == max_level || min_level == P4EST_QMAXLEVEL) {
      tempq = *owner_quad;
    }
    else {
      P4EST_ASSERT (owner_quad->level == min_level);
      tempr = *owner_quad;
      tempr.level++;
      p4est_quadrant_sibling (&tempr, &tempq,
                              p8est_edge_corners[owner_edge][0]);
    }
    type = (int8_t) (P8EST_LN_E_OFFSET + owner_edge);
    p4est_lnodes_push_binfo (&touching_procs, &all_procs, send_buf_info,
                             recv_buf_info, inode_sharers, owner_proc, rank,
                             &tempq, owner_tid, type, num_inodes);
  }
  for (i = 0; i < nodes_per_edge; i++) {
    inode = sc_array_push (inodes);
    *inode = -((p4est_locidx_t) owner_proc + 1);
  }

  sc_array_reset (&touching_procs);
  sc_array_reset (&all_procs);
}

/** p8est_lnodes_face_node_transform:
 *
 * Compute the transformation from an independent face nodes position on a face
 * to the position on the face of the touching element.
 */
static inline void
p8est_lnodes_face_node_transform (int orig_f, int f, int orientation,
                                  bool * flipj, bool * flipk, bool * swapjk)
{
  int                 ref = p8est_face_permutation_refs[f][orig_f];
  int                 set = p8est_face_permutation_sets[ref][orientation];
  int                 c0 = p8est_face_permutations[set][0];
  int                 c1 = p8est_face_permutations[set][1];
  int                 c2 = p8est_face_permutations[set][2];
  *flipj = (c1 < c0);
  *flipk = (c2 < c0);
  *swapjk = ((c0 ^ c1) == 1);
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
  size_t              zz;
  p4est_lnodes_data_t *data = Data;
  sc_array_t         *sides = &(info->sides);
  size_t              count = sides->elem_count;
  p4est_iter_face_side_t *fside;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *inode;
  sc_array_t         *inode_sharers = data->inode_sharers;
  p4est_locidx_t     *local_elem_nodes = data->local_elem_nodes;
  p4est_locidx_t     *ghost_elem_nodes = data->ghost_elem_nodes;
  p4est_locidx_t     *elnode;
  sc_array_t         *send_buf_info = data->send_buf_info;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  sc_array_t          touching_procs;
  int                *ip;
  p4est_topidx_t      tid;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_locidx_t      quadrants_offset;
  p4est_locidx_t     *qids;
  p4est_locidx_t      qid;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) inodes->elem_count;
  bool               *is_local;
  int                 f;
  p4est_locidx_t      nid;
  int                 proc;
  int                 rank = info->p4est->mpirank;
  p4est_topidx_t      owner_tid = P4EST_TOPIDX_MAX;
  int                 owner_face = -1;
  p4est_quadrant_t   *owner_quad = NULL;
  size_t              owner_side = 2;
  p4est_quadrant_t  **q;
  p4est_quadrant_t    tempq;
  int                 owner_proc = INT_MAX;
  int                 nodes_per_face = data->nodes_per_face;
  int                 nodes_per_elem = data->nodes_per_elem;
  int               **face_nodes = data->face_nodes;
  bool                is_hanging;
  int                 limit;
  int                 i, j;
#ifndef P4_TO_P8
  int                 stride;
#else
  int                 k, l;
  int                 nodes_per_edge = data->nodes_per_edge;
  bool                flipj, flipk, swapjk;
  int                 jind, kind, lind;
#endif
  p4est_locidx_t      start_node;
  bool                intra_tree = info->intra_tree;
  int8_t              type;

  p4est_lnodes_face_simple_callback (info, data);

  if (intra_tree) {
    fside = sc_array_index (sides, 0);
    owner_face = fside->face;
    owner_tid = fside->treeid;
    owner_side = 0;
    if (fside->is_hanging) {
      owner_quad = fside->is.hanging.quad[0];
      is_local = fside->is.hanging.is_local;
    }
    else {
      owner_quad = fside->is.full.quad;
      is_local = &(fside->is.full.is_local);
    }
    if (*is_local) {
      owner_proc = rank;
    }
    else {
      /* TODO: replace when the new ghost_layer API is implemented */
      owner_proc = p4est_comm_find_owner (info->p4est, owner_tid, owner_quad,
                                          rank);
      P4EST_ASSERT (owner_proc != rank);
    }
  }

  sc_array_init (&touching_procs, sizeof (int));
  for (zz = 0; zz < count; zz++) {
    fside = sc_array_index (sides, zz);
    f = fside->face;
    tid = fside->treeid;
    tree = p4est_array_index_topidx (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    is_hanging = fside->is_hanging;
    if (!is_hanging) {
      limit = 1;
      is_local = &(fside->is.full.is_local);
      qids = &(fside->is.full.quadid);
      q = &(fside->is.full.quad);
    }
    else {
      limit = P4EST_CHILDREN / 2;
      is_local = fside->is.hanging.is_local;
      qids = fside->is.hanging.quadid;
      q = fside->is.hanging.quad;
    }

    for (i = 0; i < limit; i++) {
      qid = qids[i];
      if (is_local[i]) {
        proc = rank;
        qid += quadrants_offset;
      }
      else {
        /* TODO: replace when the new ghost_layer API is implemented */
        proc = p4est_comm_find_owner (info->p4est, tid, q[i], rank);
        P4EST_ASSERT (proc != rank);
        ip = sc_array_push (&touching_procs);
        *ip = proc;
      }

      if (!intra_tree && i == 0 && tid <= owner_tid) {
        if (tid < owner_tid || p4est_quadrant_compare (q[i], owner_quad) < 0) {
          owner_proc = proc;
          owner_tid = tid;
          owner_quad = q[i];
          owner_face = f;
          owner_side = zz;
        }
      }
#ifdef P4EST_DEBUG
      else {
        P4EST_ASSERT ((zz == 0 && i == 0) || owner_tid < tid ||
                      (owner_tid == tid &&
                       p4est_quadrant_compare (owner_quad, q[i]) < 0));
      }
#endif
    }
  }

  for (zz = 0; zz < count; zz++) {
    fside = sc_array_index (sides, zz);
    f = fside->face;
    tid = fside->treeid;
    tree = p4est_array_index_topidx (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    is_hanging = fside->is_hanging;
    if (!is_hanging) {
      limit = 1;
      is_local = &(fside->is.full.is_local);
      qids = &(fside->is.full.quadid);
    }
    else {
      limit = P4EST_CHILDREN / 2;
      is_local = fside->is.hanging.is_local;
      qids = fside->is.hanging.quadid;
    }

    for (i = 0; i < limit; i++) {
      qid = qids[i];
      if (is_local[i]) {
        qid += quadrants_offset;
#ifndef P4_TO_P8
        start_node = num_inodes + ((zz == owner_side) ? 0 :
                                   info->orientation ? 0 : nodes_per_face -
                                   1);
        stride = (zz == owner_side) ? 1 : info->orientation ? 1 : -1;
        for (j = 0; j < nodes_per_face; j++, start_node += stride) {
          nid = qid * nodes_per_elem + face_nodes[f][j];
          elnode = is_local[i] ? &(local_elem_nodes[nid]) :
            &(ghost_elem_nodes[nid]);
          *elnode = start_node;
        }
#else
        if (zz == owner_side) {
          flipj = false;
          flipk = false;
          swapjk = false;
        }
        else {
          p8est_lnodes_face_node_transform (owner_face, f, info->orientation,
                                            &flipj, &flipk, &swapjk);
        }
        start_node = num_inodes;
        for (l = 0, k = 0; k < nodes_per_edge; k++) {
          for (j = 0; j < nodes_per_edge; j++, l++) {
            nid = qid * nodes_per_elem + face_nodes[f][l];
            elnode = is_local[i] ? &(local_elem_nodes[nid]) :
              &(ghost_elem_nodes[nid]);
            jind = flipj ? (nodes_per_edge - 1 - j) : j;
            kind = flipk ? (nodes_per_edge - 1 - k) : k;
            lind = swapjk ? (nodes_per_edge * jind + kind) :
              (nodes_per_edge * kind + jind);
            *elnode = start_node + lind;
          }
        }
#endif
      }
    }
  }
  sc_array_sort (&touching_procs, sc_int_compare);
  sc_array_uniq (&touching_procs, sc_int_compare);

  count = touching_procs.elem_count;
  if (count) {
    tempq = *owner_quad;
    type = (int8_t) owner_face;
    p4est_lnodes_push_binfo (NULL, &touching_procs, send_buf_info,
                             recv_buf_info, inode_sharers, owner_proc, rank,
                             &tempq, owner_tid, type, num_inodes);
  }
  for (i = 0; i < nodes_per_face; i++) {
    inode = sc_array_push (inodes);
    *inode = -((p4est_locidx_t) owner_proc + 1);
  }

  sc_array_reset (&touching_procs);
}

/* p4est_lnodes_volume_callback:
 *
 * Create independent nodes and set volume nodes to point to them.
 */
static void
p4est_lnodes_volume_callback (p4est_iter_volume_info_t * info, void *Data)
{
  p4est_lnodes_data_t *data = Data;
  p4est_tree_t       *tree = p4est_array_index_topidx (info->p4est->trees,
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
    elem_nodes[nid] = num_inodes + (p4est_locidx_t) i;
    inode = sc_array_push (inodes);
    *inode = -((p4est_locidx_t) rank + 1);
  }
}

/* p4est_lnodes_missing_proc_corner:
 *
 * We do not know which process owns the node at corner c of quadrant q.
 * Find the owner, add the node to the independent nodes array,
 * and add information describing the node to the owner's receive buffer.
 * \return the index in inodes of the new nodes.
 */
static              p4est_locidx_t
p4est_lnodes_missing_proc_corner (p4est_quadrant_t * q, p4est_topidx_t tid,
                                  int c, p4est_t * p4est,
                                  p4est_lnodes_data_t * data)
{
  int                 i;
  p4est_quadrant_t    tempq, tempr;
  p4est_quadrant_t    ownerq;
  int                 ownerc;
  p4est_topidx_t      ownertid = tid;
  p4est_topidx_t      ntid;
  int                 f;
#ifdef P4_TO_P8
  int                 e;
#endif
  sc_array_t          quads, treeids;
  size_t              zz;
  size_t              count;
  p4est_quadrant_t   *ptemp;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *inode;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) inodes->elem_count;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  p4est_lnodes_buf_info_t *binfo;
  int                 owner_proc;

  p4est_quadrant_smallest_corner_descendent (q, &tempq, c);
  ownerq = tempq;

  for (i = 0; i < P4EST_DIM; i++) {
    f = p4est_corner_faces[c][i];
    ntid = p4est_quadrant_face_neighbor_extra (&tempq, tid, f, &tempr,
                                               p4est->connectivity);
    if (ntid == -1) {
      continue;
    }
    if (ntid < ownertid) {
      ownertid = ntid;
      ownerq = tempr;
    }
    else if (ownertid == ntid && p4est_quadrant_compare (&tempr, &ownerq) < 0) {
      ownerq = tempr;
    }
  }
#ifdef P4_TO_P8
  for (i = 0; i < 3; i++) {
    e = p8est_corner_edges[c][i];
    p8est_quadrant_edge_neighbor (&tempq, e, &tempr);
    if (!p4est_quadrant_is_inside_root (&tempr)) {
      sc_array_init (&quads, sizeof (p4est_quadrant_t));
      sc_array_init (&treeids, sizeof (p4est_topidx_t));
      p8est_quadrant_edge_neighbor_extra (&tempq, tid, e,
                                          &quads, &treeids,
                                          p4est->connectivity);
      count = quads.elem_count;
      for (zz = 0; zz < count; zz++) {
        ptemp = sc_array_index (&quads, zz);
        ntid = *((p4est_topidx_t *) sc_array_index (&treeids, zz));
        if (ntid < ownertid) {
          ownertid = ntid;
          ownerq = *ptemp;
        }
        else if (ownertid == ntid &&
                 p4est_quadrant_compare (ptemp, &ownerq) < 0) {
          ownerq = *ptemp;
        }
      }
      sc_array_reset (&quads);
      sc_array_reset (&treeids);
    }
    else if (tid == ownertid && p4est_quadrant_compare (&tempr, &ownerq) < 0) {
      ownerq = tempr;
    }
  }
#endif

  p4est_quadrant_corner_neighbor (&tempq, c, &tempr);
  if (!p4est_quadrant_is_inside_root (&tempr)) {
    sc_array_init (&quads, sizeof (p4est_quadrant_t));
    sc_array_init (&treeids, sizeof (p4est_topidx_t));
    p4est_quadrant_corner_neighbor_extra (&tempq, tid, c, &quads, &treeids,
                                          p4est->connectivity);
    count = quads.elem_count;
    for (zz = 0; zz < count; zz++) {
      ptemp = sc_array_index (&quads, zz);
      ntid = *((p4est_topidx_t *) sc_array_index (&treeids, zz));
      if (ntid < ownertid) {
        ownertid = ntid;
        ownerq = *ptemp;
      }
      else if (ownertid == ntid &&
               p4est_quadrant_compare (ptemp, &ownerq) < 0) {
        ownerq = *ptemp;
      }
    }
    sc_array_reset (&quads);
    sc_array_reset (&treeids);
  }
  else if (tid == ownertid && p4est_quadrant_compare (&tempr, &ownerq) < 0) {
    ownerq = tempr;
  }

  ownerc = p4est_quadrant_child_id (&ownerq);
  owner_proc = p4est_comm_find_owner (p4est, ownertid, &ownerq,
                                      p4est->mpirank);

  P4EST_ASSERT (owner_proc != p4est->mpirank);

  binfo = sc_array_push (&(recv_buf_info[owner_proc]));
  binfo->q = ownerq;
  binfo->q.p.which_tree = ownertid;
  binfo->type = P4EST_LN_C_OFFSET + ownerc;
  binfo->first_index = num_inodes;
  binfo->send_sharers = true;
  binfo->share_offset = -1;
  binfo->share_count = -1;
  inode = sc_array_push (inodes);
  *inode = -((p4est_locidx_t) owner_proc + 1);

  return num_inodes;
}

#ifdef P4_TO_P8
/* p8est_lnodes_missing_proc_edge:
 *
 * We do not know which process owns the nodes at edge e of quadrant q.
 * Find the owner, add the nodes to the independent nodes array, and
 * add info describing the nodes to the owner's receive buffer.
 * Also, set the element nodes of the quadrants of hface that touch edge e
 * to point to the new nodes.
 */
static void
p8est_lnodes_missing_proc_edge (p4est_quadrant_t * q, p4est_topidx_t tid,
                                int e, p4est_t * p4est,
                                p4est_lnodes_data_t * data,
                                p4est_iter_face_side_t * hface)
{
  int                 i, j;
  p4est_quadrant_t    tempq[2], tempr[2];
  p4est_quadrant_t    ownerq[2], final_ownerq;
  int                 ownere;
  int                 ownerc[2];
  p4est_topidx_t      ownertid[2];
  p4est_topidx_t      ntid;
  int                 f;
  sc_array_t          quads, treeids;
  size_t              zz;
  size_t              count;
  p4est_quadrant_t   *ptemp;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *inode;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) inodes->elem_count;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  p4est_lnodes_buf_info_t *binfo;
  int                 owner_proc;
  int                 orientation;
  int                 nodes_per_edge = data->nodes_per_edge;
  int                 nodes_per_elem = data->nodes_per_elem;
  int                 c;
  p4est_locidx_t      qid;
  bool                is_local;
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree;
  p4est_locidx_t      quadrants_offset;
  p4est_locidx_t     *local_elem_nodes = data->local_elem_nodes;
  p4est_locidx_t     *ghost_elem_nodes = data->ghost_elem_nodes;
  p4est_locidx_t      nid;
  int                 stride;
  p4est_locidx_t      start_node;
  int               **edge_nodes = data->edge_nodes;

  p4est_quadrant_sibling (q, &tempq[0], p8est_edge_corners[e][0]);
  p4est_quadrant_sibling (&tempq[0], &tempq[1], p8est_edge_corners[e][1]);
  ownerq[0] = tempq[0];
  ownerq[1] = tempq[1];
  ownertid[0] = tid;
  ownertid[1] = tid;
  ownere = e;

  for (i = 0; i < 2; i++) {
    f = p8est_edge_faces[e][i];
    for (j = 0; j < 2; j++) {
      ntid = p4est_quadrant_face_neighbor_extra (&tempq[j], tid, f, &tempr[j],
                                                 p4est->connectivity);
      if (ntid == -1) {
        continue;
      }
      if (ntid < ownertid[j]) {
        ownertid[j] = ntid;
        ownerq[j] = tempr[j];
      }
      else if (ownertid[j] == ntid &&
               p4est_quadrant_compare (&tempr[j], &ownerq[j]) < 0) {
        ownerq[j] = tempr[j];
      }
    }
  }

  for (j = 0; j < 2; j++) {
    p8est_quadrant_edge_neighbor (&tempq[j], e, &tempr[j]);
    if (!p4est_quadrant_is_inside_root (&tempr[j])) {
      sc_array_init (&quads, sizeof (p4est_quadrant_t));
      sc_array_init (&treeids, sizeof (p4est_topidx_t));
      p8est_quadrant_edge_neighbor_extra (&tempq[j], tid, e, &quads, &treeids,
                                          p4est->connectivity);
      count = quads.elem_count;
      for (zz = 0; zz < count; zz++) {
        ptemp = sc_array_index (&quads, zz);
        ntid = *((p4est_topidx_t *) sc_array_index (&treeids, zz));
        if (ntid < ownertid[j]) {
          ownertid[j] = ntid;
          ownerq[j] = *ptemp;
        }
        else if (ownertid[j] == ntid &&
                 p4est_quadrant_compare (ptemp, &ownerq[j]) < 0) {
          ownerq[j] = *ptemp;
        }
      }
      sc_array_reset (&quads);
      sc_array_reset (&treeids);
    }
    else if (tid == ownertid[j] &&
             p4est_quadrant_compare (&tempr[j], &ownerq[j]) < 0) {
      ownerq[j] = tempr[j];
    }
  }

  P4EST_ASSERT (ownertid[0] == ownertid[1]);
  P4EST_ASSERT (p4est_quadrant_is_sibling (&ownerq[0], &ownerq[1]));
  ownerc[0] = p4est_quadrant_child_id (&ownerq[0]);
  ownerc[1] = p4est_quadrant_child_id (&ownerq[1]);
  ownere = p8est_child_corner_edges[ownerc[0]][ownerc[1]];
  orientation = (ownerc[0] < ownerc[1]) ? 0 : 1;
  final_ownerq = (ownerc[0] < ownerc[1]) ? ownerq[0] : ownerq[1];
  owner_proc = p4est_comm_find_owner (p4est, ownertid[0], &final_ownerq,
                                      p4est->mpirank);

  P4EST_ASSERT (owner_proc != p4est->mpirank);

  binfo = sc_array_push (&(recv_buf_info[owner_proc]));
  binfo->q = final_ownerq;
  binfo->q.p.which_tree = ownertid[0];
  binfo->type = P8EST_LN_E_OFFSET + ownere;
  binfo->first_index = num_inodes;
  binfo->send_sharers = true;
  binfo->share_offset = -1;
  binfo->share_count = -1;

  for (i = 0; i < nodes_per_edge; i++) {
    inode = sc_array_push (inodes);
    *inode = -((p4est_locidx_t) owner_proc + 1);
  }

  f = hface->face;
  tree = p4est_array_index_topidx (trees, tid);
  quadrants_offset = tree->quadrants_offset;
  for (i = 0; i < 2; i++) {
    c = p4est_quadrant_child_id (&tempq[i]);
    j = p8est_corner_face_corners[c][f];
    P4EST_ASSERT (j >= 0);
    qid = hface->is.hanging.quadid[j];
    is_local = hface->is.hanging.is_local[j];
    if (is_local) {
      qid += quadrants_offset;
    }
    start_node = num_inodes + (orientation ? nodes_per_edge - 1 : 0);
    stride = orientation ? -1 : 1;
    for (j = 0; j < nodes_per_edge; j++, start_node += stride) {
      nid = qid * nodes_per_elem + edge_nodes[e][j];
      if (is_local) {
        local_elem_nodes[nid] = start_node;
      }
      else {
        ghost_elem_nodes[nid] = start_node;
      }
    }
  }
}
#endif

/* p4est_lnodes_hface_fix:
 *
 * Coming out of the iterate loop, the element nodes that belong to the hanging
 * corner and hanging edges of the smaller side of a hanging face are not set:
 * get the correct values from the neighbors that do have those values and pass
 * them along.
 */
static void
p4est_lnodes_hface_fix (p4est_t * p4est, p4est_iter_face_side_t * hface,
                        p4est_lnodes_data_t * data)
{
  int                 i;
  p4est_topidx_t      tid = hface->treeid;
  p4est_locidx_t      qid[P4EST_CHILDREN / 2];
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree = p4est_array_index_topidx (trees, tid);
  p4est_locidx_t      quadrants_offset = tree->quadrants_offset;
  bool               *is_local = hface->is.hanging.is_local;
  int                *corner_nodes = data->corner_nodes;
  int                 nodes_per_elem = data->nodes_per_elem;
  p4est_locidx_t      nid, nidval, nid2;
  int                 f = hface->face;
  p4est_locidx_t     *local_elem_nodes = data->local_elem_nodes;
  p4est_locidx_t     *ghost_elem_nodes = data->ghost_elem_nodes;
  int                 c;
  int                 i2;
#ifdef P4_TO_P8
  int                 c2;
  int                 nodes_per_edge = data->nodes_per_edge;
  int               **edge_nodes = data->edge_nodes;
  int                 j, k;
  int                 i1;
  int                 e, e2;
#endif

  P4EST_ASSERT (hface->is_hanging);
  for (i = 0; i < P4EST_CHILDREN / 2; i++) {
    qid[i] = hface->is.hanging.quadid[i];
    if (is_local[i]) {
      qid[i] += quadrants_offset;
    }
  }

  for (i = 0; i < P4EST_CHILDREN / 2; i++) {
    i2 = P4EST_CHILDREN / 2 - 1 - i;
    if (!is_local[i2]) {
      continue;
    }
#ifndef P4_TO_P8
    c = p4est_face_corners[p4est_zface_to_rface[f]][i];
#else
    c = p4est_face_corners[f][i];
#endif
    nid = qid[i] * nodes_per_elem + corner_nodes[c];
    nid2 = qid[i2] * nodes_per_elem + corner_nodes[c];
    if (is_local[i]) {
      nidval = local_elem_nodes[nid];
      P4EST_ASSERT (nidval >= 0);
    }
    else {
      nidval = ghost_elem_nodes[nid];
      if (nidval == -1) {
        nidval = p4est_lnodes_missing_proc_corner (hface->is.hanging.quad[i],
                                                   tid, c, p4est, data);
        ghost_elem_nodes[nid] = nidval;
      }
    }
    P4EST_ASSERT (local_elem_nodes[nid2] == -1);
    local_elem_nodes[nid2] = nidval;
  }

#ifdef P4_TO_P8
  for (i = 0; i < 4; i++) {
    e = p8est_face_edges[f][i];
    e2 = p8est_face_edges[f][i ^ 1];
    for (j = 0; j < 2; j++) {
      c = p8est_edge_corners[e][j];
      c2 = p8est_edge_corners[e2][j];
      i1 = p8est_corner_face_corners[c][f];
      i2 = p8est_corner_face_corners[c2][f];
      if (!is_local[i2]) {
        continue;
      }
      for (k = 0; k < nodes_per_edge; k++) {
        nid = qid[i1] * nodes_per_elem + edge_nodes[e][k];
        nid2 = qid[i2] * nodes_per_elem + edge_nodes[e][k];
        if (is_local[i1]) {
          nidval = local_elem_nodes[nid];
          P4EST_ASSERT (nidval >= 0);
        }
        else {
          nidval = ghost_elem_nodes[nid];
          if (nidval == -1) {
            P4EST_ASSERT (k == 0);
            p8est_lnodes_missing_proc_edge (hface->is.hanging.quad[i1], tid,
                                            e, p4est, data, hface);
            nidval = ghost_elem_nodes[nid];
            P4EST_ASSERT (nidval >= 0);
          }
        }
        P4EST_ASSERT (local_elem_nodes[nid2] == -1);
        local_elem_nodes[nid2] = nidval;
      }
    }
  }
#endif
}

#ifdef P4_TO_P8
/* p8est_lnodes_hedge_fix:
 *
 * Coming out of the iterate loop, the element nodes that belong to the hanging
 * corner of the smaller side of a hanging edge are not set: get the correct
 * values from the neighbors that do have those values and pass them along.
 */
static void
p8est_lnodes_hedge_fix (p4est_t * p4est, p8est_iter_edge_side_t * hedge,
                        p4est_lnodes_data_t * data)
{
  int                 i;
  p4est_topidx_t      tid = hedge->treeid;
  p4est_locidx_t      qid[2];
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree = p4est_array_index_topidx (trees, tid);
  p4est_locidx_t      quadrants_offset = tree->quadrants_offset;
  bool               *is_local = hedge->is.hanging.is_local;
  int                *corner_nodes = data->corner_nodes;
  int                 nodes_per_elem = data->nodes_per_elem;
  p4est_locidx_t      nid, nidval, nid2;
  int                 e = hedge->edge;
  p4est_locidx_t     *local_elem_nodes = data->local_elem_nodes;
  p4est_locidx_t     *ghost_elem_nodes = data->ghost_elem_nodes;
  int                 c, c2;

  P4EST_ASSERT (hedge->is_hanging);
  for (i = 0; i < 2; i++) {
    qid[i] = hedge->is.hanging.quadid[i];
    if (is_local[i]) {
      qid[i] += quadrants_offset;
    }
  }

  for (i = 0; i < 2; i++) {
    c = p8est_edge_corners[e][i];
    c2 = p4est_face_corners[e][1 - i];
    if (!is_local[1 - i]) {
      continue;
    }
    nid = qid[i] * nodes_per_elem + corner_nodes[c];
    nid2 = qid[1 - i] * nodes_per_elem + corner_nodes[c];
    if (is_local[i]) {
      nidval = local_elem_nodes[nid];
      P4EST_ASSERT (nidval >= 0);
    }
    else {
      nidval = ghost_elem_nodes[nid];
      if (nidval == -1) {
        nidval = p4est_lnodes_missing_proc_corner (hedge->is.hanging.quad[i],
                                                   tid, c, p4est, data);
        ghost_elem_nodes[nid] = nidval;
      }
    }
    P4EST_ASSERT (local_elem_nodes[nid2] == -1);
    local_elem_nodes[nid2] = nidval;
  }
}
#endif

/* buf_info_t objects equivalence is based on the quadrants that describe them
 * and the face/edge/corner of that quadrant.  The first index of the nodes they
 * describe is not considered.
 */
static int
p4est_lnodes_binfo_compare (const void *a, const void *b)
{
  const p4est_lnodes_buf_info_t *bufa = a;
  const p4est_lnodes_buf_info_t *bufb = b;
  int                 piggy_compar =
    p4est_quadrant_compare_piggy (&(bufa->q), &(bufb->q));
  if (piggy_compar) {
    return piggy_compar;
  }
  return (bufa->type - bufb->type);
}

static              bool
p4est_lnodes_binfo_is_equal (const void *a, const void *b)
{
  const p4est_lnodes_buf_info_t *bufa = a;
  const p4est_lnodes_buf_info_t *bufb = b;
  return (p4est_quadrant_is_equal_piggy (&(bufa->q), &(bufb->q)) &&
          bufa->type == bufb->type);
}

#ifdef P4_TO_P8
/* p8est_lnodes_code_fix:
 *
 * Coming out of the iterate loop, a hanging edge that touches a hanging face
 * has its bit in the face code set to 1, whereas _decode expects it to be 0.
 */
static              int16_t
p8est_lnodes_code_fix (int16_t face_code)
{
  int                 i;

  if (!face_code) {
    return face_code;
  }
  else {
    for (i = 0; i < 3; i++) {
      if (face_code & (0x0001 << (i + 6))) {
        /* zero out touching edges */
        face_code &= ~(0x0001 << (((i + 1) % 2) + 3));
        face_code &= ~(0x0001 << (((i + 2) % 2) + 3));
      }
    }
  }

  return face_code;
}
#endif

static void
p4est_lnodes_init_data (p4est_lnodes_data_t * data, int p, p4est_t * p4est,
                        sc_array_t * ghost_layer, p4est_lnodes_t * lnodes)
{
  int                 i, j, n;
  int                 npel;
  int                 npv;
  int                 vcount;
  int                 npf;
  int                 fcount[P4EST_DIM * 2];
  int                 f;
  int                 bcount;
  int                 c;
#ifdef P4_TO_P8
  int                 e;
  int                 eshift;
  int                 k;
  int                 npe;
  int                 ecount[12];
#endif
  p4est_locidx_t      nlq = p4est->local_num_quadrants;
  p4est_locidx_t      ngq = (p4est_locidx_t) ghost_layer->elem_count;
  p4est_locidx_t      nlcdp = nlq * P4EST_CHILDREN;
  p4est_locidx_t      ngcdp = ngq * P4EST_CHILDREN;
#ifdef P4_TO_P8
  p4est_locidx_t      nledp = nlq * 12;
  p4est_locidx_t      ngedp = ngq * 12;
#endif
  p4est_locidx_t      nlen;
  p4est_locidx_t      ngen;
  int                 mpisize = p4est->mpisize;

#ifndef P4_TO_P8
  npel = data->nodes_per_elem = (p + 1) * (p + 1);
  npv = data->nodes_per_volume = (p - 1) * (p - 1);
  npf = data->nodes_per_face = p - 1;
  fcount[0] = fcount[1] = fcount[2] = fcount[3] = 0;
#else
  npel = data->nodes_per_elem = (p + 1) * (p + 1) * (p + 1);
  npv = data->nodes_per_volume = (p - 1) * (p - 1) * (p - 1);
  npf = data->nodes_per_face = (p - 1) * (p - 1);
  fcount[0] = fcount[1] = fcount[2] = fcount[3] = fcount[4] = fcount[5] = 0;
  npe = data->nodes_per_edge = p - 1;
  ecount[0] = ecount[1] = ecount[2] = ecount[3] = ecount[4] = ecount[5] = 0;
  ecount[6] = ecount[7] = ecount[8] = ecount[9] = ecount[10] = ecount[11] = 0;
#endif
  vcount = 0;
  nlen = nlq * npel;
  ngen = ngq * npel;

  data->volume_nodes = P4EST_ALLOC (int, npv);
  for (i = 0; i < P4EST_DIM * 2; i++) {
    data->face_nodes[i] = P4EST_ALLOC (int, npf);
  }
#ifdef P4_TO_P8
  for (i = 0; i < 12; i++) {
    data->edge_nodes[i] = P4EST_ALLOC (int, npf);
  }
#endif

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
          data->corner_nodes[c] = n;
          break;
        }
      }
    }
#ifdef P4_TO_P8
  }
#endif

  data->local_cdp = P4EST_ALLOC (p4est_lnodes_cdp_t, nlcdp);
  memset (data->local_cdp, -1, nlcdp * sizeof (p4est_lnodes_cdp_t));
  data->ghost_cdp = P4EST_ALLOC (p4est_lnodes_cdp_t, ngcdp);
  memset (data->ghost_cdp, -1, ngcdp * sizeof (p4est_lnodes_cdp_t));
#ifdef P4_TO_P8
  data->local_edp = P4EST_ALLOC (p8est_lnodes_edp_t, nledp);
  memset (data->local_edp, -1, nledp * sizeof (p8est_lnodes_edp_t));
  data->ghost_edp = P4EST_ALLOC (p8est_lnodes_edp_t, ngedp);
  memset (data->ghost_edp, -1, ngedp * sizeof (p8est_lnodes_edp_t));
#endif

  data->local_elem_nodes = lnodes->local_nodes;
  data->ghost_elem_nodes = P4EST_ALLOC (p4est_locidx_t, ngen);
  memset (data->ghost_elem_nodes, -1, ngen * sizeof (p4est_locidx_t));

  data->hfaces = sc_array_new (sizeof (p4est_iter_face_side_t));
#ifdef P4_TO_P8
  data->hedges = sc_array_new (sizeof (p8est_iter_edge_side_t));
#endif
  data->inodes = sc_array_new (sizeof (p4est_locidx_t));
  data->inode_sharers = sc_array_new (sizeof (int));
  data->send_buf_info = P4EST_ALLOC (sc_array_t, mpisize);
  data->recv_buf_info = P4EST_ALLOC (sc_array_t, mpisize);
  data->sorters = P4EST_ALLOC (sc_array_t, mpisize);
  for (i = 0; i < mpisize; i++) {
    sc_array_init (&(data->send_buf_info[i]),
                   sizeof (p4est_lnodes_buf_info_t));
    sc_array_init (&(data->recv_buf_info[i]),
                   sizeof (p4est_lnodes_buf_info_t));
    sc_array_init (&(data->sorters[i]), sizeof (p4est_lnodes_sorter_t));
  }

  data->face_codes = lnodes->face_code;
  data->global_offsets =
    P4EST_ALLOC_ZERO (p4est_gloidx_t, p4est->mpisize + 1);
}

static void
p4est_lnodes_reset_data (p4est_lnodes_data_t * data, p4est_t * p4est)
{
  int                 mpisize = p4est->mpisize;
  int                 i;

  P4EST_FREE (data->volume_nodes);
  for (i = 0; i < P4EST_DIM * 2; i++) {
    P4EST_FREE (data->face_nodes[i]);
  }
#ifdef P4_TO_P8
  for (i = 0; i < 12; i++) {
    P4EST_FREE (data->edge_nodes[i]);
  }
#endif

  sc_array_destroy (data->hfaces);
#ifdef P4_TO_P8
  sc_array_destroy (data->hedges);
#endif
  sc_array_destroy (data->inodes);
  sc_array_destroy (data->inode_sharers);
  for (i = 0; i < mpisize; i++) {
    sc_array_reset (&(data->send_buf_info[i]));
    sc_array_reset (&(data->recv_buf_info[i]));
    sc_array_reset (&(data->sorters[i]));
  }
  P4EST_FREE (data->send_buf_info);
  P4EST_FREE (data->recv_buf_info);
  P4EST_FREE (data->sorters);
  /* do not free face_codes: controlled by lnodes_t */
  P4EST_FREE (data->global_offsets);
}

static void
p4est_lnodes_count_nodes (p4est_lnodes_data_t * data, p4est_t * p4est,
                          sc_array_t * ghost_layer, p4est_lnodes_t * lnodes)
{
  p4est_locidx_t      num_inodes;
  p4est_locidx_t      nlq = p4est->local_num_quadrants;
  p4est_locidx_t      npel;
  p4est_locidx_t      nlen;
  p4est_locidx_t      li;
  p4est_locidx_t      inidx;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *inode;
  p4est_topidx_t     *local_en = data->local_elem_nodes;
  int                 i;
  int                 rank = p4est->mpirank;
  int                 mpisize = p4est->mpisize;
  p4est_locidx_t      count = 0;
  p4est_locidx_t     *global_num_indep =
    P4EST_ALLOC (p4est_locidx_t, mpisize);
  sc_array_t         *sorter_array = &(data->sorters[rank]);
  p4est_lnodes_sorter_t *sorter;
  p4est_locidx_t      rankval = -((p4est_locidx_t) rank + 1);

  npel = data->nodes_per_elem;
  nlen = npel * nlq;
  num_inodes = (p4est_locidx_t) inodes->elem_count;
  for (li = 0; li < nlen; li++) {
    inidx = local_en[li];
    if (inidx >= 0) {
      inode = sc_array_index (inodes, (size_t) inidx);
      if (*inode == rankval) {
        *inode = count;
        sorter = sc_array_push (sorter_array);
        sorter->local_index = count++;
        sorter->inode_index = inidx;
      }
    }
  }

  lnodes->owned_count = count;
  MPI_Allgather (&count, 1, P4EST_MPI_LOCIDX, global_num_indep, 1,
                 P4EST_MPI_LOCIDX, p4est->mpicomm);
  for (i = 0; i < mpisize; i++) {
    data->global_offsets[i + 1] = data->global_offsets[i] +
      (p4est_gloidx_t) global_num_indep[i];
  }

  P4EST_FREE (global_num_indep);
}

static              bool
p4est_lnodes_test_comm (p4est_t * p4est, p4est_lnodes_data_t * data)
{
  int                 mpisize = p4est->mpisize;
  int                 i, j;
  sc_array_t         *send;
  sc_array_t         *send_buf_info = data->send_buf_info;
  sc_array_t         *recv, *recv2;
  sc_array_t         *recv_buf;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  size_t              count, count2, zz;
  int                 mpiret;
  sc_array_t          send_requests;
  MPI_Request        *send_request;
  MPI_Status          probe_status, recv_status;
  int                 num_recv_queries = 0;
  int                *num_recv_expect = P4EST_ALLOC_ZERO (int, mpisize);
  int                 byte_count;
  size_t              elem_count;
  p4est_lnodes_buf_info_t *binfo, *binfo2, *prev;

  sc_array_init (&send_requests, sizeof (MPI_Request));
  for (i = 0; i < mpisize; i++) {
    send = &(send_buf_info[i]);
    count = send->elem_count;
    if (count > 0) {
      P4EST_ASSERT (i != p4est->mpirank);
      send_request = sc_array_push (&send_requests);
      mpiret = MPI_Isend (send->array,
                          (int) (count * sizeof (p4est_lnodes_buf_info_t)),
                          MPI_BYTE, i, P4EST_COMM_LNODES_TEST, p4est->mpicomm,
                          send_request);
      SC_CHECK_MPI (mpiret);
    }
    recv = &(recv_buf_info[i]);
    count = recv->elem_count;
    if (count) {
      P4EST_ASSERT (i != p4est->mpirank);
      num_recv_queries++;
      num_recv_expect[i]++;
    }
  }

  recv_buf = P4EST_ALLOC (sc_array_t, mpisize);
  for (i = 0; i < mpisize; i++) {
    sc_array_init (&(recv_buf[i]), sizeof (p4est_lnodes_buf_info_t));
  }
  for (i = 0; i < num_recv_queries; i++) {
    mpiret =
      MPI_Probe (MPI_ANY_SOURCE, P4EST_COMM_LNODES_TEST, p4est->mpicomm,
                 &probe_status);
    SC_CHECK_MPI (mpiret);
    j = probe_status.MPI_SOURCE;
    P4EST_ASSERT (j != p4est->mpirank && num_recv_expect[j] == 1);
    recv = &(recv_buf[j]);
    mpiret = MPI_Get_count (&probe_status, MPI_BYTE, &byte_count);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (byte_count % ((int) sizeof (p4est_lnodes_buf_info_t)) == 0);
    elem_count = ((size_t) byte_count) / sizeof (p4est_lnodes_buf_info_t);
    sc_array_resize (recv, elem_count);
    mpiret = MPI_Recv (recv->array, byte_count, MPI_BYTE, j,
                       P4EST_COMM_LNODES_TEST, p4est->mpicomm, &recv_status);
    SC_CHECK_MPI (mpiret);
    num_recv_expect[j]--;

    recv2 = &(recv_buf_info[j]);
    count2 = recv2->elem_count;
    count = 0;
    prev = NULL;
    for (zz = 0; zz < count2; zz++) {
      binfo2 = sc_array_index (recv2, zz);
      if (zz > 0 && p4est_lnodes_binfo_is_equal (prev, binfo2)) {
        continue;
      }
      binfo = sc_array_index (recv, count++);
      P4EST_ASSERT (p4est_quadrant_is_equal_piggy
                    (&(binfo->q), &(binfo2->q)));
      P4EST_ASSERT (binfo->type == binfo2->type);
      P4EST_ASSERT (binfo->send_sharers == binfo2->send_sharers);
      if (!binfo->send_sharers) {
        P4EST_ASSERT (binfo->share_count == binfo2->share_count);
      }
      prev = binfo2;
    }
    P4EST_ASSERT (count == elem_count);
  }

  if (send_requests.elem_count > 0) {
    mpiret = MPI_Waitall ((int) send_requests.elem_count,
                          (MPI_Request *) send_requests.array,
                          MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  sc_array_reset (&send_requests);
  for (i = 0; i < mpisize; i++) {
    sc_array_reset (&(recv_buf[i]));
  }

  P4EST_FREE (recv_buf);
  P4EST_FREE (num_recv_expect);
  return true;
}

static void
p4est_lnodes_pass (p4est_t * p4est, p4est_lnodes_data_t * data)
{
  int                 mpisize = p4est->mpisize;
  int                 i, j, k;
  int                 limit;
  sc_array_t         *send, *send_info;
  sc_array_t         *send_buf_info = data->send_buf_info;
  sc_array_t         *send_buf;
  sc_array_t         *recv, *recv_info;
  sc_array_t         *recv_buf;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  size_t              count, info_count, zz;
  int                 mpiret;
  sc_array_t          send_requests;
  MPI_Request        *send_request;
  MPI_Status          probe_status, recv_status;
  int                 num_recv_queries = 0;
  int                *num_recv_expect = P4EST_ALLOC_ZERO (int, mpisize);
  int                 byte_count;
  size_t              elem_count;
  p4est_lnodes_buf_info_t *binfo, *prev;
  size_t              index, prev_index;
  int                 nodes_per_face = data->nodes_per_face;
#ifdef P4_TO_P8
  int                 nodes_per_edge = data->nodes_per_edge;
#endif
  int8_t              type;
  p4est_locidx_t     *lp;
  int                *ip;
  p4est_locidx_t     *inode, *inode_prev;
  sc_array_t         *inode_sharers = data->inode_sharers;
  sc_array_t         *inodes = data->inodes;
  int                 share_proc;
  int                 share_count;
  size_t              send_count;
  sc_array_t         *sorter_array;
  p4est_lnodes_sorter_t *sorter;

  sc_array_init (&send_requests, sizeof (MPI_Request));
  send_buf = P4EST_ALLOC (sc_array_t, mpisize);
  for (i = 0; i < mpisize; i++) {
    sc_array_init (&(send_buf[i]), sizeof (p4est_locidx_t));
  }
  for (i = 0; i < mpisize; i++) {
    send_info = &(send_buf_info[i]);
    count = send_info->elem_count;
    if (count > 0) {
      P4EST_ASSERT (i != p4est->mpirank);
      send = &(send_buf[i]);
      for (zz = 0; zz < count; zz++) {
        binfo = sc_array_index (send_info, zz);
        index = (size_t) binfo->first_index;
        type = binfo->type;
        if (type >= P4EST_LN_C_OFFSET) {
          limit = 1;
        }
#ifdef P4_TO_P8
        else if (type >= P8EST_LN_E_OFFSET) {
          limit = nodes_per_edge;
        }
#endif
        else {
          P4EST_ASSERT (0 <= type && type < P4EST_DIM * 2);
          limit = nodes_per_face;
        }
        for (j = 0; j < limit; j++) {
          lp = sc_array_push (send);
          inode = sc_array_index (inodes, index++);
          P4EST_ASSERT (*inode >= 0);
          *lp = *inode;
        }
        if (binfo->send_sharers) {
          lp = sc_array_push (send);
          *lp = (p4est_locidx_t) binfo->share_count;
          P4EST_ASSERT (binfo->share_count > 0);
          index = (size_t) binfo->share_offset;
          share_count = (int) binfo->share_count;
          for (j = 0; j < share_count; j++) {
            lp = sc_array_push (send);
            share_proc = *((int *) sc_array_index (inode_sharers, index++));
            *lp = (p4est_locidx_t) share_proc;
            P4EST_ASSERT (0 <= share_proc && share_proc < mpisize);
          }
        }
      }
      send_count = send->elem_count;
      send_request = sc_array_push (&send_requests);
      mpiret = MPI_Isend (send->array,
                          (int) (send_count * sizeof (p4est_locidx_t)),
                          MPI_BYTE, i, P4EST_COMM_LNODES_PASS, p4est->mpicomm,
                          send_request);
      SC_CHECK_MPI (mpiret);
    }
    recv_info = &(recv_buf_info[i]);
    count = recv_info->elem_count;
    if (count) {
      P4EST_ASSERT (i != p4est->mpirank);
      num_recv_queries++;
      num_recv_expect[i]++;
    }
  }

  recv_buf = P4EST_ALLOC (sc_array_t, mpisize);
  for (i = 0; i < mpisize; i++) {
    sc_array_init (&(recv_buf[i]), sizeof (p4est_locidx_t));
  }
  for (i = 0; i < num_recv_queries; i++) {
    mpiret =
      MPI_Probe (MPI_ANY_SOURCE, P4EST_COMM_LNODES_PASS, p4est->mpicomm,
                 &probe_status);
    SC_CHECK_MPI (mpiret);
    j = probe_status.MPI_SOURCE;
    P4EST_ASSERT (j != p4est->mpirank && num_recv_expect[j] == 1);
    recv = &(recv_buf[j]);
    recv_info = &(recv_buf_info[j]);
    mpiret = MPI_Get_count (&probe_status, MPI_BYTE, &byte_count);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (byte_count % ((int) sizeof (p4est_locidx_t)) == 0);
    elem_count = ((size_t) byte_count) / sizeof (p4est_locidx_t);
    sc_array_resize (recv, elem_count);
    mpiret = MPI_Recv (recv->array, byte_count, MPI_BYTE, j,
                       P4EST_COMM_LNODES_PASS, p4est->mpicomm, &recv_status);
    SC_CHECK_MPI (mpiret);
    num_recv_expect[j]--;

    sorter_array = &(data->sorters[j]);
    info_count = recv_info->elem_count;
    count = 0;
    prev = NULL;
    for (zz = 0; zz < info_count; zz++) {
      binfo = sc_array_index (recv_info, zz);
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
      index = (size_t) binfo->first_index;
      if (zz > 0 && p4est_lnodes_binfo_is_equal (prev, binfo)) {
        binfo->share_offset = prev->share_offset;
        binfo->share_count = prev->share_count;
        prev_index = (size_t) prev->first_index;
        for (k = 0; k < limit; k++) {
          inode = sc_array_index (inodes, index);
          inode_prev = sc_array_index (inodes, prev_index++);
          P4EST_ASSERT (*inode == -((p4est_locidx_t) j + 1));
          P4EST_ASSERT (*inode_prev >= 0);
          *inode = *inode_prev;
          sorter = sc_array_push (sorter_array);
          sorter->local_index = *inode_prev;
          sorter->inode_index = index++;
        }
      }
      else {
        for (k = 0; k < limit; k++) {
          inode = sc_array_index (inodes, index);
          lp = sc_array_index (recv, count++);
          P4EST_ASSERT (*inode == -((p4est_locidx_t) j + 1));
          P4EST_ASSERT (*lp >= 0);
          *inode = *lp;
          sorter = sc_array_push (sorter_array);
          sorter->local_index = *lp;
          sorter->inode_index = index++;
        }
        if (binfo->send_sharers) {
          lp = sc_array_index (recv, count++);
          share_count = (int) (*lp);
          P4EST_ASSERT (share_count > 0);
          binfo->share_count = (int8_t) share_count;
          binfo->share_offset = (p4est_locidx_t) inode_sharers->elem_count;
          for (k = 0; k < share_count; k++) {
            ip = sc_array_push (inode_sharers);
            lp = sc_array_index (recv, count++);
            *ip = (int) (*lp);
            P4EST_ASSERT (0 <= *ip && *ip < mpisize);
          }
        }
        prev = binfo;
      }
    }
    P4EST_ASSERT (count == elem_count);
    sc_array_sort (sorter_array, p4est_lnodes_sorter_compare);
  }

  if (send_requests.elem_count > 0) {
    mpiret = MPI_Waitall ((int) send_requests.elem_count,
                          (MPI_Request *) send_requests.array,
                          MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  sc_array_reset (&send_requests);
  for (i = 0; i < mpisize; i++) {
    sc_array_reset (&(send_buf[i]));
    sc_array_reset (&(recv_buf[i]));
  }
  P4EST_FREE (send_buf);
  P4EST_FREE (recv_buf);
  P4EST_FREE (num_recv_expect);
}

static void
p4est_lnodes_global_and_sharers (p4est_lnodes_data_t * data,
                                 p4est_lnodes_t * lnodes, p4est_t * p4est)
{
  int                 i, j, k, l;
  int                 mpisize = p4est->mpisize;
  sc_array_t          inode_to_global;
  sc_array_t         *s_array;
  p4est_lnodes_sorter_t *sorter;
  p4est_lnodes_sorter_t *prev = NULL;
  p4est_gloidx_t     *global_offsets = data->global_offsets;
  p4est_gloidx_t     *gnodes;
  p4est_gloidx_t      offset;
  size_t              count, zz;
  p4est_locidx_t      icount = 0, gcount = 0, gcount2 = 0;
  p4est_locidx_t     *lp, li;
  p4est_locidx_t     *elnodes = lnodes->local_nodes;
  p4est_locidx_t      nlen = lnodes->num_local_elements * lnodes->vnodes;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) data->inodes->elem_count;
  p4est_locidx_t      inidx;
  int                *comm_proc;
  int                 comm_proc_count;
  sc_array_t         *inode_sharers = data->inode_sharers;
  sc_array_t         *sharers;
  p4est_lnodes_rank_t *lrank;
  sc_array_t         *binfo_array;
  p4est_lnodes_buf_info_t *binfo;
  p4est_locidx_t      share_offset;
  int                 share_count;
  p4est_locidx_t      index;
  int                 limit;
  int                 nodes_per_face = data->nodes_per_face;
#ifdef P4_TO_P8
  int                 nodes_per_edge = data->nodes_per_edge;
#endif
  int                 proc;
  int                 shareidx;
  p4est_locidx_t      gidx;
  sc_array_t         *shared_nodes;

  sc_array_init (&inode_to_global, sizeof (p4est_locidx_t));
  sc_array_resize (&inode_to_global, (size_t) num_inodes);

  for (i = 0; i < mpisize; i++) {
    s_array = &(data->sorters[i]);
    count = s_array->elem_count;
    if (s_array->elem_count == 0) {
      continue;
    }
    P4EST_ASSERT (sc_array_is_sorted (s_array, p4est_lnodes_sorter_compare));
    if (i == p4est->mpirank) {
      lnodes->owned_offset = gcount;
      P4EST_ASSERT (count == (size_t) lnodes->owned_count);
    }
    for (zz = 0; zz < count; zz++) {
      sorter = sc_array_index (s_array, zz);
      lp = sc_array_index (&inode_to_global, sorter->inode_index);
      *lp = gcount;
      if (zz == 0 || sorter->local_index > prev->local_index) {
        gcount++;
      }
      prev = sorter;
    }
    icount += count;
  }
  P4EST_ASSERT (icount == num_inodes);

  lnodes->num_indep_nodes = gcount;

  lnodes->global_nodes = gnodes = P4EST_ALLOC (p4est_gloidx_t, gcount);

  for (i = 0; i < mpisize; i++) {
    s_array = &(data->sorters[i]);
    count = s_array->elem_count;
    if (s_array->elem_count == 0) {
      continue;
    }
    offset = global_offsets[i];
    for (zz = 0; zz < count; zz++) {
      sorter = sc_array_index (s_array, zz);
      if (zz == 0 || sorter->local_index > prev->local_index) {
        gnodes[gcount2++] = offset + (p4est_gloidx_t) sorter->local_index;
        P4EST_ASSERT (gnodes[gcount2 - 1] < global_offsets[i + 1]);
      }
      prev = sorter;
    }
  }
  P4EST_ASSERT (gcount2 == gcount);

  for (li = 0; li < nlen; li++) {
    inidx = elnodes[li];
    P4EST_ASSERT (0 <= inidx && inidx < num_inodes);
    lp = sc_array_index (&inode_to_global, (size_t) inidx);
    elnodes[li] = *lp;
  }

  comm_proc = P4EST_ALLOC_ZERO (int, mpisize);
  count = inode_sharers->elem_count;
  for (zz = 0; zz < count; zz++) {
    i = *((int *) sc_array_index (inode_sharers, zz));
    comm_proc[i]++;
  }
  comm_proc_count = 0;
  lnodes->sharers = sharers = sc_array_new (sizeof (p4est_lnodes_rank_t));
  for (i = 0; i < mpisize; i++) {
    if (comm_proc[i]) {
      lrank = sc_array_push (sharers);
      lrank->rank = i;
      sc_array_init (&(lrank->shared_nodes), sizeof (p4est_locidx_t));
      comm_proc[i] = comm_proc_count++;
    }
    else {
      comm_proc[i] = -1;
    }
  }

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
        binfo = sc_array_index (binfo_array, zz);
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
        index = binfo->first_index;
        share_offset = binfo->share_offset;
        share_count = (int) binfo->share_count;
        for (k = 0; k < limit; k++) {
          gidx = *((p4est_locidx_t *) sc_array_index (&inode_to_global,
                                                      (size_t) index++));
          for (l = 0; l < share_count; l++) {
            proc = *((int *) sc_array_index (inode_sharers,
                                             (size_t) share_offset +
                                             (size_t) l));
            shareidx = comm_proc[proc];
            P4EST_ASSERT (shareidx >= 0);
            lrank = sc_array_index_int (sharers, shareidx);
            P4EST_ASSERT (lrank->rank == proc);
            lp = sc_array_push (&(lrank->shared_nodes));
            *lp = gidx;
          }
        }
      }
    }
  }

  for (i = 0; i < comm_proc_count; i++) {
    lrank = sc_array_index_int (sharers, i);
    sc_array_sort (&(lrank->shared_nodes), p4est_locidx_compare);
    proc = lrank->rank;
    lrank->shared_mine_offset = -1;
    lrank->shared_mine_count = 0;
    lrank->owned_offset = -1;
    lrank->owned_count = 0;
    shared_nodes = &(lrank->shared_nodes);
    count = shared_nodes->elem_count;
    for (zz = 0; zz < count; zz++) {
      gidx = *((p4est_locidx_t *) sc_array_index (shared_nodes, zz));
      if (lnodes->owned_offset <= gidx &&
          gidx < (lnodes->owned_offset + lnodes->owned_count)) {
        if (lrank->shared_mine_count == 0) {
          lrank->shared_mine_offset = (p4est_locidx_t) count;
        }
        lrank->shared_mine_count++;
      }
      if (global_offsets[proc] <= gnodes[gidx] &&
          gnodes[gidx] < global_offsets[proc + 1]) {
        if (lrank->owned_count == 0) {
          lrank->owned_offset = gidx;
        }
        lrank->owned_count++;
      }
    }
    if (proc == p4est->mpirank) {
      lrank->owned_count = lnodes->owned_count;
      lrank->owned_offset = lnodes->owned_offset;
    }
  }
  sc_array_reset (&inode_to_global);
  P4EST_FREE (comm_proc);
}

p4est_lnodes_t     *
p4est_lnodes_new (p4est_t * p4est, sc_array_t * ghost_layer, int degree)
{
  p4est_iter_face_t   fiter;
  p4est_iter_volume_t viter;
  p4est_iter_corner_t citer;
#ifdef P4_TO_P8
  p8est_iter_edge_t   eiter;
#endif
  p4est_lnodes_data_t data;
  p4est_iter_face_side_t *hface;
#ifdef P4_TO_P8
  p8est_iter_edge_side_t *hedge;
#endif
  p4est_locidx_t      nel;
  p4est_locidx_t      nlen;
  p4est_locidx_t      li;
#ifdef P4EST_DEBUG
  size_t              zz;
  p4est_lnodes_buf_info_t *last, *new;
#endif
  int                 i;
  int                 mpisize = p4est->mpisize;
  p4est_lnodes_t     *lnodes = P4EST_ALLOC (p4est_lnodes_t, 1);

  P4EST_ASSERT (degree >= 1);

  lnodes->degree = degree;
  lnodes->num_local_elements = nel = p4est->local_num_quadrants;
#ifndef P4_TO_P8
  lnodes->vnodes = (degree + 1) * (degree + 1);
  lnodes->face_code = P4EST_ALLOC_ZERO (int8_t, nel);
#else
  lnodes->vnodes = (degree + 1) * (degree + 1) * (degree + 1);
  lnodes->face_code = P4EST_ALLOC_ZERO (int16_t, nel);
#endif
  nlen = nel * lnodes->vnodes;
  lnodes->local_nodes = P4EST_ALLOC (p4est_locidx_t, nlen);
  memset (lnodes->local_nodes, -1, nlen * sizeof (p4est_locidx_t));

  p4est_lnodes_init_data (&data, degree, p4est, ghost_layer, lnodes);
  if (degree == 1) {
    viter = NULL;
    fiter = p4est_lnodes_face_simple_callback;
#ifdef P4_TO_P8
    eiter = p8est_lnodes_edge_simple_callback;
#endif
  }
  else {
    viter = p4est_lnodes_volume_callback;
    fiter = p4est_lnodes_face_callback;
#ifdef P4_TO_P8
    eiter = p8est_lnodes_edge_callback;
#endif
  }
  citer = p4est_lnodes_corner_callback;

  p4est_iterate (p4est, ghost_layer, &data, viter, fiter,
#ifdef P4_TO_P8
                 eiter,
#endif
                 citer);

  p4est_lnodes_count_nodes (&data, p4est, ghost_layer, lnodes);

  P4EST_FREE (data.local_cdp);
  P4EST_FREE (data.ghost_cdp);
#ifdef P4_TO_P8
  P4EST_FREE (data.local_edp);
  P4EST_FREE (data.ghost_edp);
#endif

  while (data.hfaces->elem_count) {
    hface = sc_array_pop (data.hfaces);
    p4est_lnodes_hface_fix (p4est, hface, &data);
  }
#ifdef P4_TO_P8
  while (data.hedges->elem_count) {
    hedge = sc_array_pop (data.hedges);
    p8est_lnodes_hedge_fix (p4est, hedge, &data);
  }
#endif

  P4EST_FREE (data.ghost_elem_nodes);

#ifdef P4EST_DEBUG
  for (li = 0; li < nlen; li++) {
    P4EST_ASSERT (lnodes->local_nodes[li] >= 0);
  }
#endif

  for (i = 0; i < mpisize; i++) {
    sc_array_sort (&(data.send_buf_info[i]), p4est_lnodes_binfo_compare);
    sc_array_sort (&(data.recv_buf_info[i]), p4est_lnodes_binfo_compare);
#ifdef P4EST_DEBUG
    last = NULL;
    for (zz = 0; zz < data.send_buf_info[i].elem_count; zz++) {
      new = sc_array_index (&(data.send_buf_info[i]), zz);
      if (zz > 0) {
        P4EST_ASSERT (p4est_lnodes_binfo_compare (last, new) < 0);
      }
      last = new;
    }
    last = NULL;
    for (zz = 0; zz < data.recv_buf_info[i].elem_count; zz++) {
      new = sc_array_index (&(data.recv_buf_info[i]), zz);
      if (zz > 0) {
        if (p4est_lnodes_binfo_is_equal (last, new)) {
          P4EST_VERBOSE ("Doubled receive buffer entry.\n");
        }
      }
      last = new;
    }
#endif
  }

#ifdef P4EST_MPI
  P4EST_ASSERT (p4est_lnodes_test_comm (p4est, &data));

  p4est_lnodes_pass (p4est, &data);
#endif

  p4est_lnodes_global_and_sharers (&data, lnodes, p4est);

  p4est_lnodes_reset_data (&data, p4est);

#ifdef P4_TO_P8
  for (li = 0; li < nel; li++) {
    lnodes->face_code[li] = p8est_lnodes_code_fix (lnodes->face_code[li]);
  }
#endif

  return lnodes;
}

void
p4est_lnodes_destroy (p4est_lnodes_t * lnodes)
{
  size_t              zz, count;
  p4est_lnodes_rank_t *lrank;
  P4EST_FREE (lnodes->local_nodes);
  P4EST_FREE (lnodes->global_nodes);
  P4EST_FREE (lnodes->face_code);
  count = lnodes->sharers->elem_count;
  for (zz = 0; zz < count; zz++) {
    lrank = sc_array_index (lnodes->sharers, zz);
    sc_array_reset (&(lrank->shared_nodes));
  }
  sc_array_destroy (lnodes->sharers);
  P4EST_FREE (lnodes);
}
