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
 * quadrants q1, q3, and p are present, but q0 and q2, and thus processes 0 and
 * 2, also share the nodes created on that edge.  The edp associated with q1's
 * edge will have one value set to 0, and the edp associated with q3's edge
 * will have one value set to 2.
 */
typedef struct p8est_lnodes_edp
{
  int                 face[2];
}
p8est_lnodes_edp_t;
#endif

/** inode: independent node.
 * If the node is only needed by the local process, then share_offset = -1 and
 * share_count = 0.  Otherwise these values index into the inode_sharers array
 * for a list of all processes that share it (local process included).
 */
typedef struct p4est_lnodes_inode
{
  int                 owner;
  p4est_locidx_t      share_offset;
  int8_t              share_count;
}
p4est_lnodes_inode_t;

/** buf_info: encodes/decodes the transmission of node information. */
typedef struct p4est_lnodes_buf_info
{
  p4est_quadrant_t    q;        /* p.which_tree filled */
  int8_t              type;     /* which nodes it shares */
  bool                send_sharers;     /* whether the sharers are included in the message */
  p4est_locidx_t      first_index;      /* in inodes array, first node to/from */
}
p4est_lnodes_buf_info_t;

typedef struct p4est_lnodes_iter_data
{
  p4est_lnodes_cdp_t *local_cdp;        /* num of local quads * corners per quad */
  p4est_lnodes_cdp_t *ghost_cdp;        /* num of ghost quads * corners per quad */
#ifdef P4_TO_P8
  p8est_lnodes_edp_t *local_edp;        /* num of local quads * edges per quad */
  p8est_lnodes_edp_t *ghost_edp;        /* num of ghost quads * edges per quad */
#endif
  p4est_locidx_t     *local_elem_nodes; /* number of local q's * nodes per q */
  p4est_locidx_t     *ghost_elem_nodes; /* number of ghost q's * nodes per q */
  sc_array_t         *hfaces;   /* p4est_iter_face_side_t: store hanging faces */
#ifdef P4_TO_P8
  sc_array_t         *hedges;   /* p8est_iter_edge_side_t: store hanging edges */
#endif
  sc_array_t         *inodes;   /* p4est_lnodes_inode_t */
  sc_array_t         *inode_sharers;    /* int */
  sc_array_t         *send_buf_info;    /* one for each proc: type buf_info_t */
  sc_array_t         *recv_buf_info;    /* one for each proc: type buf_info_t */
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
}
p4est_lnodes_iter_data_t;

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
p4est_lnodes_face_simple_callback (p4est_iter_face_info_t * info, void *data)
{
  int                 i;
  int                 c;
  int                 f;
  size_t              zz;
  p4est_lnodes_iter_data_t *iter_data = data;
  sc_array_t         *sides = &(info->sides);
  size_t              count = sides->elem_count;
  p4est_iter_face_side_t *fside;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_topidx_t      tid;
  p4est_locidx_t      qid;
  p4est_lnodes_cdp_t *local_cdp = iter_data->local_cdp;
  p4est_lnodes_cdp_t *ghost_cdp = iter_data->ghost_cdp;
  p4est_locidx_t      cdpid;
#ifdef P4_TO_P8
  p8est_lnodes_edp_t *local_edp = iter_data->local_edp;
  p8est_lnodes_edp_t *ghost_edp = iter_data->ghost_edp;
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
  sc_array_t         *hfaces = iter_data->hfaces;
  int                 fdir;
  int                 rank = info->p4est->mpirank;
  p4est_locidx_t      hqid[P4EST_CHILDREN / 2];
  p4est_locidx_t      quadrants_offset;
#ifndef P4_TO_P8
  int8_t             *face_codes = iter_data->face_codes;
#else
  int16_t            *face_codes = iter_data->face_codes;
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
          if (p8est_edge_faces[he[1]][0] == f) {
            eside = 0;
          }
          else {
            eside = 1;
          }
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
        if (p8est_edge_faces[e][0] == f) {
          eside = 0;
        }
        else {
          eside = 1;
        }
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
p8est_lnodes_edge_simple_callback (p8est_iter_edge_info_t * info, void *data)
{
  int                 i;
  int                 c;
  size_t              zz;
  p4est_lnodes_iter_data_t *iter_data = data;
  sc_array_t         *sides = &(info->sides);
  size_t              count = sides->elem_count;
  p8est_iter_edge_side_t *eside;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_topidx_t      tid;
  p4est_locidx_t      qid;
  p4est_lnodes_cdp_t *local_cdp = iter_data->local_cdp;
  p4est_lnodes_cdp_t *ghost_cdp = iter_data->ghost_cdp;
  p4est_locidx_t      cdpid;
  int                 e;
  bool                is_local, *h_is_local;
  int                 procs[2];
  p8est_iter_edge_side_t *hedge;
  sc_array_t         *hedges = iter_data->hedges;
  int                 edir;
  int                 rank = info->p4est->mpirank;
  p4est_locidx_t      hqid[2];
  p4est_locidx_t      quadrants_offset;
  int16_t            *face_codes = iter_data->face_codes;
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
  int                 f = p4est_corner_faces[c][dir];
#ifndef P4_TO_P8
  int                 j = p4est_face_corners[f][0] == c ? 0 : 1;
#else
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
    if (cside->treeid != ntid) {
      continue;
    }
    if (cside->corner != nc) {
      continue;
    }
    level = cside->quad->level;
    P4EST_ASSERT (nlevel <= level && level <= nlevel + 2);
    if (level == nlevel) {
      return nproc;
    }
    else {
      return -2;
    }
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
      if (cside->treeid != ntid) {
        continue;
      }
      if (cside->corner != nc) {
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
      if (cside->treeid != ntid) {
        continue;
      }
      if (cside->corner != nc) {
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
    if (eside->treeid != ntid) {
      continue;
    }
    if (eside->edge != ne) {
      continue;
    }
    if (!eside->is_hanging) {
      return nproc;
    }
    else {
      return -2;
    }
  }

  /* not reached because the face neighbor should be in the edge iterate info
   */
  SC_ABORT_NOT_REACHED ();

  return -2;
}
#endif

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
p4est_lnodes_corner_callback (p4est_iter_corner_info_t * info, void *data)
{
  int                 i;
  size_t              zz;
  p4est_lnodes_iter_data_t *iter_data = data;
  sc_array_t         *sides = &(info->sides);
  size_t              count = sides->elem_count;
  p4est_iter_corner_side_t *cside;
  sc_array_t         *inodes = iter_data->inodes;
  p4est_lnodes_inode_t *inode;
  sc_array_t         *inode_sharers = iter_data->inode_sharers;
  p4est_lnodes_cdp_t *local_cdp = iter_data->local_cdp;
  p4est_lnodes_cdp_t *ghost_cdp = iter_data->ghost_cdp;
  p4est_lnodes_cdp_t *cdp;
  p4est_locidx_t     *local_elem_nodes = iter_data->local_elem_nodes;
  p4est_locidx_t     *ghost_elem_nodes = iter_data->ghost_elem_nodes;
  p4est_locidx_t     *elnode;
  sc_array_t         *send_buf_info = iter_data->send_buf_info;
  sc_array_t         *recv_buf_info = iter_data->recv_buf_info;
  p4est_lnodes_buf_info_t *binfo;
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
  p4est_locidx_t      sharers_offset =
    (p4est_locidx_t) inode_sharers->elem_count;
  int8_t              new_count = 0;
  int                *corner_nodes = iter_data->corner_nodes;
  int                 nodes_per_elem = iter_data->nodes_per_elem;
  p4est_locidx_t      quadrants_offset;

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
      ip = sc_array_push (&touching_procs);
      *ip = proc;
      ip = sc_array_push (&all_procs);
      *ip = proc;
    }

    if (proc < owner_proc) {
      owner_proc = proc;
      owner_tid = tid;
      owner_quad = q;
      owner_corner = c;
    }
    else if (proc == owner_proc) {
      if (tid < owner_tid) {
        owner_proc = proc;
        owner_tid = tid;
        owner_quad = q;
        owner_corner = c;
      }
      else if (tid == owner_tid) {
        if (p4est_quadrant_compare (q, owner_quad) < 0) {
          owner_proc = proc;
          owner_tid = tid;
          owner_quad = q;
          owner_corner = c;
        }
      }
    }

    cdpid = qid * P4EST_CHILDREN + c;
    cdp = is_local ? &(local_cdp[cdpid]) : &(ghost_cdp[cdpid]);
    for (i = 0; i < P4EST_DIM; i++) {
      proc = cdp->face[i];
      if (proc >= 0 && proc != rank) {
        ip = sc_array_push (&all_procs);
        *ip = proc;
      }
      else if (proc == -1) {
        proc = p4est_lnodes_missing_proc_cdp_face (q, tid, c, i, info);
        if (proc >= 0 && proc != rank) {
          ip = sc_array_push (&all_procs);
          *ip = proc;
        }
      }
#ifdef P4_TO_P8
      proc = cdp->edge[i];
      if (proc >= 0 && proc != rank) {
        ip = sc_array_push (&all_procs);
        *ip = proc;
      }
      else if (proc == -1) {
        proc = p8est_lnodes_missing_proc_cdp_edge (q, tid, c, i, info);
        if (proc >= 0 && proc != rank) {
          ip = sc_array_push (&all_procs);
          *ip = proc;
        }
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
  inode->owner = owner_proc;
  count = all_procs.elem_count;
  if (count) {
    p4est_quadrant_smallest_corner_descendent (owner_quad, &tempq,
                                               owner_corner);
    inode->share_offset = sharers_offset;
    ip = sc_array_push (inode_sharers);
    *ip = rank;
    new_count++;
    for (zz = 0; zz < count; zz++, new_count++) {
      ip = sc_array_push (inode_sharers);
      proc = *((int *) sc_array_index (&all_procs, zz));
      *ip = proc;
      if (owner_proc == rank) {
        binfo = sc_array_push (&(send_buf_info[proc]));
        binfo->q = tempq;
        binfo->q.p.which_tree = owner_tid;
        binfo->type = (int8_t) (P4EST_LN_C_OFFSET + owner_corner);
        if (sc_array_bsearch (&touching_procs, &proc, sc_int_compare) >= 0) {
          binfo->send_sharers = false;
        }
        else {
          binfo->send_sharers = true;
        }
        binfo->first_index = num_inodes;
      }
      else if (proc == owner_proc) {
        binfo = sc_array_push (&(recv_buf_info[proc]));
        binfo->q = tempq;
        binfo->q.p.which_tree = owner_tid;
        binfo->type = (int8_t) (P4EST_LN_C_OFFSET + owner_corner);
        binfo->send_sharers = false;
        binfo->first_index = num_inodes;
      }
    }
    inode->share_count = new_count;
  }
  else {
    P4EST_ASSERT (owner_proc == rank);
    inode->share_offset = -1;
    inode->share_count = 0;
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
p8est_lnodes_edge_callback (p8est_iter_edge_info_t * info, void *data)
{
  int                 i;
  int                 j;
  size_t              zz;
  p4est_lnodes_iter_data_t *iter_data = data;
  sc_array_t         *sides = &(info->sides);
  size_t              count = sides->elem_count;
  p8est_iter_edge_side_t *eside;
  sc_array_t         *inodes = iter_data->inodes;
  p4est_lnodes_inode_t *inode;
  sc_array_t         *inode_sharers = iter_data->inode_sharers;
  p8est_lnodes_edp_t *local_edp = iter_data->local_edp;
  p8est_lnodes_edp_t *ghost_edp = iter_data->ghost_edp;
  p8est_lnodes_edp_t *edp;
  p4est_locidx_t     *local_elem_nodes = iter_data->local_elem_nodes;
  p4est_locidx_t     *ghost_elem_nodes = iter_data->ghost_elem_nodes;
  p4est_locidx_t     *elnode;
  sc_array_t         *send_buf_info = iter_data->send_buf_info;
  sc_array_t         *recv_buf_info = iter_data->recv_buf_info;
  p4est_lnodes_buf_info_t *binfo;
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
  p4est_locidx_t      sharers_offset =
    (p4est_locidx_t) inode_sharers->elem_count;
  int8_t              new_count = 0;
  int                 nodes_per_edge = iter_data->nodes_per_edge;
  int                 nodes_per_elem = iter_data->nodes_per_elem;
  int               **edge_nodes = iter_data->edge_nodes;
  bool                is_hanging;
  int                 limit;
  int                 stride;
  int                 orientation;
  int8_t              min_level = P4EST_QMAXLEVEL;
  int8_t              max_level = 0;
  p4est_locidx_t      start_node;

  p8est_lnodes_edge_simple_callback (info, data);

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
        ip = sc_array_push (&touching_procs);
        *ip = proc;
        ip = sc_array_push (&all_procs);
        *ip = proc;
      }

      if (proc < owner_proc) {
        owner_proc = proc;
        owner_tid = tid;
        owner_quad = q[i];
        owner_edge = e;
        owner_orientation = orientation;
      }
      else if (proc == owner_proc) {
        if (tid < owner_tid) {
          owner_proc = proc;
          owner_tid = tid;
          owner_quad = q[i];
          owner_edge = e;
          owner_orientation = orientation;
        }
        else if (tid == owner_tid) {
          if (p4est_quadrant_compare (q[i], owner_quad) < 0) {
            owner_proc = proc;
            owner_tid = tid;
            owner_quad = q[i];
            owner_edge = e;
            owner_orientation = orientation;
          }
        }
      }

      if (is_hanging) {
        edpid = qid * 12 + e;
        edp = is_local[i] ? &(local_edp[edpid]) : &(ghost_edp[edpid]);
        for (j = 0; j < 2; j++) {
          proc = edp->face[j];
          if (proc >= 0 && proc != rank) {
            ip = sc_array_push (&all_procs);
            *ip = proc;
          }
          else if (proc == -1) {
            proc = p8est_lnodes_missing_proc_edp (q[i], tid, e, j, i, info);
            if (proc >= 0 && proc != rank) {
              ip = sc_array_push (&all_procs);
              *ip = proc;
            }
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
    ip = sc_array_push (inode_sharers);
    *ip = rank;
    new_count++;
    for (zz = 0; zz < count; zz++, new_count++) {
      ip = sc_array_push (inode_sharers);
      proc = *((int *) sc_array_index (&all_procs, zz));
      *ip = proc;
      if (owner_proc == rank) {
        binfo = sc_array_push (&(send_buf_info[proc]));
        binfo->q = tempq;
        binfo->q.p.which_tree = owner_tid;
        binfo->type = (int8_t) (P8EST_LN_E_OFFSET + owner_edge);
        if (sc_array_bsearch (&touching_procs, &proc, sc_int_compare) >= 0) {
          binfo->send_sharers = false;
        }
        else {
          binfo->send_sharers = true;
        }
        binfo->first_index = num_inodes;
      }
      else if (proc == owner_proc) {
        binfo = sc_array_push (&(recv_buf_info[proc]));
        binfo->q = tempq;
        binfo->q.p.which_tree = owner_tid;
        binfo->type = (int8_t) (P8EST_LN_E_OFFSET + owner_edge);
        binfo->send_sharers = false;
        binfo->first_index = num_inodes;
      }
    }
    for (i = 0; i < nodes_per_edge; i++) {
      inode = sc_array_push (inodes);
      inode->owner = owner_proc;
      inode->share_offset = sharers_offset;
      inode->share_count = new_count;
    }
  }
  else {
    P4EST_ASSERT (owner_proc == rank);
    for (i = 0; i < nodes_per_edge; i++) {
      inode = sc_array_push (inodes);
      inode->owner = rank;
      inode->share_offset = -1;
      inode->share_count = 0;
    }
  }

  sc_array_reset (&touching_procs);
  sc_array_reset (&all_procs);
}

/** p8est_lnodes_face_node_transform:
 *
 * Compute the transformation from an independent face nodes position on a face
 * to the position on the face of the touching element.
 */
static void
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
p4est_lnodes_face_callback (p4est_iter_face_info_t * info, void *data)
{
  size_t              zz;
  p4est_lnodes_iter_data_t *iter_data = data;
  sc_array_t         *sides = &(info->sides);
  size_t              count = sides->elem_count;
  p4est_iter_face_side_t *fside;
  sc_array_t         *inodes = iter_data->inodes;
  p4est_lnodes_inode_t *inode;
  sc_array_t         *inode_sharers = iter_data->inode_sharers;
  p4est_locidx_t     *local_elem_nodes = iter_data->local_elem_nodes;
  p4est_locidx_t     *ghost_elem_nodes = iter_data->ghost_elem_nodes;
  p4est_locidx_t     *elnode;
  sc_array_t         *send_buf_info = iter_data->send_buf_info;
  sc_array_t         *recv_buf_info = iter_data->recv_buf_info;
  p4est_lnodes_buf_info_t *binfo;
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
  p4est_locidx_t      sharers_offset =
    (p4est_locidx_t) inode_sharers->elem_count;
  int8_t              new_count = 0;
  int                 nodes_per_face = iter_data->nodes_per_face;
  int                 nodes_per_elem = iter_data->nodes_per_elem;
  int               **face_nodes = iter_data->face_nodes;
  bool                is_hanging;
  int                 limit;
  int                 i, j;
#ifndef P4_TO_P8
  int                 stride;
#else
  int                 k, l;
  int                 nodes_per_edge = iter_data->nodes_per_edge;
  bool                flipj, flipk, swapjk;
  int                 jind, kind, lind;
#endif
  p4est_locidx_t      start_node;

  p4est_lnodes_face_simple_callback (info, data);

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
        ip = sc_array_push (&touching_procs);
        *ip = proc;
      }

      if (proc < owner_proc) {
        owner_proc = proc;
        owner_tid = tid;
        owner_quad = q[i];
        owner_face = f;
        owner_side = zz;
      }
      else if (proc == owner_proc) {
        if (tid < owner_tid) {
          owner_proc = proc;
          owner_tid = tid;
          owner_quad = q[i];
          owner_face = f;
          owner_side = zz;
        }
        else if (tid == owner_tid) {
          if (p4est_quadrant_compare (q[i], owner_quad) < 0) {
            owner_proc = proc;
            owner_tid = tid;
            owner_quad = q[i];
            owner_face = f;
            owner_side = zz;
          }
        }
      }
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
      }
#ifndef P4_TO_P8
      start_node = num_inodes + ((zz == owner_side) ? 0 :
                                 info->orientation ? 0 : nodes_per_face - 1);
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
  sc_array_sort (&touching_procs, sc_int_compare);
  sc_array_uniq (&touching_procs, sc_int_compare);

  count = touching_procs.elem_count;
  if (count) {
    tempq = *owner_quad;
    ip = sc_array_push (inode_sharers);
    *ip = rank;
    new_count++;
    for (zz = 0; zz < count; zz++, new_count++) {
      ip = sc_array_push (inode_sharers);
      proc = *((int *) sc_array_index (&touching_procs, zz));
      *ip = proc;
      if (owner_proc == rank) {
        binfo = sc_array_push (&(send_buf_info[proc]));
        binfo->q = tempq;
        binfo->q.p.which_tree = owner_tid;
        binfo->type = (int8_t) owner_face;
        binfo->send_sharers = false;
        binfo->first_index = num_inodes;
      }
      else if (proc == owner_proc) {
        binfo = sc_array_push (&(recv_buf_info[proc]));
        binfo->q = tempq;
        binfo->q.p.which_tree = owner_tid;
        binfo->type = (int8_t) owner_face;
        binfo->send_sharers = false;
        binfo->first_index = num_inodes;
      }
    }
    for (i = 0; i < nodes_per_face; i++) {
      inode = sc_array_push (inodes);
      inode->owner = owner_proc;
      inode->share_offset = sharers_offset;
      inode->share_count = new_count;
    }
  }
  else {
    P4EST_ASSERT (owner_proc == rank);
    for (i = 0; i < nodes_per_face; i++) {
      inode = sc_array_push (inodes);
      inode->owner = rank;
      inode->share_offset = -1;
      inode->share_count = 0;
    }
  }

  sc_array_reset (&touching_procs);
}

/* p4est_lnodes_volume_callback:
 *
 * Create independent nodes and set volume nodes to point to them.
 */
static void
p4est_lnodes_volume_callback (p4est_iter_volume_info_t * info, void *data)
{
  p4est_lnodes_iter_data_t *iter_data = data;
  p4est_tree_t       *tree = p4est_array_index_topidx (info->p4est->trees,
                                                       info->treeid);
  p4est_locidx_t      qid = info->quadid + tree->quadrants_offset;
  p4est_locidx_t     *elem_nodes = iter_data->local_elem_nodes;
  sc_array_t         *inodes = iter_data->inodes;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) inodes->elem_count;
  p4est_locidx_t      nid;
  p4est_lnodes_inode_t *inode;
  int                 nodes_per_volume = iter_data->nodes_per_volume;
  int                *volume_nodes = iter_data->volume_nodes;
  int                 nodes_per_elem = iter_data->nodes_per_elem;
  int                 i;

  for (i = 0; i < nodes_per_volume; i++) {
    nid = qid * nodes_per_elem + volume_nodes[i];
    elem_nodes[nid] = num_inodes + (p4est_locidx_t) i;
    inode = sc_array_push (inodes);
    inode->owner = info->p4est->mpirank;
    inode->share_offset = -1;
    inode->share_count = 0;
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
                                  p4est_lnodes_iter_data_t * iter_data)
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
  sc_array_t         *inodes = iter_data->inodes;
  p4est_lnodes_inode_t *inode;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) inodes->elem_count;
  sc_array_t         *recv_buf_info = iter_data->recv_buf_info;
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
    else if (ownertid == ntid) {
      if (p4est_quadrant_compare (&tempr, &ownerq) < 0) {
        ownerq = tempr;
      }
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
        else if (ownertid == ntid) {
          if (p4est_quadrant_compare (&tempr, &ownerq) < 0) {
            ownerq = *ptemp;
          }
        }
      }
      sc_array_reset (&quads);
      sc_array_reset (&treeids);
    }
    else {
      if (p4est_quadrant_compare (&tempr, &ownerq) < 0) {
        ownerq = tempr;
      }
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
      else if (ownertid == ntid) {
        if (p4est_quadrant_compare (ptemp, &ownerq) < 0) {
          ownerq = *ptemp;
        }
      }
    }
    sc_array_reset (&quads);
    sc_array_reset (&treeids);
  }
  else {
    if (p4est_quadrant_compare (&tempr, &ownerq) < 0) {
      ownerq = tempr;
    }
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
  inode = sc_array_push (inodes);
  inode->owner = owner_proc;
  inode->share_offset = -1;
  inode->share_count = -1;

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
                                p4est_lnodes_iter_data_t * iter_data,
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
  sc_array_t         *inodes = iter_data->inodes;
  p4est_lnodes_inode_t *inode;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) inodes->elem_count;
  sc_array_t         *recv_buf_info = iter_data->recv_buf_info;
  p4est_lnodes_buf_info_t *binfo;
  int                 owner_proc;
  int                 orientation;
  int                 nodes_per_edge = iter_data->nodes_per_edge;
  int                 nodes_per_elem = iter_data->nodes_per_elem;
  int                 c;
  p4est_locidx_t      qid;
  bool                is_local;
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree;
  p4est_locidx_t      quadrants_offset;
  p4est_locidx_t     *local_elem_nodes = iter_data->local_elem_nodes;
  p4est_locidx_t     *ghost_elem_nodes = iter_data->ghost_elem_nodes;
  p4est_locidx_t      nid;
  int                 stride;
  p4est_locidx_t      start_node;
  int               **edge_nodes = iter_data->edge_nodes;

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
      else if (ownertid[j] == ntid) {
        if (p4est_quadrant_compare (&tempr[j], &ownerq[j]) < 0) {
          ownerq[j] = tempr[j];
        }
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
        else if (ownertid[j] == ntid) {
          if (p4est_quadrant_compare (ptemp, &ownerq[j]) < 0) {
            ownerq[j] = *ptemp;
          }
        }
      }
      sc_array_reset (&quads);
      sc_array_reset (&treeids);
    }
    else {
      if (p4est_quadrant_compare (&tempr[j], &ownerq[j]) < 0) {
        ownerq[j] = tempr[j];
      }
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

  for (i = 0; i < nodes_per_edge; i++) {
    inode = sc_array_push (inodes);
    inode->owner = owner_proc;
    inode->share_offset = -1;
    inode->share_count = -1;
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
                        p4est_lnodes_iter_data_t * iter_data)
{
  int                 i;
  p4est_topidx_t      tid = hface->treeid;
  p4est_locidx_t      qid[P4EST_CHILDREN / 2];
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree = p4est_array_index_topidx (trees, tid);
  p4est_locidx_t      quadrants_offset = tree->quadrants_offset;
  bool               *is_local = hface->is.hanging.is_local;
  int                *corner_nodes = iter_data->corner_nodes;
  int                 nodes_per_elem = iter_data->nodes_per_elem;
  p4est_locidx_t      nid, nidval, nid2;
  int                 f = hface->face;
  p4est_locidx_t     *local_elem_nodes = iter_data->local_elem_nodes;
  p4est_locidx_t     *ghost_elem_nodes = iter_data->ghost_elem_nodes;
  int                 c, c2;
  int                 i2;
#ifdef P4_TO_P8
  int                 nodes_per_edge = iter_data->nodes_per_edge;
  int               **edge_nodes = iter_data->edge_nodes;
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
                                                   tid, c, p4est, iter_data);
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
                                            e, p4est, iter_data, hface);
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
                        p4est_lnodes_iter_data_t * iter_data)
{
  int                 i;
  p4est_topidx_t      tid = hedge->treeid;
  p4est_locidx_t      qid[2];
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree = p4est_array_index_topidx (trees, tid);
  p4est_locidx_t      quadrants_offset = tree->quadrants_offset;
  bool               *is_local = hedge->is.hanging.is_local;
  int                *corner_nodes = iter_data->corner_nodes;
  int                 nodes_per_elem = iter_data->nodes_per_elem;
  p4est_locidx_t      nid, nidval, nid2;
  int                 e = hedge->edge;
  p4est_locidx_t     *local_elem_nodes = iter_data->local_elem_nodes;
  p4est_locidx_t     *ghost_elem_nodes = iter_data->ghost_elem_nodes;
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
                                                   tid, c, p4est, iter_data);
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
p4est_lnodes_buf_info_compar (const void *a, const void *b)
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
p4est_lnodes_init_iter_data (p4est_lnodes_iter_data_t * iter_data,
                             p4est_t * p4est, sc_array_t * ghost_layer, int p)
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
  npel = iter_data->nodes_per_elem = (p + 1) * (p + 1);
  npv = iter_data->nodes_per_volume = (p - 1) * (p - 1);
  npf = iter_data->nodes_per_face = p - 1;
  fcount[0] = fcount[1] = fcount[2] = fcount[3] = 0;
#else
  npel = iter_data->nodes_per_elem = (p + 1) * (p + 1) * (p + 1);
  npv = iter_data->nodes_per_volume = (p - 1) * (p - 1) * (p - 1);
  npf = iter_data->nodes_per_face = (p - 1) * (p - 1);
  fcount[0] = fcount[1] = fcount[2] = fcount[3] = fcount[4] = fcount[5] = 0;
  npe = iter_data->nodes_per_edge = p - 1;
  ecount[0] = ecount[1] = ecount[2] = ecount[3] = ecount[4] = ecount[5] = 0;
  ecount[6] = ecount[7] = ecount[8] = ecount[9] = ecount[10] = ecount[11] = 0;
#endif
  vcount = 0;
  nlen = nlq * npel;
  ngen = ngq * npel;

  iter_data->volume_nodes = P4EST_ALLOC (int, npv);
  for (i = 0; i < P4EST_DIM * 2; i++) {
    iter_data->face_nodes[i] = P4EST_ALLOC (int, npf);
  }
#ifdef P4_TO_P8
  for (i = 0; i < 12; i++) {
    iter_data->edge_nodes[i] = P4EST_ALLOC (int, npf);
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
#endif
        switch (i == 0 ? 0 : i == p ? 1 : 2) {
        case 0:
          bcount++;
          break;
        case 1:
          f = 1;
          c |= 1;
#ifdef P4_TO_P8
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
#ifdef P4_TO_P8
        switch (k == 0 ? 0 : k == p ? 1 : 2) {
        case 0:
          f = 4;
          e <<= 1;
          bcount++;
          break;
        case 1:
          f = 5;
          c |= 4;
          e <<= 1;
          e++;
          bcount++;
          break;
        default:
          eshift = 8;
          break;
        }
#endif
        switch (bcount) {
        case 0:
          iter_data->volume_nodes[vcount++] = n;
          break;
        case 1:
          iter_data->face_nodes[f][fcount[f]++] = n;
          break;
#ifdef P4_TO_P8
        case 2:
          P4EST_ASSERT (eshift >= 0);
          e += eshift;
          iter_data->edge_nodes[e][ecount[e]++] = n;
          break;
#endif
        default:
          iter_data->corner_nodes[c] = n;
          break;
        }
      }
    }
#ifdef P4_TO_P8
  }
#endif

  iter_data->local_cdp = P4EST_ALLOC (p4est_lnodes_cdp_t, nlcdp);
  memset (iter_data->local_cdp, -1, nlcdp * sizeof (p4est_lnodes_cdp_t));
  iter_data->ghost_cdp = P4EST_ALLOC (p4est_lnodes_cdp_t, ngcdp);
  memset (iter_data->ghost_cdp, -1, ngcdp * sizeof (p4est_lnodes_cdp_t));
#ifdef P4_TO_P8
  iter_data->local_edp = P4EST_ALLOC (p8est_lnodes_edp_t, nledp);
  memset (iter_data->local_edp, -1, nledp * sizeof (p8est_lnodes_edp_t));
  iter_data->ghost_edp = P4EST_ALLOC (p8est_lnodes_edp_t, ngedp);
  memset (iter_data->ghost_edp, -1, ngedp * sizeof (p8est_lnodes_edp_t));
#endif

  iter_data->local_elem_nodes = P4EST_ALLOC (p4est_locidx_t, nlen);
  memset (iter_data->local_elem_nodes, -1, nlen * sizeof (p4est_locidx_t));
  iter_data->ghost_elem_nodes = P4EST_ALLOC (p4est_locidx_t, ngen);
  memset (iter_data->ghost_elem_nodes, -1, ngen * sizeof (p4est_locidx_t));

  iter_data->hfaces = sc_array_new (sizeof (p4est_iter_face_side_t));
#ifdef P4_TO_P8
  iter_data->hedges = sc_array_new (sizeof (p8est_iter_edge_side_t));
#endif
  iter_data->inodes = sc_array_new (sizeof (p4est_lnodes_inode_t));
  iter_data->inode_sharers = sc_array_new (sizeof (int));
  iter_data->send_buf_info = P4EST_ALLOC (sc_array_t, mpisize);
  iter_data->recv_buf_info = P4EST_ALLOC (sc_array_t, mpisize);
  for (i = 0; i < mpisize; i++) {
    sc_array_init (&(iter_data->send_buf_info[i]),
                   sizeof (p4est_lnodes_buf_info_t));
    sc_array_init (&(iter_data->recv_buf_info[i]),
                   sizeof (p4est_lnodes_buf_info_t));
  }

#ifndef P4_TO_P8
  iter_data->face_codes = P4EST_ALLOC_ZERO (int8_t, nlq);
#else
  iter_data->face_codes = P4EST_ALLOC_ZERO (int16_t, nlq);
#endif
}

static void
p4est_lnodes_reset_iter_data (p4est_lnodes_iter_data_t * data,
                              p4est_t * p4est)
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

  P4EST_FREE (data->local_cdp);
  P4EST_FREE (data->ghost_cdp);
#ifdef P4_TO_P8
  P4EST_FREE (data->local_edp);
  P4EST_FREE (data->ghost_edp);
#endif
  /* do not free *_elem_nodes: control given to lnodes_t */
  sc_array_destroy (data->hfaces);
#ifdef P4_TO_P8
  sc_array_destroy (data->hedges);
#endif
  sc_array_destroy (data->inodes);
  sc_array_destroy (data->inode_sharers);
  for (i = 0; i < mpisize; i++) {
    sc_array_reset (&(data->send_buf_info[i]));
    sc_array_reset (&(data->recv_buf_info[i]));
  }
  P4EST_FREE (data->send_buf_info);
  P4EST_FREE (data->recv_buf_info);
  /* do not free face_codes: control given to lnodes_t */
}

static p4est_locidx_t *
p4est_lnodes_order_nodes (p4est_lnodes_iter_data_t * data, p4est_t * p4est,
                          sc_array_t * ghost_layer)
{
  p4est_locidx_t      num_inodes;
  p4est_locidx_t     *inode_order;
  p4est_locidx_t      nlq = p4est->local_num_quadrants;
  p4est_locidx_t      npel;
  p4est_locidx_t      nlen;
  p4est_locidx_t      li;
  p4est_locidx_t      inidx;
  sc_array_t         *inodes = data->inodes;
  p4est_lnodes_inode_t *inode;
  p4est_topidx_t     *local_en = data->local_elem_nodes;
  int                 rank = p4est->mpirank;
  p4est_locidx_t      count = 0;

  npel = data->nodes_per_elem;
  nlen = npel * nlq;
  num_inodes = (p4est_locidx_t) inodes->elem_count;
  inode_order = P4EST_ALLOC (p4est_locidx_t, num_inodes);
  memset (inode_order, -1, num_inodes * sizeof (p4est_locidx_t));
  for (li = 0; li < nlen; li++) {
    inidx = local_en[li];
    if (inidx >= 0) {
      inode = sc_array_index (inodes, (size_t) inidx);
      if (inode->owner == rank) {
        inode_order[inidx] = count++;
      }
    }
  }

  return inode_order;
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
  p4est_lnodes_iter_data_t iter_data;
  p4est_locidx_t     *inode_owner;
  p4est_iter_face_side_t *hface;
#ifdef P4_TO_P8
  p8est_iter_edge_side_t *hedge;
#endif
#ifdef P4EST_DEBUG
  p4est_locidx_t      li;
  p4est_locidx_t      nlen;
#endif

  P4EST_ASSERT (degree >= 1);

  p4est_lnodes_init_iter_data (&iter_data, p4est, ghost_layer, degree);
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

  p4est_iterate (p4est, ghost_layer, &iter_data, viter, fiter,
#ifdef P4_TO_P8
                 eiter,
#endif
                 citer);

  inode_owner = p4est_lnodes_order_nodes (&iter_data, p4est, ghost_layer);

  while (iter_data.hfaces->elem_count) {
    hface = sc_array_pop (iter_data.hfaces);
    p4est_lnodes_hface_fix (p4est, hface, &iter_data);
  }
#ifdef P4_TO_P8
  while (iter_data.hedges->elem_count) {
    hedge = sc_array_pop (iter_data.hedges);
    p8est_lnodes_hedge_fix (p4est, hedge, &iter_data);
  }
#endif

#ifdef P4EST_DEBUG
  nlen = iter_data.nodes_per_elem * p4est->local_num_quadrants;
  for (li = 0; li < nlen; li++) {
    P4EST_ASSERT (iter_data.local_elem_nodes[li] >= 0);
  }
#endif

  P4EST_FREE (inode_owner);

  P4EST_FREE (iter_data.local_elem_nodes);
  P4EST_FREE (iter_data.ghost_elem_nodes);
  P4EST_FREE (iter_data.face_codes);
  p4est_lnodes_reset_iter_data (&iter_data, p4est);
  return NULL;
}
