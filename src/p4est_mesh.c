/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

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
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_ghost.h>
#include <p4est_mesh.h>
#endif
#include <sc_ranges.h>

#ifdef P4EST_MPI

typedef struct
{
  bool                expect_query, expect_reply;
  size_t              recv_offset;
  sc_array_t          send_first, send_second;
  sc_array_t          recv_first, recv_second;
}
p4est_node_peer_t;

#endif

#ifndef P4_TO_P8

/** Generate a neighbor of a quadrant for a given node.
 *
 * The neighbor numbering is given below.
 *
 * Neighbor numbering for q, node=0:
 *
 *      ------+------+
 *            |  q   |
 *      nnum=2|nnum=3|
 *            |      |
 *      ------+------+
 *            |      |
 *      nnum=0|nnum=1|
 *
 * Neighbor numbering for q, node=1:
 *
 *            +------+------
 *            |  q   |
 *            |nnum=2|num=3
 *            |      |
 *            +------+------
 *            |      |
 *            |nnum=0|nnum=1
 *
 * Neighbor numbering for q, node=2:
 *
 *            |      |
 *      nnum=2|nnum=3|
 *            |      |
 *      ------+------+
 *            |  q   |
 *      nnum=0|nnum=1|
 *            |      |
 *      ------+------+
 *
 * Neighbor numbering for q, node=3:
 *
 *            |      |
 *            |nnum=2|nnum=3
 *            |      |
 *            +------+------
 *            |  q   |
 *            |nnum=0|nnum=1
 *            |      |
 *            +------+------
 *
 * \param [in]  q             the quadrant whose possible node \a node neighbor
 *                            will be built.
 * \param [in]  node          the node of the quadrant \a q whose possible node
 *                            neighbor list will be built.  This is given in
 *                            pixel (Morton-) ordering.
 * \param [in]  nnum          neighbor number in the ordering described above,
 *                            if nnum==node then it is the corner neighbor.
 * \param [in]  neighor_rlev  the relative level of the neighbor compared to
 *                            the level of \a q.
 * \param [out] neighbor      the neighbor that will be filled.
 * \param [out] neighbor_node the neighbor's node which shares with \a q
 *                            the node \a node.
 */
static void
p4est_possible_node_neighbor (const p4est_quadrant_t * q, int node,
                              int nnum, int neighbor_rlev,
                              p4est_quadrant_t * neighbor, int *neighbor_node)
{
  int                 nnode;
  const int           nlevel = (int) q->level + neighbor_rlev;
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);
  const p4est_qcoord_t nh = P4EST_QUADRANT_LEN (nlevel);
  const p4est_qcoord_t qx = q->x;
  const p4est_qcoord_t qy = q->y;
  p4est_qcoord_t      cornerx, cornery;
  p4est_quadrant_t    n;
#ifdef P4EST_DEBUG
  int                 qcid;
#endif

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (-1 <= neighbor_rlev && neighbor_rlev <= 1);
  P4EST_ASSERT (0 <= nlevel && nlevel <= P4EST_QMAXLEVEL);
  P4EST_ASSERT (node + nnum != 3);

  P4EST_QUADRANT_INIT (&n);

  switch (node) {
  case 0:
    cornerx = qx;
    cornery = qy;
    break;
  case 1:
    cornerx = qx + qh;
    cornery = qy;
    break;
  case 2:
    cornerx = qx;
    cornery = qy + qh;
    break;
  case 3:
    cornerx = qx + qh;
    cornery = qy + qh;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }

#ifdef P4EST_DEBUG
  /* Check to see if it is possible to construct the neighbor */
  qcid = p4est_quadrant_child_id (q);
  P4EST_ASSERT (neighbor_rlev >= 0 || qcid == node);
#endif

  nnode = 3 - nnum;
  n.level = (int8_t) nlevel;
  switch (nnum) {
  case 0:
    n.x = cornerx - nh;
    n.y = cornery - nh;
    break;
  case 1:
    n.x = cornerx;
    n.y = cornery - nh;
    break;
  case 2:
    n.x = cornerx - nh;
    n.y = cornery;
    break;
  case 3:
    n.x = cornerx;
    n.y = cornery;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }

  *neighbor = n;
  *neighbor_node = nnode;

  P4EST_ASSERT (p4est_quadrant_is_extended (neighbor));
}

void
p4est_order_local_vertices (p4est_t * p4est,
                            bool identify_periodic,
                            p4est_locidx_t * num_uniq_local_vertices,
                            p4est_locidx_t * quadrant_to_local_vertex)
{
  const int           rank = p4est->mpirank;
  p4est_connectivity_t *conn = p4est->connectivity;
  int                 qcid, transform;
  int                 neighbor_node;
  int                 face, corner, nnum, rlev, tree_corner;
  int                 neighbor_proc;
  bool                face_contact[4];
  bool                quad_contact[4];
  p4est_topidx_t      first_local_tree = p4est->first_local_tree;
  p4est_topidx_t      last_local_tree = p4est->last_local_tree;
  p4est_topidx_t      jt, num_trees = conn->num_trees;
  p4est_locidx_t      Ntotal = 0;
  p4est_locidx_t      il, Ncells = p4est->local_num_quadrants;
  p4est_locidx_t      vertex_num;
  p4est_locidx_t      lqid;
  p4est_locidx_t      neighbor_tree;
  size_t              ctree;
  size_t              zz, numz_quadrants;
  ssize_t             lnid;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  p4est_locidx_t     *tree_offset;
  p4est_tree_t       *tree, *ntree;
  sc_array_t         *trees = p4est->trees;
  sc_array_t         *quadrants, ctransforms;
  p4est_quadrant_t    neighbor, cneighbor;
  p4est_quadrant_t   *q;
  p4est_corner_transform_t *ct;

  P4EST_ASSERT (p4est_is_valid (p4est));

  P4EST_QUADRANT_INIT (&neighbor);
  P4EST_QUADRANT_INIT (&cneighbor);

  sc_array_init (&ctransforms, sizeof (p4est_corner_transform_t));

  /* figure out the offset of each tree into the local element id */
  tree_offset = P4EST_ALLOC_ZERO (p4est_locidx_t, num_trees);
  if (first_local_tree >= 0) {
    tree_offset[first_local_tree] = 0;
    for (jt = first_local_tree; jt < last_local_tree; ++jt) {
      tree = sc_array_index (trees, jt);
      tree_offset[jt + 1] = tree_offset[jt] + tree->quadrants.elem_count;
    }
  }
  else {
    P4EST_ASSERT (first_local_tree == -1 && last_local_tree == -2);
  }

  /* Initialize vertex list to all -1.  This way we know which values
   * get set because legitimate values are >= 0.
   */
  for (il = 0; il < 4 * Ncells; ++il) {
    quadrant_to_local_vertex[il] = -1;
  }

  /* loop over all local trees to generate the connetivity list */
  for (jt = first_local_tree, vertex_num = 0, lqid = 0;
       jt <= last_local_tree; ++jt) {
    for (face = 0; face < 4; ++face) {
      face_contact[face] = (conn->tree_to_tree[4 * jt + face] != jt ||
                            (identify_periodic &&
                             (int) conn->tree_to_face[4 * jt + face] !=
                             face));
    }
    tree = sc_array_index (p4est->trees, jt);
    quadrants = &tree->quadrants;
    numz_quadrants = quadrants->elem_count;

    /* Find the neighbors of each quadrant */
    for (zz = 0; zz < numz_quadrants; ++zz, ++lqid) {
      /* this quadrant may be on the boundary with a range of processors */
      q = sc_array_index (quadrants, zz);

      /* loop over the corners of the quadrant */
      for (corner = 0; corner < 4; ++corner) {

        /* Check to see if we have a new vertex */
        if (quadrant_to_local_vertex[lqid * 4 + corner] == -1) {
          quadrant_to_local_vertex[lqid * 4 + corner] = vertex_num;

          /* loop over the possible neighbors and set the new vertex */
          for (nnum = 0; nnum < 4; ++nnum) {
            /* Don't search for the quadrant q */
            if (3 - nnum == corner)
              continue;

            qcid = p4est_quadrant_child_id (q);

            /* loop over possible neighbor sizes */
            for (rlev = -1; rlev < 2; ++rlev) {
              /* can't check for quadrants larger than the root */
              if (q->level == 0 && rlev < 0)
                continue;
              /* can't check for quadrants larger unless child id
               * and corner line up
               */
              if (qcid != corner && rlev < 0)
                continue;

              /* get possible neighbor */
              p4est_possible_node_neighbor (q, corner, nnum, rlev,
                                            &neighbor, &neighbor_node);

              if (p4est_quadrant_is_inside_root (&neighbor)) {
                /* neighbor is in the same tree */

                neighbor_proc = p4est_comm_find_owner (p4est, jt, &neighbor,
                                                       rank);

                /* Neighbor is remote so we don't number its node */
                if (neighbor_proc != rank)
                  continue;

                lnid = sc_array_bsearch (quadrants, &neighbor,
                                         p4est_quadrant_compare);
                if (lnid != -1) {
                  lnid += tree_offset[jt];
                  /* We have found a neighbor in the same tree */
                  quadrant_to_local_vertex[lnid * 4 + neighbor_node]
                    = vertex_num;

                  /* No need to check for more quadrants for this neighbor */
                  continue;
                }
              }
              else {
                /* the neighbor is in a neighboring tree or multiple
                 * if it is a neighbor across the corner of the tree
                 */

                quad_contact[0] = (neighbor.y < 0);
                quad_contact[1] = (neighbor.x >= rh);
                quad_contact[2] = (neighbor.y >= rh);
                quad_contact[3] = (neighbor.x < 0);

                if ((quad_contact[0] || quad_contact[2]) &&
                    (quad_contact[1] || quad_contact[3])) {
                  /* Neighbor is across a corner */
                  for (tree_corner = 0; tree_corner < 4; ++tree_corner) {
                    if (quad_contact[(tree_corner + 3) % 4]
                        && quad_contact[tree_corner]) {
                      break;
                    }
                  }
                  p4est_find_corner_transform (conn, jt, tree_corner,
                                               &ctransforms);
                  for (ctree = 0; ctree < ctransforms.elem_count; ++ctree) {
                    ct = sc_array_index (&ctransforms, ctree);
                    neighbor_tree = ct->ntree;

                    /* Don't use corner identification in the same tree */
                    if (!identify_periodic && neighbor_tree == jt)
                      continue;

                    cneighbor = neighbor;
                    p4est_quadrant_transform_corner (&cneighbor,
                                                     (int) ct->ncorner, true);

                    neighbor_proc = p4est_comm_find_owner (p4est,
                                                           neighbor_tree,
                                                           &cneighbor, rank);

                    /* Neighbor is remote so we don't number its node */
                    if (neighbor_proc != rank)
                      continue;

                    ntree = sc_array_index (trees, neighbor_tree);

                    lnid = sc_array_bsearch (&ntree->quadrants, &cneighbor,
                                             p4est_quadrant_compare);
                    if (lnid != -1) {
                      lnid += tree_offset[neighbor_tree];
                      neighbor_node = (int) ct->ncorner;
                      /* We have found a corner neighbor */
                      quadrant_to_local_vertex[lnid * 4 + neighbor_node]
                        = vertex_num;
                    }
                  }
                }
                else {
                  /* Neighbor is across a face */
                  for (face = 0; face < 4; ++face) {
                    if (quad_contact[face] && face_contact[face]) {
                      neighbor_tree = conn->tree_to_tree[4 * jt + face];
                      break;
                    }
                  }
                  if (face == 4) {
                    /* this quadrant ran across a face with no neighbor */
                    continue;
                  }
                  /* transform the neighbor into the other tree's
                   * coordinates
                   */
                  transform = p4est_find_face_transform (conn, jt, face);
                  p4est_quadrant_translate_face (&neighbor, face);
                  p4est_quadrant_transform_face (&neighbor, &cneighbor,
                                                 transform);

                  neighbor_proc = p4est_comm_find_owner (p4est,
                                                         neighbor_tree,
                                                         &cneighbor, rank);
                  /* Neighbor is remote so we don't number its node */
                  if (neighbor_proc != rank)
                    continue;

                  ntree = sc_array_index (trees, neighbor_tree);

                  lnid = sc_array_bsearch (&ntree->quadrants, &cneighbor,
                                           p4est_quadrant_compare);
                  if (lnid != -1) {
                    lnid += tree_offset[neighbor_tree];
                    neighbor_node = p4est_node_transform (neighbor_node,
                                                          transform);

                    /* We have found a face neighbor */
                    quadrant_to_local_vertex[lnid * 4 + neighbor_node]
                      = vertex_num;
                  }
                }
              }
            }
          }
          ++vertex_num;
        }
      }
    }
  }

  Ntotal = vertex_num;
  P4EST_FREE (tree_offset);
  sc_array_reset (&ctransforms);

  *num_uniq_local_vertices = Ntotal;
}

#endif /* !P4_TO_P8 */

/** Determine the owning tree for a node and clamp it inside the domain.
 *
 * If the node is on the boundary, assign the lowest tree to own it.
 * Clamp it inside the tree bounds if necessary.
 *
 * \param [in] p4est    The p4est to work on.
 * \param [in] treeid   Original tree index for this node.
 * \param [in] n        The node to work on.
 * \param [out] c       The clamped node in owning tree coordinates.
 *                      Its piggy data will be filled with owning tree id.
 */
static void
p4est_node_canonicalize (p4est_t * p4est, p4est_topidx_t treeid,
                         const p4est_quadrant_t * n, p4est_quadrant_t * c)
{
  int                 face_axis[P4EST_DIM];
  int                 quad_contact[2 * P4EST_DIM];
  int                 contacts, face, corner;
  size_t              ctreez;
  p4est_topidx_t      ntreeid, lowest;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_quadrant_t    tmpq, o;
#ifndef P4_TO_P8
  int                 transform;
  p4est_corner_transform_t *ct;
  sc_array_t          ctransforms, *cta;
#else
  int                 edge;
  int                 ftransform[9];
  size_t              etreez;
  p4est_topidx_t      ntreeid2;
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  p8est_corner_info_t ci;
  p8est_corner_transform_t *ct;
  sc_array_t         *eta, *cta;
#endif

  P4EST_ASSERT (treeid >= 0 && treeid < conn->num_trees);
  P4EST_ASSERT (p4est_quadrant_is_node (n, false));

  P4EST_QUADRANT_INIT (&tmpq);
  P4EST_QUADRANT_INIT (&o);

  lowest = treeid;
  p4est_node_clamp_inside (n, c);
  c->p.which_tree = -1;

  /* Check if the quadrant is inside the tree */
#ifndef P4_TO_P8
  quad_contact[0] = (n->y == 0);
  quad_contact[1] = (n->x == P4EST_ROOT_LEN);
  quad_contact[2] = (n->y == P4EST_ROOT_LEN);
  quad_contact[3] = (n->x == 0);
  face_axis[0] = quad_contact[1] || quad_contact[3];
  face_axis[1] = quad_contact[0] || quad_contact[2];
  contacts = face_axis[0] + face_axis[1];
#else
  quad_contact[0] = (n->x == 0);
  quad_contact[1] = (n->x == P4EST_ROOT_LEN);
  quad_contact[2] = (n->y == 0);
  quad_contact[3] = (n->y == P4EST_ROOT_LEN);
  quad_contact[4] = (n->z == 0);
  quad_contact[5] = (n->z == P4EST_ROOT_LEN);
  face_axis[0] = quad_contact[0] || quad_contact[1];
  face_axis[1] = quad_contact[2] || quad_contact[3];
  face_axis[2] = quad_contact[4] || quad_contact[5];
  contacts = face_axis[0] + face_axis[1] + face_axis[2];
#endif
  if (contacts == 0) {
    goto endfunction;
  }

  /* Check face neighbors */
#ifdef P4EST_DEBUG
  ntreeid = -1;
#endif
  for (face = 0; face < 2 * P4EST_DIM; ++face) {
    if (!quad_contact[face]) {
      /* The node is not touching this face */
      continue;
    }
    ntreeid = conn->tree_to_tree[2 * P4EST_DIM * treeid + face];
    if (ntreeid == treeid
        && ((int) conn->tree_to_face[2 * P4EST_DIM * treeid + face] ==
            face)) {
      /* The node touches a face with no neighbor */
      continue;
    }
    if (ntreeid > lowest) {
      /* This neighbor tree is higher, so we keep the ownership */
      continue;
    }
#ifndef P4_TO_P8
    /* Transform the node into the other tree's coordinates */
    transform = p4est_find_face_transform (conn, treeid, face);
    tmpq = *n;
    p4est_quadrant_translate_face (&tmpq, face);
    p4est_quadrant_transform_face (&tmpq, &o, transform);
#else
    ntreeid2 = p8est_find_face_transform (conn, treeid, face, ftransform);
    P4EST_ASSERT (ntreeid2 == ntreeid);
    p8est_quadrant_transform_face (n, &o, ftransform);
#endif
    if (ntreeid < lowest) {
      /* we have found a new owning tree */
      p4est_node_clamp_inside (&o, c);
      lowest = ntreeid;
    }
    else {
      P4EST_ASSERT (lowest == ntreeid);
      p4est_node_clamp_inside (&o, &tmpq);
      if (p4est_quadrant_compare (&tmpq, c) < 0) {
        /* same tree (periodic) and the new position is lower than the old */
        *c = tmpq;
      }
    }
  }
  P4EST_ASSERT (ntreeid >= 0);
  if (contacts == 1) {
    goto endfunction;
  }

#ifdef P4_TO_P8
  P4EST_ASSERT (contacts >= 2);
  eta = &ei.edge_transforms;
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
  for (edge = 0; edge < 12; ++edge) {
    if (!(quad_contact[p8est_edge_faces[edge][0]] &&
          quad_contact[p8est_edge_faces[edge][1]])) {
      continue;
    }
    p8est_find_edge_transform (conn, treeid, edge, &ei);
    for (etreez = 0; etreez < eta->elem_count; ++etreez) {
      et = sc_array_index (eta, etreez);
      ntreeid = et->ntree;
      if (ntreeid > lowest) {
        /* This neighbor tree is higher, so we keep the ownership */
        continue;
      }
      p8est_quadrant_transform_edge (n, &o, &ei, et, false);
      if (ntreeid < lowest) {
        p4est_node_clamp_inside (&o, c);
        lowest = ntreeid;
      }
      else {
        P4EST_ASSERT (lowest == ntreeid);
        p4est_node_clamp_inside (&o, &tmpq);
        if (p4est_quadrant_compare (&tmpq, c) < 0) {
          /* same tree (periodic) and the new position is lower than the old */
          *c = tmpq;
        }
      }
    }
  }
  sc_array_reset (eta);
  eta = NULL;
  et = NULL;
  if (contacts == 2) {
    goto endfunction;
  }
#endif

  P4EST_ASSERT (contacts == P4EST_DIM);
#ifndef P4_TO_P8
  cta = &ctransforms;
  sc_array_init (cta, sizeof (p4est_corner_transform_t));
#else
  cta = &ci.corner_transforms;
  sc_array_init (cta, sizeof (p8est_corner_transform_t));
#endif
  for (corner = 0; corner < P4EST_CHILDREN; ++corner) {
#ifndef P4_TO_P8
    if (!(quad_contact[(corner + 3) % 4] && quad_contact[corner])) {
      continue;
    }
    p4est_find_corner_transform (conn, treeid, corner, cta);
#else
    if (!(quad_contact[p8est_corner_faces[corner][0]] &&
          quad_contact[p8est_corner_faces[corner][1]] &&
          quad_contact[p8est_corner_faces[corner][2]])) {
      continue;
    }
    p8est_find_corner_transform (conn, treeid, corner, &ci);
#endif
    for (ctreez = 0; ctreez < cta->elem_count; ++ctreez) {
      ct = sc_array_index (cta, ctreez);
      ntreeid = ct->ntree;
      if (ntreeid > lowest) {
        /* This neighbor tree is higher, so we keep the ownership */
        continue;
      }
      o.level = P4EST_MAXLEVEL;
      p4est_quadrant_transform_corner (&o, (int) ct->ncorner, false);
      if (ntreeid < lowest) {
        p4est_node_clamp_inside (&o, c);
        lowest = ntreeid;
      }
      else {
        P4EST_ASSERT (lowest == ntreeid);
        p4est_node_clamp_inside (&o, &tmpq);
        if (p4est_quadrant_compare (&tmpq, c) < 0) {
          /* same tree (periodic) and the new position is lower than the old */
          *c = tmpq;
        }
      }
    }
  }
  sc_array_reset (cta);

endfunction:
  c->p.which_tree = lowest;
}

static              bool
p4est_nodes_foreach (void **item, const void *u)
{
  const sc_hash_array_data_t *internal_data = u;
  const p4est_locidx_t *new_node_number = internal_data->user_data;

  *item = (void *) (long) new_node_number[(long) *item];

  return true;
}

p4est_nodes_t      *
p4est_nodes_new (p4est_t * p4est, sc_array_t * ghost_layer)
{
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
#ifdef P4EST_MPI
  const int           twopeerw = 2 * p4est_num_ranges;
  int                 mpiret;
  int                 owner, prev, start;
  int                 first_peer, last_peer, nwin;
  int                 num_send_queries, num_send_nonzero, num_recv_queries;
  int                 byte_count, elem_count;
  int                 local_send_count, local_recv_count;
  int                 my_ranges[twopeerw];
  int                *procs;
  int                *all_ranges;
  int                *old_sharers, *new_sharers;
  char               *this_base;
  bool                found;
  size_t              first_size, second_size, this_size;
  size_t              num_sharers, old_position, new_position;
  p4est_qcoord_t     *xyz;
  p4est_topidx_t     *ttt;
  p4est_locidx_t     *node_number;
  p4est_node_peer_t  *peers, *peer;
  p4est_indep_t       inkey;
  sc_array_t          send_requests;
  sc_recycle_array_t *orarr, *nrarr;
  MPI_Request        *send_request;
  MPI_Status          probe_status, recv_status;
#endif
#if defined (P4EST_MPI) || defined (P4_TO_P8)
  int                 l;
#endif
  int                 k;
  int                 qcid, face;
  int                *nonlocal_ranks;
  size_t              zz, position;
  int8_t             *local_status, *quad_status;
  p4est_topidx_t      jt;
  p4est_locidx_t      il, first, second;
  p4est_locidx_t      num_local_nodes, quad_indeps[P4EST_CHILDREN];
  p4est_locidx_t      num_owned_indeps, offset_owned_indeps;
  p4est_locidx_t      num_indep_nodes, dup_indep_nodes, all_face_hangings;
  p4est_locidx_t      num_face_hangings, dup_face_hangings;
  p4est_locidx_t     *local_nodes, *quad_nodes;
  p4est_locidx_t     *new_node_number;
  p4est_tree_t       *tree;
  p4est_nodes_t      *nodes;
  p4est_quadrant_t    c, n, p;
  p4est_quadrant_t   *q, *qpp[3], *r;
  p4est_indep_t      *in;
  sc_array_t         *quadrants;
  sc_array_t         *inda, *faha;
  sc_array_t         *shared_indeps;
  sc_hash_array_t    *indep_nodes;
  sc_hash_array_t    *face_hangings;
#ifndef P4_TO_P8
  p4est_hang2_t      *fh;
#else
  int                 edge, corner;
  p4est_locidx_t      num_face_hangings_end, num_edge_hangings_begin;
  p4est_locidx_t      num_edge_hangings, dup_edge_hangings;
  p8est_hang4_t      *fh;
  p8est_hang2_t      *eh;
  sc_array_t          exist_array;
  sc_array_t         *edha;
  sc_hash_array_t    *edge_hangings;
#endif

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_ASSERT (ghost_layer != NULL);

  P4EST_QUADRANT_INIT (&c);
  P4EST_QUADRANT_INIT (&n);
  P4EST_QUADRANT_INIT (&p);
  qpp[0] = NULL;
  qpp[1] = qpp[2] = &p;

  /* initialize the node structure to return */
  nodes = P4EST_ALLOC (p4est_nodes_t, 1);
  memset (nodes, -1, sizeof (*nodes));
  faha = &nodes->face_hangings;
#ifdef P4_TO_P8
  edha = &nodes->edge_hangings;
#endif
  shared_indeps = &nodes->shared_indeps;
  sc_array_init (shared_indeps, sizeof (sc_recycle_array_t));

  /* compute number of local quadrant corners */
  nodes->num_local_quadrants = p4est->local_num_quadrants;
  num_local_nodes =             /* same type */
    P4EST_CHILDREN * nodes->num_local_quadrants;

  /* Store hanging node status:
   * 0 for independent, 1 for face hanging, 2 for edge hanging.
   */
  local_status = P4EST_ALLOC (int8_t, num_local_nodes);
  memset (local_status, -1, num_local_nodes * sizeof (*local_status));

  /* Store the local node index for each corner of the elements.
   */
  nodes->local_nodes = local_nodes =
    P4EST_ALLOC (p4est_locidx_t, num_local_nodes);
  memset (local_nodes, -1, num_local_nodes * sizeof (*local_nodes));

  indep_nodes = sc_hash_array_new (sizeof (p4est_indep_t),
                                   p4est_node_hash_piggy_fn,
                                   p4est_node_equal_piggy_fn, NULL);
#ifndef P4_TO_P8
  face_hangings = sc_hash_array_new (sizeof (p4est_hang2_t),
                                     p4est_node_hash_piggy_fn,
                                     p4est_node_equal_piggy_fn, NULL);
#else
  face_hangings = sc_hash_array_new (sizeof (p8est_hang4_t),
                                     p4est_node_hash_piggy_fn,
                                     p4est_node_equal_piggy_fn, NULL);
  edge_hangings = sc_hash_array_new (sizeof (p8est_hang2_t),
                                     p4est_node_hash_piggy_fn,
                                     p4est_node_equal_piggy_fn, NULL);
  sc_array_init (&exist_array, sizeof (int));
#endif

  /* This first loop will fill the local_status array with hanging status.
   * It will also collect all independent nodes relevant for the elements.
   */
  num_indep_nodes = dup_indep_nodes = all_face_hangings = 0;
  quad_nodes = local_nodes;
  quad_status = local_status;
  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    tree = p4est_array_index_topidx (p4est->trees, jt);
    quadrants = &tree->quadrants;

    /* determine hanging node status and collect all anchored nodes */
    for (zz = 0; zz < quadrants->elem_count;
         quad_nodes += P4EST_CHILDREN, quad_status += P4EST_CHILDREN, ++zz) {
      qpp[0] = q = sc_array_index (quadrants, zz);
      qcid = p4est_quadrant_child_id (q);
      if (q->level > 0) {
        p4est_quadrant_parent (q, &p);
      }
#ifdef P4EST_DEBUG
      else {
        P4EST_QUADRANT_INIT (&p);
      }
#endif

      /* assign independent node and face hanging node status */
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        if (k == qcid || k == P4EST_CHILDREN - 1 - qcid || q->level == 0) {
          quad_status[k] = 0;
          continue;
        }
#ifndef P4_TO_P8
        if (k == p4est_hanging_corner[qcid][0]) {
          face = p4est_hanging_face[qcid][0];
        }
        else {
          P4EST_ASSERT (k == p4est_hanging_corner[qcid][1]);
          face = p4est_hanging_face[qcid][1];
        }
#else
        face = p8est_child_corner_faces[qcid][k];
        if (face == -1) {
          P4EST_ASSERT (p8est_child_corner_edges[qcid][k] >= 0);
          continue;
        }
#endif
        p4est_quadrant_face_neighbor (&p, face, &n);
        if (p4est_quadrant_exists (p4est, ghost_layer, jt, &n, NULL)) {
          quad_status[k] = 1;
#ifdef P4_TO_P8
          for (l = 0; l < 4; ++l) {
            corner = p8est_face_corners[face][l];
            if (corner != qcid && corner != k) {
              quad_status[corner] = 2;
            }
          }
#endif
          ++all_face_hangings;
        }
        else {
          quad_status[k] = 0;
        }
      }

#ifdef P4_TO_P8
      /* assign edge hanging node status */
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        if (quad_status[k] == -1) {
          edge = p8est_child_corner_edges[qcid][k];
          P4EST_ASSERT (edge >= 0 && edge < 12);
          p8est_quadrant_edge_neighbor (&p, edge, &n);
          quad_status[k] = (int8_t)
            (p4est_quadrant_exists (p4est, ghost_layer, jt, &n,
                                    &exist_array) ? 2 : 0);
        }
      }
#endif

      /* collect all independent nodes related to the element */
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        P4EST_ASSERT (quad_status[k] >= 0 || quad_status[k] <= 2);
        p4est_quadrant_corner_node (qpp[quad_status[k]], k, &n);
        p4est_node_canonicalize (p4est, jt, &n, &c);
        r = sc_hash_array_insert_unique (indep_nodes, &c, &position);
        if (r != NULL) {
          *r = c;
          P4EST_ASSERT (num_indep_nodes == (p4est_locidx_t) position);
          ++num_indep_nodes;
        }
        else {
          ++dup_indep_nodes;
        }
        P4EST_ASSERT ((p4est_locidx_t) position < num_indep_nodes);
        quad_nodes[k] = (p4est_locidx_t) position;
      }
    }
  }
  P4EST_ASSERT (num_indep_nodes + dup_indep_nodes == num_local_nodes);
#ifdef P4_TO_P8
  sc_array_reset (&exist_array);
#endif
  inda = &indep_nodes->a;
  P4EST_ASSERT (num_indep_nodes == (p4est_locidx_t) inda->elem_count);

  /* Reorder independent nodes by their global treeid and z-order index. */
  new_node_number = P4EST_ALLOC (p4est_locidx_t, num_indep_nodes);
  for (il = 0; il < num_indep_nodes; ++il) {
    in = sc_array_index (inda, (size_t) il);
    in->pad8 = 0;               /* shared by 0 other processors so far */
    in->pad16 = (int16_t) (-1);
    in->p.piggy3.local_num = il;
  }
  sc_array_sort (inda, p4est_quadrant_compare_piggy);
  for (il = 0; il < num_indep_nodes; ++il) {
    in = sc_array_index (inda, (size_t) il);
    new_node_number[in->p.piggy3.local_num] = il;
#ifndef P4EST_MPI
    in->p.piggy3.local_num = il;
#endif
  }

  /* Re-synchronize hash array and local nodes */
  indep_nodes->internal_data.user_data = new_node_number;
  sc_hash_foreach (indep_nodes->h, p4est_nodes_foreach);
  indep_nodes->internal_data.user_data = NULL;
  for (il = 0; il < num_local_nodes; ++il) {
    P4EST_ASSERT (local_nodes[il] >= 0 && local_nodes[il] < num_indep_nodes);
    local_nodes[il] = new_node_number[local_nodes[il]];
  }
#ifndef P4EST_MPI
  num_owned_indeps = num_indep_nodes;
  offset_owned_indeps = 0;
#else
  num_owned_indeps = 0;         /* will be computed below */
  offset_owned_indeps = -1;     /* will be computed below */
#endif
  P4EST_FREE (new_node_number);

#ifdef P4EST_MPI
  /* Fill send buffers and number owned nodes. */
  first_size = P4EST_DIM * sizeof (p4est_qcoord_t) + sizeof (p4est_topidx_t);
  first_size = SC_MAX (first_size, sizeof (p4est_locidx_t));
  procs = P4EST_ALLOC_ZERO (int, (size_t) num_procs);
  all_ranges = P4EST_ALLOC (int, twopeerw * num_procs);
  peers = P4EST_ALLOC (p4est_node_peer_t, num_procs);
  sc_array_init (&send_requests, sizeof (MPI_Request));
  for (k = 0; k < num_procs; ++k) {
    peer = peers + k;
    peer->expect_query = peer->expect_reply = false;
    peer->recv_offset = 0;
    sc_array_init (&peer->send_first, first_size);
    sc_array_init (&peer->recv_first, first_size);
    sc_array_init (&peer->send_second, 1);
    sc_array_init (&peer->recv_second, 1);
  }
  first_peer = num_procs;
  last_peer = -1;
  prev = 0;
  for (il = 0; il < num_indep_nodes; ++il) {
    in = sc_array_index (inda, (size_t) il);
    owner = p4est_comm_find_owner (p4est, in->p.which_tree,
                                   (p4est_quadrant_t *) in, prev);
    if (owner != rank) {
      peer = peers + owner;
      xyz = sc_array_push (&peer->send_first);
      xyz[0] = in->x;
      xyz[1] = in->y;
#ifdef P4_TO_P8
      xyz[2] = in->z;
#endif
      ttt = (p4est_topidx_t *) (&xyz[P4EST_DIM]);
      *ttt = in->p.which_tree;
      in->p.piggy1.owner_rank = owner;
      if (first_peer == num_procs) {
        first_peer = owner;
      }
      last_peer = owner;
      ++procs[owner];
    }
    else {
      if (offset_owned_indeps == -1) {
        offset_owned_indeps = il;
      }
      in->p.piggy3.local_num = num_owned_indeps++;
    }
    P4EST_ASSERT (prev <= owner);
    prev = owner;
  }
  if (offset_owned_indeps == -1) {
    P4EST_ASSERT (num_owned_indeps == 0);
    offset_owned_indeps = 0;
  }

  /* Distribute global information about who is sending to who. */
  nwin = sc_ranges_compute (num_procs, procs, rank, first_peer, last_peer,
                            p4est_num_ranges, my_ranges);
#ifdef P4EST_STATS
  sc_ranges_statistics (SC_LP_STATISTICS, p4est->mpicomm, num_procs, procs,
                        rank, p4est_num_ranges, my_ranges);
#endif
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Allgather (my_ranges, twopeerw, MPI_INT,
                            all_ranges, twopeerw, MPI_INT, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);
  }
  P4EST_VERBOSEF ("Peer ranges %d/%d first %d last %d owned %lld/%lld\n",
                  nwin, p4est_num_ranges, first_peer, last_peer,
                  (long long) num_owned_indeps, (long long) num_indep_nodes);

  /* Send queries to the owners of the independent nodes that I share. */
  num_send_queries = num_send_nonzero = local_send_count = 0;
  for (l = 0; l < nwin; ++l) {
    for (k = my_ranges[2 * l]; k <= my_ranges[2 * l + 1]; ++k) {
      peer = peers + k;
      if (k == rank) {
        P4EST_ASSERT (peer->send_first.elem_count == 0);
        continue;
      }
      send_request = sc_array_push (&send_requests);
      this_size = peer->send_first.elem_count * first_size;
      mpiret = MPI_Isend (peer->send_first.array, (int) this_size,
                          MPI_BYTE, k, P4EST_COMM_NODES_QUERY,
                          p4est->mpicomm, send_request);
      SC_CHECK_MPI (mpiret);
      local_send_count += (int) peer->send_first.elem_count;
      ++num_send_queries;
      if (this_size > 0) {
        ++num_send_nonzero;
        peer->expect_reply = true;
      }
    }
  }

  /* Prepare to receive queries */
  num_recv_queries = local_recv_count = 0;
  for (k = 0; k < num_procs; ++k) {
    if (k == rank) {
      continue;
    }
    for (l = 0; l < p4est_num_ranges; ++l) {
      start = all_ranges[k * twopeerw + 2 * l];
      if (start == -1 || start > rank) {
        break;
      }
      if (rank <= all_ranges[k * twopeerw + 2 * l + 1]) {
        peers[k].expect_query = true;
        ++num_recv_queries;
        break;
      }
    }
  }
  P4EST_VERBOSEF ("Node queries send %d nonz %d recv %d\n",
                  num_send_queries, num_send_nonzero, num_recv_queries);

  /* Receive queries and look up the reply information */
  P4EST_QUADRANT_INIT (&inkey);
  inkey.level = P4EST_MAXLEVEL;
  for (l = 0; l < num_recv_queries; ++l) {
    mpiret = MPI_Probe (MPI_ANY_SOURCE, P4EST_COMM_NODES_QUERY,
                        p4est->mpicomm, &probe_status);
    SC_CHECK_MPI (mpiret);
    k = probe_status.MPI_SOURCE;
    peer = peers + k;
    P4EST_ASSERT (k != rank && peer->expect_query);
    mpiret = MPI_Get_count (&probe_status, MPI_BYTE, &byte_count);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (byte_count % first_size == 0);
    elem_count = byte_count / (int) first_size;
    local_recv_count += elem_count;
    sc_array_resize (&peer->recv_first, (size_t) elem_count);
    mpiret = MPI_Recv (peer->recv_first.array, byte_count, MPI_BYTE,
                       k, P4EST_COMM_NODES_QUERY,
                       p4est->mpicomm, &recv_status);
    SC_CHECK_MPI (mpiret);
    peer->expect_query = false;
    for (zz = 0; zz < peer->recv_first.elem_count; ++zz) {
      xyz = sc_array_index (&peer->recv_first, zz);
      inkey.x = xyz[0];
      inkey.y = xyz[1];
#ifdef P4_TO_P8
      inkey.z = xyz[2];
#endif
      ttt = (p4est_topidx_t *) (&xyz[P4EST_DIM]);
      inkey.p.which_tree = *ttt;
      found = sc_hash_array_lookup (indep_nodes, &inkey, &position);
      P4EST_ASSERT (found);
      P4EST_ASSERT (position >= offset_owned_indeps &&
                    position < offset_owned_indeps + num_owned_indeps);
      node_number = (p4est_locidx_t *) xyz;
      *node_number = (p4est_locidx_t) position - offset_owned_indeps;
      in = sc_array_index (inda, position);
      P4EST_ASSERT (p4est_node_equal_piggy_fn (&inkey, in, NULL));
      P4EST_ASSERT (in->pad8 >= 0);
      num_sharers = (size_t) in->pad8;
      P4EST_ASSERT (num_sharers <= shared_indeps->elem_count);
      SC_CHECK_ABORT (num_sharers < (size_t) INT8_MAX,
                      "Max independent node sharer limit exceeded");
      if (num_sharers == shared_indeps->elem_count) {
        nrarr = sc_array_push (shared_indeps);
        sc_recycle_array_init (nrarr, (num_sharers + 1) * sizeof (int));
      }
      else {
        nrarr = sc_array_index (shared_indeps, num_sharers);
      }
      new_sharers = sc_recycle_array_insert (nrarr, &new_position);
      if (num_sharers > 0) {
        P4EST_ASSERT (in->pad16 >= 0);
        old_position = (size_t) in->pad16;
        orarr = sc_array_index (shared_indeps, num_sharers - 1);
        old_sharers = sc_recycle_array_remove (orarr, old_position);
        memcpy (new_sharers, old_sharers, num_sharers * sizeof (int));
      }
      new_sharers[num_sharers] = k;
      SC_CHECK_ABORT (new_position <= (size_t) INT16_MAX,
                      "Max independent node count limit exceeded");
      ++in->pad8;
      in->pad16 = (int16_t) new_position;
    }
  }

  /* Assemble and send reply information.  This is variable size.
   * (p4est_locidx_t)      Node number in this processor's ordering
   * (int8_t)              Number of sharers (not including this processor)
   * num_sharers * (int)   The ranks of all sharers.
   */
  second_size = sizeof (p4est_locidx_t) + sizeof (int8_t);
  for (k = 0; k < num_procs; ++k) {
    peer = peers + k;
    if (peer->recv_first.elem_count == 0) {
      continue;
    }
    for (zz = 0; zz < peer->recv_first.elem_count; ++zz) {
      node_number = sc_array_index (&peer->recv_first, zz);
      in =
        sc_array_index (inda, (size_t) (*node_number + offset_owned_indeps));
      P4EST_ASSERT (p4est_quadrant_is_node ((p4est_quadrant_t *) in, true));
      P4EST_ASSERT (in->pad8 >= 0);
      num_sharers = (size_t) in->pad8;
      P4EST_ASSERT (num_sharers <= shared_indeps->elem_count);
      this_size = second_size + num_sharers * sizeof (int);
      this_base = sc_array_push_bytes (&peer->send_second, this_size);
      *(p4est_locidx_t *) this_base = *node_number;
      *(int8_t *) (this_base + sizeof (p4est_locidx_t)) = in->pad8;
      if (num_sharers > 0) {
        P4EST_ASSERT (in->pad16 >= 0);
        nrarr = sc_array_index (shared_indeps, num_sharers - 1);
        new_sharers = sc_array_index (&nrarr->a, (size_t) in->pad16);
        memcpy (this_base + second_size, new_sharers,
                num_sharers * sizeof (int));
      }
    }
    send_request = sc_array_push (&send_requests);
    mpiret = MPI_Isend (peer->send_second.array,
                        (int) peer->send_second.elem_count,
                        MPI_BYTE, k, P4EST_COMM_NODES_REPLY,
                        p4est->mpicomm, send_request);
    SC_CHECK_MPI (mpiret);
    sc_array_reset (&peer->recv_first);
  }
#endif /* P4EST_MPI */

  /* This second loop will collect and assign all hanging nodes. */
  num_face_hangings = dup_face_hangings = 0;    /* still unknown */
#ifdef P4_TO_P8
  num_edge_hangings = dup_edge_hangings = 0;    /* still unknown */
  num_edge_hangings_begin = num_indep_nodes + all_face_hangings;
#endif
  quad_nodes = local_nodes;
  quad_status = local_status;
  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    tree = p4est_array_index_topidx (p4est->trees, jt);
    quadrants = &tree->quadrants;

    /* collect all face and edge hanging nodes */
    for (zz = 0; zz < quadrants->elem_count;
         quad_nodes += P4EST_CHILDREN, quad_status += P4EST_CHILDREN, ++zz) {
      q = sc_array_index (quadrants, zz);
      qcid = p4est_quadrant_child_id (q);

      /* create hanging nodes and assign related independent nodes */
      memcpy (quad_indeps, quad_nodes, P4EST_CHILDREN * sizeof (*quad_nodes));
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        if (quad_status[k] == 1) {
          P4EST_ASSERT (qcid != k && quad_indeps[qcid] != quad_indeps[k]);
#ifndef P4_TO_P8
          P4EST_ASSERT (k == p4est_hanging_corner[qcid][0] ||
                        k == p4est_hanging_corner[qcid][1]);
#else
          P4EST_ASSERT (p8est_child_corner_faces[qcid][k] >= 0);
#endif
          p4est_quadrant_corner_node (q, k, &n);
          p4est_node_canonicalize (p4est, jt, &n, &c);
          r = sc_hash_array_insert_unique (face_hangings, &c, &position);
          if (r != NULL) {
            *r = c;
            P4EST_ASSERT (num_face_hangings == (p4est_locidx_t) position);
#ifndef P4_TO_P8
            fh = (p4est_hang2_t *) r;
            first = quad_indeps[qcid];
            second = quad_indeps[k];
            if (first < second) {
              fh->p.piggy.depends[0] = first;
              fh->p.piggy.depends[1] = second;
            }
            else {
              fh->p.piggy.depends[0] = second;
              fh->p.piggy.depends[1] = first;
            }
#else
            fh = (p8est_hang4_t *) r;
            fh->p.piggy.depends[0] = quad_indeps[qcid];
            fh->p.piggy.depends[1] = quad_indeps[k];
            fh->p.piggy.depends[2] = -1;
            fh->p.piggy.depends[3] = -1;
            face = p8est_child_corner_faces[qcid][k];
            for (l = 0; l < 4; ++l) {
              corner = p8est_face_corners[face][l];
              if (corner != qcid && corner != k) {
                if (fh->p.piggy.depends[2] == -1) {
                  fh->p.piggy.depends[2] = quad_indeps[corner];
                }
                else {
                  P4EST_ASSERT (fh->p.piggy.depends[3] == -1);
                  fh->p.piggy.depends[3] = quad_indeps[corner];
                }
              }
            }
            qsort (fh->p.piggy.depends,
                   4, sizeof (p4est_locidx_t), p4est_locidx_compare);
#endif
            ++num_face_hangings;
          }
          else {
            ++dup_face_hangings;
          }
          quad_nodes[k] =       /* same type */
            num_indep_nodes + (p4est_locidx_t) position;
        }
#ifdef P4_TO_P8
        else if (quad_status[k] == 2) {
          P4EST_ASSERT (qcid != k && quad_indeps[qcid] != quad_indeps[k]);
          P4EST_ASSERT (p8est_child_corner_edges[qcid][k] >= 0);
          p4est_quadrant_corner_node (q, k, &n);
          p4est_node_canonicalize (p4est, jt, &n, &c);
          r = sc_hash_array_insert_unique (edge_hangings, &c, &position);
          if (r != NULL) {
            *r = c;
            P4EST_ASSERT (num_edge_hangings == (p4est_locidx_t) position);
            eh = (p8est_hang2_t *) r;
            first = quad_indeps[qcid];
            second = quad_indeps[k];
            if (first < second) {
              eh->p.piggy.depends[0] = first;
              eh->p.piggy.depends[1] = second;
            }
            else {
              eh->p.piggy.depends[0] = second;
              eh->p.piggy.depends[1] = first;
            }
            ++num_edge_hangings;
          }
          else {
            ++dup_edge_hangings;
          }
          quad_nodes[k] =       /* same type */
            num_edge_hangings_begin + (p4est_locidx_t) position;
        }
#endif
      }
    }
  }
  P4EST_ASSERT (num_face_hangings + dup_face_hangings == all_face_hangings);
  P4EST_FREE (local_status);
  sc_hash_array_rip (face_hangings, faha);
  P4EST_ASSERT (num_face_hangings == (p4est_locidx_t) faha->elem_count);
#ifdef P4_TO_P8
  sc_hash_array_rip (edge_hangings, edha);
  P4EST_ASSERT (num_edge_hangings == (p4est_locidx_t) edha->elem_count);

  /* Correct the offsets of edge hanging nodes */
  num_face_hangings_end = num_indep_nodes + num_face_hangings;
  for (il = 0; il < num_local_nodes; ++il) {
    if (local_nodes[il] >= num_edge_hangings_begin) {
      local_nodes[il] -= dup_face_hangings;
      P4EST_ASSERT (local_nodes[il] >= num_face_hangings_end);
    }
    else {
      P4EST_ASSERT (local_nodes[il] >= 0 &&
                    local_nodes[il] < num_face_hangings_end);
    }
  }
#endif

  /* Allocate remaining output data structures */
  nodes->num_owned_indeps = num_owned_indeps;
  nodes->offset_owned_indeps = offset_owned_indeps;
  sc_hash_array_rip (indep_nodes, inda = &nodes->indep_nodes);
  nonlocal_ranks = nodes->nonlocal_ranks =
    P4EST_ALLOC (int, num_indep_nodes - num_owned_indeps + 1);
  nonlocal_ranks[num_indep_nodes - num_owned_indeps] = -1;
  nodes->global_owned_indeps = P4EST_ALLOC (p4est_locidx_t, num_procs);
  nodes->global_owned_indeps[rank] = num_owned_indeps;
  indep_nodes = NULL;

#ifdef P4EST_MPI
  /* Receive the replies. */
  for (l = 0; l < num_send_nonzero; ++l) {
    mpiret = MPI_Probe (MPI_ANY_SOURCE, P4EST_COMM_NODES_REPLY,
                        p4est->mpicomm, &probe_status);
    SC_CHECK_MPI (mpiret);
    k = probe_status.MPI_SOURCE;
    peer = peers + k;
    P4EST_ASSERT (k != rank && peer->expect_reply);
    mpiret = MPI_Get_count (&probe_status, MPI_BYTE, &byte_count);
    SC_CHECK_MPI (mpiret);
    sc_array_resize (&peer->recv_second, byte_count);
    mpiret = MPI_Recv (peer->recv_second.array, byte_count, MPI_BYTE,
                       k, P4EST_COMM_NODES_REPLY,
                       p4est->mpicomm, &recv_status);
    SC_CHECK_MPI (mpiret);
    peer->expect_reply = false;
  }

  /* Convert the receive buffers into the output data structures. */
  for (il = 0; il < num_indep_nodes; ++il) {
    if (il == offset_owned_indeps) {
      il += num_owned_indeps;
      if (il == num_indep_nodes) {
        break;
      }
    }
    in = sc_array_index (inda, (size_t) il);
    k = in->p.piggy1.owner_rank;
    P4EST_ASSERT (k >= 0 && k != rank && k < num_procs);
    *nonlocal_ranks++ = k;
    peer = peers + k;
    P4EST_ASSERT (peer->recv_offset + second_size <=
                  peer->recv_second.elem_count);
    this_base = sc_array_index (&peer->recv_second, peer->recv_offset);
    in->p.piggy3.local_num = *(p4est_locidx_t *) this_base;
    num_sharers = (size_t) *(int8_t *) (this_base + sizeof (p4est_locidx_t));
    P4EST_ASSERT (num_sharers > 0);
    this_size = second_size + num_sharers * sizeof (int);
    P4EST_ASSERT (peer->recv_offset + this_size <=
                  peer->recv_second.elem_count);
    if (shared_indeps->elem_count < num_sharers) {
      position = shared_indeps->elem_count;
      sc_array_resize (shared_indeps, num_sharers);
      for (zz = position; zz < num_sharers; ++zz) {
        nrarr = sc_array_index (shared_indeps, zz);
        sc_recycle_array_init (nrarr, (zz + 1) * sizeof (int));
      }
    }
    else {
      nrarr = sc_array_index (shared_indeps, num_sharers - 1);
    }
    new_sharers = sc_recycle_array_insert (nrarr, &new_position);
    memcpy (new_sharers, this_base + second_size, num_sharers * sizeof (int));
    for (zz = 0; zz < num_sharers; ++zz) {
      if (new_sharers[zz] == rank) {
        new_sharers[zz] = k;
        break;
      }
    }
    P4EST_ASSERT (zz < num_sharers);
    SC_CHECK_ABORT (new_position <= (size_t) INT16_MAX,
                    "Max independent node count limit exceeded");
    in->pad8 = (int8_t) num_sharers;
    in->pad16 = (int16_t) new_position;
    peer->recv_offset += this_size;
  }

  /* Wait and close all send requests. */
  if (send_requests.elem_count > 0) {
    mpiret = MPI_Waitall ((int) send_requests.elem_count,
                          (MPI_Request *) send_requests.array,
                          MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }

  /* Clean up allocated communications memory. */
  sc_array_reset (&send_requests);
  for (k = 0; k < num_procs; ++k) {
    peer = peers + k;
    P4EST_ASSERT (peer->recv_offset == peer->recv_second.elem_count);
    sc_array_reset (&peer->send_first);
    /* peer->recv_first has been reset above */
    sc_array_reset (&peer->send_second);
    sc_array_reset (&peer->recv_second);
  }
  P4EST_FREE (peers);
  P4EST_FREE (all_ranges);
  P4EST_FREE (procs);

  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Allgather (&num_owned_indeps, 1, P4EST_MPI_LOCIDX,
                            nodes->global_owned_indeps, 1, P4EST_MPI_LOCIDX,
                            p4est->mpicomm);
    SC_CHECK_MPI (mpiret);
  }
#endif /* P4EST_MPI */

  /* Print some statistics and clean up. */
  P4EST_VERBOSEF ("Collected %lld independent nodes with %lld duplicates\n",
                  (long long) num_indep_nodes, (long long) dup_indep_nodes);
  P4EST_VERBOSEF ("Collected %lld face hangings with %lld duplicates\n",
                  (long long) num_face_hangings,
                  (long long) dup_face_hangings);
#ifdef P4_TO_P8
  P4EST_VERBOSEF ("Collected %lld edge hangings with %lld duplicates\n",
                  (long long) num_edge_hangings,
                  (long long) dup_edge_hangings);
#endif
#ifdef P4EST_MPI
  P4EST_VERBOSEF ("Owned nodes %lld/%lld max shared count %llu\n",
                  (long long) num_owned_indeps,
                  (long long) num_indep_nodes,
                  (unsigned long long) shared_indeps->elem_count);
#endif

  P4EST_ASSERT (*nonlocal_ranks == -1);
  P4EST_ASSERT (p4est_nodes_is_valid (p4est, nodes));

  return nodes;
}

void
p4est_nodes_destroy (p4est_nodes_t * nodes)
{
  size_t              zz;
  sc_recycle_array_t *rarr;

  sc_array_reset (&nodes->indep_nodes);
  sc_array_reset (&nodes->face_hangings);
#ifdef P4_TO_P8
  sc_array_reset (&nodes->edge_hangings);
#endif
  P4EST_FREE (nodes->local_nodes);

  for (zz = 0; zz < nodes->shared_indeps.elem_count; ++zz) {
    rarr = sc_array_index (&nodes->shared_indeps, zz);
    sc_recycle_array_reset (rarr);
  }
  sc_array_reset (&nodes->shared_indeps);
  P4EST_FREE (nodes->nonlocal_ranks);
  P4EST_FREE (nodes->global_owned_indeps);

  P4EST_FREE (nodes);
}

bool
p4est_nodes_is_valid (p4est_t * p4est, p4est_nodes_t * nodes)
{
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  int                 k, prev, owner, pshare, ocount;
  int                *sharers, *sorted;
  bool                failed;
  size_t              zz, position, sharez, num_sharers, max_sharers;
  p4est_topidx_t      otree, ntree;
  p4est_locidx_t      il, num_indep_nodes, local_num;
  p4est_locidx_t      num_owned_indeps, offset_owned_indeps, end_owned_indeps;
  p4est_indep_t      *in;
  sc_recycle_array_t *rarr;

  failed = false;
  sorted = P4EST_ALLOC (int, INT8_MAX);

  max_sharers = nodes->shared_indeps.elem_count;
  for (zz = 0; zz < max_sharers; ++zz) {
    rarr = sc_array_index (&nodes->shared_indeps, zz);
    num_sharers = zz + 1;
    for (position = 0; position < rarr->a.elem_count; ++position) {
      sharers = sc_array_index (&rarr->a, position);
      memcpy (sorted, sharers, num_sharers * sizeof (int));
      qsort (sorted, num_sharers, sizeof (int), sc_int_compare);
      prev = -1;
      for (sharez = 0; sharez < num_sharers; ++sharez) {
        k = sorted[sharez];
        if (prev >= k || k == rank || k >= num_procs) {
          P4EST_NOTICE ("p4est nodes invalid sharers 1\n");
          failed = true;
          goto failtest;
        }
        prev = k;
      }
    }
  }

  num_indep_nodes = (p4est_locidx_t) nodes->indep_nodes.elem_count;
  num_owned_indeps = nodes->num_owned_indeps;
  offset_owned_indeps = nodes->offset_owned_indeps;
  end_owned_indeps = offset_owned_indeps + num_owned_indeps;
  k = prev = 0;
  otree = 0;
  for (il = 0; il < num_indep_nodes; ++il) {
    in = sc_array_index (&nodes->indep_nodes, (size_t) il);
    ntree = in->p.piggy3.which_tree;
    local_num = in->p.piggy3.local_num;
    if (ntree < otree || ntree >= p4est->connectivity->num_trees) {
      P4EST_NOTICE ("p4est nodes invalid tree\n");
      failed = true;
      goto failtest;
    }
    otree = ntree;
    if (in->pad8 < 0 || (num_sharers = (size_t) in->pad8) > max_sharers) {
      P4EST_NOTICE ("p4est nodes invalid sharer count\n");
      failed = true;
      goto failtest;
    }
    if (il < offset_owned_indeps || il >= end_owned_indeps) {
      owner = nodes->nonlocal_ranks[k++];
      if (owner < prev || owner == rank || owner >= num_procs) {
        P4EST_NOTICE ("p4est nodes invalid owner\n");
        failed = true;
        goto failtest;
      }
      if (local_num < 0 || local_num >= nodes->global_owned_indeps[owner]) {
        P4EST_NOTICE ("p4est nodes invalid non-owned index\n");
        failed = true;
        goto failtest;
      }
      prev = owner;
      if (num_sharers < 1) {
        P4EST_NOTICE ("p4est nodes invalid non-owned sharing 1\n");
        failed = true;
        goto failtest;
      }
    }
    else {
      if (local_num != il - offset_owned_indeps) {
        P4EST_NOTICE ("p4est nodes invalid owned index\n");
        failed = true;
        goto failtest;
      }
      owner = prev = rank;
    }
    ocount = 0;
    if (num_sharers > 0) {
      rarr = sc_array_index (&nodes->shared_indeps, num_sharers - 1);
      if (in->pad16 < 0 ||
          (position = (size_t) in->pad16) >= rarr->a.elem_count) {
        P4EST_NOTICE ("p4est nodes invalid sharer position\n");
        failed = true;
        goto failtest;
      }
      sharers = sc_array_index (&rarr->a, position);
      memcpy (sorted, sharers, num_sharers * sizeof (int));
      qsort (sorted, num_sharers, sizeof (int), sc_int_compare);
      pshare = -1;
      for (zz = 0; zz < num_sharers; ++zz) {
        if (sorted[zz] <= pshare || sorted[zz] == rank) {
          P4EST_NOTICE ("p4est nodes invalid sharers 2\n");
          failed = true;
          goto failtest;
        }
        if (sorted[zz] == owner) {
          ++ocount;
        }
        pshare = sorted[zz];
      }
    }
    if (owner != rank && ocount != 1) {
      P4EST_NOTICE ("p4est nodes invalid non-owned sharing 2\n");
      failed = true;
      goto failtest;
    }
  }
  if (nodes->nonlocal_ranks [k] != -1) {
    P4EST_NOTICE ("p4est nodes invalid safeguard\n");
    failed = true;
    goto failtest;
  }

  /* TODO: Test hanging nodes and local corners. */

failtest:
  P4EST_FREE (sorted);

  return !p4est_comm_sync_flag (p4est, failed, MPI_BOR);
}

p4est_neighborhood_t *
p4est_neighborhood_new (p4est_t * p4est)
{
  bool                success;
  p4est_topidx_t      local_num_trees, flt, nt;
  p4est_locidx_t      local_num_quadrants, lsum;
  p4est_tree_t       *tree;
  p4est_nodes_t      *nodes;
  p4est_neighborhood_t *nhood;
  sc_array_t          ghost_layer;

  P4EST_ASSERT (p4est_is_valid (p4est));

  sc_array_init (&ghost_layer, sizeof (p4est_quadrant_t));
  success = p4est_build_ghost_layer (p4est, &ghost_layer);
  P4EST_ASSERT (success);

  nodes = p4est_nodes_new (p4est, &ghost_layer);
  p4est_nodes_destroy (nodes);
  sc_array_reset (&ghost_layer);

  if (p4est->first_local_tree < 0) {
    flt = 0;
    local_num_trees = 0;
  }
  else {
    flt = p4est->first_local_tree;
    local_num_trees = p4est->last_local_tree - flt + 1; /* type ok */
  }
  local_num_quadrants = p4est->local_num_quadrants;

  nhood = P4EST_ALLOC (p4est_neighborhood_t, 1);
  nhood->cumulative_count = P4EST_ALLOC (p4est_locidx_t, local_num_trees + 1);
  nhood->element_offsets =
    P4EST_ALLOC (p4est_locidx_t, local_num_quadrants + 1);
  nhood->local_neighbors = sc_array_new (sizeof (p4est_locidx_t));

  lsum = 0;
  for (nt = 0; nt < local_num_trees; ++nt) {
    nhood->cumulative_count[nt] = lsum;
    tree = p4est_array_index_topidx (p4est->trees, flt + nt);   /* type ok */
    lsum += (p4est_locidx_t) tree->quadrants.elem_count;        /* type ok */
  }
  P4EST_ASSERT (lsum == local_num_quadrants);
  nhood->cumulative_count[nt] = lsum;

  return nhood;
}

void
p4est_neighborhood_destroy (p4est_neighborhood_t * nhood)
{
  P4EST_FREE (nhood->cumulative_count);
  P4EST_FREE (nhood->element_offsets);
  sc_array_destroy (nhood->local_neighbors);

  P4EST_FREE (nhood);
}

/* EOF p4est_mesh.c */
