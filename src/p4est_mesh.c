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
 *                      Its user data will be filled with owning tree id.
 * \return              Returns the owning tree id.
 */
static              p4est_topidx_t
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

  p4est_node_clamp_inside (n, c);
  c->p.which_tree = lowest = treeid;

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
    return lowest;
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
      c->p.which_tree = lowest = ntreeid;
    }
    else {
      P4EST_ASSERT (lowest == ntreeid);
      p4est_node_clamp_inside (&o, &tmpq);
      if (p4est_quadrant_compare (&tmpq, c) < 0) {
        /* same tree (periodic) and the new position is lower than the old */
        *c = tmpq;
        c->p.which_tree = lowest;
      }
    }
  }
  P4EST_ASSERT (ntreeid >= 0);
  if (contacts == 1) {
    return lowest;
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
        c->p.which_tree = lowest = ntreeid;
      }
      else {
        P4EST_ASSERT (lowest == ntreeid);
        p4est_node_clamp_inside (&o, &tmpq);
        if (p4est_quadrant_compare (&tmpq, c) < 0) {
          /* same tree (periodic) and the new position is lower than the old */
          *c = tmpq;
          c->p.which_tree = lowest;
        }
      }
    }
  }
  sc_array_reset (eta);
  eta = NULL;
  et = NULL;
  if (contacts == 2) {
    return lowest;
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
      p4est_quadrant_transform_corner (&o, ct->ncorner, false);
      if (ntreeid < lowest) {
        p4est_node_clamp_inside (&o, c);
        c->p.which_tree = lowest = ntreeid;
      }
      else {
        P4EST_ASSERT (lowest == ntreeid);
        p4est_node_clamp_inside (&o, &tmpq);
        if (p4est_quadrant_compare (&tmpq, c) < 0) {
          /* same tree (periodic) and the new position is lower than the old */
          *c = tmpq;
          c->p.which_tree = lowest;
        }
      }
    }
  }
  sc_array_reset (cta);

  return lowest;
}

void
p4est_collect_nodes (p4est_t * p4est, sc_array_t * ghost_layer)
{
  int                 k;
  int                 qcid;
  int                 face;
  size_t              zz;
  p4est_qcoord_t      hanging[P4EST_CHILDREN];
  p4est_topidx_t      jt;
  p4est_locidx_t      num_local_nodes;
  p4est_locidx_t     *local_nodes, *quad_nodes;
  p4est_tree_t       *tree;
  p4est_quadrant_t    p, n, c;
  p4est_quadrant_t   *q, *qp[2];
  sc_array_t         *quadrants;
  sc_array_t          exist_array;
#ifdef P4_TO_P8
  int                 l, edge, corner;
#ifdef P4EST_DEBUG
  int                 check;
#endif
#endif

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_ASSERT (ghost_layer != NULL);

  sc_array_init (&exist_array, sizeof (int));
  P4EST_QUADRANT_INIT (&p);
  P4EST_QUADRANT_INIT (&n);
  P4EST_QUADRANT_INIT (&c);
  qp[0] = NULL;
  qp[1] = &p;
  num_local_nodes = P4EST_CHILDREN * p4est->local_num_quadrants;
  local_nodes = P4EST_ALLOC (p4est_locidx_t, num_local_nodes);

  quad_nodes = local_nodes;
  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    tree = p4est_array_index_topidx (p4est->trees, jt);
    quadrants = &tree->quadrants;

    for (zz = 0; zz < quadrants->elem_count;
         quad_nodes += P4EST_CHILDREN, ++zz) {
      qp[0] = q = sc_array_index (quadrants, zz);
      qcid = p4est_quadrant_child_id (q);
      if (q->level > 0) {
        p4est_quadrant_parent (q, &p);
      }
#ifdef P4EST_DEBUG
      else {
        P4EST_QUADRANT_INIT (&p);
      }
#endif

      /* establish hanging node status of all corners */
      hanging[0] = hanging[1] = hanging[2] = hanging[3] = -1;
#ifdef P4_TO_P8
      hanging[4] = hanging[5] = hanging[6] = hanging[7] = -1;
#endif
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        if (k == qcid || k == P4EST_CHILDREN - 1 - qcid || q->level == 0) {
          hanging[k] = 0;
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
          hanging[k] = 1;
#ifdef P4_TO_P8
#ifdef P4EST_DEBUG
          check = 0;
#endif
          for (l = 0; l < 4; ++l) {
            corner = p8est_face_corners[face][l];
            if (corner != qcid && corner != k) {
              hanging[corner] = 1;
#ifdef P4EST_DEBUG
              ++check;
#endif
            }
          }
          P4EST_ASSERT (check == 2);
#endif
        }
        else {
          hanging[k] = 0;
        }
      }
#ifdef P4_TO_P8
#ifdef P4EST_DEBUG
      check = 0;
#endif
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        if (hanging[k] == -1) {
          edge = p8est_child_corner_edges[qcid][k];
          P4EST_ASSERT (edge >= 0 && edge < 12);
          p8est_quadrant_edge_neighbor (&p, edge, &n);
          hanging[k] = (int) p4est_quadrant_exists
            (p4est, ghost_layer, jt, &n, &exist_array);
#ifdef P4EST_DEBUG
          ++check;
#endif
        }
      }
      P4EST_ASSERT (check <= 3);
#endif

      /* assign all independent and hanging nodes to the element */
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        P4EST_ASSERT (hanging[k] == 0 || hanging[k] == 1);
        p4est_quadrant_corner_node (qp[hanging[k]], k, &n);
        p4est_node_canonicalize (p4est, jt, &n, &c);

      }
    }
  }

  P4EST_FREE (local_nodes);
  sc_array_reset (&exist_array);
}

p4est_neighborhood_t *
p4est_neighborhood_new (p4est_t * p4est)
{
  bool                success;
  p4est_topidx_t      local_num_trees, flt, nt;
  p4est_locidx_t      local_num_quadrants, lsum;
  p4est_tree_t       *tree;
  p4est_neighborhood_t *nhood;
  sc_array_t          ghost_layer;

  P4EST_ASSERT (p4est_is_valid (p4est));

  sc_array_init (&ghost_layer, sizeof (p4est_quadrant_t));
  success = p4est_build_ghost_layer (p4est, &ghost_layer);
  P4EST_ASSERT (success);

  p4est_collect_nodes (p4est, &ghost_layer);
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
