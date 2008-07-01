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

#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_mesh.h>

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

void
p4est_possible_node_neighbor (const p4est_quadrant_t * q, int node,
                              int nnum, int neighbor_rlev,
                              p4est_quadrant_t * neighbor, int *neighbor_node)
{
  int                 nnode;
  const int           nlevel = (int) q->level + neighbor_rlev;
  const p4est_qcoord_t qh =
    (p4est_qcoord_t) (1 << (P4EST_MAXLEVEL - q->level));
  const p4est_qcoord_t nh = (p4est_qcoord_t) (1 << (P4EST_MAXLEVEL - nlevel));
  const p4est_qcoord_t qx = q->x;
  const p4est_qcoord_t qy = q->y;
  p4est_qcoord_t      cornerx, cornery;
  p4est_quadrant_t    n;
#ifdef P4EST_HAVE_DEBUG
  int                 qcid;
#endif

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (-1 <= neighbor_rlev && neighbor_rlev <= 1);
  P4EST_ASSERT (0 <= nlevel && nlevel <= P4EST_MAXLEVEL);

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

#ifdef P4EST_HAVE_DEBUG
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

/* EOF p4est_mesh.c */
