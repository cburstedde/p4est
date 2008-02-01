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

#include <p4est_base.h>
#include <p4est_mesh.h>
#include <p4est_algorithms.h>
#include <p4est_communication.h>

void
p4est_order_local_vertices (p4est_t * p4est,
                            int32_t * num_uniq_local_vertices,
                            int32_t * quadrant_to_local_vertex)
{
  const int           rank = p4est->mpirank;
  int                 qcid, transform;
  int                 zcorner, neighbor_node;
  int                 face, corner, nnum, rlev, tree_corner;
  int                 face_contact[4];
  int                 quad_contact[4];
  int32_t             Ntotal = 0;
  int32_t             Ncells = p4est->local_num_quadrants;
  int32_t             first_local_tree = p4est->first_local_tree;
  int32_t             last_local_tree = p4est->last_local_tree;
  p4est_array_t      *trees = p4est->trees;
  p4est_connectivity_t *conn = p4est->connectivity;
  int32_t             num_trees = conn->num_trees;
  int32_t             i, j;
  int32_t             vertex_num;
  int32_t             rh = (1 << P4EST_MAXLEVEL);
  int32_t             num_quadrants;
  int32_t             lqid;
  int32_t             neighbor_proc, neighbor_tree, ctree;
  int32_t            *tree_offset;
  ssize_t             lnid;
  p4est_tree_t       *tree, *ntree;
  p4est_array_t      *quadrants, corner_info;
  p4est_quadrant_t    neighbor, cneighbor;
  p4est_quadrant_t   *q;
  p4est_corner_info_t *ci;

  P4EST_ASSERT (p4est_is_valid (p4est));

  P4EST_QUADRANT_INIT (&neighbor);
  P4EST_QUADRANT_INIT (&cneighbor);

  p4est_array_init (&corner_info, sizeof (p4est_corner_info_t));

  /* figure out the offset of each tree into the local element id */
  tree_offset = P4EST_ALLOC_ZERO (int32_t, num_trees);
  P4EST_CHECK_ALLOC (tree_offset);
  if (first_local_tree >= 0) {
    tree_offset[first_local_tree] = 0;
    for (j = first_local_tree; j < last_local_tree; ++j) {
      tree = p4est_array_index (trees, j);
      tree_offset[j + 1] = tree_offset[j] + tree->quadrants.elem_count;
    }
  }
  else {
    P4EST_ASSERT (first_local_tree == -1 && last_local_tree == -2);
  }

  /* Initialize vertex list to all -1.  This way we know which values
   * get set because legitimate values are >= 0.
   */
  for (i = 0; i < 4 * Ncells; ++i) {
    quadrant_to_local_vertex[i] = -1;
  }

  /* loop over all local trees to generate the connetivity list */
  for (j = first_local_tree, vertex_num = 0, lqid = 0;
       j <= last_local_tree; ++j) {
    for (face = 0; face < 4; ++face) {
      face_contact[face] = (conn->tree_to_tree[4 * j + face] != j ||
                            conn->tree_to_face[4 * j + face] != face);
    }
    tree = p4est_array_index (p4est->trees, j);
    quadrants = &tree->quadrants;
    num_quadrants = quadrants->elem_count;

    /* Find the neighbors of each quadrant */
    for (i = 0; i < num_quadrants; ++i, ++lqid) {
      /* this quadrant may be on the boundary with a range of processors */
      q = p4est_array_index (quadrants, i);

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

              if (p4est_quadrant_is_inside (&neighbor)) {
                /* neighbor is in the same tree */

                neighbor_proc = p4est_comm_find_owner (p4est, j, &neighbor,
                                                       rank);

                /* Neighbor is remote so we don't number its node */
                if (neighbor_proc != rank)
                  continue;

                lnid = p4est_array_bsearch (quadrants, &neighbor,
                                            p4est_quadrant_compare);
                if (lnid != -1) {
                  lnid += tree_offset[j];
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
                  p4est_find_corner_info (conn, j, tree_corner, &corner_info);
                  for (ctree = 0; ctree < corner_info.elem_count; ++ctree) {
                    ci = p4est_array_index (&corner_info, ctree);
                    neighbor_tree = ci->ntree;
                    zcorner = p4est_corner_to_zorder[ci->ncorner];
                    cneighbor = neighbor;
                    p4est_quadrant_corner (&cneighbor, zcorner, 1);

                    neighbor_proc = p4est_comm_find_owner (p4est,
                                                           neighbor_tree,
                                                           &cneighbor, rank);

                    /* Neighbor is remote so we don't number its node */
                    if (neighbor_proc != rank)
                      continue;

                    ntree = p4est_array_index (trees, neighbor_tree);

                    lnid = p4est_array_bsearch (&ntree->quadrants, &cneighbor,
                                                p4est_quadrant_compare);
                    if (lnid != -1) {
                      lnid += tree_offset[neighbor_tree];
                      neighbor_node = zcorner;
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
                      neighbor_tree = conn->tree_to_tree[4 * j + face];
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
                  transform = p4est_find_face_transform (conn, j, face);
                  p4est_quadrant_translate (&neighbor, face);
                  p4est_quadrant_transform (&neighbor, &cneighbor, transform);

                  neighbor_proc = p4est_comm_find_owner (p4est,
                                                         neighbor_tree,
                                                         &cneighbor, rank);
                  /* Neighbor is remote so we don't number its node */
                  if (neighbor_proc != rank)
                    continue;

                  ntree = p4est_array_index (trees, neighbor_tree);

                  lnid = p4est_array_bsearch (&ntree->quadrants, &cneighbor,
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
  p4est_array_reset (&corner_info);

  *num_uniq_local_vertices = Ntotal;
}

void
p4est_possible_node_neighbor (const p4est_quadrant_t * q, int node,
                              int nnum, int neighbor_rlev,
                              p4est_quadrant_t * neighbor, int *neighbor_node)
{
  int                 nnode;
  const int           nlevel = q->level + neighbor_rlev;
  const int32_t       qh = (1 << (P4EST_MAXLEVEL - q->level));
  const int32_t       nh = (1 << (P4EST_MAXLEVEL - nlevel));
  const int32_t       qx = q->x;
  const int32_t       qy = q->y;
  int32_t             cornerx, cornery;
  p4est_quadrant_t    n;

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
    P4EST_ASSERT_NOT_REACHED ();
    break;
  }

#ifdef P4EST_HAVE_DEBUG
  /* Check to see if it is possible to construct the neighbor */
  int                 qcid = p4est_quadrant_child_id (q);
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
    P4EST_ASSERT_NOT_REACHED ();
    break;
  }

  *neighbor = n;
  *neighbor_node = nnode;

  P4EST_ASSERT (p4est_quadrant_is_extended (neighbor));
}

/* EOF p4est_mesh.c */
