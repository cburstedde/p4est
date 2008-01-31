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
  but WITHOUT ANY WARRANTY; wconst int8_t        p4est_corner_to_zorder[5] = { 0, 1, 3, 2, 4 };ithout even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <p4est_mesh.h>
#include <p4est_base.h>
#include <p4est_algorithms.h>
#include <p4est_communication.h>

static const int8_t fully_owned_flag = 0x01;
static const int8_t any_face_flag = 0x02;

void
p4est_order_local_vertices (p4est_t * p4est,
                            int32_t *
                            num_uniq_local_vertices,
                            int32_t * quadrant_to_local_vertex)
{
  const int           rank = p4est->mpirank;
  int32_t             Ntotal = 0;
  int32_t             Ncells = p4est->local_num_quadrants;
  int32_t             first_local_tree = p4est->first_local_tree;
  int32_t             last_local_tree = p4est->last_local_tree;
  p4est_array_t      *trees = p4est->trees;
  p4est_connectivity_t *conn = p4est->connectivity;
  int32_t             num_trees = conn->num_trees;
  int32_t             next_tree;
  int32_t             i, j;
  int32_t             vertex_num;
  int32_t             qh;
  int32_t             rh = (1 << P4EST_MAXLEVEL);
  int32_t             transform;
  int32_t             num_quadrants;
  int32_t             any_face, face_contact[4];
  int32_t             quad_contact[4];
  int32_t             found_neighbors;
  int32_t             lqid, lnid;
  int32_t             neighbor_proc, neighbor_tree, ctree;
  int32_t            *tree_offset;
  int8_t             *tree_flags;
  int8_t              face, corner, nnum, rlev, tree_corner, zcorner;
  int8_t              qcid, neighbor_node;
  p4est_tree_t       *tree, *ntree;
  p4est_array_t      *quadrants, corner_info;
  p4est_quadrant_t    neighbor;
  p4est_quadrant_t   *q;
  p4est_quadrant_t    mylow, nextlow;
  p4est_corner_info_t *ci;

  P4EST_ASSERT (p4est_is_valid (p4est));

  P4EST_QUADRANT_INIT (&mylow);
  P4EST_QUADRANT_INIT (&nextlow);
  P4EST_QUADRANT_INIT (&neighbor);

  P4EST_ASSERT (p4est->global_first_indices[3 * rank + 0] ==
                first_local_tree);
  mylow.x = p4est->global_first_indices[3 * rank + 1];
  mylow.y = p4est->global_first_indices[3 * rank + 2];
  mylow.level = P4EST_MAXLEVEL;
  next_tree = p4est->global_first_indices[3 * (rank + 1) + 0];
  P4EST_ASSERT (next_tree == last_local_tree
                || next_tree == last_local_tree + 1);
  nextlow.x = p4est->global_first_indices[3 * (rank + 1) + 1];
  nextlow.y = p4est->global_first_indices[3 * (rank + 1) + 2];
  nextlow.level = P4EST_MAXLEVEL;

  p4est_array_init (&corner_info, sizeof (p4est_corner_info_t));

  /* figure out the offset of each tree into the local element id */
  tree_offset = P4EST_ALLOC_ZERO (int32_t, num_trees);
  P4EST_CHECK_ALLOC (tree_offset);
  if (num_trees > 0) {
    tree_offset[first_local_tree] = 0;
    for (j = first_local_tree; j < last_local_tree; ++j) {
      tree = p4est_array_index (trees, j);
      tree_offset[j + 1] = tree_offset[j] + tree->quadrants.elem_count;
    }
  }

  /* tree status flags (max 8 per tree) */
  tree_flags = P4EST_ALLOC (int8_t, num_trees);
  P4EST_CHECK_ALLOC (tree_flags);
  for (i = 0; i < num_trees; ++i) {
    tree_flags[i] = 0x00;
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
    any_face = 0;
    for (face = 0; face < 4; ++face) {
      face_contact[face] = (conn->tree_to_tree[4 * j + face] != j);
      any_face = any_face || face_contact[face];
    }
    if (any_face) {
      tree_flags[j] |= any_face_flag;
    }
    tree = p4est_array_index (p4est->trees, j);
    quadrants = &tree->quadrants;
    num_quadrants = quadrants->elem_count;

    /* Find the neighbors of each quadrant */
    for (i = 0; i < num_quadrants; ++i, ++lqid) {
      /* this quadrant may be on the boundary with a range of processors */
      q = p4est_array_index (quadrants, i);
      qh = (1 << (P4EST_MAXLEVEL - q->level));

      /* loop over the corners of the quadrant */
      for (corner = 0; corner < 4; ++corner) {

        /* Check to see if we have a new vertex */
        if (quadrant_to_local_vertex[lqid * 4 + corner] == -1) {
          quadrant_to_local_vertex[lqid * 4 + corner] = vertex_num;

          /* loop over the possible neighbors and set the new vertex */
          for (nnum = 0; nnum < 4; ++nnum) {
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
              p4est_possible_node_neigbors (q, corner, nnum, rlev,
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
                /* the neightbor is in a neighboring tree or multiple
                 * if it is a neighbor across the corner of the tree
                 */
              }
            }
          }
        }
        ++vertex_num;
      }
    }
  }

  Ntotal = vertex_num;
  P4EST_FREE (tree_offset);
  P4EST_FREE (tree_flags);

#if 0
  /* ----- Cut Here ---- */

  Ntotal = p4est->local_num_quadrants * 4;

  for (i = 0; i < Ncells; ++i) {
    quadrant_to_local_vertex[4 * i + 0] = 4 * i + 0;
    quadrant_to_local_vertex[4 * i + 1] = 4 * i + 1;
    quadrant_to_local_vertex[4 * i + 2] = 4 * i + 2;
    quadrant_to_local_vertex[4 * i + 3] = 4 * i + 3;
  }

  /* ----- To Cut Here ---- */
#endif

  *num_uniq_local_vertices = Ntotal;
}

void
p4est_possible_node_neigbors (p4est_quadrant_t * q,
                              int32_t node,
                              int8_t nnum,
                              int8_t neighbor_rlev,
                              p4est_quadrant_t * neighbor,
                              int8_t * neighbor_node)
{
  p4est_quadrant_t    n;
  int8_t              nnode;
  int32_t             qh = (1 << (P4EST_MAXLEVEL - q->level));
  int32_t             nh =
    (1 << (P4EST_MAXLEVEL - (q->level + neighbor_rlev)));
  int32_t             qx = q->x;
  int32_t             qy = q->y;
  int32_t             cornerx, cornery;

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
  }

#ifdef P4EST_HAVE_DEBUG
  /* Check to see if it is possible to construct the neighbor */
  int8_t              qcid = p4est_quadrant_child_id (q);
  P4EST_ASSERT (neighbor_rlev >= 0 || qcid == node);
#endif

  nnode = (int8_t) (3 - nnum);
  n.level = (int8_t) (q->level + neighbor_rlev);
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
  }

  *neighbor = n;
  *neighbor_node = nnode;

  P4EST_ASSERT (p4est_quadrant_is_extended (neighbor));
}

/* EOF p4est_mesh.c */
