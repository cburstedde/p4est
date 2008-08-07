/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2008 Carsten Burstedde, Lucas Wilcox.

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

#ifndef P8EST_MESH_H
#define P8EST_MESH_H

#include <p8est.h>

/** Store an independent node.
 * Keep this in sync with the p8est_t data structure.
 */
typedef struct p8est_indep
{
  p4est_qcoord_t      x, y, z;
  int8_t              level, pad8;
  int16_t             pad16;
  union p4est_indep_data
  {
    void               *unused;
    p4est_topidx_t      which_tree;
    struct
    {
      p4est_topidx_t      which_tree;
      int                 owner_rank;
    }
    piggy1;
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_topidx_t      from_tree;
    }
    piggy_unused2;
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_locidx_t      local_num;
    }
    piggy3;
  }
  p;
}
p8est_indep_t;

/** Store a hanging node that depends on two independent nodes.
 * Keep this in sync with the p8est_t data structure.
 */
typedef struct p8est_hang2
{
  p4est_qcoord_t      x, y, z;
  int8_t              level, pad8;
  int16_t             pad16;
  union p8est_hang2_data
  {
    void               *unused;
    p4est_topidx_t      which_tree;
    struct
    {
      p4est_topidx_t      which_tree;
      int                 owner_rank;
    }
    piggy_unused1;
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_topidx_t      from_tree;
    }
    piggy_unused2;
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_locidx_t      depends[2];
    }
    piggy;
  }
  p;
}
p8est_hang2_t;

/** Store a hanging node that depends on four independent nodes.
 * Keep this in sync with the p8est_t data structure.
 */
typedef struct p8est_hang4
{
  p4est_qcoord_t      x, y, z;
  int8_t              level, pad8;
  int16_t             pad16;
  union p8est_hang4_data
  {
    void               *unused;
    p4est_topidx_t      which_tree;
    struct
    {
      p4est_topidx_t      which_tree;
      int                 owner_rank;
    }
    piggy_unused1;
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_topidx_t      from_tree;
    }
    piggy_unused2;
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_locidx_t      depends[4];
    }
    piggy;
  }
  p;
}
p8est_hang4_t;

/** This structure holds complete parallel node information.
 *
 * All nodes are canonicalized and store their tree id in piggy3.which_tree.
 * Their index in their owner's ordering is stored in piggy3.local_num.
 *
 * Canonicalized nodes are unique and either independent, face or edge hanging.
 * The local_nodes table is of dimension 8 * num_local_quadrants
 * and encodes the node indexes for all corners of all quadrants.  Let
 * ni := indep_nodes.elem_count,
 * fi := face_hangings.elem_count,
 * ei := edge_hangings.elem_count.
 * If for l := local_nodes[k]
 * l >= 0 && l < ni: l indexes into indep_nodes.
 * l >= ni && l < ni + fi: l - ni indexes into face_hangings.
 * l >= ni + fi && l < ni + fi + ei: l - ni - fi indexes into edge_hangings.
 * No other values for l are permitted.
 *
 * The array shared_indeps holds lists of node sharers (not including rank).
 * The entry shared_indeps[i] is of type sc_recycle_array_t
 * and holds the list of nodes with i + 1 sharers.
 * For each independent node, its member pad8 holds the number of sharers
 * and its member pad16 holds the position in the assigned recycle array
 * if this number fits into an int16_t.  If this limit is exceeded, the
 * array shared_offsets is filled with these positions as one p4est_locidx_t
 * per independent node, and all pad16 members are set to -1.  To recognize
 * the latter situation you can check for shared_offsets != NULL.
 *
 * Each processor owns num_owned_indeps of the stored independent nodes.
 * The first independent owned node is at index offset_owned_indeps.
 * The table nonlocal_ranks contains the ranks of all stored non-owned nodes.
 * The table global_owned_nodes holds the number of owned nodes for each rank.
 */
typedef struct p8est_nodes
{
  p4est_locidx_t      num_local_quadrants;
  p4est_locidx_t      num_owned_indeps, num_owned_shared;
  p4est_locidx_t      offset_owned_indeps;
  sc_array_t          indep_nodes;
  sc_array_t          face_hangings;
  sc_array_t          edge_hangings;
  p4est_locidx_t     *local_nodes;
  sc_array_t          shared_indeps;
  p4est_locidx_t     *shared_offsets;
  int                *nonlocal_ranks;
  p4est_locidx_t     *global_owned_indeps;
}
p8est_nodes_t;

/** This mesh structure holds complete neighborhood information.
 * cumulative_count[i]   is the sum of local quadrants in the
 *                       local trees 0..i-1. i == local_num_trees is allowed.
 * element_offsets[i]    is the offset into local_neighbors for local
 *                       element i. i == local_num_quadrants is allowed
 *                       and contains the number of stored local neighbors.
 */
typedef struct
{
  p4est_locidx_t     *cumulative_count;
  p4est_locidx_t     *element_offsets;
  sc_array_t         *local_neighbors;
}
p8est_neighborhood_t;

/** Create node information. */
p8est_nodes_t      *p8est_nodes_new (p8est_t * p8est,
                                     sc_array_t * ghost_layer);

/** Destroy node information. */
void                p8est_nodes_destroy (p8est_nodes_t * nodes);

/** Check node information for internal consistency. */
bool                p8est_nodes_is_valid (p8est_t * p8est,
                                          p8est_nodes_t * nodes);

/** Create neighborhood information.
 */
p8est_neighborhood_t *p8est_neighborhood_new (p8est_t * p8est);

/** Destroy neighborhood information.
 */
void                p8est_neighborhood_destroy (p8est_neighborhood_t * nhood);

#endif /* P8EST_MESH_H */
