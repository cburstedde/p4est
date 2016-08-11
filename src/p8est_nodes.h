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

#ifndef P8EST_NODES_H
#define P8EST_NODES_H

#include <p8est.h>
#include <p8est_ghost.h>

SC_EXTERN_C_BEGIN;

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
      p4est_locidx_t      local_num;
    }
    piggy_unused3;
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
      p4est_locidx_t      local_num;
    }
    piggy_unused3;
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
 * Nodes are unique and either independent, face or edge hanging.
 * Independent nodes store their owner's tree id in piggy3.which_tree.
 * The index in their owner's ordering is stored in piggy3.local_num.
 * Hanging nodes store their owner's tree id in piggy.which_tree.
 * The numbers of their associated independent nodes are in piggy.depends[].
 *
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
 * The table global_owned_indeps holds the number of owned nodes for each rank.
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

/** Create node information.
 * \param [in] ghost    Ghost layer.  If this is NULL, then only
 *                      tree- and processor-local nodes will be matched
 *                      and all others duplicated, all nodes will be
 *                      counted as independent with no sharers, and
 *                      nodes->global_owned_indeps will be NULL;
 *                      this also works for a corner-unbalanced forest,
 *                      but nodes may not be numbered uniquely in this case.
 */
p8est_nodes_t      *p8est_nodes_new (p8est_t * p8est, p8est_ghost_t * ghost);

/** Destroy node information. */
void                p8est_nodes_destroy (p8est_nodes_t * nodes);

/** Check node information for internal consistency. */
int                 p8est_nodes_is_valid (p8est_t * p8est,
                                          p8est_nodes_t * nodes);

SC_EXTERN_C_END;

#endif /* !P8EST_NODES_H */
