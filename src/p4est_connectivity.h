/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007 Carsten Burstedde, Lucas Wilcox.

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

#ifndef __P4EST_CONNECTIVITY_H__
#define __P4EST_CONNECTIVITY_H__

#include <stdint.h>

typedef struct p4est_connectivity
{
  int32_t             num_trees;
  int32_t             num_vertices;
  int32_t            *tree_to_vertex;   /* allocated [0][0]..[0][3]..
                                           [num_trees-1][0]..[num_trees-1][3]

                                           Note: This is stored in
                                           right-hand-rule order.
                                         */
  int32_t            *tree_to_tree;     /* allocated [0][0]..[0][3]..
                                           [num_trees-1][0]..[num_trees-1][3] */
  int8_t             *tree_to_face;     /* allocated [0][0]..[0][3]..
                                           [num_trees-1][0]..[num_trees-1][3] */
  double             *vertices; /* allocated [0][0]..[0][2]..
                                   [num_vertices-1][0]..
                                   [num_vertices-1][2] */

}
p4est_connectivity_t;

/** Contains integers 0..7 denoting the type of inter-tree transformation.
 * The first 4 transformations are rotations about 0, -90, 180, 90.
 * The second 4 transformations are mirrors along axis 0, 45, 90, 135.
 * The indices are my_face, neighbor_face, orientation.
 * The orientation index is 0 for same, 1 for opposing sense of rotation.
 */
extern const int p4est_transform_table[4][4][2];

/** Allocate a connectivity structure
 * \param [in] num_trees    Number of trees in the forest.
 * \param [in] num_vertices Number of total vertices.
 */
p4est_connectivity_t *p4est_connectivity_new (int32_t num_trees,
                                              int32_t num_vertices);

/** Destroy a connectivity structure
 */
void                p4est_connectivity_destroy (p4est_connectivity_t *
                                                connectivity);

/** Create a connectivity structure for thu unit square
 */
p4est_connectivity_t *p4est_connectivity_new_unitsquare (void);

/** Create a connectivity structure for a three element mesh at a corner.
 *
 * This is used for testing multi tree forests.
 */
p4est_connectivity_t *p4est_connectivity_new_corner (void);

#endif /* !__P4EST_CONNECTIVITY_H__ */
