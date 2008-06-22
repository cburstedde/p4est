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

#ifndef P8EST_CONNECTIVITY_H
#define P8EST_CONNECTIVITY_H

#include <p4est_base.h>

/** This structure holds the 3D inter-tree connectivity information.
 * Identification of separate faces and corners is possible.
 *
 * The arrays tree_to_* are stored in z ordering.
 * For vertices the order is 000 001 010 011 100 101 110 111.
 * For faces the order is -x +x -y +y -z +z.
 * They are allocated [0][0]..[0][8 or 6, resp.]..[num_trees-1][0]...
 *
 * The values for tree_to_face are 0..5 for equal orientation
 * and 6..11 for opposite orientation.
 *
 * The vertex coordinates are stored in the array vertices, allocated
 * [0][0]..[0][2]..[num_vertices-1][0]..[num_vertices-1][2].
 *
 * The arrays vertex_to_* store a variable number of entries per vertex.
 * For vertex v these are at position [vtt_offset[v]]..[vtt_offset[v+1]-1].
 * Their number for vertex v is vtt_offset[v+1] - vtt_offset[v].
 * The size of the vertex_to_* arrays is vtt_offset[num_vertices].
 */
typedef struct p8est_connectivity
{
  p4est_topidx_t      num_trees;
  p4est_topidx_t      num_vertices;

  p4est_topidx_t     *tree_to_vertex;
  p4est_topidx_t     *tree_to_tree;
  int8_t             *tree_to_face;

  double             *vertices;

  p4est_topidx_t     *vtt_offset;
  p4est_topidx_t     *vertex_to_tree;
  p4est_topidx_t     *vertex_to_vertex;
}
p8est_connectivity_t;

typedef struct
{
  p4est_topidx_t      ntree;
  int8_t              ncorner;
}
p8est_corner_info_t;

#if 0                           /* this was 2D stuff may be useful later */
/** Contains integers 0..7 denoting the type of inter-tree transformation.
 * The first 4 transformations are rotations about 0, -90, 180, 90.
 * The second 4 transformations are mirrors along axis 0, 45, 90, 135.
 * The indices are my_face, neighbor_face, orientation.
 * The orientation index is 0 for same, 1 for opposing sense of rotation.
 */
extern const int    p8est_transform_table[4][4][2];
#endif

/** Store the vertex numbers 0..7 for each tree face. */
extern const int    p8est_face_vertices[6][4];

/** Store the face numbers 0..12 for each tree face. */
extern const int    p8est_face_edges[6][4];

/** Store only the 8 out of 24 possible permutations that occur. */
extern const int    p8est_face_permutations[8][4];

/** Store the 3 occurring sets of 4 permutations per face. */
extern const int    p8est_face_permutation_sets[3][4];

/** For each face combination store the permutation set.
 * The order is [my_face][neighbor_face]
 */
extern const int    p8est_face_permutation_refs[6][6];

/** Store the vertex numbers 0..8 for each tree edge. */
extern const int    p8est_edge_vertices[12][2];

/** Allocate a connectivity structure
 * \param [in] num_trees    Number of trees in the forest.
 * \param [in] num_vertices Number of total vertices.
 * \param [in] num_vtt      Number of total trees in vertex_to_tree array.
 * \param [in] alloc_vertices   Boolean flag for vertex memory allocation.
 */
p8est_connectivity_t *p8est_connectivity_new (p4est_topidx_t num_trees,
                                              p4est_topidx_t num_vertices,
                                              p4est_topidx_t num_vtt,
                                              bool alloc_vertices);

/** Destroy a connectivity structure
 */
void                p8est_connectivity_destroy (p8est_connectivity_t *
                                                connectivity);

/** Examine a connectivity structure
 * \return  Returns 1 if structure is valid, 0 otherwise.
 */
bool                p8est_connectivity_is_valid (p8est_connectivity_t *
                                                 connectivity);

/** Create a connectivity structure for the unit cube.
 */
p8est_connectivity_t *p8est_connectivity_new_unitcube (void);

/** Create a connectivity structure for a mostly periodic unit cube.
 * The left and right faces are identified, and bottom and top opposite.
 * Front and back are not identified.
 */
p8est_connectivity_t *p8est_connectivity_new_periodic (void);

/** Create a connectivity structure that contains a few cubes.
 * These are rotated against each other to stress the topology routines.
 */
p8est_connectivity_t *p8est_connectivity_new_rotcubes (void);

#if 0
/** Fills an array with information about corner neighbors.
 * \param [in,out]  corner_info  Array of p8est_corner_info_t members.
 */
void                p8est_find_corner_info (p8est_connectivity_t *
                                            connectivity,
                                            p4est_topidx_t itree, int icorner,
                                            sc_array_t * corner_info);
#endif

/** Fills arrays encoding the axis combinations for a face transform.
 * \param [out] my_axis      The coordinate axis sequence of the origin face.
 * \param [out] target_axis  The coordinate axis sequence of the target face.
 * \param [out] edge_reverse Edge reverse flag for axes 0, 1; face code for 2.
 * \return   Returns the face neighbor tree if it exists, -1 otherwise.
 */
p4est_topidx_t      p8est_find_face_transform (p8est_connectivity_t *
                                               connectivity,
                                               p4est_topidx_t itree,
                                               int iface, int my_axis[3],
                                               int target_axis[3],
                                               int edge_reverse[3]);

#endif /* !P8EST_CONNECTIVITY_H */
