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
 * They are allocated [0][0]..[0][N-1]..[num_trees-1][0]..[num_trees-1][N-1].
 * where N is 6 for tree and face, 8 for vertex, 12 for edge.
 *
 * The values for tree_to_face are in 0..23
 * where ttf % 6 gives the face number and ttf / 6 the orientation code.
 *
 * The vertex coordinates are stored in the array vertices, allocated
 * [0][0]..[0][2]..[num_vertices-1][0]..[num_vertices-1][2].
 *
 * The edges are only stored and counted when they are relevant for balance.
 * Otherwise the tree_to_edge entry must be -1 and this edge will be ignored.
 *
 * The arrays edge_to_* store a variable number of entries per edge.
 * For edge e these are at position [ett_offset[e]]..[ett_offset[e+1]-1].
 * Their number for edge e is ett_offset[e+1] - ett_offset[e].
 * The size of the edge_to_* arrays is num_ett = ett_offset[num_edges].
 * The edge_to_edge array holds values in 0..23, where the lower 12 indicate
 * one orientation and the higher 12 the opposite orientation.
 *
 * The arrays vertex_to_* store a variable number of entries per vertex.
 * For vertex v these are at position [vtt_offset[v]]..[vtt_offset[v+1]-1].
 * Their number for vertex v is vtt_offset[v+1] - vtt_offset[v].
 * The size of the vertex_to_* arrays is num_vtt = vtt_offset[num_vertices].
 */
typedef struct p8est_connectivity
{
  p4est_topidx_t      num_vertices;
  p4est_topidx_t      num_trees;
  p4est_topidx_t      num_edges;
  p4est_topidx_t      num_corners;

  double             *vertices;
  p4est_topidx_t     *tree_to_vertex;

  p4est_topidx_t     *tree_to_tree;
  int8_t             *tree_to_face;

  p4est_topidx_t     *tree_to_edge;
  p4est_topidx_t     *ett_offset;
  p4est_topidx_t     *edge_to_tree;
  int8_t             *edge_to_edge;

  p4est_topidx_t     *tree_to_corner;
  p4est_topidx_t     *ctt_offset;
  p4est_topidx_t     *corner_to_tree;
  int8_t             *corner_to_corner;
}
p8est_connectivity_t;

typedef struct
{
  p4est_topidx_t      ntree;
  int8_t              nedge, naxis[3], nflip, corners;
}
p8est_edge_transform_t;

typedef struct
{
  int8_t              iedge, iflip;
  sc_array_t          edge_transforms;
}
p8est_edge_info_t;

#if 0
typedef struct
{
  p4est_topidx_t      ntree;
  int8_t              ncorner;
}
p8est_corner_info_t;
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
 * \param [in] num_vertices   Number of total vertices (i.e. geometric points).
 * \param [in] num_trees      Number of trees in the forest.
 * \param [in] num_edges      Number of balance-relevant identified edges.
 * \param [in] num_ett        Number of total trees in edge_to_tree array.
 * \param [in] num_corners    Number of balance-relevant identified corners.
 * \param [in] num_ctt        Number of total trees in corner_to_tree array.
 * \param [in] alloc_vertices   Specify if coordinate array should be alloced.
 */
p8est_connectivity_t *p8est_connectivity_new (p4est_topidx_t num_vertices,
                                              p4est_topidx_t num_trees,
                                              p4est_topidx_t num_edges,
                                              p4est_topidx_t num_ett,
                                              p4est_topidx_t num_corners,
                                              p4est_topidx_t num_ctt,
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
 * The left and right faces are identified, and bottom and top rotated.
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

/** Fills an array with information about edge neighbors.
 * \param [in] itree    The number of the originating tree.
 * \param [in] iedge    The number of the originating edge.
 * \param [in,out] ei   A p8est_edge_info_t structure with initialized array.
 */
void                p8est_find_edge_transform (p8est_connectivity_t *
                                               connectivity,
                                               p4est_topidx_t itree,
                                               int iedge,
                                               p8est_edge_info_t * ei);

#endif /* !P8EST_CONNECTIVITY_H */
