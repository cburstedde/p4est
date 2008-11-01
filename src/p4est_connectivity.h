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

#ifndef P4EST_CONNECTIVITY_H
#define P4EST_CONNECTIVITY_H

#ifdef P4_TO_P8
#error "Including a p4est header with P4_TO_P8 defined"
#endif

#include <p4est_base.h>

SC_EXTERN_C_BEGIN;

/* Increase this number whenever the on-disk format for
 * p4est_connectivity, p4est, or any other 2D data structure changes.
 * The format for reading and writing must be the same.
 */
#define P4EST_ONDISK_FORMAT 0x2000003

/** This structure holds the 2D inter-tree connectivity information.
 * Identification of separate faces and corners is possible.
 *
 * The arrays tree_to_* are stored in right-hand rule ordering.
 * They are allocated [0][0]..[0][3]..[num_trees-1][0]..[num_trees-1][3].
 *
 * The values for tree_to_face are 0..3 for equal orientation
 * and 4..7 for opposite orientation, both in right-hand rule.
 *
 * The vertex coordinates are stored in the array vertices, allocated
 * [0][0]..[0][2]..[num_vertices-1][0]..[num_vertices-1][2].
 *
 * The arrays vertex_to_* store a variable number of entries per vertex.
 * For vertex v these are at position [vtt_offset[v]]..[vtt_offset[v+1]-1].
 * Their number for vertex v is vtt_offset[v+1] - vtt_offset[v].
 * The size of the vertex_to_* arrays is vtt_offset[num_vertices].
 */
typedef struct p4est_connectivity
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
p4est_connectivity_t;

typedef struct
{
  p4est_topidx_t      ntree;
  int8_t              ncorner;  /* this corner is in z-order */
}
p4est_corner_transform_t;

/** Mappings between right-hand rule and z-ordering. */
extern const int    p4est_corner_to_zorder[5];
extern const int    p4est_zface_to_rface[4];

/** Contains integers 0..7 denoting the type of inter-tree transformation.
 * The first 4 transformations are rotations about 0, -90, 180, 90.
 * The second 4 transformations are mirrors along axis 0, 45, 90, 135.
 * The indices are my_face, neighbor_face, orientation.
 * The orientation index is 0 for same, 1 for opposing sense of rotation.
 */
extern const int    p4est_transform_table[4][4][2];

/** Store the face numbers in the face neighbor's system. */
extern const int    p4est_face_dual[4];

/** Store the hanging face number in the big neighbor of a small quadrant. */
extern const int    p4est_face_child_hang[4][4];

/** Hanging corners and faces indexed by child id, two each. */
extern const int    p4est_hanging_corner[4][2];
extern const int    p4est_hanging_face[4][2];

/** Allocate a connectivity structure
 * \param [in] num_trees    Number of trees in the forest.
 * \param [in] num_vertices Number of total vertices.
 * \param [in] num_vtt      Number of total trees in vertex_to_tree array.
 * \param [in] alloc_vertices   Flag for vertex coordinate allocation.
 */
p4est_connectivity_t *p4est_connectivity_new (p4est_topidx_t num_trees,
                                              p4est_topidx_t num_vertices,
                                              p4est_topidx_t num_vtt,
                                              bool alloc_vxyz);

/** Destroy a connectivity structure.
 */
void                p4est_connectivity_destroy (p4est_connectivity_t *
                                                connectivity);

/** Examine a connectivity structure.
 * \return          Returns true if structure is valid, false otherwise.
 */
bool                p4est_connectivity_is_valid (p4est_connectivity_t *
                                                 connectivity);

/** Check two connectivity structures for equality.
 * \return          Returns true if structures are equal, false otherwise.
 */
bool                p4est_connectivity_is_equal (p4est_connectivity_t * conn1,
                                                 p4est_connectivity_t *
                                                 conn2);

/** Save a connectivity structure to disk.
 * \param [in] filename         Name of the file to write.
 * \param [in] connectivity     Valid connectivity structure.
 * \note            Aborts on file errors.
 */
void                p4est_connectivity_save (const char *filename,
                                             p4est_connectivity_t *
                                             connectivity);

/** Load a connectivity structure from disk.
 * \param [in] filename Name of the file to read.
 * \param [out] length  Size in bytes of connectivity on disk or NULL.
 * \return              Returns a valid connectivity structure.
 * \note                Aborts on file errors or invalid data.
 */
p4est_connectivity_t *p4est_connectivity_load (const char *filename,
                                               long *length);

/** Create a connectivity structure for the unit square
 */
p4est_connectivity_t *p4est_connectivity_new_unitsquare (void);

/** Create a connectivity structure for a three tree mesh at a corner.
 */
p4est_connectivity_t *p4est_connectivity_new_corner (void);

/** Create a connectivity structure for a five tree moebius band.
 */
p4est_connectivity_t *p4est_connectivity_new_moebius (void);

/** Create a connectivity structure for a six tree star.
 */
p4est_connectivity_t *p4est_connectivity_new_star (void);

/** Create a connectivity structure for an all-periodic unit square.
 * The left and right faces are identified, and bottom and top opposite.
 */
p4est_connectivity_t *p4est_connectivity_new_periodic (void);

/** Returns the transformation number from a tree to a neighbor tree.
 * \param [int] iface   Face id in right-hand-rule order.
 * \return   Returns -1 if there is no neighbor at that face, or 0..7.
 */
int                 p4est_find_face_transform (p4est_connectivity_t *
                                               connectivity,
                                               p4est_topidx_t itree,
                                               int iface);

/** Fills an array with information about corner neighbors.
 * \param [int] icorner          Corner id in right-hand-rule order.
 * \param [in,out] ctransforms   Array of p4est_corner_transform_t members.
 */
void                p4est_find_corner_transform (p4est_connectivity_t *
                                                 connectivity,
                                                 p4est_topidx_t itree,
                                                 int icorner,
                                                 sc_array_t * ctransforms);

SC_EXTERN_C_END;

#endif /* !P4EST_CONNECTIVITY_H */
