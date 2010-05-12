/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
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

#ifndef P4EST_CONNECTIVITY_H
#define P4EST_CONNECTIVITY_H

#ifdef P4_TO_P8
#error "Including a p4est header with P4_TO_P8 defined"
#endif

#include <p4est_base.h>

SC_EXTERN_C_BEGIN;

/* spatial dimension */
#define P4EST_DIM 2
#define P4EST_FACES (2 * P4EST_DIM)
#define P4EST_CHILDREN 4
#define P4EST_HALF (P4EST_CHILDREN / 2)

/* size of face transformation encoding */
#define P4EST_FTRANSFORM 9

/* p4est identification string */
#define P4EST_STRING "p4est"

/* Increase this number whenever the on-disk format for
 * p4est_connectivity, p4est, or any other 2D data structure changes.
 * The format for reading and writing must be the same.
 */
#define P4EST_ONDISK_FORMAT 0x2000007

/** This structure holds the 2D inter-tree connectivity information.
 * Identification of arbitrary faces and corners is possible.
 *
 * The arrays tree_to_* are stored in z ordering.
 * For corners the order wrt. yx is 00 01 10 11.
 * For faces the order is -x +x -y +y.
 * They are allocated [0][0]..[0][3]..[num_trees-1][0]..[num_trees-1][3].
 *
 * The values for tree_to_face are 0..7
 * where ttf % 4 gives the face number and ttf / 4 the face orientation code.
 * The orientation is 0 for edges that are aligned in z-order,
 * and 1 for edges that are running opposite in z-order.
 *
 * It is valid to specify num_vertices as 0.
 * In this case vertices and tree_to_vertex are set to NULL.
 * Otherwise the vertex coordinates are stored in the array vertices as
 * [0][0]..[0][2]..[num_vertices-1][0]..[num_vertices-1][2].
 *
 * The corners are only stored when they connect trees.
 * Otherwise the tree_to_corner entry must be -1 and this corner is ignored.
 * If num_corners == 0, tree_to_corner and corner_to_* arrays are set to NULL.
 *
 * The arrays corner_to_* store a variable number of entries per corner.
 * For corner c these are at position [ctt_offset[c]]..[ctt_offset[c+1]-1].
 * Their number for corner c is ctt_offset[c+1] - ctt_offset[c].
 * The size of the corner_to_* arrays is num_ctt = ctt_offset[num_corners].
 */
typedef struct p4est_connectivity
{
  p4est_topidx_t      num_vertices;
  p4est_topidx_t      num_trees;
  p4est_topidx_t      num_corners;

  double             *vertices;
  p4est_topidx_t     *tree_to_vertex;

  p4est_topidx_t     *tree_to_tree;
  int8_t             *tree_to_face;

  p4est_topidx_t     *tree_to_corner;
  p4est_topidx_t     *ctt_offset;
  p4est_topidx_t     *corner_to_tree;
  int8_t             *corner_to_corner;
}
p4est_connectivity_t;

typedef struct
{
  p4est_topidx_t      ntree;
  int8_t              ncorner;
}
p4est_corner_transform_t;

typedef struct
{
  p4est_topidx_t      icorner;
  sc_array_t          corner_transforms;
}
p4est_corner_info_t;

/** Store the corner numbers 0..4 for each tree face. */
extern const int    p4est_face_corners[4][2];

/** Store the face numbers in the face neighbor's system. */
extern const int    p4est_face_dual[4];

/** Store the face numbers 0..3 for each tree corner. */
extern const int    p4est_corner_faces[4][2];

/** Store the face corner numbers for the faces touching a tree corner. */
extern const int    p4est_corner_face_corners[4][4];

/** Store the faces for each child and corner, can be -1. */
extern const int    p4est_child_corner_faces[4][4];

/** Allocate a connectivity structure.
 * \param [in] num_vertices   Number of total vertices (i.e. geometric points).
 * \param [in] num_trees      Number of trees in the forest.
 * \param [in] num_corners    Number of tree-connecting corners.
 * \param [in] num_ctt        Number of total trees in corner_to_tree array.
 * \return                    A connectivity structure with allocated arrays.
 */
p4est_connectivity_t *p4est_connectivity_new (p4est_topidx_t num_vertices,
                                              p4est_topidx_t num_trees,
                                              p4est_topidx_t num_corners,
                                              p4est_topidx_t num_ctt);

/** Allocate a connectivity structure and populate from constants.
 * \param [in] num_vertices   Number of total vertices (i.e. geometric points).
 * \param [in] num_trees      Number of trees in the forest.
 * \param [in] num_corners    Number of tree-connecting corners.
 * \param [in] coff           Corner-to-tree offsets (num_corners + 1 values).
 * \return                    The connectivity is checked for validity.
 */
p4est_connectivity_t *p4est_connectivity_new_copy (p4est_topidx_t
                                                   num_vertices,
                                                   p4est_topidx_t num_trees,
                                                   p4est_topidx_t num_corners,
                                                   const double *vertices,
                                                   const p4est_topidx_t * ttv,
                                                   const p4est_topidx_t * ttt,
                                                   const int8_t * ttf,
                                                   const p4est_topidx_t * ttc,
                                                   const p4est_topidx_t *
                                                   coff,
                                                   const p4est_topidx_t * ctt,
                                                   const int8_t * ctc);

/** Destroy a connectivity structure.
 */
void                p4est_connectivity_destroy (p4est_connectivity_t *
                                                connectivity);

/** Examine a connectivity structure.
 * \return          Returns true if structure is valid, false otherwise.
 */
int                 p4est_connectivity_is_valid (p4est_connectivity_t *
                                                 connectivity);

/** Check two connectivity structures for equality.
 * \return          Returns true if structures are equal, false otherwise.
 */
int                 p4est_connectivity_is_equal (p4est_connectivity_t * conn1,
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

/** Create a connectivity structure for the unit square.
 */
p4est_connectivity_t *p4est_connectivity_new_unitsquare (void);

/** Create a connectivity structure for an all-periodic unit square.
 */
p4est_connectivity_t *p4est_connectivity_new_periodic (void);

/** Create a connectivity structure for a periodic unit square.
 * The left and right faces are identified, and bottom and top opposite.
 */
p4est_connectivity_t *p4est_connectivity_new_rotwrap (void);

/** Create a connectivity structure for a three-tree mesh around a corner.
 */
p4est_connectivity_t *p4est_connectivity_new_corner (void);

/** Create a connectivity structure for a five-tree moebius band.
 */
p4est_connectivity_t *p4est_connectivity_new_moebius (void);

/** Create a connectivity structure for a six-tree star.
 */
p4est_connectivity_t *p4est_connectivity_new_star (void);

/** Fills arrays encoding the axis combinations for a face transform.
 * \param [in]  itree       The number of the originating tree.
 * \param [in]  iface       The number of the originating face.
 * \param [out] ftransform  This array holds 9 integers.
 *              [0,2]       The coordinate axis sequence of the origin face.
 *              [3,5]       The coordinate axis sequence of the target face.
 *              [6,8]       Edge reverse flag for axis 0; face code for 1.
 *              [1,4,7]     0 (unused for compatibility with 3D).
 * \return                  The face neighbor tree if it exists, -1 otherwise.
 */
p4est_topidx_t      p4est_find_face_transform (p4est_connectivity_t *
                                               connectivity,
                                               p4est_topidx_t itree,
                                               int iface, int ftransform[]);

/** Fills an array with information about corner neighbors.
 * \param [in] itree    The number of the originating tree.
 * \param [in] icorner  The number of the originating corner.
 * \param [in,out] ci   A p4est_corner_info_t structure with initialized array.
 */
void                p4est_find_corner_transform (p4est_connectivity_t *
                                                 connectivity,
                                                 p4est_topidx_t itree,
                                                 int icorner,
                                                 p4est_corner_info_t * ci);

/** Return a pointer to a p4est_corner_transform_t array element. */
/*@unused@*/
static inline p4est_corner_transform_t *
p4est_corner_array_index (sc_array_t * array, size_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p4est_corner_transform_t));
  P4EST_ASSERT (it < array->elem_count);

  return
    (p4est_corner_transform_t *) (array->array +
                                  sizeof (p4est_corner_transform_t) * it);
}

SC_EXTERN_C_END;

#endif /* !P4EST_CONNECTIVITY_H */
