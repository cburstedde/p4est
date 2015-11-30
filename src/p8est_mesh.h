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

/** \file p8est_mesh.h
 *
 * forest topology in a conventional mesh format
 *
 * \ingroup p8est
 */

#ifndef P8EST_MESH_H
#define P8EST_MESH_H

#include <p8est_ghost.h>

SC_EXTERN_C_BEGIN;

/** This structure contains complete mesh information on the forest.
 * It stores the locally relevant neighborhood, that is, all locally owned
 * quadrants and one layer of adjacent ghost quadrants and their owners.
 *
 * For each local quadrant, its tree number is stored in quad_to_tree. The
 * quad_to_tree array is NULL by default and can be enabled using
 * \ref p8est_mesh_new_ext.
 * For each ghost quadrant, its owner rank is stored in ghost_to_proc.
 * For each level, an array of local quadrant numbers is stored in quad_level.
 * The quad_level array is NULL by default and can be enabled using
 * \ref p8est_mesh_new_ext.
 *
 * The quad_to_quad list stores one value for each local quadrant's face.
 * This value is in 0..local_num_quadrants-1 for local quadrants, or in
 * local_num_quadrants + (0..ghost_num_quadrants-1) for ghost quadrants.
 * The quad_to_face list has equally many entries which are either:
 * 1. A value of v = 0..23 indicates one same-size neighbor.
 *    This value is decoded as v = r * 6 + nf, where nf = 0..5 is the
 *    neighbor's connecting face number and r = 0..3 is the relative
 *    orientation of the neighbor's face, see p8est_connectivity.h.
 * 2. A value of v = 24..119 indicates a double-size neighbor.
 *    This value is decoded as v = 24 + h * 24 + r * 6 + nf, where
 *    r and nf are as above and h = 0..3 is the number of the subface.
 *    TODO: Define what perspective is used to define h?
 *          We may prefer it to be the perspective of the neighbor,
 *          as is the usual convention for quad_to_face.
 *          Our own perspective can be derived from our child id.
 *          On the other hand, getting the child id requires quadrant access,
 *          which we may want to avoid in most cases.
 * 3. A value of v = -24..-1 indicates four half-size neighbors.
 *    In this case the corresponding quad_to_quad index points into the
 *    quad_to_half array that stores four quadrant numbers per index,
 *    and the orientation of the smaller faces follows from 24 + v.
 *    The entries of quad_to_half encode between local and ghost quadrant
 *    in the same way as the quad_to_quad values described above.
 *    TODO: Define exactly in which sequence the four small neighbors
 *          are stored in the current version of the code.
 *          We may subsequently consider reordering them.
 *
 * A quadrant on the boundary of the forest sees itself and its face number.
 *
 * The quad_to_edge list stores edge neighbors that are not face neighbors.
 * On the inside of a tree, there are one or two of those depending on size.
 * Between trees, there can be any number of same- or different sized neighbors.
 * For same-tree same/double-size edge neighbors, we record their number in
 * quad_to_edge by the same convention as described for quad_to_quad above.
 *
 * For half-size neighbors and all inter-tree neighbors, the quad_to_edge value
 * is in
 *    local_num_quadrants + local_num_ghosts + [0 .. local_num_edges - 1]
 * where the offset by local quadrants and ghosts is implicitly subtracted.
 * It indexes into edge_offset, which encodes a group of edge neighbors.
 * Each member of a group may be one same/double-size quadrant or two half-size
 * quadrants; this is determined by the value of the edge_edge field as follows.
 * 1. A value of e = 0..23 indicates one same-size neighbor.
 *    This value is encoded as e = r * 12 + ne, where ne = 0..11 is the
 *    neighbor's connecting edge number and r = 0..1 indicates an edge flip.
 * 2. A value of e = 24..71 indicates a double-size neighbor.
 *    This value is decoded as e = 24 + h * 24 + r * 12 + ne, where
 *    r and ne are as above and h = 0..1 is the number of the subedge.
 *    TODO: How do we interpret h?
 *          Is it what the respective neighbor sees the small quadrant as?
 * 3. A value of e = -24..-1 indicates two half-size neighbors.
 *    In this case the corresponding edge_to_quad index points into the
 *    quad_to_hedge array that stores two quadrant numbers per index,
 *    and the orientation of the smaller edges follows from 24 + e.
 *    The entries of quad_to_hedge encode between local and ghost quadrant
 *    in the same way as the quad_to_quad values described above.
 *    TODO: In what sequence are these neighbors stored in quad_to_hedge?
 *          Compare this is the same convention as with quad_to_half.
 *
 * TODO: Idea to remove quad_to_hedge in favor of storing two quadrant ids
 *       in the offset structure.  This would save some memory.
 *
 * TODO: With the current code test_mesh fails in debug mode (assertion).
 *
 * The quad_to_corner list stores corner neighbors that are not face or edge
 * neighbors.  On the inside of a tree, there is precisely one such neighbor
 * per corner.  In this case, its index is encoded as described above for
 * quad_to_quad.  The neighbor's matching corner number is always diagonally
 * opposite, that is, corner number ^ 7.
 *
 * On the inside of an inter-tree face, we have precisely one corner neighbor.
 * If a corner is across an inter-tree edge or corner, then the number of
 * corner neighbors may be any non-negative number.  In all inter-tree cases,
 * the quad_to_corner value is in
 *    local_num_quadrants + local_num_ghosts + [0 .. local_num_corners - 1]
 * where the offset by local quadrants and ghosts is implicitly subtracted.
 * It indexes into corner_offset, which encodes a group of corner neighbors.
 * Each group contains the quadrant numbers encoded as usual for quad_to_quad
 * in corner_quad, and the corner number from the neighbor as corner_corner.
 *
 * Intra-tree corners and inter-tree face and corner corners are implemented.
 * Inter-tree-edge corners are NOT IMPLEMENTED and are assigned the value -2.
 * Currently we do NOT exclude inter-tree-corners of inter-tree-edge neighbors.
 * Corners with no diagonal neighbor at all are assigned the value -1.
 */
typedef struct
{
  p4est_locidx_t      local_num_quadrants;
  p4est_locidx_t      ghost_num_quadrants;

  p4est_topidx_t     *quad_to_tree;     /**< tree index for each local quad,
                                             NULL by default */
  int                *ghost_to_proc;    /**< processor for each ghost quad */

  p4est_locidx_t     *quad_to_quad;     /**< one index for each of the 6 faces */
  int8_t             *quad_to_face;     /**< encodes orientation/2:1 status */
  sc_array_t         *quad_to_half;     /**< stores half-size neighbors */
  sc_array_t         *quad_level;       /**< stores lists of per-level quads,
                                             NULL by default */

  /* These members are NULL if connect_t is less than P4EST_CONNECT_EDGE */
  /* CAUTION: not yet implemented */
  p4est_locidx_t      local_num_edges;  /**< half-size and tree-boundary edges */
  p4est_locidx_t     *quad_to_edge;     /**< 12 indices for each local quad */
  /* TODO: use these containers at all and use them properly */
  sc_array_t         *edge_offset;      /**< local_num_edges + 1 entries */
  sc_array_t         *edge_quad;        /**< edge_offset indexes into this */
  sc_array_t         *edge_edge;        /**< and this one too (type int8_t) */

  /* These members are NULL if the connect_t is not P4EST_CONNECT_CORNER */
  /* CAUTION: tree-boundary corners do not exclude tree-boundary edge corners */
  p4est_locidx_t      local_num_corners;        /* tree-boundary corners */
  p4est_locidx_t     *quad_to_corner;   /* 8 indices for each local quad */
  sc_array_t         *corner_offset;    /* local_num_corners + 1 entries */
  sc_array_t         *corner_quad;      /* corner_offset indexes into this */
  sc_array_t         *corner_corner;    /* and this one too (type int8_t) */
}
p8est_mesh_t;

/** This structure can be used as the status of a face neighbor iterator.
  * It always contains the face and subface of the neighbor to be processed.
  */
typedef struct
{
  /* forest information */
  p8est_t            *p4est;
  p8est_ghost_t      *ghost;
  p8est_mesh_t       *mesh;

  /* quadrant information */
  p4est_topidx_t      which_tree;
  p4est_locidx_t      quadrant_id;      /* tree-local quadrant index */
  p4est_locidx_t      quadrant_code;    /* 6 * (quadrant_id + tree_offset) */

  /* neighbor information */
  int                 face;     /* Face number in 0..5. */
  int                 subface;  /* Hanging neighbor number in 0..3. */

  /* internal information */
  p4est_locidx_t      current_qtq;
}
p8est_mesh_face_neighbor_t;

/** Calculate the memory usage of the mesh structure.
 * \param [in] mesh     Mesh structure.
 * \return              Memory used in bytes.
 */
size_t              p8est_mesh_memory_used (p8est_mesh_t * mesh);

/** Create a p8est_mesh structure.
 * This function does not populate the quad_to_tree and quad_level fields.
 * To populate them, use \ref p8est_mesh_new_ext.
 * \param [in] p8est    A forest that is fully 2:1 balanced.
 * \param [in] ghost    The ghost layer created from the provided p4est.
 * \param [in] btype    Determines the highest codimension of neighbors.
 * \return              A fully allocated mesh structure.
 */
p8est_mesh_t       *p8est_mesh_new (p8est_t * p8est, p8est_ghost_t * ghost,
                                    p8est_connect_type_t btype);

/** Destroy a p8est_mesh structure.
 * \param [in] mesh     Mesh structure previously created by p8est_mesh_new.
 */
void                p8est_mesh_destroy (p8est_mesh_t * mesh);

/** Find a quadrant based on its cumulative number in the local forest.
 * If the quad_to_tree field of the mesh structure exists, this is O(1).
 * Otherwise, we perform a binary search over the processor-local trees.
 *
 * \param [in]  p8est           Forest to be worked with.
 * \param [in]  mesh            A mesh derived from the forest.
 * \param [in]  cumulative_id   Cumulative index over all trees of quadrant.
 *                              Must refer to a local (non-ghost) quadrant.
 * \param [in,out] which_tree   If not NULL, the input value can be -1
 *                              or an initial guess for the quadrant's tree
 *                              and output is the tree of returned quadrant.
 * \param [out] quadrant_id     If not NULL, the number of quadrant in tree.
 * \return                      The identified quadrant.
 */
p8est_quadrant_t   *p8est_mesh_quadrant_cumulative (p8est_t * p8est,
                                                    p8est_mesh_t * mesh,
                                                    p4est_locidx_t
                                                    cumulative_id,
                                                    p4est_topidx_t *
                                                    which_tree,
                                                    p4est_locidx_t *
                                                    quadrant_id);

/** Initialize a mesh neighbor iterator by quadrant index.
 * \param [out] mfn         A p8est_mesh_face_neighbor_t to be initialized.
 * \param [in]  which_tree  Tree of quadrant whose neighbors are looped over.
 * \param [in]  quadrant_id Index relative to which_tree of quadrant.
 */
void                p8est_mesh_face_neighbor_init2 (p8est_mesh_face_neighbor_t
                                                    * mfn, p8est_t * p8est,
                                                    p8est_ghost_t * ghost,
                                                    p8est_mesh_t * mesh,
                                                    p4est_topidx_t which_tree,
                                                    p4est_locidx_t
                                                    quadrant_id);

/** Initialize a mesh neighbor iterator by quadrant pointer.
 * \param [out] mfn         A p8est_mesh_face_neighbor_t to be initialized.
 * \param [in]  which_tree  Tree of quadrant whose neighbors are looped over.
 * \param [in]  quadrant    Pointer to quadrant contained in which_tree.
 */
void                p8est_mesh_face_neighbor_init (p8est_mesh_face_neighbor_t
                                                   * mfn, p8est_t * p8est,
                                                   p8est_ghost_t * ghost,
                                                   p8est_mesh_t * mesh,
                                                   p4est_topidx_t which_tree,
                                                   p8est_quadrant_t *
                                                   quadrant);

/** Move the iterator forward to loop around neighbors of the quadrant.
 * \param [in,out] mfn      Internal status of the iterator.
 * \param [out]    ntree    If not NULL, the tree number of the neighbor.
 * \param [out]    nquad    If not NULL, the quadrant number within tree.
 *                          For ghosts instead the number in ghost layer.
 * \param [out]    nface    If not NULL, neighbor's face as in p8est_mesh_t.
 * \param [out]    nrank    If not NULL, the owner process of the neighbor.
 * \return                  Either a real quadrant or one from the ghost layer.
 *                          Returns NULL when the iterator is done.
 */
p8est_quadrant_t   *p8est_mesh_face_neighbor_next (p8est_mesh_face_neighbor_t
                                                   * mfn,
                                                   p4est_topidx_t * ntree,
                                                   p4est_locidx_t * nquad,
                                                   int *nface, int *nrank);

/** Get the user data for the current face neighbor.
 * \param [in]     mfn           Internal status of the iterator.
 * \param [in]     ghost_data    Data for the ghost quadrants that has been
 *                               synchronized with p4est_ghost_exchange_data.
 * \return                       A pointer to the user data for the current
 *                               neighbor.
 */
void               *p8est_mesh_face_neighbor_data (p8est_mesh_face_neighbor_t
                                                   * mfn, void *ghost_data);

SC_EXTERN_C_END;

#endif /* !P8EST_MESH_H */
