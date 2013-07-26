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

#ifndef P8EST_MESH_H
#define P8EST_MESH_H

#include <p8est_ghost.h>

SC_EXTERN_C_BEGIN;

/** This structure contains complete mesh information on the forest.
 * It stores the locally relevant neighborhood, that is, all locally owned
 * quadrants and one layer of adjacent ghost quadrants and their owners.
 *
 * All vertices of the locally owned quadrants are stored in xyz triples.
 * For each local quadrant, its eight corners point into the vertices array.
 * Some vertices are duplicated and no effort is made to uniquify them.
 * For each ghost quadrant, its owner rank is stored in ghost_to_proc,
 * and its number in it's owners range of local quadrants in ghost_to_index.
 *
 * The quad_to_quad list stores one value for each local quadrant's face.
 * This value is in 0..local_num_quadrants-1 for local quadrants, or in
 * local_num_quadrants + (0..ghost_num_quadrants-1) for ghost quadrants.
 * The quad_to_face list has equally many entries which are either:
 * 1. A value of v = 0..23 indicates one same-size neighbor.
 *    This value is decoded as v = r * 6 + nf, where nf = 0..5 is the
 *    neigbbor's connecting face number and r = 0..3 is the relative
 *    orientation of the neighbor's face, see p8est_connectivity.h.
 * 2. A value of v = 24..119 indicates a double-size neighbor.
 *    This value is decoded as v = 24 + h * 24 + r * 6 + nf, where
 *    r and nf are as above and h = 0..3 is the number of the subface.
 * 3. A value of v = -24..-1 indicates four half-size neighbors.
 *    In this case the corresponding quad_to_quad index points into the
 *    quad_to_half array which stores four quadrant numbers per index,
 *    and the orientation of the smaller faces follows from 24 + v.
 *    The entries of quad_to_half encode between local and ghost quadrant
 *    in the same way as the quad_to_quad values described above.
 * A quadrant on the boundary of the forest sees itself and its face number.
 */
typedef struct
{
  p4est_locidx_t      local_num_vertices;
  p4est_locidx_t      local_num_quadrants;
  p4est_locidx_t      ghost_num_quadrants;

  p4est_topidx_t     *quad_to_tree;     /* Tree index for each quad */

  /* Vertex information for each local quadrant */
  double             *vertices;
  p4est_locidx_t     *quad_to_vertex;   /* 8 indices for each local quad */
  int                *ghost_to_proc;    /* 1 integer for each ghost quad */
  p4est_locidx_t     *ghost_to_index;   /* 1 remote index for each ghost */

  p4est_locidx_t     *quad_to_quad;     /* 1 index for each of the 6 faces */
  int8_t             *quad_to_face;     /* encodes orientation/2:1 status */
  sc_array_t         *quad_to_half;     /* stores half-size neigbors */
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
}
p8est_mesh_face_neighbor_t;

/** Calculate the memory usage of the mesh structure.
 * \param [in] mesh     Mesh structure.
 * \return              Memory used in bytes.
 */
size_t              p8est_mesh_memory_used (p8est_mesh_t * mesh);

/** Create a p8est_mesh structure.
 * The vertex information will be filled if p8est->connectivity contains
 * vertices.  Currently only face neighborhood information is stored.
 * \param [in] p8est    A forest that is fully 2:1 balanced.
 * \param [in] ghost    The ghost layer created from the provided p4est.
 * \param [in] btype    Currently ignored, only face neighbors are stored.
 * \return              A fully allocated mesh structure.
 */
p8est_mesh_t       *p8est_mesh_new (p8est_t * p8est, p8est_ghost_t * ghost,
                                    p8est_connect_type_t btype);

/** Destroy a p8est_mesh structure.
 * \param [in] mesh     Mesh structure previously created by p8est_mesh_new.
 */
void                p8est_mesh_destroy (p8est_mesh_t * mesh);

/** Find a quadrant based on its cumulative number in the local forest.
 * \param [in]  p8est           Forest to be worked with.
 * \param [in]  cumulative_id   Cumulative index over all trees of quadrant.
 * \param [in,out] which_tree   If not NULL, the input value can be -1
 *                              or an initial guess for the quadrant's tree
 *                              and output is the tree of returned quadrant.
 * \param [out] quadrant_id     If not NULL, the number of quadrant in tree.
 * \return                      The identified quadrant.
 */
p8est_quadrant_t   *p8est_mesh_quadrant_cumulative (p8est_t * p8est,
                                                    p4est_locidx_t
                                                    cumulative_id,
                                                    p4est_topidx_t
                                                    * which_tree,
                                                    p4est_locidx_t
                                                    * quadrant_id);

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
                                                   p8est_quadrant_t
                                                   * quadrant);

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

SC_EXTERN_C_END;

#endif /* !P8EST_MESH_H */
