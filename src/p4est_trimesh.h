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

/** \file p4est_trimesh.h
 *
 * Generate a conforming triangle mesh from a 2:1 balanced p4est.
 * This mesh is represented by augmenting the \ref p4est_lnodes structure.
 */

#ifndef P4EST_TRIMESH_H
#define P4EST_TRIMESH_H

#include <p4est_lnodes.h>

SC_EXTERN_C_BEGIN;

#if 0
typedef struct p4est_tnode_t
{
  p4est_qcoord_t      x;
  p4est_qcoord_t      y;
  p4est_topidx_t      which_tree;
}
p4est_tnode_t;
#endif

/** Lookup table structure defining a conforming triangle mesh.
 *
 * The \a lnodes member encodes the process-relavent corners and faces.
 * Triangle-shaped element and corner entities are always included.
 * It can be created with or without including faces as mesh entities.
 * The members of \a lnodes are reinterpreted; cf. \ref p4est_lnodes.h :
 *  - degree is set to 0.
 *  - vnodes is the maxium number of nodes per element, 9 (corners only)
 *    or 25 (with faces).  Each element gets this amount of memory in the
 *    \a element_nodes member.  Unused positions are set to -1.
 *    The position of the corner nodes wrt. the element are as follows:
 *
 *        [ center node, 4 x hanging face midpoint as corner node
 *          in face order (or -1), 4 x corner node in corner order ].
 *
 *    If face midpoint nodes are requested, their positions are appended:
 *
 *        [ 4 x face connecting quadrant center with hanging face midpoint
 *          in order (or -1 if face not hanging), 4 x face connecting center
 *          with a corner in order, 8 x midpoint of half faces:
 *          2 each for face f in order ].
 *
 *    If a face is not hanging, the first half face midpoint position serves
 *    to store the full face midpoint and the second becomes -1 instead.
 *  - face_code as defined in \ref p4est_lnodes.h encodes hanging neighbors.
 *    Each valid face_code determines one possible triangle layout.
 *  - According to the node layout, the nodes of the elemnt are encoded.
 */
typedef struct p4est_trimesh
{
#if 0
  sc_array_t         *onodes;   /**< owned nodes: p4est_tnode_t */
  sc_array_t         *snodes;   /**< shared nodes: p4est_tnode_t */
#endif
  p4est_lnodes_t     *lnodes;   /**< Mesh metadata; \see p4est_trimesh_t. */
  char               *nflags;   /**< One byte for each local node.  The
                                     low 2 bits encode the codimension.
                                     The 4-th lowest bit says boundary. */
}
p4est_trimesh_t;

/** Generate a conforming triangle mesh from a 2:1 balance forest.
 * \param [in] p4est    Valid forest after 2:1 (at least face) balance.
 * \param [in] ghost    Ghost layer created form \b p4est.  May be NULL.
 * \param [in] with_faces   If false, only number triangles and mesh corners.
 *                          Otherwise, include each triangle face as entity.
 * \return              Valid conforming triangle mesh structure.
 */
p4est_trimesh_t    *p4est_trimesh_new (p4est_t * p4est,
                                       p4est_ghost_t * ghost, int with_faces);

/** Free the memory in a conforming triangle mesh structure.
 * \param [in] trimesh      Memory is deallocated.  Do not use after return.
 */
void                p4est_trimesh_destroy (p4est_trimesh_t * trimesh);

SC_EXTERN_C_END;

#endif /* !P4EST_TRIMESH_H */
