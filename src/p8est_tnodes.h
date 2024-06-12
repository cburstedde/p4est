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

/** \file p8est_tnodes.h
 *
 * Generate a conforming tetrahedron mesh from a 2:1 balanced p8est.
 * This mesh is represented by augmenting the \ref p8est_lnodes structure.
 */

#ifndef P8EST_TNODES_H
#define P8EST_TNODES_H

#include <p8est_lnodes.h>

SC_EXTERN_C_BEGIN;

/** Type for the code of each cube simplex point. */
typedef int8_t      p8est_tnodes_eind_code_t;

/** Type for the simplex node index into a cubic lattice. */
typedef int16_t     p8est_tnodes_eindex_t;

/** Type for the simplex sort key. */
typedef int16_t     p8est_tnodes_simplex_key_t;

/** Type for a simplex of the elementary refinement forest. */
typedef struct p8est_tnodes_simplex p8est_tnodes_simplex_t;

/** Memory for a simplex of the elementary refinement forest. */
struct p8est_tnodes_simplex
{
  p8est_tnodes_simplex_t *parent;   /**< Pointer to parent simplex. */
  p8est_tnodes_eindex_t nodes[4];   /**< Indices of corner nodes. */
  p8est_tnodes_eindex_t lemnode;    /**< Longest edge midpoint. */
  int8_t              ledge[2];     /**< Corners of longest edge. */
  int8_t              sibid;        /**< Number among its siblings. */
  int8_t              index;        /**< Sequence number in array. */
  int8_t              level;        /**< Depth in elementary forest. */
  p8est_tnodes_simplex_key_t key;   /**< Sort key in forest. */
};

typedef struct p8est_tnodes_context
{
  /** A node code for each elementary node number.
   * The code is an 8-bit number where the four high bits
   * contain the node's codimension, i. e. 0 for the volume,
   * and the low bits the number of the cube entity it touches.
   */
  p8est_tnodes_eind_code_t *eind_code;

  /** The elementary forest contains a simplicial refinement.
   * The root simplices cover the reference cube, and we include
   * all nodes subdivided by longest edge bisection down to
   * a similar refinement of half the size.
   * The type is \ref p8est_tnodes_simplex_t.
   */
  sc_array_t         *eforest;
}
p8est_tnodes_context_t;

/* Produce a context for simplicial element subdivision. */
p8est_tnodes_context_t *p4est_tnodes_context_new (void);

/* Free a context for simplicial element subdivision. */
void                p8est_tnodes_context_destroy (p8est_tnodes_context_t *
                                                  econ);

/* End the new code to produce a systematic recursive subdivision. */

/** Integer type to store the bits of an element configuration. */
typedef uint32_t    p8est_tnodes_config_t;

/** Lookup table structure defining a conforming tetrahedral mesh.
 *
 * The \a lnodes member encodes process-relavent corners, edges and faces.
 * Tetrahedron-shaped volume and corner entities are always included.
 * Can be created with or without including faces and/or edges as entities.
 * The members of \a lnodes are reinterpreted; cf. \ref p8est_lnodes.h :
 *  - degree is set to 0.
 *  - vnodes is the maxium number of nodes per element.
 */
typedef struct p8est_tnodes
{
  int                 full_style;       /**< Full style subdivision? */
  int                 with_faces;       /**< Include tetrahedron faces? */
  int                 with_edges;       /**< Include tetrahedron edges? */
  p8est_tnodes_config_t *configuration; /**< One entry per element. */
  p4est_gloidx_t      global_toffset;   /**< Global tetrahedron offset
                                             for the current process. */
  p4est_locidx_t     *local_tcount;     /**< Tetrahedron count per process
                                             (has mpisize entries). */
  p4est_locidx_t     *local_toffset;    /**< Tetrahedron offsets per local
                                             element and one beyond. */

  p8est_lnodes_t     *lnodes;   /**< Element and tetrahedron node data. */
  int                 lnodes_owned;     /**< Boolean: ownership of \a lnodes. */
  sc_array_t         *simplex_lnodes;   /**< Local nodes of local simplices.
                                             Each array entry holds 4 node
                                             indices of type p4est_locidx_t. */

  struct p8est_tnodes_private *pri;  /**< Private member not to access. */
}
p8est_tnodes_t;

/** Generate a conforming tetrahedron mesh from a 2:1 balance forest.
 * \param [in] p4est    Valid forest after 2:1 (at least face) balance.
 * \param [in] ghost    Ghost layer created from \b p4est.  Even with MPI,
 *                      it may be NULL to number the nodes purely locally.
 *                      In this case, nodes on a parallel boundary will be
 *                      considered as local for each touching process.
 *                      No shared nodes will be created.
 * \param [in] full_style   Half or full subdivision for unrefined elements.
 * \param [in] with_faces   If true, include each face of the tetrahedral
 *                          mesh as a node, otherwise ignore all faces.
 * \param [in] with_edges   If true, include each edge of the tetrahedral
 *                          mesh as a node, otherwise ignore all edges.
 * \return              Valid conforming tetrahedron mesh structure.
 */
p8est_tnodes_t     *p8est_tnodes_new (p8est_t * p4est,
                                      p8est_ghost_t * ghost, int full_style,
                                      int with_faces, int with_edges);

/** Generate a conforming tetrahedron mesh from a Q2 nodes structure.
 * \param [in] lnodes                   Valid nodes structure of degree 2.
 * \param [in] lnodes_take_ownership    Boolean: we will own \c lnodes.
 * \param [in] construction_flags       Currently must be 0.
 * \return                              Valid conforming tetrahedron mesh.
 *                     Each tetrahedron is strictly contained in one element
 *                     of the p8est hexahedral mesh underlying \c lnodes.
 *                     Each element contains from 4 to 48 tetrahedra.
 *                     The tetrahedra are right-handed with respect to the
 *                     tree coordinate system containing their element.
 */
p8est_tnodes_t     *p8est_tnodes_new_Q2 (p8est_lnodes_t * lnodes,
                                         int lnodes_take_ownership,
                                         int construction_flags);

/** Free the memory in a conforming tetrahedron mesh structure.
 * \param [in] tnodes      Memory is deallocated.  Do not use after return.
 */
void                p8est_tnodes_destroy (p8est_tnodes_t * tnodes);

SC_EXTERN_C_END;

#endif /* !P8EST_TNODES_H */
