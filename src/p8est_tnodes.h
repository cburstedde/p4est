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
 * This mesh is based on a given \ref p8est_lnodes structure.
 * Additional lookup tables define the simplex mesh.
 *
 * \ingroup p8est
 */

#ifndef P8EST_TNODES_H
#define P8EST_TNODES_H

#include <p8est_lnodes.h>

SC_EXTERN_C_BEGIN;

/** Number of corners of a simplex. */
#define P8EST_TNODES_SIMPLEX_CORNERS 4

/** Number of coarse simplices in a cube. */
#define P8EST_TNODES_CUBE_SIMPLICES 6

/** Lookup table structure defining a conforming tetrahedral mesh. */
typedef struct p8est_tnodes
{
  int                 mpisize;          /**< Number of parallel processes. */
  int                 Qdegree;          /**< Degree of original lnodes. */
  int                 Pdegree;          /**< Degree of simplex space. */

  /* Counts of simplices */
  p4est_gloidx_t      global_toffset;   /**< Global tetrahedron offset
                                             for the current process. */
  p4est_gloidx_t      global_tcount;    /**< Global tetrahedron count. */
  p4est_locidx_t     *local_tcount;     /**< Tetrahedron count per process
                                             (has mpisize entries). */

  /** Offsets into local tetrahedra per element and one beyond.
   * The number of local elements is the array length of \c element_bits. */
  p4est_locidx_t     *local_element_offset;

  sc_array_t         *simplex_level;    /**< Simplex level l as an int8_t,
                                             reference volume 2**{-l} / 6. */
  sc_array_t         *simplices;        /**< Vertex indices of local
                                             simplices.  Each array entry
                                             holds 4 int8_t variables. */
  /** For each local element, one or eight bytes of flag bits.
   * For degree 1, there is one byte per local element storing 6 bits.
   * For degree 2, there are eight bytes per local element of this kind.
   * A bit is set if the corresponding elementary simplex exists.
   * Simplices may be omitted at a hanging face or edge.
   */
  sc_array_t         *element_bits;
}
p8est_tnodes_t;

/** For every simplex number the normal direction of its outside face.
 * The normal is for every k = 2i + j, with i in [0, 3), j in [0, 2).
 * The edge direction is i, and j indexes the normals of the two
 * faces touching an edge of this direction in ascending order.
 * The result is a normal direction in [0, 3).
 */
extern const int    p8est_tnodes_face_normal[6];

/** The number of the congruent simplex that starts at the antipode.
 * There is exactly one simplex that starts at the antipode corner
 * for every one starting at the anchor corner, in reverse order.
 */
extern const int    p8est_tnodes_simplex_reverse[6];

/** The simplex that shares the face direction and swaps the edge.
 */
extern const int    p8est_tnodes_simplex_edgeswap[6];

/** Compute the number of the parent simplex that contains this one.
 * \param [in] p    Child id of parent element in [0, 8).
 * \param [in] c    Child id of child element in [0, 8).
 * \param [in] k    Number of simplex in [0, 6) within child.
 * return           Number of simplex in [0, 6) within parent.
 */
int                 p4est_tnodes_simplex_parent (int p, int c, int k);

/** Verify that a given child simplex is contained in a parent
 * \param [in] p    Number of parent element in [0, 8).
 * \param [in] kp   Number of simplex in [0, 6) within parent.
 * \param [in] c    Number of child element in [0, 8).
 * \param [in] kc   Number of simplex in [0, 6) within child.
 * \return          True if contained, false if not.
 */
int                 p8est_tnodes_simplex_parent_is_valid
  (int p, int kp, int c, int kc);

/** Based on a face code and the parent's id, calculate simplex count.
 * \param [in] fc   Valid face code as defined in \ref p8est_lnodes.h.
 * \param [in] pc   Child id of quadrant's parent.
 * \return          Element simplex count between 2 and 6 inclusive.
 */
int                 p8est_tnodes_quadrant_Q1_simplices
  (p8est_lnodes_code_t fc, int pc);

/** Based on the face code of an element, calculate contained simplices.
 * \param [in] fc   Valid face code as defined in \ref p8est_lnodes.h.
 * \return          Element simplex count between 24 and 48 inclusive.
 */
int                 p8est_tnodes_quadrant_Q2_simplices
  (p8est_lnodes_code_t fc);

/** Generate a conforming tetrahedral mesh from a Q1 lnodes structure.
 * \param [in] p4est    Forest underlying the mesh.
 *                      It must not contain any root-level elements:
 *                      to avoid this, it should be refined a priori.
 * \param [in] lnodes   Valid node structure of degree 1.
 *                      Must be derived from the \c p4est.
 * \return              Valid conforming tetrahedral mesh.
 *                      Each tetrahedron overlaps one or more elements.  It
 *                      is assigned to exactly one of their owner processes.
 *                      The tetrahedra are right-handed with respect to the
 *                      tree coordinate system containing their element.
 */
p8est_tnodes_t     *p8est_tnodes_new_Q1_P1 (p8est_t *p4est,
                                            p8est_lnodes_t *lnodes);

/** Generate a conforming tetrahedral mesh from a Q2 lnodes structure.
 * \param [in] p4est    Forest underlying the mesh.
 *                      It must not contain any elements at P4EST_QMAXLEVEL.
 *                      To avoid this, it should not be refined that deep.
 * \param [in] lnodes   Valid node structure of degree 2.
 *                      Must be derived from the \c p4est.
 * \return              Valid conforming tetrahedral mesh.
 *                      Each tetrahedron is contained in exactly one process.
 *                      The tetrahedra are right-handed with respect to the
 *                      tree coordinate system containing their element.
 */
p8est_tnodes_t     *p8est_tnodes_new_Q2_P1 (p8est_t *p4est,
                                            p8est_lnodes_t *lnodes);

/** Generate a conforming tetrahedron mesh from a Q2 lnodes structure.
 *
 * This function uses a method distinct from \ref p8est_tnodes_new_Q2_P1
 * but produces an identical result.  It is useful for verification.
 *
 * \param [in] p4est    Forest underlying the mesh.
 *                      It must not contain any elements at P4EST_QMAXLEVEL.
 *                      To avoid this, it should not be refined that deep.
 * \param [in] lnodes   Valid node structure of degree 2.
 *                      Must be derived from the \c p4est.
 * \return              Valid conforming tetrahedral mesh.
 *                      Each tetrahedron is contained in exactly one process.
 *                      The tetrahedra are right-handed with respect to the
 *                      tree coordinate system containing their element.
 */
p8est_tnodes_t     *p8est_tnodes_new_Q2_P1_ind (p8est_t *p4est,
                                                p8est_lnodes_t *lnodes);

/** Calculate memory allocated in a tnodes structure.
 * \param [in] tnodes   Valid tnodes structure.
 * \return              Total memory allocation in bytes.
 */
size_t              p8est_tnodes_memory_used (p8est_tnodes_t *tnodes);

/** Free the memory in a conforming tetrahedron mesh structure.
 * \param [in] tnodes      Memory is deallocated.  Do not use after return.
 */
void                p8est_tnodes_destroy (p8est_tnodes_t * tnodes);

SC_EXTERN_C_END;

#endif /* !P8EST_TNODES_H */
