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
  int                 Qdegree;          /**< Degree of original lnodes. */
  int                 Pdegree;          /**< Degree of simplex space. */

  /* Counts of simplices */
  p4est_gloidx_t      global_toffset;   /**< Global tetrahedron offset
                                             for the current process. */
  p4est_gloidx_t      global_tcount;    /**< Global tetrahedron count. */
  p4est_locidx_t     *local_tcount;     /**< Tetrahedron count per process
                                             (has mpisize entries). */

  /** Offsets into local triangles per element and one beyond. */
  p4est_locidx_t     *local_element_offset;

  /** Offsets into local triangles, zero indexed from local_first_tree
   * to local_last_tree + 1 inclusive.  Length 1 on empty processes. */
  p4est_topidx_t     *local_tree_offset;

  sc_array_t         *simplex_level;    /**< Simplex refinement level l,
                                             reference volume 2**{-l} / 6. */
  sc_array_t         *simplices;        /**< Vertex indices of local
                                             simplices.  Each array entry
                                             holds 4 int8_t variables. */
  /** For each element, one or eight bytes of flag bits.
   * For degree 1, there is one byte per local element storing 6 bits.
   * For degree 2, there are eight bytes per local element of this kind.
   * A bit is set if the corresponding elementary simplex exists.
   * Simplices may be omitted at a hanging face or edge.
   */
  sc_array_t         *element_bits;
}
p8est_tnodes_t;

/** Compute the number of the parent simplex that contains this one.
 * \param [in] p    Number of parent element in [0, 4).
 * \param [in] c    Number of child element in [0, 8).
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

/** Generate a conforming tetrahedron mesh from a Q2 lnodes structure.
 * Obsolete code that provides calls for generating node coordinates.
 * \param [in] p4est                    Forest underlying the mesh.
 * \param [in] lnodes                   Valid node structure of degree 2.
 *                                      Must be derived from the \c p8est.
 * \return                              Valid conforming tetrahedron mesh.
 *                     Each tetrahedron is strictly contained in one element
 *                     of the p8est hexahedral mesh underlying \c lnodes.
 *                     Each element contains from 4 to 48 tetrahedra.
 *                     The tetrahedra are right-handed with respect to the
 *                     tree coordinate system containing their element.
 */
p8est_tnodes_t     *p8est_tnodes_new_Q2_P1_exp (p8est_t *p4est,
                                                p8est_lnodes_t *lnodes);

/** Generate a conforming triangle mesh from a Q1 lnodes structure.
 * \param [in] p4est    Forest underlying the mesh.
 *                      It must not contain any root-level elements:
 *                      to avoid this, it should be refined a priori.
 * \param [in] lnodes   Valid node structure of degree 1.
 *                      Must be derived from the \c p4est.
 * \return              Valid conforming triangle mesh.
 *                      Some fields are ignored in view of eventual removal.
 *                      Each triangle overlaps one or more elements.  It is
 *                      assigned to exactly one of their owner processes.
 *                      The triangles are right-handed with respect to the tree
 *                      coordinate system containing their element.
 */
p8est_tnodes_t     *p8est_tnodes_new_Q1_P1 (p8est_t *p4est,
                                            p8est_lnodes_t *lnodes);

/** Generate a conforming triangle mesh from a Q2 lnodes structure.
 * \param [in] p4est    Forest underlying the mesh.
 *                      It must not contain any elements at P4EST_QMAXLEVEL.
 *                      To avoid this, it should not be refined that deep.
 * \param [in] lnodes   Valid node structure of degree 2.
 *                      Must be derived from the \c p4est.
 * \return              Valid conforming triangle mesh.
 *                      Some fields are ignored in view of eventual removal.
 *                      Each triangle is contained in exactly one processes.
 *                      The triangles are right-handed with respect to the
 *                      tree coordinate system containing their element.
 */
p8est_tnodes_t     *p8est_tnodes_new_Q2_P1 (p8est_t *p4est,
                                            p8est_lnodes_t *lnodes);

/** Free the memory in a conforming tetrahedron mesh structure.
 * \param [in] tnodes      Memory is deallocated.  Do not use after return.
 */
void                p8est_tnodes_destroy (p8est_tnodes_t * tnodes);

SC_EXTERN_C_END;

#endif /* !P8EST_TNODES_H */
