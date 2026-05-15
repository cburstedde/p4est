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

/** \file p4est_tnodes.h
 *
 * Generate a conforming triangle mesh from a 2:1 balanced p4est.
 * This mesh is based on a given \ref p4est_lnodes structure.
 * Additional lookup tables define the simplex mesh.
 *
 * \ingroup p4est
 */

#ifndef P4EST_TNODES_H
#define P4EST_TNODES_H

#include <p4est_lnodes.h>

SC_EXTERN_C_BEGIN;

/** Number of corners of a simplex. */
#define P4EST_TNODES_SIMPLEX_CORNERS 3

/** Number of coarse simplices in a cube. */
#define P4EST_TNODES_CUBE_SIMPLICES 2

/** Lookup table structure defining a conforming triangle mesh.
 *
 * Trying to conform to latest status of paper repository:
 *
 *     7c96f3bbefad364e3fa657272bca757d13d82a82
 *     d43b2e54f939b186ef765c65638fde2fe792aa55
 *     6bc25f04355eef8d73ec53bdcf6f5915a5748559
 *     711e76748721665bdebb3d5f0bfd53dbd1702a8e
 *
 * In the meantime, we have added a Q2 recursive bisection construction
 * that works in both 2D and 3D and appears to be functional.
 * Rewriting the paper from scratch in its own repository:
 *
 *     6db206b36bcbc3f602bfce39ef03756b9cc2c845
 */
typedef struct p4est_tnodes
{
  int                 mpisize;          /**< Number of parallel processes. */
  int                 Qdegree;          /**< Degree of original lnodes. */
  int                 Pdegree;          /**< Degree of simplex space. */

  /* Counts of simplices */
  p4est_gloidx_t      global_toffset;   /**< Global triangle offset
                                             for the current process. */
  p4est_gloidx_t      global_tcount;    /**< Global triangle count. */
  p4est_locidx_t     *local_tcount;     /**< Triangle count for each process
                                             (has mpisize entries). */

  /** Offsets into local triangles per element and one beyond.
   * The number of local elements is the array length of \c element_bits. */
  p4est_locidx_t     *local_element_offset;

  /** Offsets into local triangles, zero indexed from local_first_tree
   * to local_last_tree + 1 inclusive.  Length 1 on empty processes. */
  p4est_topidx_t     *local_tree_offset;

  sc_array_t         *simplex_level;    /**< Simplex refinement level l,
                                             reference volume 2**{-l} / 2. */
  sc_array_t         *simplices;        /**< Vertex indices of local
                                             simplices.  Each array entry
                                             holds 3 int8_t variables. */
  /** For each local element, one or four bytes of flag bits.
   * For degree 1, there is one byte per local element storing 2 bits.
   * For degree 2, there are four bytes per local element of this kind.
   * A bit is set if the corresponding elementary simplex exists.
   * Simplices may be omitted at a hanging face.
   */
  sc_array_t         *element_bits;

}
p4est_tnodes_t;

/** The number of the congruent simplex that starts at the antipode.
 * There is exactly one simplex that starts at the antipode corner
 * for every one starting at the anchor corner, in reverse order.
 */
extern const int    p4est_tnodes_simplex_reverse[2];

/** Compute the number of the parent simplex that contains this one.
 * \param [in] p    Child id of parent element in [0, 4).
 * \param [in] c    Child id of child element in [0, 4).
 * \param [in] k    Number of simplex in [0, 2) within child.
 * return           Number of simplex in [0, 2) within parent.
 */
int                 p4est_tnodes_simplex_parent (int p, int c, int k);

/** Verify that a given child simplex is contained in a parent
 * \param [in] p    Number of parent element in [0, 4).
 * \param [in] kp   Number of simplex in [0, 2) within parent.
 * \param [in] c    Number of child element in [0, 4).
 * \param [in] kc   Number of simplex in [0, 2) within child.
 * \return          True if contained, false if not.
 */
int                 p4est_tnodes_simplex_parent_is_valid
  (int p, int kp, int c, int kc);

/** Based on the face code of an element, calculate contained simplices.
 * \param [in] fc   Valid face code as defined in \ref p4est_lnodes.h.
 * \return          Element simplex count between 6 and 8 inclusive.
 */
int                 p4est_tnodes_quadrant_Q2_simplices
  (p4est_lnodes_code_t fc);

/** Generate a conforming triangle mesh from a Q2 lnodes structure.
 * Obsolete code that provides calls for generating node coordinates.
 * \param [in] p4est                    Forest underlying the mesh.
 * \param [in] lnodes                   Valid node structure of degree 2.
 *                                      Must be derived from the \c p4est.
 * \return                              Valid conforming triangle mesh.
 *                     Each triangle is strictly contained in one element
 *                     of the p4est quadrilateral mesh underlying \c lnodes.
 *                     Each element contains between 4 and 8 triangles.
 *                     The triangles are right-handed with respect to the
 *                     tree coordinate system containing their element.
 */
p4est_tnodes_t     *p4est_tnodes_new_Q2_P1_exp (p4est_t *p4est,
                                                p4est_lnodes_t *lnodes);

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
p4est_tnodes_t     *p4est_tnodes_new_Q1_P1 (p4est_t *p4est,
                                            p4est_lnodes_t *lnodes);

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
p4est_tnodes_t     *p4est_tnodes_new_Q2_P1 (p4est_t *p4est,
                                            p4est_lnodes_t *lnodes);

/** Calculate memory allocated in a tnodes structure.
 * \param [in] tnodes   Valid tnodes structure.
 * \return              Total memory allocation in bytes.
 */
size_t              p4est_tnodes_memory_used (p4est_tnodes_t *tnodes);

/** Free the memory in a conforming triangle mesh structure.
 * \param [in] tnodes      Memory is deallocated.  Do not use after return.
 */
void                p4est_tnodes_destroy (p4est_tnodes_t * tnodes);

SC_EXTERN_C_END;

#endif /* !P4EST_TNODES_H */
