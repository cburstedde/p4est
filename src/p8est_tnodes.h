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

#include <p8est_geometry.h>

SC_EXTERN_C_BEGIN;

/** Flag values for tnodes construction. */
typedef enum        p8est_tnodes_flags
{
  /** The default flags have no bits set. */
  P8EST_TNODES_FLAGS_NONE = 0,

  /** Generate geometric coordinates for nodes on the tree boundary.
   * Since the \ref p8est_connectivity may be periodic, the same lnode
   * entry (see \ref p8est_lnodes) may be referenced from more than one
   * coordinate location.  If periodicity is not expected, and the geometry
   * is continuous across tree boundaries, this flag is not needed.
   * Otherwise, setting it disambiguates the coordinates between multiple
   * instances for the same lnode entry.  This enables for example the
   * visualization of the periodic unit square as a factual square, or the
   * setup of geometries that are artificially mapped for display. */
  P8EST_TNODES_COORDS_SEPARATE = 0x01,

  /** Allocate and fill local_element_level array in \ref p8est_tnodes. */
  P8EST_TNODES_ELEMENT_LEVEL = 0x02,

  /** Allocate and fill local_simplex_level array in \ref p8est_tnodes. */
  P8EST_TNODES_SIMPLEX_LEVEL = 0x04,

  /** Allocate and fill both level arrays in \ref p8est_tnodes. */
  P8EST_TNODES_STORE_LEVELS = 0x06,

  /** All flags have all bits set. */
  P8EST_TNODES_FLAGS_ALL = -1
}
p8est_tnodes_flags_t;

/** Lookup table structure defining a conforming tetrahedral mesh.
 */
typedef struct p8est_tnodes
{
  p4est_gloidx_t      global_toffset;   /**< Global tetrahedron offset
                                             for the current process. */
  p4est_gloidx_t      global_tcount;    /**< Global tetrahedron count. */
  p4est_locidx_t     *local_tcount;     /**< Tetrahedron count per process
                                             (has mpisize entries). */

  /** Offsets into local triangles per element and one beyond. */
  p4est_locidx_t     *local_element_offset;
  /** The level of a p8est element applies to all simplices within.
   * Depending on the simplex construction, other elements may
   * overlap the same simplex, but no more than this element.
   * The tree reference volume of this element is 8**{-l}.
   * Array may be NULL if the information is not provided. */
  int8_t             *local_element_level;

  int                 local_first_child;        /**< First child id on
                                             process, or -1 if empty. */
  p4est_topidx_t      local_first_tree; /**< First local tree on process,
                                             -1 if process has no elements. */
  p4est_topidx_t      local_last_tree;  /**< Last local tree on process,
                                             -2 if process has no elements. */
  /** Offsets into local triangles, zero indexed from local_first_tree
   * to local_last_tree + 1 inclusive.  Length 1 on empty processes. */
  p4est_topidx_t     *local_tree_offset;

  sc_array_t         *simplex_level;    /**< Simplex refinement level l,
                                             reference volume 2**{-l} / 6.
                                             May be NULL if not provided. */
  sc_array_t         *simplices;        /**< Vertex indices of local
                                             simplices.  Each array entry
                                             holds 4 \ref p4est_locidx_t. */
  sc_array_t         *coord_to_lnode;   /**< This pointer may be NULL, in
                                             which case \c simplices indexes
                                             into both the local nodes from
                                             \ref p4est_lnodes and the \c
                                             coordinates below.  Otherwise,
                                             the simplex array indexes into \c
                                             coordinates, and this array maps
                                             a coordinate to its local node. */
  sc_array_t         *coordinates;      /**< Each entry is a double 3-tuple. */

  /* deprecated members below */
  p8est_lnodes_t     *lnodes;   /**< Element and tetrahedron node data. */
  int                 lnodes_owned;     /**< Boolean: ownership of \a lnodes. */
}
p8est_tnodes_t;

/** Generate a conforming tetrahedron mesh from a Q2 nodes structure.
 * \param [in] p8est                    Forest underlying the mesh.
 * \param [in] lnodes                   Valid node structure of degree 2.
 *                                      Must be derived from the \c p8est.
 * \param [in] geom                     If NULL, we create tree relative
 *                                      reference coordinates in [0, 1]^3.
 *                                      Otherwise we apply \c geom.
 *                                      Any geometry might also be passed
 *                                      to the VTK output routine, but
 *                                      shall not given in both places.
 * \param [in] construction_flags       Currently must be 0.
 * \return                              Valid conforming tetrahedron mesh.
 *                     Each tetrahedron is strictly contained in one element
 *                     of the p8est hexahedral mesh underlying \c lnodes.
 *                     Each element contains from 4 to 48 tetrahedra.
 *                     The tetrahedra are right-handed with respect to the
 *                     tree coordinate system containing their element.
 */
p8est_tnodes_t     *p8est_tnodes_new_Q2_P1 (p8est_t *p8est,
                                            p8est_lnodes_t *lnodes,
                                            p8est_geometry_t *geom,
                                            int construction_flags);

/** Free the memory in a conforming tetrahedron mesh structure.
 * \param [in] tnodes      Memory is deallocated.  Do not use after return.
 */
void                p8est_tnodes_destroy (p8est_tnodes_t * tnodes);

SC_EXTERN_C_END;

#endif /* !P8EST_TNODES_H */
