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
 * This mesh is derived from a \ref p4est_lnodes structure.
 * We create additional lookup tables to define the mesh.
 *
 * \ingroup p4est
 */

#ifndef P4EST_TNODES_H
#define P4EST_TNODES_H

#include <p4est_geometry.h>

SC_EXTERN_C_BEGIN;

/** Flag values for tnodes construction. */
typedef enum        p4est_tnodes_flags
{
  /** The default flags have no bits set. */
  P4EST_TNODES_FLAGS_NONE = 0,

  /** Generate geometric coordinates for nodes on the tree boundary.
   * Since the \ref p4est_connectivity may be periodic, the same lnode
   * entry (see \ref p4est_lnodes) may be referenced from more than one
   * coordinate location.  If periodicity is not expected, this flag is not
   * needed.  Otherwise, setting it disambiguates the coordinates between
   * multiple instances for the same lnode entry.  This enables for example
   * the visualization of the periodic unit square as a factual square. */
  P4EST_TNODES_COORDS_SEPARATE = 0x01,

  /** Allocate and fill local_element_level array in \ref p4est_tnodes. */
  P4EST_TNODES_ELEMENT_LEVEL = 0x02,

  /** Allocate and fill local_simplex_level array in \ref p4est_tnodes. */
  P4EST_TNODES_SIMPLEX_LEVEL = 0x04,

  /** Allocate and fill both level arrays in \ref p4est_tnodes. */
  P4EST_TNODES_STORE_LEVELS = 0x06,

  /** All flags have all bits set. */
  P4EST_TNODES_FLAGS_ALL = -1
}
p4est_tnodes_flags_t;

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
 */
typedef struct p4est_tnodes
{
  p4est_gloidx_t      global_toffset;   /**< Global triangle offset
                                             for the current process. */
  p4est_gloidx_t      global_tcount;    /**< Global triangle count. */
  p4est_locidx_t     *local_tcount;     /**< Triangle count for each process
                                             (has mpisize entries). */

  /** Offsets into local triangles per element and one beyond. */
  p4est_locidx_t     *local_element_offset;
  /** The level of a p4est element applies to all simplices within.
   * Depending on the simplex construction, other elements may
   * overlap the same simplex, but no more than this element.
   * The tree reference volume of this element is 4**{-l}.
   * Array may be NULL if the information is not provided. */
  int8_t             *local_element_level;

  p4est_topidx_t      local_first_tree; /**< First local tree on process,
                                             -1 if process has no elements. */
  p4est_topidx_t      local_last_tree;  /**< Last local tree on process,
                                             -2 if process has no elements. */
  /** Offsets into local triangles, zero indexed from local_first_tree
   * to local_last_tree + 1 inclusive.  Length 1 on empty processes. */
  p4est_topidx_t     *local_tree_offset;

  sc_array_t         *simplex_level;    /**< Simplex refinement level l,
                                             reference volume 2**{-l} / 2.
                                             May be NULL if not provided. */
  sc_array_t         *simplices;        /**< Vertex indices of local
                                             simplices.  Each array entry
                                             holds 3 \ref p4est_locidx_t. */
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
  p4est_lnodes_t     *lnodes;   /**< Element and triangle node data. */
  int                 lnodes_owned;     /**< Boolean: ownership of \a lnodes. */
}
p4est_tnodes_t;

/** Generate a conforming triangle mesh from a Q2 nodes structure.
 * \param [in] p4est                    Forest underlying the mesh.
 * \param [in] lnodes                   Valid node structure of degree 2.
 *                                      Must be derived from the \c p4est.
 * \param [in] geom                     If NULL, we create tree relative
 *                                      reference coordinates in [0, 1]^2.
 *                                      Otherwise we apply \c geom.
 *                                      Any geometry might also be passed
 *                                      to the VTK output routine, but
 *                                      shall not given in both places.
 * \param [in] construction_flags       Currently must be 0.
 * \return                              Valid conforming triangle mesh.
 *                     Each triangle is strictly contained in one element
 *                     of the p4est quadrilateral mesh underlying \c lnodes.
 *                     Each element contains between 4 and 8 triangles.
 *                     The triangles are right-handed with respect to the
 *                     tree coordinate system containing their element.
 */
p4est_tnodes_t     *p4est_tnodes_new_Q2_P1 (p4est_t *p4est,
                                            p4est_lnodes_t *lnodes,
                                            p4est_geometry_t *geom,
                                            int construction_flags);

/** Free the memory in a conforming triangle mesh structure.
 * \param [in] tnodes      Memory is deallocated.  Do not use after return.
 */
void                p4est_tnodes_destroy (p4est_tnodes_t * tnodes);

SC_EXTERN_C_END;

#endif /* !P4EST_TNODES_H */
