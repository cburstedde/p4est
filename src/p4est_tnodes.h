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
 * This mesh is represented by augmenting the \ref p4est_lnodes structure.
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
  P4EST_TNODES_COORDS_SEPARATE = 0x01
}
p4est_tnodes_flags_t;

/** Integer type to store the bits of an element configuration. */
typedef uint8_t     p4est_tnodes_config_t;

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
 * The below configuration convention does not apply to it.
 *
 * The \a lnodes member encodes the process-relavent corners and faces.
 * Triangle-shaped volume and corner entities are always included.
 * It can be created with or without including faces as mesh entities.
 * The members of \a lnodes are reinterpreted; cf. \ref p4est_lnodes.h :
 *  - degree is set to 0.
 *  - vnodes is the maxium number of nodes per element, 9 (corners only)
 *    or 25 (with faces).  Each element gets this amount of memory in the
 *    \a element_nodes member.  Unused positions are set to -1.
 *    The position of the nodes wrt. the element are as follows:
 *
 *          y        3
 *          +----------------+
 *          |  2 23  8 24  3 |
 *          | 15 11 22 12 18 |
 *        0 |  5 14  4 17  6 | 1
 *          | 13  9 20 10 16 |
 *          |  0 19  7 21  1 |
 *          +----------------+-> x
                     2
 *
 *    The nodes 0--3 are always triangle corner nodes.
 *    The nodes 9--24 are always triangle face nodes.
 *    The nodes 4--8 may be either depending on the configuration.
 *    The face numbers are displayed on the outside for completeness.
 *
 * There are 16 configurations for splitting an element into triangles.
 * Each configuration is encoded by one on bit for each split face
 * counted from the right (set the least significant bit for face 0).
 * Configuration 0 may have two additional bits set: bit 5 if there is
 * a full split of the element into four triangles, and no bit 5 when
 * there is a half split.  The half split sets bit 4 for child 1 and 2.
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
  /** The level of an element applies to all simplices within. */
  int8_t             *local_element_level;

  p4est_topidx_t      local_first_tree; /**< First local tree on process,
                                             -1 if process has no elements. */
  p4est_topidx_t      local_last_tree;  /**< Last local tree on process,
                                             -2 if process has no elements. */
  /** Offsets into local triangles, zero indexed from local_first_tree
   * to local_last_tree + 1 inclusive.  Length 1 on empty processes. */
  p4est_topidx_t     *local_tree_offset;

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
  p4est_tnodes_config_t *configuration; /**< One entry per element. */
  p4est_lnodes_t     *lnodes;   /**< Element and triangle node data. */
  int                 lnodes_owned;     /**< Boolean: ownership of \a lnodes. */
  struct p4est_tnodes_private *pri;     /**< Private member not to access. */
}
p4est_tnodes_t;

/** The iterator state to go through the triangles in a \ref p4est_tnodes_t.
 * The traversal is process-local and not collective, all members are local.
 */
typedef struct p4est_tnodes_iter
{
  /* context members */
  struct p4est_tnodes_iter_private *pri;   /**< Private member not to access. */
  p4est_t            *p4est;        /**< The forest backing the mesh. */
  p4est_tnodes_t     *tnodes;       /**< The triangle mesh structure. */

  /* define current triangle */
  p4est_topidx_t      which_tree;   /**< The current tree number. */
  p4est_locidx_t      which_quad;   /**< The local quadrant number is
                                         relative to process, not tree. */
  p4est_quadrant_t   *quadrant;     /**< The current local quadrant. */
  p4est_locidx_t      triangle;     /**< The current local triangle. */

  /* properties of current triangle */
  int                 corner_nodes[3];      /**< Element node number in 0..8. */
  int                 face_nodes[3];        /**< Element node number in 4..24.
                                                 If no faces nodes are stored,
                                                 these are all set to -1. */
}
p4est_tnodes_iter_t;

/** There are 16 elementary triangles in a quadrant.
 * We list them in order of ascending configurations.
 * We list the corner nodes before the face nodes.
 * The triangles begin with the lowest numbered quadrant node
 * and proceed right-handed.  The faces run right-handed, too,
 * where the first face touches the first and second corner.
 */
extern const int    p4est_tnodes_triangle_nodes[16][6];

/** For each distinct configuration, the number of corner and face
 * nodes and then the number of triangles in an element.
 * They are indexed by running number and then by codimension
 * in the sequence of corner, face, volume.
 */
extern const int    p4est_tnodes_lookup_counts[6][3];

/** For each configuration, lookup a distinct combination of number
 * of corner and face nodes in \ref p4est_tnodes_lookup_counts.
 * The configurations are indexed bitwise by the first five bits
 * from 0 to 16 inclusive plus one more, where configuration 0's
 * three subconfigurations have indices 0, 16, 17, and the other
 * configurations are indexed with their true numbers 1--15.
 */
extern const int    p4est_tnodes_config_lookup[18];

/** For each configuration the list of corner nodes padded with -1. */
extern const int    p4est_tnodes_config_corners[18][9];

/** For each configuration the list of face nodes padded with -1. */
extern const int    p4est_tnodes_config_faces[18][16];

/** For each configuration the list of triangles in the quadrant. */
extern const int    p4est_tnodes_config_triangles[18][8];

/** Generate a conforming triangle mesh from a 2:1 balance forest.
 * \param [in] p4est    Valid forest after 2:1 (at least face) balance.
 * \param [in] ghost    Ghost layer created from \b p4est.  Even with MPI,
 *                      it may be NULL to number the nodes purely locally.
 *                      In this case, nodes on a parallel boundary will be
 *                      considered as local for each touching process.
 *                      No shared nodes will be created.
 * \param [in] full_style   Half or full subdivision for unrefined elements.
 * \param [in] with_faces   If false, only number triangles and corner nodes.
 *                          Otherwise, include each triangle face as a node.
 * \return              Valid conforming triangle mesh structure.
 */
p4est_tnodes_t     *p4est_tnodes_new (p4est_t * p4est,
                                      p4est_ghost_t * ghost,
                                      int full_style, int with_faces);

/** Generate a conforming triangle mesh from a Q2 nodes structure.
 * \param [in] p4est                    Forest underlying the mesh.
 * \param [in] geom                     If NULL, we create tree relative
 *                                      reference coordinates in [0, 1]^2.
 *                                      Otherwise we apply \c geom.
 *                                      Any geometry should either be passed
 *                                      here, or to the VTK output routine,
 *                                      but not given in both places.
 * \param [in] lnodes                   Valid node structure of degree 2.
 *                                      Must be derived from \c p4est.
 * \param [in] lnodes_take_ownership    Boolean: we will own \c lnodes.
 * \param [in] construction_flags       Currently must be 0.
 * \return                              Valid conforming triangle mesh.
 *                     Each triangle is strictly contained in one element
 *                     of the p4est quadrilateral mesh underlying \c lnodes.
 *                     Each element contains between 4 and 8 triangles.
 *                     The triangles are right-handed with respect to the
 *                     tree coordinate system containing their element.
 */
p4est_tnodes_t     *p4est_tnodes_new_Q2_P1 (p4est_t *p4est,
                                            p4est_geometry_t *geom,
                                            p4est_lnodes_t *lnodes,
                                            int lnodes_take_ownership,
                                            int construction_flags);

/** Free the memory in a conforming triangle mesh structure.
 * \param [in] tnodes      Memory is deallocated.  Do not use after return.
 */
void                p4est_tnodes_destroy (p4est_tnodes_t * tnodes);

/** Create an iterator through the triangles in a tnodes structure.
 * The iterator may be used in a for loop as follows:
 *
 *     for (iter = p4est_tnodes_iter_new (p4est, tnodes);
 *          iter != NULL; p4est_tnodes_iter_next (&iter))
 *
 * \param [in] p4est    The forest is needed to access its elements,
 *                      which contain the triangles to iterate through.
 * \param [in] tnodes   Valid tnodes structure created from the \a p4est.
 * \return              Iterator pointing to the first triangle in order
 *                      or NULL if the local process has no triangles.
 */
p4est_tnodes_iter_t *p4est_tnodes_iter_new (p4est_t * p4est,
                                            p4est_tnodes_t * tnodes);

/** Advance to next triangle in a \ref p4est_tnodes_iter_t iterator.
 * This function must no longer be called on a NULL iterator.
 * If it is called on the last triangle, the iterator becomes NULL.
 * \param [in, out] piter       This pointer must not be NULL.
 *                              It must point to an iterator that is
 *                              also not NULL.  On output it becomes NULL
 *                              when called on the last triangle.
 *                              Otherwise its state advances to the next.
 */
void                p4est_tnodes_iter_next (p4est_tnodes_iter_t ** piter);

SC_EXTERN_C_END;

#endif /* !P4EST_TNODES_H */
