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

#include <p4est_lnodes.h>

SC_EXTERN_C_BEGIN;

/** Lookup table structure defining a conforming triangle mesh.
 *
 * Trying to conform to latest status of paper repository:
 *
 *     6bc25f04355eef8d73ec53bdcf6f5915a5748559
 *     711e76748721665bdebb3d5f0bfd53dbd1702a8e
 *
 * The \a lnodes member encodes the process-relavent corners and faces.
 * Triangle-shaped element and corner entities are always included.
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
 *    The nodes 4--8 may be either.
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
  int                 full_style;       /**< Full style subdivision? */
  int                 with_faces;       /**< Include triangle faces? */
  uint8_t            *configuration;    /**< One entry per element. */
  p4est_locidx_t     *local_toffset;    /**< Triangle offsets per local
                                             element and one beyond. */
  p4est_gloidx_t     *global_toffset;   /**< Global triangle offsets.
                                             Has mpisize + 1 entries. */
  p4est_lnodes_t     *lnodes;   /**< Element and triangle node data. */
}
p4est_tnodes_t;

/** For each configuration, the number of corner and face nodes.
 * The configurations are indexed bitwise by the first five bits
 * from 0 to 16 inclusive plus one more, where configuration 0's
 * three subconfigurations have indices 0, 16, 17, and the other
 * configurations are indexed with their true numbers 1--15.
 */
extern const int p4est_tnodes_config_count[18][2];

/** For each configuration the list of corner nodes padded with -1. */
extern const int p4est_tnodes_config_corners[18][9];

/** For each configuration the list of face nodes padded with -1. */
extern const int p4est_tnodes_config_faces[18][16];

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

/** Free the memory in a conforming triangle mesh structure.
 * \param [in] tnodes      Memory is deallocated.  Do not use after return.
 */
void                p4est_tnodes_destroy (p4est_tnodes_t * tnodes);

SC_EXTERN_C_END;

#endif /* !P4EST_TNODES_H */
