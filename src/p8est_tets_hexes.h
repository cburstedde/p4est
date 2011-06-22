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

#ifndef P8EST_TETS_HEXES_H
#define P8EST_TETS_HEXES_H

#include <p4est_base.h>

typedef struct p8est_tetgen
{
  /** The node array contains a triplet of double coordinates per node. */
  sc_array_t         *nodes;

  /** The tet array contains a quartet of p4est_topidx_t nodes per tet. */
  sc_array_t         *tets;

  /** The element_attributes array can contain one int attribute per tet. */
  sc_array_t         *tet_attributes;
}
p8est_tetgen_t;

/** Read nodes from a tetgen .node file.
 * \param [in] nodefile     Name of file in tetgen .node format.
 * \return                  An array with three double coordinates per node,
 *                          or NULL on file error.
 */
sc_array_t         *p8est_tetgen_read_node (const char *nodefile);

/** Read tetrahedra from a tetgen .ele file.
 * \param [in] elefile          Name of file in tetgen .ele format.
 * \param [in] num_nodes        If nonnegative, (exclusive) upper node number.
 * \param [in,out] attributes   If not NULL, an array will be created
 *                              if the .ele file contains attributes.
 * \return                      An array with four p4est_topidx_t nodes
 *                              per tet, or NULL on file error.
 */
sc_array_t         *p8est_tetgen_read_ele (const char *elefile,
                                           p4est_topidx_t num_nodes,
                                           sc_array_t ** attributes);

/** Read element and node information from a tetgen base name.
 * The names for element and node files are derived from base name by suffix.
 * \param [in] tetgenbasename   Base name for tetgen files (without suffix).
 * \return                      A populated p8est_tetgen_t structure
 *                              or NULL on file error.
 */
p8est_tetgen_t     *p8est_tetgen_read (const char *tetgenbasename);

/** Destroy all memory associated with a p8est_tetgen_t structure.
 * \param [in] ptg          Allocated p8est_tetgen_t structure.
 */
void                p8est_tetgen_destroy (p8est_tetgen_t * ptg);

/** Change all left-handed tetrahedra to right-handed ones.
 * \param [in,out] ptg      Structure with node and tet information.
 * \return                  The number of tets that were flipped.
 */
p4est_topidx_t      p8est_tetgen_make_righthanded (p8est_tetgen_t * ptg);

#endif /* !P8EST_TETS_HEXES */
