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

#ifndef P4EST_PLEX_H
#define P4EST_PLEX_H

#include <p4est.h>

SC_EXTERN_C_BEGIN;

/** Create the data necessary to create a PETsc DMPLEX representation of a
 * forest.  The forest must be at least face balanced (see p4est_balance()).
 * See test/test_plex2.c for example usage.
 *
 * All arrays should be initialized to hold sizeof (p4est_locidx_t), except
 * for \a out_remotes, which should be initialized to hold
 * (2 * sizeof (p4est_locidx_t)).
 *
 * \param[in]     p4est                 the forest
 * \param[in]     ctype                 the type of adjacency for the overlap
 * \param[in]     overlap               the number of layers of overlap (zero
 *                                      is acceptable)
 * \param[out]    first_local_quad      the local quadrants are assigned
 *                                      contiguous plex indices, starting with
 *                                      this index
 * \param[in,out] out_points_per_dim    filled with argument for
 *                                      DMPlexCreateFromDAG()
 * \param[in,out] out_cone_sizes        filled with argument for
 *                                      DMPlexCreateFromDAG()
 * \param[in,out] out_cones             filled with argument for
 *                                      DMPlexCreateFromDAG()
 * \param[in,out] out_cone_orientations filled with argument for
 *                                      DMPlexCreateFromDAG()
 * \param[in,out] out_vertex_coords     filled with argument for
 *                                      DMPlexCreateFromDAG()
 * \param[in,out] out_children          filled with argument for
 *                                      DMPlexSetTree()
 * \param[in,out] out_parents           filled with argument for
 *                                      DMPlexSetTree()
 * \param[in,out] out_childids          filled with argument for
 *                                      DMPlexSetTree()
 * \param[in,out] out_leaves            filled with argument for
 *                                      PetscSFSetGraph()
 * \param[in,out] out_remotes           filled with argument for
 *                                      PetscSFSetGraph()
 */
void                p4est_get_plex_data (p4est_t * p4est,
                                         p4est_connect_type_t ctype,
                                         int overlap,
                                         p4est_locidx_t * first_local_quad,
                                         sc_array_t * out_points_per_dim,
                                         sc_array_t * out_cone_sizes,
                                         sc_array_t * out_cones,
                                         sc_array_t * out_cone_orientations,
                                         sc_array_t * out_vertex_coords,
                                         sc_array_t * out_children,
                                         sc_array_t * out_parents,
                                         sc_array_t * out_childids,
                                         sc_array_t * out_leaves,
                                         sc_array_t * out_remotes);

SC_EXTERN_C_END;
#endif /* P4EST_PLEX_H */
