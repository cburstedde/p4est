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

#ifndef P8EST_IO_H
#define P8EST_IO_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** Extract processor local quadrants' x y z level data.
 * Optionally extracts the quadrant data as well into a separate array.
 * \param [in] p8est    The forest is not modified.
 * \param [in,out] data If not NULL, pointer to a pointer that will be set
 *                      to a newly allocated array with per-quadrant data.
 *                      Must be NULL if p4est->data_size == 0.
 * \return              An array of type p8est_qcoord_t that contains
 *                      x y z level for each quadrant on this processor.
 *                      The tree information is not extracted.
 */
sc_array_t         *p8est_deflate_quadrants (p8est_t * p8est,
                                             sc_array_t ** data);

/** Create a new p4est based on serialized data.
 * Its revision counter is set to zero.
 * See p8est.h and p8est_communication.h for more information on parameters.
 * \param [in] mpicomm       A valid MPI communicator.
 * \param [in] connectivity  This is the connectivity information that
 *                           the forest is built with.  Note that p4est
 *                           does not take ownership of the memory.
 * \param [in] global_first_quadrant First global quadrant on each proc and
 *                           one beyond.  Copied into global_first_quadrant.
 *                           Local count on rank is gfq[rank + 1] - gfq[rank].
 * \param [in] pertree       The cumulative quadrant counts per tree.
 * \param [in] quadrants     Array as returned by p8est_deflate_quadrants.
 * \param [in] data          Array as from p8est_deflate_quadrants or NULL.
 *                           The elem_size of this array informs data_size.
 *                           Its elem_count equals the number of local quads.
 * \param [in] user_pointer  Assign to the user_pointer member of the p4est.
 * \return              The newly created p4est with a zero revision counter.
 */
p8est_t            *p8est_inflate (sc_MPI_Comm mpicomm,
                                   p8est_connectivity_t * connectivity,
                                   const p4est_gloidx_t *
                                   global_first_quadrant,
                                   const p4est_gloidx_t * pertree,
                                   sc_array_t * quadrants, sc_array_t * data,
                                   void *user_pointer);

SC_EXTERN_C_END;

#endif /* !P8EST_IO_H */
