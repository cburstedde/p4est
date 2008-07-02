/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef P8EST_TRILINEAR_H
#define P8EST_TRILINEAR_H

#include <p8est.h>

/**
 * p8est_trilinear_neighbor_t
 *
 * The face neighbors are ordered in -x, +x, -y, +y, -z, and +z directions.
 *
 * For each direction, the neighbor(s) may be:
 *
 * 1. Out of the domain.
 *
 *      face_neighbor_eid[direction][0] = -1.
 *      face_neighbor_eid[direction][1--3] are undefined.
 *
 * 2. As large or twice as large as the current element:
 *
 *      face_neighbor_eid[direction][0] = index
 *
 *    where ((index >= 0) && (index < local_elem_num)) if the neighbor is
 *    LOCAL, or
 *    ((index >= local_elem_num) &&
 *    (index < (local_elem_num + phantom_elem_num))) if the neighor is
 *    REMOTE.
 *
 *      face_neighbor_eid[direction][1--3] are undefined.
 *
 * 3. Half as large as the current element:
 *
 *      face_neighbor_eid[direction][i] = index
 *
 *    where index is defined as above. Note that in this case all four
 *    neighbors (half as large) must exist.
 */

typedef struct p8est_trilinear_neighbor
{
  p4est_locidx_t      face_neighbor_eid[6][4];
}
p8est_trilinear_neighbor_t;

typedef struct p8est_trilinear_phantom
{
  p4est_qcoord_t      lx, ly, lz;
  p4est_qcoord_t      size;
  int                 owner_procid;     /* remote processor id */
  p4est_locidx_t      reid;     /* remote processor element index */
}
p8est_trilinear_phantom_t;

typedef struct p8est_trilinear_neighborhood
{
  p4est_locidx_t      phantom_elem_num;
  p8est_trilinear_phantom_t *phantom_elem_table;
  p8est_trilinear_neighbor_t *local_elem_neighbor_table;
}
p8est_trilinear_neighborhood_t;

p8est_trilinear_neighborhood_t *p8est_trilinear_neighborhood_new (p8est_t *
                                                                  p8est);

/* *INDENT-OFF* */
void
p8est_trilinear_neighborhood_destroy (p8est_trilinear_neighborhood_t *
                                      neighborhood);
/* *INDENT-ON* */

#endif /* !P8EST_TRILINEAR_H */
