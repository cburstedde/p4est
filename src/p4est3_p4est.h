/*
  This file is part of p4est, version 3
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2019 individual authors
  Originally written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

/** \file p4est3_p4est.h
 *
 * Wrap standard p4est objects for use in the version 3 implementation.
 *
 * We may use the standard quadrants to construct a version 3 forest.
 * To this end, call \ref p4est3_quadrant_vtable_p4est to populate a virtual
 * table suitable for passing it to \ref p4est3_set_quadrant_vtable.
 *
 * We may use \ref p4est_connectivity_t and \ref p4est_t objects in version 3
 * by constructing the respective connectivity and forest objects from
 * a virtual table.
 * This file provides the associated convenience constructors.
 *
 * \ingroup p4est3
 */

#ifndef P4EST3_P4EST_H
#define P4EST3_P4EST_H

#include <p4est.h>
#include <p4est3.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/** Return pointer to quadrant virtual table for the standard version */
sc3_error_t        *p4est3_quadrant_vtable_p4est
  (const p4est3_quadrant_vtable_t ** qvt);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_P4EST_H */
