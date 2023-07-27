/*
  This file is part of p4est, version 3.
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

/** \file p4est3.h
 *
 * Main interface file to construct and interact with a version 3 forest.
 * The forest may be constructed piece by piece using the setter functions.
 * It may also be defined by populating and passing a forest virtual table.
 * The latter approach allows third-party objects to pass as a legal forest.
 *
 * \ingroup p4est3
 */

#ifndef P4EST3_H
#define P4EST3_H

#include <sc3_mpi.h>
#include <p4est3_quadrant_vtable.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_H */
