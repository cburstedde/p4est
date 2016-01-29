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

#ifndef P8EST_SEARCH_BUILD_H
#define P8EST_SEARCH_BUILD_H

#include <p8est_search.h>

/** \file p8est_search_build.h
 * Create a new p8est object by running p8est_search_local.
 * This allows to create a heavily coarsened forest in one pass.
 */

typedef struct
{
  p8est_t            *from;      /**< Existing forest used as template. */
  p8est_t            *p4est;     /**< New forest being built. */
}
p8est_search_build_t;

p8est_search_build_t *p8est_search_build_new (p8est_t * p4est);

p8est_t            *p8est_search_build_complete (p8est_search_build_t *
                                                 build);

#endif /* ! P8EST_SEARCH_BUILD_H */
