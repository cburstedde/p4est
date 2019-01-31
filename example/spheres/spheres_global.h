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

/*
 * This is an internal header file.
 * It is included from a c file that may be preceded by p4est_to_p8est.h.
 */

#ifndef P4EST_SPHERES_GLOBAL_H
#define P4EST_SPHERES_GLOBAL_H

#include <sc_statistics.h>
#ifndef P4_TO_P8
#include <p4est.h>
#else
#include <p8est.h>
#endif /* P4_TO_P8 */

SC_EXTERN_C_BEGIN;

typedef struct qu_data
{
}
qu_data_t;

typedef struct spheres_global
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize, mpirank;

  int                 minlevel;
  int                 maxlevel;
  int                 scaling;

  double              vdensity;

  const char         *prefix;

  p4est_connectivity_t *conn;
  p4est_t            *p4est;
}
spheres_global_t;

SC_EXTERN_C_END;

#endif /* !P4EST_SPHERES_GLOBAL_H */
