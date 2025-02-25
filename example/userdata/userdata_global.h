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
 * This is an internal header that must not be installed.
 * It may be included from 2D code or alternatively
 * from 3D code after including p4est_to_p8est.h.
 */

#ifndef P4EST_USERDATA_GLOBAL_H
#define P4EST_USERDATA_GLOBAL_H

#include <sc_options.h>
#ifndef P4_TO_P8
#include <p4est_geometry.h>
#else
#include <p8est_geometry.h>
#define p4est_userdata_run              p8est_userdata_run
#define p4est_userdata_global           p8est_userdata_global
#define p4est_userdata_global_t         p8est_userdata_global_t
#endif

/* we maintain application data in a global structure and pass it around */
typedef struct p4est_userdata_global
{
  /* global objects */
  sc_MPI_Comm         mpicomm;
  sc_options_t       *options;
  int                 help;
  int                 maxlevel;
  int                 novtk;
  int                 noint;
  int                 noext;
  const char         *configuration;
  p4est_geometry_t   *geom;
  p4est_connectivity_t *conn;
  p4est_t            *p4est;

  /* temporary work storage */
  int                 in_balance;
  p4est_locidx_t      qcount;
  p4est_locidx_t      bcount;
  sc_array_t         *qarray;
}
p4est_userdata_global_t;

/* this function contains the core of the demonstration program */
int                 p4est_userdata_run (p4est_userdata_global_t *g);

#endif /* !P4EST_USERDATA_GLOBAL_H */
