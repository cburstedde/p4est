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
 * Usage: p4est_overlap
 */

#ifndef P4_TO_P8
#include <p4est.h>
#else
#include <p8est.h>
#endif

typedef struct overlap_producer
{
  p4est_connectivity_t *proconn;
  p4est_t            *pro4est;
}
overlap_producer_t;

typedef struct overlap_consumer
{
  p4est_connectivity_t *conconn;
  p4est_t            *con4est;
}
overlap_consumer_t;

typedef struct overlap_global
{
  overlap_producer_t  pro;
  overlap_consumer_t  con;
}
overlap_global_t;

static void
overlap_apps_init (overlap_global_t * g)
{
}

static void
overlap_apps_reset (overlap_global_t * g)
{
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  overlap_global_t global, *g = &global;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpicomm = sc_MPI_COMM_WORLD;
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  overlap_apps_init (g);

  overlap_apps_reset (g);

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
