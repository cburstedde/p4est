/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2011 The University of Texas System
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

#ifndef P4_TO_P8
#include <p4est.h>
#else
#include <p8est.h>
#endif

int
main (int argc, char **argv)
{
  MPI_Comm            mpicomm;
  int                 mpiret;
  int                 size, rank;
  p4est_connectivity_t *conn;
  int                 k = 5;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpicomm, &size);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

#ifndef P4_TO_P8
  conn = p4est_connectivity_new_brick (16, 17, 0, 1);
#else
  conn = p8est_connectivity_new_brick (7, 8, 9, 0, 1, 1);
#endif

  p4est_connectivity_reorder (mpicomm, k, conn, P4EST_CONNECT_FULL);

  SC_CHECK_ABORT (p4est_connectivity_is_valid (conn), "Bad reordering");

  p4est_connectivity_destroy (conn);

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
