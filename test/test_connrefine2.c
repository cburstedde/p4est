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

#ifndef P4_TO_P8
#include <p4est_connrefine.h>
#include <p4est_vtk.h>
#else
#include <p8est_connrefine.h>
#include <p8est_vtk.h>
#endif

int
main (int argc, char **argv)
{
  int                 mpiret;
  p4est_connectivity_t *conn_in, *conn_out;
  p4est_t            *p4est;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

#ifndef P4_TO_P8
  conn_in = p4est_connectivity_new_cubed ();
#else
  conn_in = p8est_connectivity_new_rotcubes ();
#endif
  conn_out = p4est_connectivity_refine (conn_in, 5);
  p4est_connectivity_destroy (conn_in);

  p4est = p4est_new (sc_MPI_COMM_WORLD, conn_out, 0, NULL, NULL);
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_test_connrefine");

  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn_out);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
