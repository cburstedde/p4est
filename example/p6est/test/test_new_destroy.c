/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2014 The University of Texas System
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

#include <p6est.h>
#include <p6est_vtk.h>

char               *TEST_USER_POINTER;

void
init_fn (p6est_t * p6est, p4est_topidx_t which_tree,
         p4est_quadrant_t * col, p2est_quadrant_t * layer)
{
  SC_CHECK_ABORT (p6est->user_pointer == TEST_USER_POINTER,
                  "user_pointer corruption\n");
}

int
main (int argc, char **argv)
{
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
  p4est_connectivity_t *conn4;
  p6est_connectivity_t *conn;
  p6est_t            *p6est;
  double              height[3] = { 0., 0., 1. };
  int                 mpiret;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  SC_CHECK_ABORTF (argc == 2,
                   "Usage:\n%s NAME\n"
                   "  NAME=<corner|cubed|disk|periodic|rotwrap|star|unit>\n",
                   argv[0]);

  conn4 = p4est_connectivity_new_byname (argv[1]);
  conn = p6est_connectivity_new (conn4, NULL, height);

  p4est_connectivity_destroy (conn4);

  p6est = p6est_new (mpicomm, conn, 4, init_fn, TEST_USER_POINTER);
  p6est_destroy (p6est);

  p6est = p6est_new_ext (mpicomm, conn, 0, 1, 2, 1, 3, init_fn,
                         TEST_USER_POINTER);

  p6est_vtk_write_file (p6est, "p6est_test_new_destroy");

  p6est_destroy (p6est);

  p6est_connectivity_destroy (conn);

  /* exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
