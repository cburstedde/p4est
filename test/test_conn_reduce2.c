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
#include <p4est_connectivity.h>
#else
#include <p8est_connectivity.h>
#endif

static void
test_reduce (p4est_connectivity_t * conn, const char *which)
{
  SC_GLOBAL_INFOF ("Testing standard connectivity %s\n", which);
  SC_CHECK_ABORTF (p4est_connectivity_is_valid (conn),
                   "Invalid connectivity %s before reduce", which);
  p4est_connectivity_reduce (conn);
  SC_CHECK_ABORTF (p4est_connectivity_is_valid (conn),
                   "Invalid connectivity %s after reduce", which);
  SC_GLOBAL_INFOF ("Testing completion for reduced connectivity %s\n", which);
  p4est_connectivity_complete (conn);
  SC_CHECK_ABORTF (p4est_connectivity_is_valid (conn),
                   "Invalid connectivity %s after completion", which);
  p4est_connectivity_destroy (conn);
}

int
main (int argc, char *argv[])
{
  int                 mpiret;

  /* initialize MPI and p4est internals */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

#ifndef P4_TO_P8
  test_reduce (p4est_connectivity_new_unitsquare (), "unitsquare");
  test_reduce (p4est_connectivity_new_periodic (), "periodic");
  test_reduce (p4est_connectivity_new_rotwrap (), "rotwrap");
  test_reduce (p4est_connectivity_new_corner (), "corner");
  test_reduce (p4est_connectivity_new_pillow (), "pillow");
  test_reduce (p4est_connectivity_new_moebius (), "moebius");
  test_reduce (p4est_connectivity_new_star (), "star");
  test_reduce (p4est_connectivity_new_cubed (), "cubed");
  test_reduce (p4est_connectivity_new_disk (), "disk");
  test_reduce (p4est_connectivity_new_brick (3, 2, 0, 0), "brick00");
  test_reduce (p4est_connectivity_new_brick (3, 2, 0, 1), "brick01");
  test_reduce (p4est_connectivity_new_brick (3, 2, 1, 0), "brick10");
  test_reduce (p4est_connectivity_new_brick (3, 2, 1, 1), "brick11");
#else
  test_reduce (p8est_connectivity_new_unitcube (), "unitcube");
  test_reduce (p8est_connectivity_new_periodic (), "periodic");
  test_reduce (p8est_connectivity_new_rotwrap (), "rotwrap");
  test_reduce (p8est_connectivity_new_twocubes (), "twocubes");
  test_reduce (p8est_connectivity_new_twowrap (), "twowrap");
  test_reduce (p8est_connectivity_new_rotcubes (), "rotcubes");
  test_reduce (p8est_connectivity_new_brick (4, 3, 2, 0, 0, 0), "brick000");
  test_reduce (p8est_connectivity_new_brick (4, 3, 2, 0, 0, 1), "brick001");
  test_reduce (p8est_connectivity_new_brick (4, 3, 2, 0, 1, 0), "brick010");
  test_reduce (p8est_connectivity_new_brick (4, 3, 2, 0, 1, 1), "brick011");
  test_reduce (p8est_connectivity_new_brick (4, 3, 2, 1, 0, 0), "brick100");
  test_reduce (p8est_connectivity_new_brick (4, 3, 2, 1, 0, 1), "brick101");
  test_reduce (p8est_connectivity_new_brick (4, 3, 2, 1, 1, 1), "brick111");
#endif
  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
