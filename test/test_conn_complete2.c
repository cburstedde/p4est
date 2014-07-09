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
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#else
#include <p8est_connectivity.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#endif

static void
test_the_p4est (p4est_connectivity_t * conn, int N)
{
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;
  p4est_lnodes_t     *lnodes;

  p4est = p4est_new (sc_MPI_COMM_WORLD, conn, 0, NULL, NULL);
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  lnodes = p4est_lnodes_new (p4est, ghost, N);
  p4est_lnodes_destroy (lnodes);
  p4est_ghost_destroy (ghost);
  p4est_destroy (p4est);
}

static void
test_complete (p4est_connectivity_t * conn, const char *which, int test_p4est)
{
  SC_GLOBAL_INFOF ("Testing standard connectivity %s\n", which);
  SC_CHECK_ABORTF (p4est_connectivity_is_valid (conn),
                   "Invalid connectivity %s before completion", which);
  if (0 && test_p4est) {
    test_the_p4est (conn, 3);
  }

  SC_GLOBAL_INFOF ("Testing completion for connectivity %s\n", which);
  p4est_connectivity_complete (conn);
  SC_CHECK_ABORTF (p4est_connectivity_is_valid (conn),
                   "Invalid connectivity %s after completion", which);
  if (test_p4est) {
    test_the_p4est (conn, 3);
  }

  p4est_connectivity_destroy (conn);
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

#ifndef P4_TO_P8
  test_complete (p4est_connectivity_new_unitsquare (), "unitsquare", 1);
  test_complete (p4est_connectivity_new_periodic (), "2D periodic", 0);
  test_complete (p4est_connectivity_new_rotwrap (), "rotwrap", 0);
  test_complete (p4est_connectivity_new_corner (), "corner", 1);
  test_complete (p4est_connectivity_new_moebius (), "moebius", 1);
  test_complete (p4est_connectivity_new_star (), "star", 1);
  test_complete (p4est_connectivity_new_brick (3, 18, 0, 1),
                 "2D periodic brick", 0);
  test_complete (p4est_connectivity_new_brick (3, 18, 0, 0), "2D brick", 1);
#else
  test_complete (p8est_connectivity_new_unitcube (), "unitcube", 1);
  test_complete (p8est_connectivity_new_periodic (), "3D periodic", 0);
  test_complete (p8est_connectivity_new_rotwrap (), "rotwrap", 0);
  test_complete (p8est_connectivity_new_twowrap (), "twowrap", 1);
  test_complete (p8est_connectivity_new_twocubes (), "twocubes", 1);
  test_complete (p8est_connectivity_new_rotcubes (), "rotcubes", 1);
  test_complete (p8est_connectivity_new_brick (3, 2, 8, 1, 0, 1),
                 "3D periodic brick", 0);
  test_complete (p8est_connectivity_new_brick (3, 2, 8, 0, 0, 0),
                 "3D brick", 1);
#endif

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
