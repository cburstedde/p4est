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

#ifndef P4_TO_P8
#include <p4est.h>
#include <p4est_extended.h>
#include <p4est_nodes.h>
#else
#include <p8est.h>
#include <p8est_extended.h>
#include <p8est_nodes.h>
#endif

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  if (quadrant->level >= P4EST_QMAXLEVEL) {
    return 0;
  }
  if (quadrant->x == 0 && quadrant->y == 0 &&
#ifdef P4_TO_P8
      quadrant->z == 0 &&
#endif
      1) {
    return 1;
  }

  return 0;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  p4est_ghost_t      *ghost;
  p4est_nodes_t      *nodes1, *nodes2;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  /* create connectivity and forest structures */
#ifndef P4_TO_P8
  connectivity = p4est_connectivity_new_unitsquare ();
#else
  connectivity = p8est_connectivity_new_unitcube ();
#endif
  /* specifying min_quadrants is not partition independent but fun */
  p4est = p4est_new_ext (mpicomm, connectivity, 15, 0, 0, 0, NULL, NULL);

  /* refine and balance to make the number of elements interesting */
  p4est_refine (p4est, 1, refine_fn, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  /* create the ghost structure */
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

  /* create the nodes structure */
  P4EST_GLOBAL_INFO ("Making nodes with ghosts\n");
  nodes1 = p4est_nodes_new (p4est, ghost);
  P4EST_GLOBAL_INFO ("Making nodes without ghosts\n");
  nodes2 = p4est_nodes_new (p4est, NULL);

  /* clean up and exit */
  p4est_nodes_destroy (nodes1);
  p4est_nodes_destroy (nodes2);
  p4est_ghost_destroy (ghost);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
