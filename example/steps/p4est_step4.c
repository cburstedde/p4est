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

/** \file p4est_step4.c
 *
 * This 2D example program solves the Poisson equation using finite elements.
 */

/* p4est has two separate interfaces for 2D and 3D, p4est*.h and p8est*.h.
 * Most API functions are available for both dimensions.  The header file
 * p4est_to_p8est.h #define's the 2D names to the 3D names such that most code
 * only needs to be written once.  In this example, we rely on this. */
#ifndef P4_TO_P8
#include <p4est_vtk.h>
#else
#include <p8est_vtk.h>
#endif

/** Callback function to decide on refinement.
 *
 * This function is called for every processor-local quadrant in order; its
 * return value is understood as a boolean refinement flag.  We refine around a
 * h = 1/8 block with left front lower corner (5/8, 2/8, 6/8).
 */
static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  /* Compute the integer coordinate extent of a quadrant of length 2^(-3). */
  const p4est_qcoord_t eighth = P4EST_QUADRANT_LEN (3);

  /* Compute the length of the current quadrant in integer coordinates. */
  const p4est_qcoord_t length = P4EST_QUADRANT_LEN (quadrant->level);

  /* Refine if the quadrant intersects the block in question. */
  return ((quadrant->x + length > 5 * eighth && quadrant->x < 6 * eighth) &&
          (quadrant->y + length > 2 * eighth && quadrant->y < 3 * eighth) &&
#ifdef P4_TO_P8
          (quadrant->z + length > 6 * eighth && quadrant->z < 7 * eighth) &&
#endif
          1);
}

/** The main function of the step4 example program.
 *
 * It creates a connectivity and forest, refines it, and solves the Poisson
 * equation with piecewise linear finite elements.
 */
int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 startlevel, endlevel, level;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;

  /* Initialize MPI; see sc_mpi.h.
   * If configure --enable-mpi is given these are true MPI calls.
   * Else these are dummy functions that simulate a single-processor run. */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  /* These functions are optional.  If called they store the MPI rank as a
   * static variable so subsequent global p4est log messages are only issued
   * from processor zero.  Here we turn off most of the logging; see sc.h. */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_PRODUCTION);
  P4EST_GLOBAL_PRODUCTIONF
    ("This is the p4est %dD demo example/steps/%s_step4\n",
     P4EST_DIM, P4EST_STRING);

  /* Create a forest that consists of just one quadtree/octree.
   * This file is compiled for both 2D and 3D: the macro P4_TO_P8 can be
   * checked to execute dimension-dependent code. */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_unitsquare ();
#else
  conn = p8est_connectivity_new_unitcube ();
#endif

  /* Create a forest that is not refined; it consists of the root octant.
   * The p4est_new_ext function can take a startlevel for a load-balanced
   * initial uniform refinement.  Here we refine adaptively instead. */
  startlevel = 0;
  p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

  /* Refine the forest iteratively, load balancing at each iteration.
   * This is important when starting with an unrefined forest */
  endlevel = 5;
  for (level = startlevel; level < endlevel; ++level) {
    p4est_refine (p4est, 0, refine_fn, NULL);
    /* Refinement has lead to up to 8x more elements; redistribute them. */
    p4est_partition (p4est, 0, NULL);
  }
  if (startlevel < endlevel) {
    /* For finite elements this corner balance is not strictly required.
     * We call balance only once since it's more expensive than coarsen/refine
     * and partition. */
    p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
    /* We repartition with coarsening in mind to allow for
     * partition-independent a-posteriori adaptation. */
    p4est_partition (p4est, 1, NULL);
  }

  /* Write the forest to disk for visualization, one file per processor. */
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_step4");

  /* Destroy the p4est and the connectivity structure. */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  /* Verify that allocations internal to p4est and sc do not leak memory.
   * This should be called if sc_init () has been called earlier. */
  sc_finalize ();

  /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
