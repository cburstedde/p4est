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

/** \file p4est_step2.c
 *
 * This 2D example program refines a domain given by an ABAQUS .inp file.
 */

/* p4est has two separate interfaces for 2D and 3D, p4est*.h and p8est*.h.
 * Most API functions are available for both dimensions.  The header file
 * p4est_to_p8est.h #define's the 2D names to the 3D names such that most code
 * only needs to be written once.  In this example, we rely on this
 * to compile both 2D and 3D from the same source file. */
#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_vtk.h>
#endif

/** We're not using p4est->user_pointer here but take a shortcut.
 */
static int          refine_level = 0;

/** Callback function to decide on refinement.
 *
 * Refinement and coarsening is controlled by callback functions.
 * This function is called for every processor-local quadrant in order; its
 * return value is understood as a boolean refinement flag.
 *
 * Here we use uniform refinement.  Note that this function is not suitable for
 * recursive refinement and must be used in an iterative fashion.
 */
static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  return 1;
}

/** The main function of the step2 example program.
 *
 * It creates a connectivity from an ABAQUS .inp file and forest, refines it,
 * and writes a VTK file.
 */
int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 balance;
  int                 level;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  const char         *filename;

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
    ("This is the p4est %dD demo example/steps/%s_step2\n",
     P4EST_DIM, P4EST_STRING);

  /* Get the .inp file name from the list of arguments along with an optional
   * level of refinement.  */
  if (argc != 2 && argc != 3) {
    SC_GLOBAL_LERRORF ("Usage: %s <inp file name> [level of refinement]\n",
                       argv[0]);
    sc_abort ();
  }
  filename = argv[1];
  if (argc == 3)
    refine_level = atoi (argv[2]);

  /* Create a forest from the inp file with name filename  */
  conn = p4est_connectivity_read_inp (filename);
  if (conn == NULL) {
    P4EST_LERRORF ("Failed to read a valid connectivity from %s\n", filename);
    sc_abort ();
  }

#ifdef P4EST_WITH_METIS
  /* Use metis (if p4est is compiled with the flag '--with-metis') to
   * reorder the connectivity for better partitioning of the forest
   * across processors.
   */
  p4est_connectivity_reorder (mpicomm, 0, conn, P4EST_CONNECT_FACE);
#endif /* P4EST_WITH_METIS */

  /* Create a forest that is not refined; it consists of the root octant. */
  p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

  /* Refine the forest iteratively, load balancing at each iteration.
   * This is important when starting with an unrefined forest */
  for (level = 0; level < refine_level; ++level) {
    p4est_refine (p4est, 0, refine_fn, NULL);
    /* Refinement has lead to up to 8x more elements; redistribute them. */
    p4est_partition (p4est, 0, NULL);
  }

  /* If we call the 2:1 balance we ensure that neighbors do not differ in size
   * by more than a factor of 2.  This can optionally include diagonal
   * neighbors across edges or corners as well; see p4est.h.
   *
   * Note that this balance step is not strictly necessary since we are using
   * uniform refinement but may be required for other types of refinement.
   */
  balance = 1;
  if (balance) {
    p4est_balance (p4est, P4EST_CONNECT_FACE, NULL);
    p4est_partition (p4est, 0, NULL);
  }

  /* Write the forest to disk for visualization, one file per processor. */
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_step2");

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
