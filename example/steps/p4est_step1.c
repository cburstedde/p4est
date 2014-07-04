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

/** \file p4est_step1.c
 *
 * This 2D example program refines a domain based on given image data.
 * The image file hw32.h has been created with the GIMP and is compiled in.
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
#include "hw32.h"

/** The resolution of the image data in powers of two. */
#define P4EST_STEP1_PATTERN_LEVEL 5
/** The dimension of the image data. */
#define P4EST_STEP1_PATTERN_LENGTH (1 << P4EST_STEP1_PATTERN_LEVEL)
static const int    plv = P4EST_STEP1_PATTERN_LEVEL;    /**< Shortcut */
static const int    ple = P4EST_STEP1_PATTERN_LENGTH;   /**< Shortcut */
#ifdef P4_TO_P8
static const p4est_qcoord_t eighth = P4EST_QUADRANT_LEN (3);
#endif

/** Callback function to decide on refinement.
 *
 * Refinement and coarsening is controlled by callback functions.
 * This function is called for every processor-local quadrant in order; its
 * return value is understood as a boolean refinement flag.
 * In this example we use the image file hw32.h to determine the refinement.
 */
static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 tilelen;
  int                 offsi, offsj;
  int                 i, j;
  const char         *d;
  unsigned char       p[3];

  /* The connectivity chosen in main () only consists of one tree. */
  P4EST_ASSERT (which_tree == 0);

  /* We do not want to refine deeper than a given maximum level. */
  if (quadrant->level > plv) {
    return 0;
  }
#ifdef P4_TO_P8
  /* In 3D we extrude the 2D image in the z direction between [3/8, 5/8]. */
  if (quadrant->level >= 3 &&
      (quadrant->z < 3 * eighth || quadrant->z >= 5 * eighth)) {
    return 0;
  }
#endif

  /* We read the image data and refine wherever the color value is dark.
   * We can then visualize the output and highlight level > PATTERN_LEVEL. */
  tilelen = 1 << (plv - quadrant->level);       /* Pixel size of quadrant */
  offsi = quadrant->x / P4EST_QUADRANT_LEN (plv);       /* Pixel x offset */
  offsj = quadrant->y / P4EST_QUADRANT_LEN (plv);       /* Pixel y offset */
  P4EST_ASSERT (offsi >= 0 && offsj >= 0);
  for (j = 0; j < tilelen; ++j) {
    P4EST_ASSERT (offsj + j < ple);
    for (i = 0; i < tilelen; ++i) {
      P4EST_ASSERT (offsi + i < ple);
      d =
        hw32_header_data + 4 * (ple * (ple - 1 - (offsj + j)) + (offsi + i));
      HW32_HEADER_PIXEL (d, p);
      P4EST_ASSERT (p[0] == p[1] && p[1] == p[2]);      /* Grayscale image */
      if (p[0] < 128) {
        return 1;
      }
    }
  }
  return 0;
}

/** The main function of the step1 example program.
 *
 * It creates a connectivity and forest, refines it, and writes a VTK file.
 */
int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 recursive, partforcoarsen, balance;
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
    ("This is the p4est %dD demo example/steps/%s_step1\n",
     P4EST_DIM, P4EST_STRING);

  /* Create a forest that consists of just one quadtree/octree.
   * This file is compiled for both 2D and 3D: the macro P4_TO_P8 can be
   * checked to execute dimension-dependent code. */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_unitsquare ();
#else
  conn = p8est_connectivity_new_unitcube ();
#endif

  /* Create a forest that is not refined; it consists of the root octant. */
  p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

  /* Refine the forest recursively in parallel.
   * Since refinement does not change the partition boundary, this call
   * must not create an overly large number of quadrants.  A numerical
   * application would call p4est_refine non-recursively in a loop,
   * repartitioning in each iteration.
   * The P4EST_ASSERT macro only activates with --enable-debug.
   * We check against the data dimensions in example/steps/hw32.h. */
  P4EST_ASSERT (P4EST_STEP1_PATTERN_LENGTH == width);
  P4EST_ASSERT (P4EST_STEP1_PATTERN_LENGTH == height);
  recursive = 1;
  p4est_refine (p4est, recursive, refine_fn, NULL);

  /* Partition: The quadrants are redistributed for equal element count.  The
   * partition can optionally be modified such that a family of octants, which
   * are possibly ready for coarsening, are never split between processors. */
  partforcoarsen = 0;
  p4est_partition (p4est, partforcoarsen, NULL);

  /* If we call the 2:1 balance we ensure that neighbors do not differ in size
   * by more than a factor of 2.  This can optionally include diagonal
   * neighbors across edges or corners as well; see p4est.h. */
  balance = 1;
  if (balance) {
    p4est_balance (p4est, P4EST_CONNECT_FACE, NULL);
    p4est_partition (p4est, partforcoarsen, NULL);
  }

  /* Write the forest to disk for visualization, one file per processor. */
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_step1");

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
