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
 * This example program demonstrates how to manage application data.
 *
 *   p4est_userdata <OPTIONS> [<configuration> [<level>]]
 *
 * with possible configurations (default is "unit"):
 *   o unit          Refinement on the unit square.
 *   o periodic      Unit square with all-periodic boundary conditions.
 *   o brick         Refinement on a 2x3 rectangle of quadtrees.
 *   o disk2d        Refinement on a spherical disk of five trees.
 *   o corner        Refinement on a non-planar hexagon made of three trees.
 *   o moebius       Refinement on a 5-tree Moebius band embedded in 3D.
 *   o icosahedron   Refinement on the icosahedron sphere with geometry.
 *
 * The maximum refinement level may be appended (default is 4).
 *
 * The following options are recognized:
 *   --help          Display a usage and help message and exit successfully.
 *   --level         The level may alternatively be specified as an option.
 *
 * Invalid options or arguments result in an error message and exit status.
 */

/* This file is used to compile both the 2D and the 3D code, separately. */
#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#endif

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpicomm;

  /* initialize MPI subsystem */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* initialize p4est parallel logging by querying for the rank below */
  mpicomm = sc_MPI_COMM_WORLD;

  /*
   * The options 1, 1 catch signals and set an abort handler.
   * This is NOT what you want if (a) you are using p4est purely as a
   * library or (b) your main programs refer to a different framework for
   * such lowlevel functionality.  In these cases please use 0, 0.
   *
   * The log level SC_LP_APPLICATION is rather conservative.  If you'd
   * prefer total silence under normal operating conditions you may use
   * SC_LP_ERROR.  p4est does not trigger the error priority by itself, but
   * it can be used by the application developer, for example to issue log
   * messages on any usage or file I/O errors detected.
   */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_APPLICATION);

  /*
   * The setting SC_LP_APPLICATION will log levels from SC_LP_PRODUCTION
   * upwards.  Thus, if your program should print performance metrics, for
   * example, that can be accomplished (for information available on rank 0)
   * using the logging macros P4EST_PRODUCTION and P4EST_PRODUCTIONF.
   */
  p4est_init (NULL, SC_LP_APPLICATION);

  /* TO DO: write demo code */

  /* check memory balance and clean up internal registrations */
  sc_finalize ();

  /* release the MPI subsytem */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  /* this standard value returns success to the calling shell */
  return EXIT_SUCCESS;
}
