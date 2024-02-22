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

/* This program demonstrates how to initialize the p4est library and print
messages to the console. */

#include <p4est.h>              /* Include the p4est library header for parallel adaptive
                                   mesh refinement. */

int
main (int argc, char **argv)
{
  /* Initialize the MPI communicator to the default world communicator, which
     includes all MPI processes. */
  sc_MPI_Comm         mpicomm = sc_MPI_COMM_WORLD;

  /* Initialize the MPI environment with arguments from the command line. */
  int                 mpiret = sc_MPI_Init (&argc, &argv);

  /* Check the return status of MPI initialization and abort if it failed. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the SC library with 1-byte alignment and default log priority */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  /* Initialize the p4est library with default logging priority. */
  p4est_init (NULL, SC_LP_DEFAULT);

  /* Print a global production-level message saying "Hello World!". */
  P4EST_GLOBAL_PRODUCTION ("Hello World!\n");

  /* Print a production-level message from the current MPI process. */
  P4EST_PRODUCTION ("Hello World from the parallel process!\n");

  /* Finalize the MPI environment and clean up all MPI resources. */
  mpiret = sc_MPI_Finalize ();

  /* Check the return status of MPI finalization and abort if an error occurred. */
  SC_CHECK_MPI (mpiret);

  /* Return 0 from main to indicate that the program has finished successfully. */
  return 0;
}
