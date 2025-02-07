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

/* This program demonstrates how to create a new p4est structure of a unit
square domain, refine it with initial level 3, and write the forest to a VTK
file for visualization purposes. */

#include <p4est_extended.h>
#include <p4est_vtk.h>

/* also please run p4estindent on all main programs */

int
main (int argc, char **argv)
{
  /* Define the minimum level of refinement, used for the initial mesh. */
  const int           minlevel = 3;

  /* Declare the MPI communicator and initialize the MPI environment */
  sc_MPI_Comm         mpicomm = sc_MPI_COMM_WORLD;
  int                 mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* Create a new connectivity for a unit square domain */
  p4est_connectivity_t *conn;
  conn = p4est_connectivity_new_unitsquare ();

  /* Create a new p4est structure (forest of quadtrees) with the created connectivity */
  p4est_t            *p4est;

  /* Create the structure by using p4est_new_ext, which is a more general form
     of p4est_new.
     Major parameters:
     minlevel: The forest is refined at most to this level. If minlevel is negative
     or 0, then it has no effect.
     fill_uniform: If true, the forest is filled with a uniform mesh instead of
     the coarsest possible one. For the sake of reproducibility, set this to 1.
     min_quadrants: Minimum initial quadrants per processor. For reproducibility,
     set this to 0 so that the mesh is independent of the number of processes. */
  p4est = p4est_new_ext (mpicomm, conn, 0, minlevel, 1, 0, NULL, NULL);

  /* Write the forest structure to a VTK file for visualization purposes */
  /* The filename will be prefixed with "p4est_a_unitsquare" */
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_a_unitsquare");

  /* Destroy the p4est structure to free memory */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  /* Finalize the MPI environment and check for errors */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;

}
