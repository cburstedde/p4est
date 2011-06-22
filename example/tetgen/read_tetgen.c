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

#include <p8est_tets_hexes.h>

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 num_flips;
  const char         *argbasename;
  p8est_tetgen_t     *ptg;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  if (argc != 2) {
    SC_GLOBAL_LERRORF ("Usage: %s <tetgen file base name>\n", argv[0]);
    sc_abort ();
  }
  argbasename = argv[1];

  /* read tetgen nodes and tetrahedra from files */
  ptg = p8est_tetgen_read (argbasename);
  SC_CHECK_ABORTF (ptg != NULL, "Failed to read tetgen %s", argbasename);
  P4EST_GLOBAL_STATISTICSF ("Read %d nodes and %d tets %s attributes\n",
                            (int) ptg->nodes->elem_count / 3,
                            (int) ptg->tets->elem_count / 4,
                            ptg->tet_attributes != NULL ? "with" : "without");

  /* flip orientation to right-handed */
  num_flips = p8est_tetgen_make_righthanded (ptg);
  P4EST_GLOBAL_STATISTICSF ("Performed %d orientation flip(s)\n", num_flips);

  /* clean up */
  p8est_tetgen_destroy (ptg);

  sc_finalize ();
  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
