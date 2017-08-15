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
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_io.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_io.h>
#endif
#include <sc_options.h>
#include <sc_statistics.h>

#ifndef P4_TO_P8
#define P4EST_CONN_SUFFIX "p4c"
#define P4EST_FOREST_SUFFIX "p4p"
#else
#define P4EST_CONN_SUFFIX "p8c"
#define P4EST_FOREST_SUFFIX "p8p"
#endif

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpirank;
  int                 first_arg;
  int                 autopartition;
  int                 broadcasthead;
  const char         *filename;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  sc_options_t       *opt;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* initialize libsc and p4est */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* handle command line options */
  opt = sc_options_new (argv[0]);
  sc_options_add_string (opt, 'f', "file", &filename,
                         NULL, "p4est data file to load");
  sc_options_add_bool (opt, 'a', "autopartition", &autopartition,
                       0, "Create a uniform partition when reading");
  sc_options_add_bool (opt, 'b', "broadcasthead", &broadcasthead,
                       0, "Broadcast header information when reading");
  first_arg = sc_options_parse (p4est_package_id, SC_LP_INFO,
                                opt, argc, argv);
  SC_CHECK_ABORT (first_arg == argc, "Option error: no arguments required");
  if (sc_is_root ()) {
    sc_options_print_summary (p4est_package_id, SC_LP_PRODUCTION, opt);
  }

  /* only proceed if the filename is set */
  if (filename != NULL) {
    /* load the file */
    p4est = p4est_load_ext (filename, mpicomm, 0, 0,
                            autopartition, broadcasthead,
                            NULL, &connectivity);

    /* and destroy the forest read */
    p4est_destroy (p4est);
    p4est_connectivity_destroy (connectivity);
  }

  /* clean up and exit */
  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
