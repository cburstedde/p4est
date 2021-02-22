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

#include <p4est_base.h>

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  int                 num_failed_tests;
  int                 version_major, version_minor;
  const char         *version;
  char                version_tmp[32];
  const char         *unknown = "UNKNOWN";

  /* standard initialization */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* check all functions related to version numbers of p4est */
  num_failed_tests = 0;
  version = p4est_version ();
  P4EST_GLOBAL_LDEBUGF ("Full p4est version: %s\n", version);

  if (!strncmp (version, unknown, strlen (unknown))) {
    P4EST_GLOBAL_VERBOSE ("Version is unknown\n");
    version_major = p4est_version_major ();
    version_minor = p4est_version_minor ();
    if (version_major != 0 || version_minor != 0) {
      P4EST_GLOBAL_VERBOSE ("Unknown version of p4est not zero\n");
      num_failed_tests++;
    }
  }
  else {
    version_major = p4est_version_major ();
    P4EST_GLOBAL_LDEBUGF ("Major p4est version: %d\n", version_major);
    snprintf (version_tmp, 32, "%d", version_major);
    if (strncmp (version, version_tmp, strlen (version_tmp))) {
      P4EST_GLOBAL_VERBOSE ("Test failure for major version of p4est\n");
      num_failed_tests++;
    }

    version_minor = p4est_version_minor ();
    P4EST_GLOBAL_LDEBUGF ("Minor p4est version: %d\n", version_minor);
    snprintf (version_tmp, 32, "%d.%d", version_major, version_minor);
    if (strncmp (version, version_tmp, strlen (version_tmp))) {
      P4EST_GLOBAL_VERBOSE ("Test failure for minor version of p4est\n");
      num_failed_tests++;
    }
  }

  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return num_failed_tests ? EXIT_FAILURE : EXIT_SUCCESS;
}
