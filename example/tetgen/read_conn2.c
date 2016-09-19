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
#include <p4est_connectivity.h>
#else
#include <p8est_connectivity.h>
#endif

int
main (int argc, char **argv)
{
  p4est_connectivity_t *conn;
  size_t              bytes;

  if (argc != 2) {
    char               *cp, *bn;

    cp = strdup (argv[0]);
    bn = basename (cp);
    fprintf (stderr, "Usage: %s <connectivity file>\n", bn);
    free (cp);
    exit (1);
  }

  conn = p4est_connectivity_load (argv[1], &bytes);
  if (conn == NULL) {
    P4EST_PRODUCTIONF ("Failed to load connectivity file \"%s\"\n", argv[1]);
  }
  else {
    P4EST_STATISTICSF ("Loaded \"%s\" bytes %lld\n",
                       argv[1], (long long) bytes);
    P4EST_STATISTICSF ("Loaded %lld trees, %lld vertices, %lld corners\n",
                       (long long) conn->num_trees,
                       (long long) conn->num_vertices,
                       (long long) conn->num_corners);
    p4est_connectivity_destroy (conn);
  }

  return 0;
}
