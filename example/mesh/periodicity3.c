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

#include <p4est_to_p8est.h>
#ifdef P4_TO_P8
#include <p8est_connectivity.h>

int
main (int argc, char **argv)
{
  int                 i;
  size_t              bytes;
  const char         *filename = "conndebug.p8c";
  p4est_connectivity_t *conn;

  conn = p4est_connectivity_load (filename, &bytes);
  if (conn == NULL) {
    P4EST_LERRORF ("Could not read file %s\n", filename);
    return 1;
  }
  P4EST_INFOF ("Read %d bytes\n", (int) bytes);

  p4est_connectivity_complete (conn);

  for (i = 0; i < 4; ++i) {
    P4EST_VERBOSEF ("Z face %d\n", i);
    p4est_connectivity_join_faces (conn, 2 * i + 0, 2 * i + 1, 4, 5, 0);
  }
  /* p4est_connectivity_join_faces (conn, 0, 2, 3, 2, 0); */

  p4est_connectivity_destroy (conn);

  return 0;
}

#endif /* P4_TO_P8 */
