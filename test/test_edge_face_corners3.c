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

/*
 * Purpose of this program is to verify the edge face corner connections
 * that are computed in the static function p8est_compute_edge_face_corners.
 */

#include <p8est.h>
#include <p4est_to_p8est.h>

static int
p8est_compute_edge_face_corners (int edge, int face, int corners[])
{
  int                 nfound = 0;
  int                 fc, c0, c1, cx;

  P4EST_ASSERT (0 <= edge && edge < 12);
  P4EST_ASSERT (0 <= face && face < 6);

  c0 = p8est_edge_corners[edge][0];
  c1 = p8est_edge_corners[edge][1];

  for (fc = 0; fc < 4; ++fc) {
    cx = p8est_face_corners[face][fc];
    if (c0 == cx) {
      corners[0] = fc;
      ++nfound;
      continue;
    }
    if (c1 == cx) {
      corners[1] = fc;
      ++nfound;
      continue;
    }
  }
  P4EST_ASSERT (nfound <= 2);

  return nfound == 2;
}

int
main (int argc, char **argv)
{
  int                 edge, face, cs[2];
  int                 success;

  sc_init (sc_MPI_COMM_NULL, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  for (edge = 0; edge < 12; ++edge) {
    for (face = 0; face < 6; ++face) {
      cs[0] = cs[1] = 0;        /* to avoid compiler warning */
      success = p8est_compute_edge_face_corners (edge, face, cs);

      if (!success) {
        P4EST_LDEBUGF ("Nothing for %d %d\n", edge, face);
        SC_CHECK_ABORT (p8est_edge_face_corners[edge][face][0] == -1 &&
                        p8est_edge_face_corners[edge][face][1] == -1,
                        "Invalid nonexisting connection");
      }
      else {
        P4EST_LDEBUGF ("Results for %d %d are %d %d\n",
                       edge, face, cs[0], cs[1]);
        SC_CHECK_ABORT (p8est_edge_face_corners[edge][face][0] == cs[0] &&
                        p8est_edge_face_corners[edge][face][1] == cs[1],
                        "Invalid existing connection");
      }
    }
  }

  sc_finalize ();

  return 0;
}
