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
#include <p4est_bits.h>
#include <p4est_connectivity.h>
#else
#include <p8est_bits.h>
#include <p8est_connectivity.h>
#endif

static int
test_node_coordinates (p4est_quadrant_t *r, p4est_qcoord_t coords[])
{
  P4EST_ASSERT (r != NULL);
  P4EST_ASSERT (p4est_quadrant_is_node (r, 0));
  P4EST_ASSERT (coords != NULL);
  return (r->x == coords[0] &&
          r->y == coords[1] &&
#ifdef P4_TO_P8
          r->z == coords[2] &&
#endif
          1);
}

static void
test_connectivity (sc_MPI_Comm mpicomm, p4est_connectivity_t *conn)
{
  int                 mpiret;
  int                 size, rank;
  int                 face, corner;
#ifdef P4_TO_P8
  int                 edge;
#endif
  p4est_topidx_t      num_trees, tt, nt;
  p4est_qcoord_t      coords[P4EST_DIM];
  p4est_qcoord_t      coords_out[P4EST_DIM];
  p4est_quadrant_t    root, *q, node, *r;

  mpiret = sc_MPI_Comm_size (mpicomm, &size);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
  num_trees = conn->num_trees;

  if (rank != 0) {
    /* the connectivity is the same on every rank */
    return;
  }

  /* initialize quadrant storage */
  p4est_quadrant_root (q = &root);
  r = &node;

  /* loop over all trees in the connectivity */
  for (tt = 0; tt < num_trees; ++tt) {
    P4EST_INFOF ("Going through tree %ld\n", (long) tt);

    /* verify volume midpoint */
    p4est_quadrant_volume_coordinates (q, coords);
    p4est_connectivity_coordinates_canonicalize
      (conn, tt, coords, &nt, coords_out);
    SC_CHECK_ABORT (nt == tt, "Mysterious volume tree");
    SC_CHECK_ABORT (!p4est_coordinates_compare (coords, coords_out),
                    "Mysterious volume coordinates");

    for (face = 0; face < P4EST_FACES; ++face) {
      /* verify face midpoints */
      p4est_quadrant_face_coordinates (q, face, coords);
      p4est_connectivity_coordinates_canonicalize
        (conn, tt, coords, &nt, coords_out);
      SC_CHECK_ABORT (nt <= tt, "Mysterious face tree");
      SC_CHECK_ABORT (nt < tt ||
                      p4est_coordinates_compare (coords, coords_out) >= 0,
                      "Mysterious face coordinates");
    }

#ifdef P4_TO_P8
    for (edge = 0; edge < P8EST_EDGES; ++edge) {
      /* verify edge midpoints */
      p8est_quadrant_edge_coordinates (q, edge, coords);
      p4est_connectivity_coordinates_canonicalize
        (conn, tt, coords, &nt, coords_out);
      SC_CHECK_ABORT (nt <= tt, "Mysterious edge tree");
      SC_CHECK_ABORT (nt < tt ||
                      p4est_coordinates_compare (coords, coords_out) >= 0,
                      "Mysterious edge coordinates");
    }
#endif

    for (corner = 0; corner < P4EST_CHILDREN; ++corner) {
      /* verify tree corners */
      p4est_quadrant_corner_node (q, corner, r);
      p4est_quadrant_corner_coordinates (q, corner, coords);
      SC_CHECK_ABORT (test_node_coordinates (r, coords), "Node coordinates");
      p4est_connectivity_coordinates_canonicalize
        (conn, tt, coords, &nt, coords_out);
      SC_CHECK_ABORT (nt <= tt, "Mysterious corner tree");
      SC_CHECK_ABORT (nt < tt ||
                      p4est_coordinates_compare (coords, coords_out) >= 0,
                      "Mysterious corner coordinates");
    }
  }
}

static void
test_coordinates (sc_MPI_Comm mpicomm)
{
  p4est_connectivity_t *conn;

#ifndef P4_TO_P8
  conn = p4est_connectivity_new_corner ();
#else
  conn = p8est_connectivity_new_rotcubes ();
#endif

  test_connectivity (mpicomm, conn);

  p4est_connectivity_destroy (conn);
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;

  /* initialize MPI environment */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  /* establish parallel logging */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* proceed with tests */
  test_coordinates (mpicomm);

  /* clean up and exit */
  sc_finalize ();

  /* finalize MPI environment */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  /* by default exit cleanly */
  return EXIT_SUCCESS;
}
