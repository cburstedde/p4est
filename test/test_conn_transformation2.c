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
#else /* !P4_TO_P8 */
#include <p8est_connectivity.h>
#endif /* !P4_TO_P8 */

/** Checks that orientation is properly set, i.e. face corner 0 of the
 * face with the lower face index is touching face corner \a
 * orientation of the face with the higher face index.
 * \param[in] conn        p4est_connectivity structure two trees
 * \param[in] l_face      left face index
 * \param[in] r_face      right face index
 * \param[in] orientation the orientation that has been set
 */
static int
test_conn_transformation_check_orientation (p4est_connectivity_t * conn,
                                            int l_face, int r_face,
                                            int orientation)
{
  int                 neighboring_face_corner, corner_index;
  int                 lowerFaceIndex, higherFaceIndex;
  if (l_face <= r_face) {
    lowerFaceIndex = l_face;
    higherFaceIndex = r_face;
  }
  else {
    lowerFaceIndex = r_face;
    higherFaceIndex = l_face;
  }

  corner_index = p4est_face_corners[lowerFaceIndex][0];
  neighboring_face_corner =
    p4est_connectivity_face_neighbor_corner_orientation (corner_index,
                                                         lowerFaceIndex,
                                                         higherFaceIndex,
                                                         orientation);

  neighboring_face_corner =
    p4est_corner_face_corners[neighboring_face_corner][higherFaceIndex];
  P4EST_ASSERT (neighboring_face_corner == orientation);

  return 0;
}

/** Checks for each face corner if the corner indices match on both
 * sides.
 * Let face corner fci, corresponding to corner index ci, be
 * adjacent to face corner fcj, corresponding to corner index cj. We
 * test if fci is seen from fcj and vice versa.
 * \param [in] conn        p4est_connectivity structure two trees
 * \param [in] l_face      left face index
 * \param [in] r_face      right face index
 * \param [in] orientation the orientation that has been set
 */
static int
test_conn_transformation_check_face_corners (p4est_connectivity_t * conn,
                                             int l_face, int r_face,
                                             int orientation)
{
  int                 c0, c1, cx;
  int                 i;
  int                 lowerFaceIndex, higherFaceIndex;

  /* swap face indices if necessary */
  if (l_face <= r_face) {
    lowerFaceIndex = l_face;
    higherFaceIndex = r_face;
  }
  else {
    lowerFaceIndex = r_face;
    higherFaceIndex = l_face;
  }

  /* verify bijectivity of transformation */
  for (i = 0; i < P4EST_HALF; ++i) {
    c0 = p4est_face_corners[lowerFaceIndex][i];
    c1 =
      p4est_connectivity_face_neighbor_corner_orientation (c0, lowerFaceIndex,
                                                           higherFaceIndex,
                                                           orientation);
    cx =
      p4est_connectivity_face_neighbor_corner_orientation (c1,
                                                           higherFaceIndex,
                                                           lowerFaceIndex,
                                                           orientation);

    P4EST_ASSERT (c0 == cx);
  }

  return 0;
}

#ifdef P4_TO_P8
/** Checks for each face edge if the edge indices match on both
 * sides.
 * Let face edge fei, corresponding to edge index ei, be
 * adjacent to face edge fej, corresponding to edge index ej. We
 * test if fei is seen from fej and vice versa.
 * \param [in] conn        p4est_connectivity structure two trees
 * \param [in] l_face      left face index
 * \param [in] r_face      right face index
 * \param [in] orientation the orientation that has been set
 */
static int
test_conn_transformation_check_face_edges (p4est_connectivity_t * conn,
                                           int l_face, int r_face,
                                           int orientation)
{
  int                 e0, e1, ex;
  int                 i;
  int                 lowerFaceIndex, higherFaceIndex;

  /* swap face indices if necessary */
  if (l_face <= r_face) {
    lowerFaceIndex = l_face;
    higherFaceIndex = r_face;
  }
  else {
    lowerFaceIndex = r_face;
    higherFaceIndex = l_face;
  }

  /* verify bijectivity of transformation */
  for (i = 0; i < P4EST_HALF; ++i) {
    e0 = p8est_face_edges[lowerFaceIndex][i];
    e1 =
      p8est_connectivity_face_neighbor_edge_orientation (e0, lowerFaceIndex,
                                                         higherFaceIndex,
                                                         orientation);
    ex =
      p8est_connectivity_face_neighbor_edge_orientation (e1, higherFaceIndex,
                                                         lowerFaceIndex,
                                                         orientation);

    P4EST_ASSERT (e0 == ex);
  }

  return 0;
}
#endif /* P4_TO_P8 */

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  p4est_connectivity_t *conn;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* initialize libsc and p4est */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  conn = 0;

  int                 i, j;
  int                 k;

  for (i = 0; i < P4EST_FACES; ++i) {   /* set l_face */
    for (j = 0; j < P4EST_FACES; ++j) { /* set r_face */
      for (k = 0; k < P4EST_HALF; ++k) {        /* set orientation */
        P4EST_ASSERT (conn == NULL);

        /* create connectivity structure */
        conn = p4est_connectivity_new_twotrees (i, j, k);

        test_conn_transformation_check_orientation (conn, i, j, k);
        test_conn_transformation_check_face_corners (conn, i, j, k);
#ifdef P4_TO_P8
        test_conn_transformation_check_face_edges (conn, i, j, k);
#endif /* P4_TO_P8 */

        p4est_connectivity_destroy (conn);
        conn = 0;
      }
    }
  }

  /* exit */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
