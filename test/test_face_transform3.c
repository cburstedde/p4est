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
 * Purpose of this program is to verify the face transformation codes
 * that are computed in a fast but cryptic way in p8est_find_face_transform.
 */

#include <p8est.h>
#include <p4est_to_p8est.h>

int
main (int argc, char **argv)
{
  int                 my_face, target_face, orientation;
  int                 face_ref, face_perm;
  int                 low[2], high[2], swap;
  int                 i, reverse;
  int                 ft[9], gt[9];
  int                *my_axis = &ft[0];
  int                *target_axis = &ft[3];
  int                *edge_reverse = &ft[6];

  sc_init (sc_MPI_COMM_NULL, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  for (my_face = 0; my_face < 2 * P4EST_DIM; ++my_face) {
    for (target_face = 0; target_face < 2 * P4EST_DIM; ++target_face) {
      for (orientation = 0; orientation < 4; ++orientation) {

        /* find if my edges 0 and 2 are parallel to the x, y, or z-axis */
        my_axis[0] = p8est_face_edges[my_face][0] / 4;
        my_axis[1] = p8est_face_edges[my_face][2] / 4;
        target_axis[0] = target_axis[1] = -1;
        edge_reverse[0] = edge_reverse[1] = 0;

        /* find matching target vertices */
        face_ref = p8est_face_permutation_refs[my_face][target_face];
        face_perm = p8est_face_permutation_sets[face_ref][orientation];
        low[0] = low[1] =
          p8est_face_corners[target_face][p8est_face_permutations[face_perm]
                                          [0]];
        high[0] =
          p8est_face_corners[target_face][p8est_face_permutations[face_perm]
                                          [1]];
        high[1] =
          p8est_face_corners[target_face][p8est_face_permutations[face_perm]
                                          [2]];
        if (low[0] > high[0]) {
          swap = low[0];
          low[0] = high[0];
          high[0] = swap;
          edge_reverse[0] = 1;
        }
        if (low[1] > high[1]) {
          swap = low[1];
          low[1] = high[1];
          high[1] = swap;
          edge_reverse[1] = 1;
        }

        /* find matching target edges */
        for (i = 0; i < 12; ++i) {
          if (low[0] == p8est_edge_corners[i][0] &&
              high[0] == p8est_edge_corners[i][1]) {
            P4EST_ASSERT (target_axis[0] == -1);
            target_axis[0] = i / 4;
#ifndef P4EST_ENABLE_DEBUG
            if (target_axis[1] >= 0)
              break;
#endif
          }
          else if (low[1] == p8est_edge_corners[i][0] &&
                   high[1] == p8est_edge_corners[i][1]) {
            P4EST_ASSERT (target_axis[1] == -1);
            target_axis[1] = i / 4;
#ifndef P4EST_ENABLE_DEBUG
            if (target_axis[0] >= 0)
              break;
#endif
          }
        }

        /* find what axis is normal to the faces */
        my_axis[2] = my_face / 2;
        target_axis[2] = target_face / 2;
        edge_reverse[2] = 2 * (my_face % 2) + target_face % 2;

#ifdef P4EST_ENABLE_DEBUG
        for (i = 0; i < 3; ++i) {
          P4EST_ASSERT (0 <= my_axis[i] && my_axis[i] < 3);
          P4EST_ASSERT (0 <= target_axis[i] && target_axis[i] < 3);
        }
        P4EST_ASSERT (my_axis[0] != my_axis[1] &&
                      my_axis[0] != my_axis[2] && my_axis[1] != my_axis[2]);
        P4EST_ASSERT (target_axis[0] != target_axis[1] &&
                      target_axis[0] != target_axis[2] &&
                      target_axis[1] != target_axis[2]);
#endif

        /* output the results */
        P4EST_LDEBUGF
          ("Results for %d %d %d are %d %d %d %d %d %d %d %d %d\n",
           my_face, target_face, orientation, ft[0], ft[1], ft[2],
           ft[3], ft[4], ft[5], ft[6], ft[7], ft[8]);

        /* compute the transformation code in a faster way and compare */
        gt[0] = my_face < 2 ? 1 : 0;
        gt[1] = my_face < 4 ? 2 : 1;
        gt[2] = my_face / 2;
        reverse =
          p8est_face_permutation_refs[0][my_face] ^
          p8est_face_permutation_refs[0][target_face] ^
          (orientation == 0 || orientation == 3);
        gt[3 + reverse] = target_face < 2 ? 1 : 0;
        gt[3 + !reverse] = target_face < 4 ? 2 : 1;
        gt[5] = target_face / 2;
        reverse = p8est_face_permutation_refs[my_face][target_face] == 1;
        gt[6 + reverse] = orientation % 2;
        gt[6 + !reverse] = orientation / 2;
        gt[8] = 2 * (my_face % 2) + target_face % 2;

        /* ensure that both computations yield the same result */
        SC_CHECK_ABORT (!memcmp (ft, gt, 9 * sizeof (int)), "Mismatch");
      }
    }
  }

  sc_finalize ();

  return 0;
}
