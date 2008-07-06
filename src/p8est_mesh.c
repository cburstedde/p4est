/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2008 Carsten Burstedde, Lucas Wilcox.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <p8est_mesh.h>
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>

#include <p4est_to_p8est.h>

/** Get the small face neighbors of \a q.
 *
 * Gets the small face neighbors, which are half of the size assuming the
 * 2-1 constant.
 *
 * The order of the \a n[i] is given in the Morton ordering.
 *
 * \param [in]  q      The quadrant whose face neighbors will be constructed.
 * \param [in]  face   The face across which to generate the neighbors.
 * \param [out] n[0]..n[3] Filled with the four smaller face neighbors.
 * \param [out] nur[0]..nur[3] If not NULL, filled with smallest quadrants
 *                     that fit in the upper right corners of \a n.
 */
static void
p4est_quadrant_get_half_face_neighbors (const p4est_quadrant_t * q,
                                        int face, p4est_quadrant_t n[],
                                        p4est_quadrant_t nur[])
{
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);
  const p4est_qcoord_t qh_2 = P4EST_QUADRANT_LEN (q->level + 1);
  int                 i;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (q->level < P4EST_MAXLEVEL);
  P4EST_ASSERT (0 <= face && face < 2 * P4EST_DIM);

  n[0].x = q->x + ((face == 0) ? -qh_2 : (face == 1) ? qh : 0);
  n[0].y = q->y + ((face == 2) ? -qh_2 : (face == 3) ? qh : 0);
  n[0].z = q->z + ((face == 4) ? -qh_2 : (face == 5) ? qh : 0);
  switch (face / 2) {
  case 0:
    for (i = 1; i < 4; ++i) {
      n[i].x = n[0].x;
      n[i].y = n[0].y + (i & 0x01) * qh_2;
      n[i].z = n[0].z + ((i & 0x02) / 2) * qh_2;
    }
    break;
  case 1:
    for (i = 1; i < 4; ++i) {
      n[i].x = n[0].x + (i & 0x01) * qh_2;
      n[i].y = n[0].y;
      n[i].z = n[0].z + ((i & 0x02) / 2) * qh_2;
    }
    break;
  case 2:
    for (i = 1; i < 4; ++i) {
      n[i].x = n[0].x + (i & 0x01) * qh_2;
      n[i].y = n[0].y + ((i & 0x02) / 2) * qh_2;
      n[i].z = n[0].z;
    }
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }
  for (i = 0; i < 4; ++i) {
    n[i].level = (int8_t) (q->level + 1);
    P4EST_ASSERT (p4est_quadrant_is_extended (&n[i]));
  }

  if (nur != NULL) {
    const p4est_qcoord_t th = P4EST_QUADRANT_LEN (P4EST_MAXLEVEL);

    for (i = 0; i < 4; ++i) {
      nur[i].x = n[i].x + qh_2 - th;
      nur[i].y = n[i].y + qh_2 - th;
      nur[i].z = n[i].z + qh_2 - th;
      nur[i].level = P4EST_MAXLEVEL;
      P4EST_ASSERT (p4est_quadrant_is_extended (&nur[i]));
    }
  }
}

/** Get all possible face neighbors of \a q.
 *
 * Gets the face neighbors, possible assuming the 2-1 constraint.
 * If the larger quadrant doesn't exist than it is returned
 * as initialized by P4EST_QUADRANT_INIT.
 *
 * The order of \a n[0] through \a n[3] are given in Morton ordering.
 *
 * \param [in]  q      The quadrant whose face neighbors will be constructed.
 * \param [in]  face   The face across which to generate the neighbors.
 * \param [out] n[0]..n[3] Filled with the smaller possible face neighbors,
 *                     which are half of the size if they exist
 *                     or initialized to P4EST_QUADRANT_INIT.
 * \param [out] n[4]   Filled with the face neighbor, which is the same size.
 * \param [out] n[5]   Filled with the face neighbor, which is twice the size
 *                     if it exists or initialized to P4EST_QUADRANT_INIT.
 */
static void
p4est_quadrant_get_possible_face_neighbors (const p4est_quadrant_t * q,
                                            int face, p4est_quadrant_t n[])
{
  const int           qcid = p4est_quadrant_child_id (q);
  p4est_quadrant_t   *r = &n[5];

  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (q->level == P4EST_MAXLEVEL) {
    P4EST_QUADRANT_INIT (&n[0]);
    P4EST_QUADRANT_INIT (&n[1]);
  }
  else {
    p4est_quadrant_get_half_face_neighbors (q, face, n, NULL);
  }

  p4est_quadrant_face_neighbor (q, face, &n[4]);

  /* Check to see if the larger neighbor exists */
  if (((qcid >> (face / 2)) & 0x01) != (face & 0x01) || q->level == 0) {
    P4EST_QUADRANT_INIT (r);
  }
  else {
    p4est_quadrant_parent (q, r);
    p4est_quadrant_face_neighbor (r, face, r);
  }
}

#include "p4est_mesh.c"

p4est_neighborhood_t *
p4est_neighborhood_new (p4est_t * p4est)
{
  p4est_topidx_t      local_num_trees, flt, nt;
  p4est_locidx_t      local_num_quadrants, lsum;
  p4est_tree_t       *tree;
  p4est_neighborhood_t *nhood;

  P4EST_ASSERT (p4est_is_valid (p4est));

  if (p4est->first_local_tree < 0) {
    flt = 0;
    local_num_trees = 0;
  }
  else {
    flt = p4est->first_local_tree;
    local_num_trees = p4est->last_local_tree - flt + 1; /* type ok */
  }
  local_num_quadrants = p4est->local_num_quadrants;

  nhood = P4EST_ALLOC (p4est_neighborhood_t, 1);
  nhood->cumulative_count = P4EST_ALLOC (p4est_locidx_t, local_num_trees + 1);
  nhood->element_offsets =
    P4EST_ALLOC (p4est_locidx_t, local_num_quadrants + 1);
  nhood->local_neighbors = sc_array_new (sizeof (p4est_locidx_t));

  lsum = 0;
  for (nt = 0; nt < local_num_trees; ++nt) {
    nhood->cumulative_count[nt] = lsum;
    tree = p4est_array_index_topidx (p4est->trees, flt + nt);   /* type ok */
    lsum += (p4est_locidx_t) tree->quadrants.elem_count;        /* type ok */
  }
  P4EST_ASSERT (lsum == local_num_quadrants);
  nhood->cumulative_count[nt] = lsum;

  return nhood;
}

void
p4est_neighborhood_destroy (p4est_neighborhood_t * nhood)
{
  P4EST_FREE (nhood->cumulative_count);
  P4EST_FREE (nhood->element_offsets);
  sc_array_destroy (nhood->local_neighbors);

  P4EST_FREE (nhood);
}

/* EOF p8est_mesh.c */
