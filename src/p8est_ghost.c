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
#include <p8est_bits.h>
#include <p8est_communication.h>

#ifdef P4EST_ENABLE_MPI

/** Gets the procids of the owners of \a q.
 *
 * For a quadrant across the edge of a tree has possibly multiple
 * trees in which it lives, and thus multiple procs.
 *
 * \param [in]     p4est    The forest in which to search for \a q.
 * \param [in]     treeid   The tree id to which \a q belongs.
 * \param [in]     edge     The edge of the tree \a q is across from.
 * \param [in]     q        The quadrant that is being searched for.
 * \param [in,out] qprocs   Starts as an initialized array and ends with
 *                          the list of processors that \a q belongs too.
 * \param [out]    nurgood  If not NULL, check if the smallest quadrant
 *                          in the upper right corner of q after transform
 *                          has the same owner.
 */
static void
p8est_quadrant_find_tree_edge_owners (p4est_t * p4est,
                                      p4est_topidx_t treeid,
                                      int edge,
                                      const p4est_quadrant_t * q,
                                      sc_array_t * q_procs, int *nurgood)
{
  const int           rank = p4est->mpirank;
  int                *proc, nurproc;
  size_t              etree;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_quadrant_t    eq;
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  sc_array_t         *eta;

  P4EST_ASSERT (p8est_quadrant_is_outside_edge (q));

  P4EST_QUADRANT_INIT (&eq);

  /* Find all edges that are not myself or from a face neighbor */
  eta = &ei.edge_transforms;
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
  p8est_find_edge_transform (conn, treeid, edge, &ei);

  sc_array_resize (q_procs, 0);
  if (nurgood != NULL) {
    *nurgood = 1;
    if (q->level == P4EST_QMAXLEVEL)
      nurgood = NULL;
  }

  for (etree = 0; etree < eta->elem_count; ++etree) {
    et = p8est_edge_array_index (eta, etree);

    p8est_quadrant_transform_edge (q, &eq, &ei, et, 1);

    proc = (int *) sc_array_push (q_procs);
    *proc = p4est_comm_find_owner (p4est, et->ntree, &eq, rank);

    if (nurgood != NULL) {
      p4est_quadrant_last_descendant (&eq, &eq, P4EST_QMAXLEVEL);
      nurproc = p4est_comm_find_owner (p4est, et->ntree, &eq, *proc);
      *nurgood = *nurgood && (nurproc == *proc);
    }
  }

  sc_array_reset (eta);
}

#endif /* P4EST_ENABLE_MPI */

/** Get the small edge neighbors of \a q.
 *
 * Gets the two small edge neighbors, which are half of the size assuming
 * the 2-1 constant.
 *
 * The order of the \a n[i] is given in the Morton ordering.
 *
 * \param [in]  q      The quadrant whose face neighbors will be constructed.
 * \param [in]  edge   The edge across which to generate the neighbors.
 * \param [out] n[0],n[1] Filled with the two smaller face neighbors.
 * \param [out] nur[0],nur[1] If not NULL, filled with smallest quadrants
 *                     that fit in the upper right corners of \a n.
 */
static void
p8est_quadrant_get_half_edge_neighbors (const p4est_quadrant_t * q,
                                        int edge, p4est_quadrant_t n[],
                                        p4est_quadrant_t nur[])
{
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);
  const p4est_qcoord_t qh_2 = P4EST_QUADRANT_LEN (q->level + 1);

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (q->level < P4EST_QMAXLEVEL);
  P4EST_ASSERT (0 <= edge && edge < 12);

  switch (edge / 4) {
  case 0:
    n[0].x = n[1].x = q->x;
    n[0].y = n[1].y = q->y + (!(edge & 0x01) ? -qh_2 : qh);
    n[0].z = n[1].z = q->z + (!(edge & 0x02) ? -qh_2 : qh);
    n[1].x += qh_2;
    break;
  case 1:
    n[0].x = n[1].x = q->x + (!(edge & 0x01) ? -qh_2 : qh);
    n[0].y = n[1].y = q->y;
    n[0].z = n[1].z = q->z + (!(edge & 0x02) ? -qh_2 : qh);
    n[1].y += qh_2;
    break;
  case 2:
    n[0].x = n[1].x = q->x + (!(edge & 0x01) ? -qh_2 : qh);
    n[0].y = n[1].y = q->y + (!(edge & 0x02) ? -qh_2 : qh);
    n[0].z = n[1].z = q->z;
    n[1].z += qh_2;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }
  n[0].level = n[1].level = (int8_t) (q->level + 1);
  P4EST_ASSERT (p4est_quadrant_is_extended (&n[0]));
  P4EST_ASSERT (p4est_quadrant_is_extended (&n[1]));

  if (nur != NULL) {
    const p4est_qcoord_t dh = qh_2 - P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);

    nur[0].x = n[0].x + dh;
    nur[0].y = n[0].y + dh;
    nur[0].z = n[0].z + dh;
    nur[0].level = P4EST_QMAXLEVEL;
    P4EST_ASSERT (p4est_quadrant_is_extended (&nur[0]));
    nur[1].x = n[1].x + dh;
    nur[1].y = n[1].y + dh;
    nur[1].z = n[1].z + dh;
    nur[1].level = P4EST_QMAXLEVEL;
    P4EST_ASSERT (p4est_quadrant_is_extended (&nur[1]));
  }
}

/** Get all possible edge neighbors of \a q.
 *
 * Gets the edge neighbors, possible assuming the 2-1 constraint.
 * If the larger quadrant doesn't exist than it is returned
 * as initialized by P4EST_QUADRANT_INIT.
 *
 * The order of \a n[0], \a n[1] is given in Morton ordering.
 *
 * \param [in]  q      The quadrant whose edge neighbors will be constructed.
 * \param [in]  edge   The edge across which to generate the neighbors.
 * \param [out] n[0],n[1] Filled with the smaller possible edge neighbors,
 *                     which are half of the size if they exist
 *                     or initialized to P4EST_QUADRANT_INIT.
 * \param [out] n[2]   Filled with the edge neighbor, which is the same size.
 * \param [out] n[3]   Filled with the edge neighbor, which is twice the size
 *                     if it exists or initialized to P4EST_QUADRANT_INIT.
 */
static void
p8est_quadrant_get_possible_edge_neighbors (const p4est_quadrant_t * q,
                                            int edge, p4est_quadrant_t n[])
{
  const int           qcid = p4est_quadrant_child_id (q);
  p4est_quadrant_t   *r = &n[3];

  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (q->level == P4EST_QMAXLEVEL) {
    P4EST_QUADRANT_INIT (&n[0]);
    P4EST_QUADRANT_INIT (&n[1]);
  }
  else {
    p8est_quadrant_get_half_edge_neighbors (q, edge, n, NULL);
  }

  p8est_quadrant_edge_neighbor (q, edge, &n[2]);

  /* Check to see if the larger neighbor exists */
  if ((qcid != p8est_edge_corners[edge][0] &&
       qcid != p8est_edge_corners[edge][1]) || q->level == 0) {
    P4EST_QUADRANT_INIT (r);
  }
  else {
    p4est_quadrant_parent (q, r);
    p8est_quadrant_edge_neighbor (r, edge, r);
  }
}

/** Checks if a quadrant's edge is on the boundary of the forest.
 *
 * This means that the quadrant's tree doesn't have any non face neighbors.
 *
 * \param [in] p4est  The forest in which to search for \a q
 * \param [in] treeid The tree id for which \a q belongs.
 * \param [in] q      The quadrant that is in question.
 * \param [in] edge   The edge of quadrant that is in question.
 *
 * \return true if the quadrant's edge is on the boundary of the forest and
 *         false otherwise.
 */
static int
p8est_quadrant_on_edge_boundary (p4est_t * p4est, p4est_topidx_t treeid,
                                 int edge, const p4est_quadrant_t * q)
{
  int                 face;
  int                 on_boundary;
  p4est_quadrant_t    q2;
  p4est_connectivity_t *conn = p4est->connectivity;
  p8est_edge_info_t   ei;
  sc_array_t         *eta;

  P4EST_ASSERT (0 <= edge && edge < 12);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (p8est_quadrant_touches_edge (q, edge, 1)) {
    eta = &ei.edge_transforms;
    sc_array_init (eta, sizeof (p8est_edge_transform_t));
    p8est_find_edge_transform (conn, treeid, edge, &ei);

    on_boundary = (eta->elem_count == 0);
    sc_array_reset (eta);

    return on_boundary;
  }

  P4EST_QUADRANT_INIT (&q2);
  p8est_quadrant_edge_neighbor (q, edge, &q2);
  P4EST_ASSERT (!p4est_quadrant_is_outside_corner (&q2));
  P4EST_ASSERT (!p8est_quadrant_is_outside_edge (&q2));
  if (q2.x < 0) {
    face = 0;
  }
  else if (q2.x >= P4EST_ROOT_LEN) {
    face = 1;
  }
  else if (q2.y < 0) {
    face = 2;
  }
  else if (q2.y >= P4EST_ROOT_LEN) {
    face = 3;
  }
  else if (q2.z < 0) {
    face = 4;
  }
  else if (q2.z >= P4EST_ROOT_LEN) {
    face = 5;
  }
  else {
    return 0;
  }

  return
    (conn->tree_to_tree[P4EST_FACES * treeid + face] == treeid &&
     (int) conn->tree_to_face[P4EST_FACES * treeid + face] == face);
}

#include "p4est_ghost.c"
