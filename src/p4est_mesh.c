/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

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

#ifndef P4_TO_P8

#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_mesh.h>
#endif /* !P4_TO_P8 */

/* *INDENT-OFF* */
#ifndef P4_TO_P8
static const int    p4est_hanging_corner[4][2] = {
  { 1, 2 },
  { 0, 3 },
  { 0, 3 },
  { 1, 2 }};
static const int    p4est_hanging_face[4][2] = {
  { 0, 3 },
  { 0, 1 },
  { 3, 2 },
  { 1, 2 }};
#else
static const int    p8est_child_edge_face[8][12] = {
  { -1,  4,  2, -1, -1,  4,  0, -1, -1,  2,  0, -1 },
  { -1,  4,  2, -1,  4, -1, -1,  1,  2, -1, -1,  1 },
  {  4, -1, -1,  3, -1,  4,  0, -1,  0, -1, -1,  3 },
  {  4, -1, -1,  3,  4, -1, -1,  1, -1,  1,  3, -1 },
  {  2, -1, -1,  5,  0, -1, -1,  5, -1,  2,  0, -1 },
  {  2, -1, -1,  5, -1,  1,  5, -1,  2, -1, -1,  1 },
  { -1,  3,  5, -1,  0, -1, -1,  5,  0, -1, -1,  3 },
  { -1,  3,  5, -1, -1,  1,  5, -1, -1,  1,  3, -1 }};
static const int    p8est_child_corner_face[8][8] = {
  { -1, -1, -1,  4, -1,  2,  0, -1 },
  { -1, -1,  4, -1,  2, -1, -1,  1 },
  { -1,  4, -1, -1,  0, -1, -1,  3 },
  {  4, -1, -1, -1, -1,  1,  3, -1 },
  { -1,  2,  0, -1, -1, -1, -1,  5 },
  {  2, -1, -1,  1, -1, -1,  5, -1 },
  {  0, -1, -1,  3, -1,  5, -1, -1 },
  { -1,  1,  3, -1,  5, -1, -1, -1 }};
static const int    p8est_child_corner_edge[8][8] = {
  { -1,  0,  4, -1,  8, -1, -1, -1 },
  {  0, -1, -1,  5, -1,  9, -1, -1 },
  {  4, -1, -1,  1, -1, -1, 10, -1 },
  { -1,  5,  1, -1, -1, -1, -1, 11 },
  {  8, -1, -1, -1, -1,  2,  6, -1 },
  { -1,  9, -1, -1,  2, -1, -1,  7 },
  { -1, -1, 10, -1,  6, -1, -1,  3 },
  { -1, -1, -1, 11, -1,  7,  3, -1 }};
#endif
/* *INDENT-ON* */

#ifdef P4EST_MPI

/** Gets the procid of the owner of \a q.
 *
 * \param [in] p4est  The forest in which to search for \a q
 * \param [in] treeid The tree id for which \a q belongs.
 * \param [in] q      The quadrant that is being searched for.
 *
 * \return Procid of the owner of \a q
 *                or -1 if the quadrant lies outside of the mesh.
 *
 * \warning Does not work for tree edge or corner neighbors.
 */
static int
p4est_quadrant_find_owner (p4est_t * p4est, p4est_topidx_t treeid,
                           const p4est_quadrant_t * q)
{
  const int           rank = p4est->mpirank;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  p4est_connectivity_t *conn = p4est->connectivity;
  bool                quad_contact[2 * P4EST_DIM];
  int                 face;
  p4est_topidx_t      ntreeid;
  p4est_quadrant_t    nq;
#ifndef P4_TO_P8
  int                 transform;
  p4est_quadrant_t    tmpq;
#else
  int                 ftransform[9];
  p4est_topidx_t      ntreeid2;
#endif

  if (p4est_quadrant_is_inside_root (q)) {
    return p4est_comm_find_owner (p4est, treeid, q, rank);
  }

  P4EST_QUADRANT_INIT (&nq);

  /* We are outside of the unit tree */
#ifndef P4_TO_P8
  quad_contact[0] = (q->y < 0);
  quad_contact[1] = (q->x >= rh);
  quad_contact[2] = (q->y >= rh);
  quad_contact[3] = (q->x < 0);

  /* Make sure we are not a tree corner */
  P4EST_ASSERT (!((quad_contact[0] || quad_contact[2]) &&
                  (quad_contact[1] || quad_contact[3])));
#else
  quad_contact[0] = (q->x < 0);
  quad_contact[1] = (q->x >= rh);
  quad_contact[2] = (q->y < 0);
  quad_contact[3] = (q->y >= rh);
  quad_contact[4] = (q->z < 0);
  quad_contact[5] = (q->z >= rh);
  P4EST_ASSERT (((quad_contact[0] || quad_contact[1]) ? 1 : 0) +
                ((quad_contact[2] || quad_contact[3]) ? 1 : 0) +
                ((quad_contact[4] || quad_contact[5]) ? 1 : 0) == 1);
#endif

  ntreeid = -1;
  for (face = 0; face < 2 * P4EST_DIM; ++face) {
    if (quad_contact[face]) {
      ntreeid = conn->tree_to_tree[2 * P4EST_DIM * treeid + face];
      if (ntreeid == treeid
          && ((int) conn->tree_to_face[2 * P4EST_DIM * treeid + face] ==
              face)) {
        /* This quadrant goes across a face with no neighbor */
        return -1;
      }
      break;
    }
  }
  P4EST_ASSERT (face < 2 * P4EST_DIM && ntreeid >= 0);

#ifndef P4_TO_P8
  transform = p4est_find_face_transform (conn, treeid, face);
  tmpq = *q;
  p4est_quadrant_translate_face (&tmpq, face);
  p4est_quadrant_transform_face (&tmpq, &nq, transform);
#else
  ntreeid2 = p8est_find_face_transform (conn, treeid, face, ftransform);
  P4EST_ASSERT (ntreeid2 == ntreeid);
  p8est_quadrant_transform_face (q, &nq, ftransform);
#endif

  return p4est_comm_find_owner (p4est, ntreeid, &nq, rank);
}

/** Gets the procids of the owners of \a q.
 *
 * For a quadrant across the corner of a tree has possibly multiple
 * trees in which it lives, and thus multiple procs.
 *
 * \param [in]     p4est      The forest in which to search for \a q.
 * \param [in]     treeid     The tree id to which \a q belongs.
 * \param [in]     treecorner The r-corner of the tree \a q is across from.
 * \param [in]     q          The quadrant that is being searched for.
 * \param [in,out] qprocs     Starts as an initialized array and ends with
 *                            the list of processors that \a q belongs too.
 * \param [out]    nurgood    If not NULL, check if the smallest quadrant
 *                            in the upper right corner of q after transform
 *                            has the same owner.
 */
static void
p4est_quadrant_find_tree_corner_owners (p4est_t * p4est,
                                        p4est_topidx_t treeid,
                                        int treecorner,
                                        const p4est_quadrant_t * q,
                                        sc_array_t * q_procs, bool * nurgood)
{
  const int           rank = p4est->mpirank;
  int                *proc, nurproc;
  size_t              ctree;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_quadrant_t    cq;
#ifndef P4_TO_P8
  sc_array_t          ctransforms, *cta;
  p4est_corner_transform_t *ct;
#else
  p8est_corner_info_t ci;
  p8est_corner_transform_t *ct;
  sc_array_t         *cta;
#endif

  P4EST_ASSERT (p4est_quadrant_is_outside_corner (q));

  P4EST_QUADRANT_INIT (&cq);

  /* Find all corners that are not myself or from a face neighbor */
#ifndef P4_TO_P8
  cta = &ctransforms;
  sc_array_init (cta, sizeof (p4est_corner_transform_t));
  p4est_find_corner_transform (conn, treeid, treecorner, cta);
#else
  cta = &ci.corner_transforms;
  sc_array_init (cta, sizeof (p8est_corner_transform_t));
  p8est_find_corner_transform (conn, treeid, treecorner, &ci);
#endif

  sc_array_resize (q_procs, 0);
  if (nurgood != NULL) {
    *nurgood = true;
    if (q->level == P4EST_QMAXLEVEL)
      nurgood = NULL;
  }

  for (ctree = 0; ctree < cta->elem_count; ++ctree) {
    ct = sc_array_index (cta, ctree);

    cq = *q;
    p4est_quadrant_transform_corner (&cq, (int) ct->ncorner, true);

    proc = sc_array_push (q_procs);
    *proc = p4est_comm_find_owner (p4est, ct->ntree, &cq, rank);

    if (nurgood != NULL) {
      p4est_quadrant_last_descendent (&cq, &cq, P4EST_QMAXLEVEL);
      nurproc = p4est_comm_find_owner (p4est, ct->ntree, &cq, *proc);
      *nurgood = *nurgood && (nurproc == *proc);
    }
  }

  sc_array_reset (cta);
}

#endif /* P4EST_MPI */

/** Checks if quadrant exists in the local forest or the ghost layer.
 *
 * For quadrants across tree corners/edges it checks if the quadrant exists
 * in any of the corner/edge neighbors.
 *
 * \param [in]  p4est        The forest in which to search for \a q
 * \param [in]  ghost_layer  The ghost layer in which to search for \a q
 * \param [in]  treeid       The tree id for which \a q belongs.
 * \param [in]  q            The quadrant that is being searched for.
 * \param [out] exists_arr   Filled for tree corner/edge cases.  An entry
 *                           for each corner/edge neighbor is set to 1 if
 *                           it exists in the local forest or ghost_layer.
 *
 * \return true if the quadrant exists in the local forest or in the
 *              ghost_layer, and false if doesn't exist in either
 */
static              bool
p4est_quadrant_exists (p4est_t * p4est, sc_array_t * ghost_layer,
                       p4est_topidx_t treeid, const p4est_quadrant_t * q,
                       sc_array_t * exists_arr)
{
  const int           rank = p4est->mpirank;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  int                 qproc, face;
  int                *pexists;
  bool                quad_contact[2 * P4EST_DIM];
  bool                quad_corner;
  bool                exists;
  size_t              ctreeidz;
  ssize_t             lnid;
  p4est_topidx_t      tqtreeid;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_tree_t       *tree = p4est_array_index_topidx (p4est->trees, treeid);
  p4est_tree_t       *tqtree;
  p4est_quadrant_t    tq, non_existent;
  sc_array_t         *quadrants = &tree->quadrants;
  sc_array_t         *ta;
#ifndef P4_TO_P8
  int                 transform, tree_corner;
  p4est_quadrant_t    tempq;
  p4est_corner_transform_t *ct;
  sc_array_t          ctransforms;
#else
  int                 edge, corner;
  int                 ftransform[9];
  bool                face_axis[3];
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  p8est_corner_info_t ci;
  p8est_corner_transform_t *ct;
#endif
#ifdef P4EST_DEBUG
  p4est_quadrant_t   *q2;
#endif

  if (exists_arr != NULL) {
    P4EST_ASSERT (exists_arr->elem_size == sizeof (int));
    sc_array_resize (exists_arr, 0);
  }
  P4EST_QUADRANT_INIT (&non_existent);
  ta = NULL;

  if (non_existent.level == q->level) {
    return false;
  }

  /* q is in the unit domain */
  if (p4est_quadrant_is_inside_root (q)) {
    qproc = p4est_comm_find_owner (p4est, treeid, q, rank);
    if (qproc == rank) {
      lnid = sc_array_bsearch (quadrants, q, p4est_quadrant_compare);
    }
    else {
      /* off processor so search in the ghost layer */
      tq = *q;
      tq.p.piggy1.which_tree = treeid;
      lnid = sc_array_bsearch (ghost_layer, &tq,
                               p4est_quadrant_compare_piggy);
      P4EST_ASSERT (lnid == -1 ||
                    (q2 = sc_array_index_ssize_t (ghost_layer, lnid),
                     q2->p.piggy1.owner_rank == qproc));
    }
    return (lnid != -1);
  }

  /* q is in a neighboring tree */
#ifndef P4_TO_P8
  quad_contact[0] = (q->y < 0);
  quad_contact[1] = (q->x >= rh);
  quad_contact[2] = (q->y >= rh);
  quad_contact[3] = (q->x < 0);
  quad_corner =
    (quad_contact[0] || quad_contact[2]) &&
    (quad_contact[1] || quad_contact[3]);
#else
  quad_contact[0] = (q->x < 0);
  quad_contact[1] = (q->x >= rh);
  quad_contact[2] = (q->y < 0);
  quad_contact[3] = (q->y >= rh);
  quad_contact[4] = (q->z < 0);
  quad_contact[5] = (q->z >= rh);
  face_axis[0] = quad_contact[0] || quad_contact[1];
  face_axis[1] = quad_contact[2] || quad_contact[3];
  face_axis[2] = quad_contact[4] || quad_contact[5];
  quad_corner = false;
  face = edge = corner = -1;
  P4EST_ASSERT (face_axis[0] || face_axis[1] || face_axis[2]);
  if (!face_axis[1] && !face_axis[2]) {
    face = 0 + quad_contact[1];
  }
  else if (!face_axis[0] && !face_axis[2]) {
    face = 2 + quad_contact[3];
  }
  else if (!face_axis[0] && !face_axis[1]) {
    face = 4 + quad_contact[5];
  }
  else if (!face_axis[0]) {
    edge = 0 + 2 * quad_contact[5] + quad_contact[3];
    quad_corner = true;
  }
  else if (!face_axis[1]) {
    edge = 4 + 2 * quad_contact[5] + quad_contact[1];
    quad_corner = true;
  }
  else if (!face_axis[2]) {
    edge = 8 + 2 * quad_contact[3] + quad_contact[1];
    quad_corner = true;
  }
  else {
    corner = 4 * quad_contact[5] + 2 * quad_contact[3] + quad_contact[1];
    quad_corner = true;
  }
#endif
  if (quad_corner) {
    /* Neighbor is across a tree edge or corner */
    P4EST_ASSERT (exists_arr != NULL);

#ifndef P4_TO_P8
    for (tree_corner = 0; tree_corner < 4; ++tree_corner) {
      if (quad_contact[(tree_corner + 3) % 4]
          && quad_contact[tree_corner]) {
        break;
      }
    }
    ta = &ctransforms;
    sc_array_init (ta, sizeof (p4est_corner_transform_t));
    p4est_find_corner_transform (conn, treeid, tree_corner, ta);
#else
    P4EST_ASSERT (face == -1 &&
                  ((edge >= 0 && corner == -1) ||
                   (edge == -1 && corner >= 0)));
    if (edge >= 0) {
      ta = &ei.edge_transforms;
      sc_array_init (ta, sizeof (p8est_edge_transform_t));
      p8est_find_edge_transform (conn, treeid, edge, &ei);
    }
    else {
      ta = &ci.corner_transforms;
      sc_array_init (ta, sizeof (p8est_corner_transform_t));
      p8est_find_corner_transform (conn, treeid, corner, &ci);
    }
#endif
    sc_array_resize (exists_arr, ta->elem_count);

    exists = false;
    for (ctreeidz = 0; ctreeidz < ta->elem_count; ++ctreeidz) {
#ifndef P4_TO_P8
      ct = sc_array_index (ta, ctreeidz);
      tqtreeid = ct->ntree;

      tq = *q;
      p4est_quadrant_transform_corner (&tq, (int) ct->ncorner, true);
#else
      if (edge >= 0) {
        et = sc_array_index (ta, ctreeidz);
        tqtreeid = et->ntree;
        p8est_quadrant_transform_edge (q, &tq, &ei, et, true);
      }
      else {
        ct = sc_array_index (ta, ctreeidz);
        tqtreeid = ct->ntree;

        tq = *q;
        p4est_quadrant_transform_corner (&tq, (int) ct->ncorner, true);
      }
      et = NULL;
      ct = NULL;
#endif

      qproc = p4est_comm_find_owner (p4est, tqtreeid, &tq, rank);

      if (qproc == rank) {
        tqtree = p4est_array_index_topidx (p4est->trees, tqtreeid);
        lnid = sc_array_bsearch (&tqtree->quadrants, &tq,
                                 p4est_quadrant_compare);
      }
      else {
        tq.p.piggy1.which_tree = tqtreeid;
        lnid = sc_array_bsearch (ghost_layer, &tq,
                                 p4est_quadrant_compare_piggy);
        P4EST_ASSERT (lnid == -1 ||
                      (q2 = sc_array_index_ssize_t (ghost_layer, lnid),
                       q2->p.piggy1.owner_rank == qproc));
      }

      /* add the existence value */
      pexists = sc_array_index (exists_arr, ctreeidz);
      *pexists = (lnid != -1);
      exists = exists || *pexists;
    }

    sc_array_reset (ta);
    return exists;
  }
  else {
#ifndef P4_TO_P8
    /* Neighbor is across a tree face */
    tqtreeid = -1;
    for (face = 0; face < 2 * P4EST_DIM; ++face) {
      if (quad_contact[face]) {
        tqtreeid = conn->tree_to_tree[2 * P4EST_DIM * treeid + face];
        if (tqtreeid == treeid
            && ((int) conn->tree_to_face[2 * P4EST_DIM * treeid + face] ==
                face)) {
          return false;
        }
        break;
      }
    }
    P4EST_ASSERT (face < 2 * P4EST_DIM && tqtreeid >= 0);

    /* transform the neighbor into the other tree's
     * coordinates
     */
    tempq = *q;
    transform = p4est_find_face_transform (conn, treeid, face);
    p4est_quadrant_translate_face (&tempq, face);
    p4est_quadrant_transform_face (&tempq, &tq, transform);
#else
    P4EST_ASSERT (face >= 0 && edge == -1 && corner == -1);
    P4EST_ASSERT (quad_contact[face]);

    tqtreeid = p8est_find_face_transform (conn, treeid, face, ftransform);
    if (tqtreeid == -1) {
      /* there is no tree neighbor across this face */
      return false;
    }
    p8est_quadrant_transform_face (q, &tq, ftransform);
#endif

    /* find owner of the transformed quadrant */
    qproc = p4est_comm_find_owner (p4est, tqtreeid, &tq, rank);
    if (qproc == rank) {
      tqtree = p4est_array_index_topidx (p4est->trees, tqtreeid);
      lnid = sc_array_bsearch (&tqtree->quadrants, &tq,
                               p4est_quadrant_compare);
    }
    else {
      /* off processor so search in the ghost layer */
      tq.p.piggy1.which_tree = tqtreeid;
      lnid = sc_array_bsearch (ghost_layer, &tq,
                               p4est_quadrant_compare_piggy);
      P4EST_ASSERT (lnid == -1 ||
                    (q2 = sc_array_index_ssize_t (ghost_layer, lnid),
                     q2->p.piggy1.owner_rank == qproc));
    }
    return (lnid != -1);
  }
}

/** Checks if a quadrant's face is on the boundary of the forest.
 *
 * \param [in] p4est  The forest in which to search for \a q
 * \param [in] treeid The tree id for which \a q belongs.
 * \param [in] q      The quadrant or node that is in question.
 * \param [in] face   The face of quadrant that is in question.
 *
 * \return true if the quadrant's face is on the boundary of the forest and
 *         false otherwise.
 */
static              bool
p4est_quadrant_on_face_boundary (p4est_t * p4est, p4est_topidx_t treeid,
                                 int face, const p4est_quadrant_t * q)
{
  p4est_qcoord_t      dh;
  p4est_connectivity_t *conn = p4est->connectivity;
#ifdef P4_TO_P8
  p4est_qcoord_t      xyz;
#endif

  P4EST_ASSERT (0 <= face && face < 2 * P4EST_DIM);
  if (p4est_quadrant_is_node (q, false)) {
    dh = P4EST_ROOT_LEN;
  }
  else {
    P4EST_ASSERT (p4est_quadrant_is_valid (q));
    dh = P4EST_LAST_OFFSET (q->level);
  }

  if (conn->tree_to_tree[2 * P4EST_DIM * treeid + face] != treeid ||
      (int) conn->tree_to_face[2 * P4EST_DIM * treeid + face] != face) {
    return false;
  }

#ifndef P4_TO_P8
  switch (face) {
  case 0:
    return q->y == 0;
  case 1:
    return q->x == dh;
  case 2:
    return q->y == dh;
  case 3:
    return q->x == 0;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }
#else
  switch (face / 2) {
  case 0:
    xyz = q->x;
    break;
  case 1:
    xyz = q->y;
    break;
  case 2:
    xyz = q->z;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }
  return xyz == ((face & 0x01) ? dh : 0);
#endif
}

/** Get the smallest corner neighbor of \a q.
 *
 * Gets the smallest corner neighbor, which is half of the size assuming the
 * 2-1 constaint.
 *
 * \param [in]  q      The quadrant whose corner neighbor will be constructed.
 * \param [in]  corner The z-corner across which to generate the neighbor.
 * \param [out] n0     Filled with the smallest corner neighbor, which is
 *                     half of the size assuming the 2-1 constaint.
 * \param [out] n0ur   If not NULL, it is filled with smallest quadrant
 *                     that fits in the upper right corner of \a n0.
 */
static void
p4est_quadrant_get_half_corner_neighbor (const p4est_quadrant_t * q,
                                         int corner,
                                         p4est_quadrant_t * n0,
                                         p4est_quadrant_t * n0ur)
{
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);
  const p4est_qcoord_t qh_2 = P4EST_QUADRANT_LEN (q->level + 1);

  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (q->level < P4EST_QMAXLEVEL);

  n0->x = q->x + ((corner & 0x01) ? qh : -qh_2);
  n0->y = q->y + ((corner & 0x02) ? qh : -qh_2);
#ifdef P4_TO_P8
  n0->z = q->z + ((corner & 0x04) ? qh : -qh_2);
#endif
  n0->level = (int8_t) (q->level + 1);
  P4EST_ASSERT (p4est_quadrant_is_extended (n0));

  if (n0ur != NULL) {
    const p4est_qcoord_t dh = qh_2 - P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);

    n0ur->x = n0->x + dh;
    n0ur->y = n0->y + dh;
#ifdef P4_TO_P8
    n0ur->z = n0->z + dh;
#endif
    n0ur->level = P4EST_QMAXLEVEL;
    P4EST_ASSERT (p4est_quadrant_is_extended (n0ur));
  }
}

/** Get the possible corner neighbors of \a q.
 *
 * Gets the corner face neighbors, possible assuming the 2-1 constant.
 * If the larger quadrant doesn't exist than it returned as initialized by
 * P4EST_QUADRANT_INIT.
 *
 * \param [in]  q        The quadrant whose face neighbors will be constructed.
 * \param [in]  corner   The corner across which to generate the neighbors.
 *                       The corner is given in pixel order. So
 *                          2            3
 *                           +----------+
 *                           |          |
 *                           |          |
 *                           |          |
 *                           |          |
 *                           |          |
 *                           +----------+
 *                          0            1
 * \param [out] n[0]   Filled with the possible corner neighbor, which is
 *                     half of the size.
 * \param [out] n[1]   Filled with the face neighbor, which is the same size.
 * \param [out] n[2]   Filled with the face neighbor, which is twice the size
 *                     if it exists or initialized to P4EST_QUADRANT_INIT.
 *
 */
static void
p4est_quadrant_get_possible_corner_neighbors (const p4est_quadrant_t * q,
                                              int corner,
                                              p4est_quadrant_t n[])
{
  const int           qcid = p4est_quadrant_child_id (q);
  p4est_quadrant_t   *r = &n[2];

  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (q->level == P4EST_QMAXLEVEL) {
    P4EST_QUADRANT_INIT (&n[0]);
  }
  else {
    p4est_quadrant_get_half_corner_neighbor (q, corner, &n[0], NULL);
  }

  p4est_quadrant_corner_neighbor (q, corner, &n[1]);

  /* Check to see if the larger neighbor exists */
  if ((corner != qcid) || (q->level == 0)) {
    P4EST_QUADRANT_INIT (r);
  }
  else {
    p4est_quadrant_parent (q, r);
    p4est_quadrant_corner_neighbor (r, corner, r);
  }
}

/** Checks if a quadrant's corner is on the boundary of the forest.
 *
 * This means that the quadrant's tree doesn't have any non face neighbors.
 *
 * \param [in] p4est  The forest in which to search for \a q
 * \param [in] treeid The tree id for which \a q belongs.
 * \param [in] q      The quadrant that is in question.
 * \param [in] corner The z-corner of quadrant that is in question.
 *
 * \return true if the quadrant's corner is on the boundary of the forest and
 *         false otherwise.
 */
static              bool
p4est_quadrant_on_corner_boundary (p4est_t * p4est, p4est_topidx_t treeid,
                                   int corner, const p4est_quadrant_t * q)
{
  int                 face;
  bool                on_boundary = false;
  p4est_quadrant_t    q2;
  p4est_connectivity_t *conn = p4est->connectivity;
#ifndef P4_TO_P8
  int                 rcorner = p4est_corner_to_zorder[corner];
  sc_array_t          ctransforms, *cta;
#else
  int                 edge;
  p8est_edge_info_t   ei;
  p8est_corner_info_t ci;
  sc_array_t         *eta, *cta;
#endif

  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);
  P4EST_ASSERT (p4est_quadrant_is_inside_root (q));

  if (p4est_quadrant_touches_corner (q, corner)) {
#ifndef P4_TO_P8
    cta = &ctransforms;
    sc_array_init (cta, sizeof (p4est_corner_transform_t));
    p4est_find_corner_transform (conn, treeid, rcorner, cta);
#else
    cta = &ci.corner_transforms;
    sc_array_init (cta, sizeof (p8est_corner_transform_t));
    p8est_find_corner_transform (conn, treeid, corner, &ci);
#endif

    on_boundary = (cta->elem_count == 0);
    sc_array_reset (cta);

    return on_boundary;
  }

  P4EST_QUADRANT_INIT (&q2);
  p4est_quadrant_corner_neighbor (q, corner, &q2);
  P4EST_ASSERT (!p4est_quadrant_is_outside_corner (&q2));

#ifdef P4_TO_P8
  if (p8est_quadrant_is_outside_edge_extra (&q2, &edge)) {
    eta = &ei.edge_transforms;
    sc_array_init (eta, sizeof (p8est_edge_transform_t));
    p8est_find_edge_transform (conn, treeid, edge, &ei);

    on_boundary = (eta->elem_count == 0);
    sc_array_reset (eta);

    return on_boundary;
  }
  P4EST_ASSERT (!p8est_quadrant_is_outside_edge (&q2));
#endif

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
#ifdef P4_TO_P8
  else if (q2.z < 0) {
    face = 4;
  }
  else if (q2.z >= P4EST_ROOT_LEN) {
    face = 5;
  }
#endif
  else {
    return false;
  }

#ifndef P4_TO_P8
  face = p4est_zface_to_rface[face];
#endif

  return
    (conn->tree_to_tree[2 * P4EST_DIM * treeid + face] == treeid &&
     (int) conn->tree_to_face[2 * P4EST_DIM * treeid + face] == face);
}

#ifndef P4_TO_P8

/** Get the smallest face neighbors of \a q.
 *
 * Gets the smallest face neighbors, which are half of the size assuming the
 * 2-1 constant.
 *
 * The order of \a n0 and \a n1 are given in the morton ordering.
 *
 * \param [in]  q      The quadrant whose face neighbors will be constructed.
 * \param [in]  face   The face across which to generate the neighbors.  The
 *                     face is given in right hand order. So
 *                                2
 *                           +----------+
 *                           |          |
 *                           |          |
 *                          3|          |1
 *                           |          |
 *                           |          |
 *                           +----------+
 *                                0
 * \param [out] n[2]   Filled with the two possible face neighbors, which are
 *                     half of the size assuming the 2-1 constaint.
 * \param [out] nur[2] If not NULL, filled with smallest quadrants that fit
 *                     in the upper right corners of \a n.
 */
static void
p4est_quadrant_get_half_face_neighbors (const p4est_quadrant_t * q,
                                        int face,
                                        p4est_quadrant_t n[],
                                        p4est_quadrant_t nur[])
{
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);
  const p4est_qcoord_t qh_2 = P4EST_QUADRANT_LEN (q->level + 1);

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (q->level < P4EST_QMAXLEVEL);

  n[0].level = (int8_t) (q->level + 1);
  n[1].level = (int8_t) (q->level + 1);

  switch (face) {
  case 0:
    n[0].x = q->x;
    n[0].y = n[1].y = q->y - qh_2;
    n[1].x = n[0].x + qh_2;
    break;
  case 1:
    n[0].x = n[1].x = q->x + qh;
    n[0].y = q->y;
    n[1].y = n[0].y + qh_2;
    break;
  case 2:
    n[0].x = q->x;
    n[0].y = n[1].y = q->y + qh;
    n[1].x = n[0].x + qh_2;
    break;
  case 3:
    n[0].x = n[1].x = q->x - qh_2;
    n[0].y = q->y;
    n[1].y = n[0].y + qh_2;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }
  P4EST_ASSERT (p4est_quadrant_is_extended (&n[0]));
  P4EST_ASSERT (p4est_quadrant_is_extended (&n[1]));

  if (nur != NULL) {
    const p4est_qcoord_t dh = qh_2 - P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);

    nur[0].x = n[0].x + dh;
    nur[0].y = n[0].y + dh;
    nur[0].level = P4EST_QMAXLEVEL;
    P4EST_ASSERT (p4est_quadrant_is_extended (&nur[0]));
    nur[1].x = n[1].x + dh;
    nur[1].y = n[1].y + dh;
    nur[1].level = P4EST_QMAXLEVEL;
    P4EST_ASSERT (p4est_quadrant_is_extended (&nur[1]));
  }
}

/** Get the possible face neighbors of \a q.
 *
 * Gets the all face neighbors, possible assuming the 2-1 constraint.
 * If the larger quadrant doesn't exist than it is returned
 * as initialized by P4EST_QUADRANT_INIT.
 *
 * The order of \a n0 and \a n1 are given in the morton ordering.
 *
 * \param [in]  q      The quadrant whose face neighbors will be constructed.
 * \param [in]  face   The face across which to generate the neighbors.  The
 *                     face is given in right hand order. So
 *                                2
 *                           +----------+
 *                           |          |
 *                           |          |
 *                          3|          |1
 *                           |          |
 *                           |          |
 *                           +----------+
 *                                0
 * \param [out] n[0]   Filled with the first possible face neighbor, which is
 *                     half of the size if it exists or initialized to
 *                     P4EST_QUADRANT_INIT.
 * \param [out] n[1]   Filled with the second possible face neighbor, which is
 *                     half of the size if it exists or initialized to
 *                     P4EST_QUADRANT_INIT.
 * \param [out] n[2]   Filled with the face neighbor, which is the same size.
 * \param [out] n[3]   Filled with the face neighbor, which is twice the size
 *                     if it exists or initialized to P4EST_QUADRANT_INIT.
 */
static void
p4est_quadrant_get_possible_face_neighbors (const p4est_quadrant_t * q,
                                            int face, p4est_quadrant_t n[])
{
  const int           qcid = p4est_quadrant_child_id (q);
  const int           rqcid = p4est_corner_to_zorder[qcid];
  p4est_quadrant_t   *r = &n[3];

  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (q->level == P4EST_QMAXLEVEL) {
    P4EST_QUADRANT_INIT (&n[0]);
    P4EST_QUADRANT_INIT (&n[1]);
  }
  else {
    p4est_quadrant_get_half_face_neighbors (q, face, n, NULL);
  }

  p4est_quadrant_face_neighbor (q, face, &n[2]);

  /* Check to see if the larger element exists */
  if (((face != rqcid) && (face != ((rqcid + 3) % 4))) || (q->level == 0)) {
    P4EST_QUADRANT_INIT (r);
  }
  else {
    p4est_quadrant_parent (q, r);
    p4est_quadrant_face_neighbor (r, face, r);
  }
}

#endif /* !P4_TO_P8 */

#ifdef P4EST_MPI

/** This adds a quadrant to the end of a buffer.
 *
 * It crams the tree id into the user_data field of the quadrant in
 * the buffer and only adds the quadrant to the end of the buffer if
 * it is unique.
 *
 * \param [in,out] buf    \a q is added to the end if it is not alread there.
 * \param [in,out] q      the quadrant to be added.  The \c user_data field
 *                        is filled with \a treeid.
 * \param [in]     treeid the tree id of \a q.
 *
 */
static void
p4est_add_ghost_to_buf (sc_array_t * buf, p4est_topidx_t treeid,
                        int owner_rank, const p4est_quadrant_t * q)
{
  p4est_quadrant_t   *qold, *qnew;

  P4EST_ASSERT (treeid >= 0);
  P4EST_ASSERT (owner_rank >= 0);

  /* Check to see if the quadrant already is last in the array */
  if (buf->elem_count > 0) {
    qold = sc_array_index (buf, buf->elem_count - 1);
    if (treeid == qold->p.piggy1.which_tree &&
        p4est_quadrant_compare (q, qold) == 0) {
      return;
    }
  }

  qnew = sc_array_push (buf);
  *qnew = *q;

  /* Cram the tree id and the MPI rank into the user_data pointer */
  qnew->p.piggy1.which_tree = treeid;
  qnew->p.piggy1.owner_rank = owner_rank;
}

#endif /* P4EST_MPI */

bool
p4est_is_balanced (p4est_t * p4est)
{
  int                 zero = 0;
  int                 face, corner;
  int                 i, qcid;
  int                *pe0, *pe1, *pe2;
  bool                failed;
  bool                e0, e1, e0b, e1b, e2, e3;
  bool                bigger_face[2 * P4EST_DIM];
  size_t              cez;
  p4est_topidx_t      nt;
  p4est_topidx_t      first_local_tree = p4est->first_local_tree;
  p4est_topidx_t      last_local_tree = p4est->last_local_tree;
  p4est_locidx_t      li;
  p4est_locidx_t      num_quadrants;
  p4est_quadrant_t   *q;
  p4est_quadrant_t    n[P4EST_CHILDREN / 2 + 2];
  p4est_tree_t       *tree;
  sc_array_t          ghost_layer;
  sc_array_t         *quadrants;
  sc_array_t          e0_a, e1_a, e2_a;
#ifdef P4_TO_P8
  int                 edge;
  int                *pe3;
  bool                bigger_edge[12];
  size_t              big_count[12];
  sc_array_t          e3_a;
#endif

  sc_array_init (&ghost_layer, sizeof (p4est_quadrant_t));
  if (!p4est_build_ghost_layer (p4est, &ghost_layer)) {
    P4EST_NOTICE ("Ghost layer could not be built\n");
    return false;
  }

  for (i = 0; i < P4EST_CHILDREN / 2 + 2; ++i) {
    P4EST_QUADRANT_INIT (&n[i]);
  }

  failed = false;
  sc_array_init (&e0_a, sizeof (int));
  sc_array_init (&e1_a, sizeof (int));
  sc_array_init (&e2_a, sizeof (int));
#ifdef P4_TO_P8
  sc_array_init (&e3_a, sizeof (int));
#endif

  /* loop over all local trees */
  for (nt = first_local_tree; nt <= last_local_tree; ++nt) {
    tree = p4est_array_index_topidx (p4est->trees, nt);
    quadrants = &tree->quadrants;
    num_quadrants = (p4est_locidx_t) quadrants->elem_count;

    /* Find the neighboring processors of each quadrant */
    for (li = 0; li < num_quadrants; ++li) {
      q = sc_array_index (quadrants, (size_t) li);
      qcid = p4est_quadrant_child_id (q);

      /* Find face neighbors */
      for (face = 0; face < 2 * P4EST_DIM; ++face) {
        bigger_face[face] = false;

        /* If q is at a boundary then it is automatically balanced */
        if (p4est_quadrant_on_face_boundary (p4est, nt, face, q)) {
          continue;
        }

        /* Do more expensive face balance checks */
        p4est_quadrant_get_possible_face_neighbors (q, face, n);
        e0 = p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[0], NULL);
        e1 = p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[1], NULL);
#ifndef P4_TO_P8
        e0b = e1b = e0;
        i = 2;
#else
        e0b = p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[2], NULL);
        e1b = p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[3], NULL);
        i = 4;
#endif
        if (e0 != e1 || e0 != e0b || e0 != e1b) {
          P4EST_NOTICE ("Contradicting small face neighbors\n");
          failed = true;
          goto failtest;
        }
        e2 = p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[i], NULL);
        e3 = p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[i + 1], NULL);
        if ((int) e0 + (int) e2 + (int) e3 != 1) {
          P4EST_NOTICE ("Face balance quadrant mismatch\n");
          failed = true;
          goto failtest;
        }
        bigger_face[face] = e3;
      }

#ifdef P4_TO_P8
      /* Find edge neighbors */
      for (edge = 0; edge < 12; ++edge) {
        bigger_edge[edge] = false;
        big_count[edge] = 0;

        /* If q is at a boundary then it is automatically balanced */
        if (p8est_quadrant_on_edge_boundary (p4est, nt, edge, q)) {
          continue;
        }

        /* Do more expensive edge balance checks */
        p8est_quadrant_get_possible_edge_neighbors (q, edge, n);
        e0 = p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[0], &e0_a);
        e1 = p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[1], &e1_a);
        if (e0 != e1 || e0_a.elem_count != e1_a.elem_count) {
          P4EST_NOTICE ("Contradicting small edge neighbors\n");
          failed = true;
          goto failtest;
        }
        e2 = p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[2], &e2_a);
        e3 = p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[3], &e3_a);
        P4EST_ASSERT (((e0_a.elem_count == 0 && q->level == P4EST_QMAXLEVEL)
                       || e0_a.elem_count == e2_a.elem_count)
                      && ((e3_a.elem_count == 0 && (!e3 || q->level == 0))
                          || e3_a.elem_count == e2_a.elem_count));

        face = p8est_child_edge_face[qcid][edge];
        if (face >= 0 && bigger_face[face]) {
          P4EST_ASSERT (e2_a.elem_count == 0);
          if (e0 || e2 || e3) {
            P4EST_NOTICE ("Invalid edges across hanging face\n");
            failed = true;
            goto failtest;
          }
        }
        else {
          if (!e0 && !e2 && !e3) {
            P4EST_NOTICE ("Edge balance missing quadrants\n");
            failed = true;
            goto failtest;
          }
          if (e2_a.elem_count == 0 && (int) e0 + (int) e2 + (int) e3 != 1) {
            P4EST_NOTICE ("Edge balance duplicate quadrants\n");
            failed = true;
            goto failtest;
          }
          for (cez = 0; cez < e2_a.elem_count; ++cez) {
            pe0 = (e0_a.elem_count > 0) ? sc_array_index (&e0_a, cez) : &zero;
            pe1 = (e1_a.elem_count > 0) ? sc_array_index (&e1_a, cez) : &zero;
            pe2 = sc_array_index (&e2_a, cez);
            pe3 = (e3_a.elem_count > 0) ? sc_array_index (&e3_a, cez) : &zero;
            if (*pe0 + *pe2 + *pe3 != 1 || *pe0 != *pe1) {
              P4EST_NOTICE ("Edge balance quadrant mismatch\n");
              failed = true;
              goto failtest;
            }
          }
        }
        bigger_edge[edge] = e3;
        big_count[edge] = e3_a.elem_count;
      }
#endif

      /* Find corner neighbors, corner is in z-order here */
      for (corner = 0; corner < P4EST_CHILDREN; ++corner) {

        /* If q is at a boundary then it is automatically balanced */
        if (p4est_quadrant_on_corner_boundary (p4est, nt, corner, q)) {
          continue;
        }

        /* Do more expensive corner balance checks */
        p4est_quadrant_get_possible_corner_neighbors (q, corner, n);
        e0 = p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[0], &e0_a);
        e1 = p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[1], &e1_a);
        e2 = p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[2], &e2_a);
        P4EST_ASSERT (((e0_a.elem_count == 0 && q->level == P4EST_QMAXLEVEL)
                       || e0_a.elem_count == e1_a.elem_count)
                      && ((e2_a.elem_count == 0 && (!e2 || q->level == 0))
                          || e2_a.elem_count == e1_a.elem_count));

#ifndef P4_TO_P8
        if ((corner == p4est_hanging_corner[qcid][0] &&
             bigger_face[p4est_hanging_face[qcid][0]]) ||
            (corner == p4est_hanging_corner[qcid][1] &&
             bigger_face[p4est_hanging_face[qcid][1]])) {
          P4EST_ASSERT (e1_a.elem_count == 0);
          if (e0 || e1 || e2) {
            P4EST_NOTICE ("Invalid corners across hanging face\n");
            failed = true;
            goto failtest;
          }
        }
        else if (!e0 && !e1 && !e2) {
          P4EST_NOTICE ("Corner balance missing quadrants\n");
          failed = true;
          goto failtest;
        }
        else if (e1_a.elem_count == 0 && (int) e0 + (int) e1 + (int) e2 != 1) {
          P4EST_NOTICE ("Corner balance duplicate quadrants\n");
          failed = true;
          goto failtest;
        }
#else
        face = p8est_child_corner_face[qcid][corner];
        edge = p8est_child_corner_edge[qcid][corner];
        if (face >= 0 && bigger_face[face]) {
          P4EST_ASSERT (e1_a.elem_count == 0);
          if (e0 || e1 || e2) {
            P4EST_NOTICE ("Invalid corners across hanging face\n");
            failed = true;
            goto failtest;
          }
        }
        else if (e1_a.elem_count == 0) {
          if (edge >= 0 && bigger_edge[edge]) {
            P4EST_ASSERT (big_count[edge] == 0);
            if (e0 || e1 || e2) {
              P4EST_NOTICE ("Invalid corners across hanging edge\n");
              failed = true;
              goto failtest;
            }
          }
          else if ((int) e0 + (int) e1 + (int) e2 != 1) {
            P4EST_NOTICE ("Corner balance quadrant mismatch\n");
            failed = true;
            goto failtest;
          }
        }
#endif
        else {
#ifdef P4_TO_P8
          e3 = false;
          if (edge >= 0 && bigger_edge[edge]) {
            P4EST_ASSERT (big_count[edge] == e1_a.elem_count);

            /* recreate the edge neighbor information */
            p4est_quadrant_parent (q, &n[3]);
            p8est_quadrant_edge_neighbor (&n[3], edge, &n[3]);
            e3 =
              p4est_quadrant_exists (p4est, &ghost_layer, nt, &n[3], &e3_a);
            P4EST_ASSERT (e3 && big_count[edge] == e3_a.elem_count);
          }
#endif
          for (cez = 0; cez < e1_a.elem_count; ++cez) {
            pe0 = (e0_a.elem_count > 0) ? sc_array_index (&e0_a, cez) : &zero;
            pe1 = sc_array_index (&e1_a, cez);
            pe2 = (e2_a.elem_count > 0) ? sc_array_index (&e2_a, cez) : &zero;
#ifdef P4_TO_P8
            if (e3) {
              pe3 = sc_array_index (&e3_a, cez);
              if (*pe3) {
                if (*pe0 || *pe1 || *pe2) {
                  P4EST_NOTICE ("Invalid corners across hanging edge\n");
                  failed = true;
                  goto failtest;
                }
                continue;
              }
            }
#endif
            if (*pe0 + *pe1 + *pe2 != 1) {
              P4EST_NOTICE ("Corner balance quadrant mismatch\n");
              failed = true;
              goto failtest;
            }
          }
        }
      }
    }
  }

failtest:
  sc_array_reset (&ghost_layer);
  sc_array_reset (&e0_a);
  sc_array_reset (&e1_a);
  sc_array_reset (&e2_a);
#ifdef P4_TO_P8
  sc_array_reset (&e3_a);
#endif

  return !p4est_comm_sync_flag (p4est, failed, MPI_BOR);
}

bool
p4est_build_ghost_layer (p4est_t * p4est, sc_array_t * ghost_layer)
{
#ifdef P4EST_MPI
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  int                 face, corner, zcorner;
  int                 i, ncheck, ncount;
  int                 n0_proc, n0ur_proc, n1_proc;
  int                 num_peers, peer, peer_proc;
  int                 mpiret;
  bool                maxed, failed;
  bool                full_tree[2], tree_contact[2 * P4EST_DIM];
  bool                urg[P4EST_DIM - 1];
  size_t              pz;
  p4est_topidx_t      nt;
  p4est_topidx_t      first_local_tree = p4est->first_local_tree;
  p4est_topidx_t      last_local_tree = p4est->last_local_tree;
  p4est_locidx_t      li;
  p4est_locidx_t      num_quadrants;
  p4est_locidx_t      num_ghosts, ghost_offset, skipped;
  p4est_locidx_t     *send_counts, *recv_counts;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  p4est_quadrant_t    n[P4EST_CHILDREN / 2], nur[P4EST_CHILDREN / 2];
  sc_array_t          send_bufs;
  sc_array_t          procs[P4EST_DIM - 1];
  sc_array_t         *buf;
  sc_array_t         *quadrants;
  MPI_Comm            comm = p4est->mpicomm;
  MPI_Request        *recv_request, *send_request;
  MPI_Status         *recv_status, *send_status;
  MPI_Request        *recv_load_request, *send_load_request;
  MPI_Status         *recv_load_status, *send_load_status;
#ifdef P4_TO_P8
  int                 edge;
  int                 n1ur_proc;
#endif
#ifdef P4EST_DEBUG
  p4est_quadrant_t   *q2;
#endif

  P4EST_ASSERT (ghost_layer->elem_size == sizeof (p4est_quadrant_t));

  for (i = 0; i < P4EST_CHILDREN / 2; ++i) {
    P4EST_QUADRANT_INIT (&n[i]);
    P4EST_QUADRANT_INIT (&nur[i]);
  }

  failed = false;
  for (i = 0; i < P4EST_DIM - 1; ++i) {
    sc_array_init (&procs[i], sizeof (int));
  }
  skipped = 0;

  /* allocate empty send buffers */
  sc_array_init (&send_bufs, sizeof (sc_array_t));
  sc_array_resize (&send_bufs, (size_t) num_procs);
  for (i = 0; i < num_procs; ++i) {
    buf = sc_array_index_int (&send_bufs, i);
    sc_array_init (buf, sizeof (p4est_quadrant_t));
  }

  /* loop over all local trees */
  for (nt = first_local_tree; nt <= last_local_tree; ++nt) {
    tree = p4est_array_index_topidx (p4est->trees, nt);
    quadrants = &tree->quadrants;
    num_quadrants = (p4est_locidx_t) quadrants->elem_count;
    p4est_comm_tree_info (p4est, nt, full_tree, tree_contact, NULL, NULL);

    /* Find the smaller neighboring processors of each quadrant */
    for (li = 0; li < num_quadrants; ++li) {
      q = sc_array_index (quadrants, li);

      if (p4est_comm_neighborhood_owned
          (p4est, nt, full_tree, tree_contact, q)) {
        /* The 3x3 neighborhood of q is owned by this processor */
        ++skipped;
        continue;
      }

      /* Find smaller face neighbors */
      for (face = 0; face < 2 * P4EST_DIM; ++face) {
        if (q->level == P4EST_QMAXLEVEL) {
          p4est_quadrant_face_neighbor (q, face, &n[0]);
          ncheck = 0;
          ncount = 1;
        }
        else {
          p4est_quadrant_get_half_face_neighbors (q, face, n, nur);
          ncheck = ncount = P4EST_CHILDREN / 2;
        }

        n1_proc = -1;
        for (i = 0; i < ncount; ++i) {
          n0_proc = p4est_quadrant_find_owner (p4est, nt, &n[i]);
          if (i < ncheck) {
            /* Note that we will always check this
             * because it prevents deadlocks
             */
            n0ur_proc = p4est_quadrant_find_owner (p4est, nt, &nur[i]);
            if (n0_proc != n0ur_proc) {
              P4EST_NOTICE ("Small face owner inconsistency\n");
              failed = true;
              goto failtest;
            }
          }

          if (n0_proc != rank && n0_proc >= 0 && n0_proc != n1_proc) {
            buf = sc_array_index_int (&send_bufs, n0_proc);
            p4est_add_ghost_to_buf (buf, nt, rank, q);
            n1_proc = n0_proc;
          }
        }
      }

#ifdef P4_TO_P8
      /* Find smaller edge neighbors */
      for (edge = 0; edge < 12; ++edge) {
        if (q->level == P4EST_QMAXLEVEL) {
          p8est_quadrant_edge_neighbor (q, edge, &n[0]);
          maxed = true;
        }
        else {
          p8est_quadrant_get_half_edge_neighbors (q, edge, n, nur);
          maxed = false;
        }

        /* Check to see if we are a tree edge neighbor */
        P4EST_ASSERT (!p4est_quadrant_is_outside_corner (&n[0]));
        if (p8est_quadrant_is_outside_edge (&n[0])) {
          p8est_quadrant_find_tree_edge_owners (p4est, nt, edge,
                                                &n[0], &procs[0], &urg[0]);
          if (!maxed) {
            p8est_quadrant_find_tree_edge_owners (p4est, nt, edge,
                                                  &n[1], &procs[1], &urg[1]);
            P4EST_ASSERT (procs[0].elem_count == procs[1].elem_count);

            if (!urg[0] || !urg[1]) {
              P4EST_NOTICE ("Tree edge owner inconsistency\n");
              failed = true;
              goto failtest;
            }
          }

          /* Then we have to loop over multiple neighbors */
          for (pz = 0; pz < procs[0].elem_count; ++pz) {
            n0_proc = *((int *) sc_array_index (&procs[0], pz));

            if (n0_proc != rank) {
              buf = sc_array_index_int (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, rank, q);
            }

            if (!maxed) {
              n1_proc = *((int *) sc_array_index (&procs[1], pz));

              if (n1_proc != n0_proc && n1_proc != rank) {
                buf = sc_array_index_int (&send_bufs, n1_proc);
                p4est_add_ghost_to_buf (buf, nt, rank, q);
              }
            }
          }
        }
        else {
          /* We are not at a tree edge so we only have two neighbors
           * either inside the tree or across a face
           */
          n0_proc = n1_proc = p4est_quadrant_find_owner (p4est, nt, &n[0]);
          if (!maxed) {
            n1_proc = p4est_quadrant_find_owner (p4est, nt, &n[1]);
            n0ur_proc = p4est_quadrant_find_owner (p4est, nt, &nur[0]);
            n1ur_proc = p4est_quadrant_find_owner (p4est, nt, &nur[1]);

            /* Note that we will always check this
             * because it prevents deadlocks
             */
            if (n0_proc != n0ur_proc || n1_proc != n1ur_proc) {
              P4EST_NOTICE ("Small edge owner inconsistency\n");
              failed = true;
              goto failtest;
            }
          }

          if (n0_proc != rank && n0_proc >= 0) {
            buf = sc_array_index_int (&send_bufs, n0_proc);
            p4est_add_ghost_to_buf (buf, nt, rank, q);
          }

          if (n1_proc != n0_proc && n1_proc != rank && n1_proc >= 0) {
            buf = sc_array_index_int (&send_bufs, n1_proc);
            p4est_add_ghost_to_buf (buf, nt, rank, q);
          }
        }
      }
#endif

      /* Find smaller corner neighbors */
      for (corner = 0; corner < P4EST_CHILDREN; ++corner) {
#ifndef P4_TO_P8
        zcorner = p4est_corner_to_zorder[corner];
#else
        zcorner = corner;
#endif
        if (q->level == P4EST_QMAXLEVEL) {
          p4est_quadrant_corner_neighbor (q, zcorner, &n[0]);
          maxed = true;
        }
        else {
          p4est_quadrant_get_half_corner_neighbor (q, zcorner, &n[0],
                                                   &nur[0]);
          maxed = false;
        }

        /* Check to see if we are a tree corner neighbor */
        if (p4est_quadrant_is_outside_corner (&n[0])) {
          /* Then we have to loop over multiple corner neighbors */
          p4est_quadrant_find_tree_corner_owners (p4est, nt, corner,
                                                  &n[0], &procs[0], &urg[0]);
          if (!urg[0]) {
            P4EST_NOTICE ("Tree corner owner inconsistency\n");
            failed = true;
            goto failtest;
          }

          for (pz = 0; pz < procs[0].elem_count; ++pz) {
            n0_proc = *((int *) sc_array_index (&procs[0], pz));

            if (n0_proc != rank) {
              buf = sc_array_index_int (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, rank, q);
            }
          }
        }
#ifdef P4_TO_P8
        /* Check to see if we are a tree edge neighbor */
        else if (p8est_quadrant_is_outside_edge_extra (&n[0], &edge)) {
          p8est_quadrant_find_tree_edge_owners (p4est, nt, edge,
                                                &n[0], &procs[0], &urg[0]);
          if (!urg[0]) {
            P4EST_NOTICE ("Tree corner/edge owner inconsistency\n");
            failed = true;
            goto failtest;
          }

          /* Then we have to loop over multiple edge neighbors */
          for (pz = 0; pz < procs[0].elem_count; ++pz) {
            n0_proc = *((int *) sc_array_index (&procs[0], pz));

            if (n0_proc != rank) {
              buf = sc_array_index_int (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, rank, q);
            }
          }
        }
#endif
        else {
          /* We are not at a tree edge or corner so
           * we only have one corner neighbor
           */
          n0_proc = p4est_quadrant_find_owner (p4est, nt, &n[0]);
          if (!maxed) {
            n0ur_proc = p4est_quadrant_find_owner (p4est, nt, &nur[0]);

            /* Note that we will always check this
             * because it prevents deadlocks
             */
            if (n0_proc != n0ur_proc) {
              P4EST_NOTICE ("Small corner owner inconsistency\n");
              failed = true;
              goto failtest;
            }
          }

          if (n0_proc != rank && n0_proc >= 0) {
            buf = sc_array_index_int (&send_bufs, n0_proc);
            p4est_add_ghost_to_buf (buf, nt, rank, q);
          }
        }
      }
    }
  }

failtest:
  if (p4est_comm_sync_flag (p4est, failed, MPI_BOR)) {
    for (i = 0; i < num_procs; ++i) {
      buf = sc_array_index_int (&send_bufs, i);
      sc_array_reset (buf);
    }
    sc_array_reset (&send_bufs);
    for (i = 0; i < P4EST_DIM - 1; ++i) {
      sc_array_reset (&procs[i]);
    }
    sc_array_reset (ghost_layer);

    return false;
  }

  /* Count the number of peers that I send to and receive from */
  for (i = 0, num_peers = 0; i < num_procs; ++i) {
    buf = sc_array_index_int (&send_bufs, i);
    if (buf->elem_count > 0)
      ++num_peers;
  }

  recv_request = P4EST_ALLOC (MPI_Request, 2 * num_peers);
  recv_status = P4EST_ALLOC (MPI_Status, 2 * num_peers);

  send_request = P4EST_ALLOC (MPI_Request, 2 * num_peers);
  send_status = P4EST_ALLOC (MPI_Status, 2 * num_peers);

  recv_counts = P4EST_ALLOC (p4est_locidx_t, 2 * num_peers);
  send_counts = recv_counts + num_peers;

  recv_load_request = recv_request + num_peers;
  recv_load_status = recv_status + num_peers;

  send_load_request = send_request + num_peers;
  send_load_status = send_status + num_peers;

  /* Post receives for the counts of ghosts to be received */
  for (i = 0, peer = 0; i < num_procs; ++i) {
    buf = sc_array_index_int (&send_bufs, i);
    if (buf->elem_count > 0) {
      peer_proc = i;
      P4EST_ASSERT (peer_proc != rank);
      P4EST_LDEBUGF ("ghost layer post count receive from %d\n", peer_proc);
      mpiret = MPI_Irecv (recv_counts + peer, 1, P4EST_MPI_LOCIDX,
                          peer_proc, P4EST_COMM_GHOST_COUNT, comm,
                          recv_request + peer);
      SC_CHECK_MPI (mpiret);
      ++peer;
    }
  }

  /* Send the counts of ghosts that are going to be sent */
  for (i = 0, peer = 0; i < num_procs; ++i) {
    buf = sc_array_index_int (&send_bufs, i);
    if (buf->elem_count > 0) {
      peer_proc = i;
      send_counts[peer] = (p4est_locidx_t) buf->elem_count;
      P4EST_LDEBUGF ("ghost layer post count send %lld to %d\n",
                     (long long) send_counts[peer], peer_proc);
      mpiret = MPI_Isend (send_counts + peer, 1, P4EST_MPI_LOCIDX,
                          peer_proc, P4EST_COMM_GHOST_COUNT,
                          comm, send_request + peer);
      SC_CHECK_MPI (mpiret);
      ++peer;
    }
  }

  /* Wait for the counts */
  if (num_peers > 0) {
    mpiret = MPI_Waitall (num_peers, recv_request, recv_status);
    SC_CHECK_MPI (mpiret);

    mpiret = MPI_Waitall (num_peers, send_request, send_status);
    SC_CHECK_MPI (mpiret);
  }

#ifdef P4EST_DEBUG
  for (i = 0; i < num_peers; ++i) {
    P4EST_ASSERT (recv_request[i] == MPI_REQUEST_NULL);
  }
  for (i = 0; i < num_peers; ++i) {
    P4EST_ASSERT (send_request[i] == MPI_REQUEST_NULL);
  }
#endif

  /* Count ghosts */
  for (peer = 0, num_ghosts = 0; peer < num_peers; ++peer) {
    P4EST_ASSERT (recv_counts[peer] > 0);
    num_ghosts += recv_counts[peer];    /* same type */
  }
  P4EST_VERBOSEF ("Total quadrants skipped %lld ghosts to receive %lld\n",
                  (long long) skipped, (long long) num_ghosts);

  /* Allocate space for the ghosts */
  sc_array_resize (ghost_layer, (size_t) num_ghosts);

  /* Post receives for the ghosts */
  for (i = 0, peer = 0, ghost_offset = 0; i < num_procs; ++i) {
    buf = sc_array_index_int (&send_bufs, i);
    if (buf->elem_count > 0) {
      peer_proc = i;
      P4EST_LDEBUGF
        ("ghost layer post ghost receive %lld quadrants from %d\n",
         (long long) recv_counts[peer], peer_proc);
      mpiret =
        MPI_Irecv (ghost_layer->array +
                   ghost_offset * sizeof (p4est_quadrant_t),
                   (int) (recv_counts[peer] * sizeof (p4est_quadrant_t)),
                   MPI_BYTE, peer_proc, P4EST_COMM_GHOST_LOAD, comm,
                   recv_load_request + peer);
      SC_CHECK_MPI (mpiret);
      ghost_offset += recv_counts[peer];        /* same type */
      ++peer;
    }
  }
  P4EST_ASSERT (ghost_offset == num_ghosts);

  /* Send the ghosts */
  for (i = 0, peer = 0; i < num_procs; ++i) {
    buf = sc_array_index_int (&send_bufs, i);
    if (buf->elem_count > 0) {
      peer_proc = i;
      P4EST_ASSERT ((p4est_locidx_t) buf->elem_count == send_counts[peer]);
      P4EST_LDEBUGF ("ghost layer post ghost send %lld quadrants to %d\n",
                     (long long) send_counts[peer], peer_proc);
      mpiret =
        MPI_Isend (buf->array,
                   (int) (send_counts[peer] * sizeof (p4est_quadrant_t)),
                   MPI_BYTE, peer_proc, P4EST_COMM_GHOST_LOAD, comm,
                   send_load_request + peer);
      SC_CHECK_MPI (mpiret);
      ++peer;
    }
  }

  /* Wait for everything */
  if (num_peers > 0) {
    mpiret = MPI_Waitall (num_peers, recv_load_request, recv_load_status);
    SC_CHECK_MPI (mpiret);

    mpiret = MPI_Waitall (num_peers, send_load_request, send_load_status);
    SC_CHECK_MPI (mpiret);
  }

  /* Clean up */
  P4EST_FREE (recv_counts);

#ifdef P4EST_DEBUG
  for (i = 0; i < num_peers; ++i) {
    P4EST_ASSERT (recv_load_request[i] == MPI_REQUEST_NULL);
  }
  for (i = 0; i < num_peers; ++i) {
    P4EST_ASSERT (send_load_request[i] == MPI_REQUEST_NULL);
  }
  q2 = NULL;
  for (li = 0; li < num_ghosts; ++li) {
    q = sc_array_index (ghost_layer, (size_t) li);
    P4EST_ASSERT (p4est_quadrant_is_valid (q));
    P4EST_ASSERT (q->p.piggy1.which_tree >= 0 &&
                  q->p.piggy1.which_tree < p4est->connectivity->num_trees);
    P4EST_ASSERT (q->p.piggy1.owner_rank >= 0 &&
                  q->p.piggy1.owner_rank < num_procs);
    P4EST_ASSERT (q->p.piggy1.owner_rank != rank);
    if (q2 != NULL) {
      P4EST_ASSERT (p4est_quadrant_compare_piggy (q2, q) < 0);
      P4EST_ASSERT (q2->p.piggy1.owner_rank <= q->p.piggy1.owner_rank);
    }
    q2 = q;
  }
#endif

  P4EST_FREE (recv_request);
  P4EST_FREE (recv_status);
  P4EST_FREE (send_request);
  P4EST_FREE (send_status);

  for (i = 0; i < num_procs; ++i) {
    buf = sc_array_index_int (&send_bufs, i);
    sc_array_reset (buf);
  }
  sc_array_reset (&send_bufs);
  for (i = 0; i < P4EST_DIM - 1; ++i) {
    sc_array_reset (&procs[i]);
  }

#else
  /* If we are not running with mpi then we don't need to do anything */
  sc_array_reset (ghost_layer);
#endif

  return true;
}

#ifndef P4_TO_P8

/** Generate a neighbor of a quadrant for a given node.
 *
 * The neighbor numbering is given below.
 *
 * Neighbor numbering for q, node=0:
 *
 *      ------+------+
 *            |  q   |
 *      nnum=2|nnum=3|
 *            |      |
 *      ------+------+
 *            |      |
 *      nnum=0|nnum=1|
 *
 * Neighbor numbering for q, node=1:
 *
 *            +------+------
 *            |  q   |
 *            |nnum=2|num=3
 *            |      |
 *            +------+------
 *            |      |
 *            |nnum=0|nnum=1
 *
 * Neighbor numbering for q, node=2:
 *
 *            |      |
 *      nnum=2|nnum=3|
 *            |      |
 *      ------+------+
 *            |  q   |
 *      nnum=0|nnum=1|
 *            |      |
 *      ------+------+
 *
 * Neighbor numbering for q, node=3:
 *
 *            |      |
 *            |nnum=2|nnum=3
 *            |      |
 *            +------+------
 *            |  q   |
 *            |nnum=0|nnum=1
 *            |      |
 *            +------+------
 *
 * \param [in]  q             the quadrant whose possible node \a node neighbor
 *                            will be built.
 * \param [in]  node          the node of the quadrant \a q whose possible node
 *                            neighbor list will be built.  This is given in
 *                            pixel (Morton-) ordering.
 * \param [in]  nnum          neighbor number in the ordering described above,
 *                            if nnum==node then it is the corner neighbor.
 * \param [in]  neighor_rlev  the relative level of the neighbor compared to
 *                            the level of \a q.
 * \param [out] neighbor      the neighbor that will be filled.
 * \param [out] neighbor_node the neighbor's node which shares with \a q
 *                            the node \a node.
 */
static void
p4est_possible_node_neighbor (const p4est_quadrant_t * q, int node,
                              int nnum, int neighbor_rlev,
                              p4est_quadrant_t * neighbor, int *neighbor_node)
{
  int                 nnode;
  const int           nlevel = (int) q->level + neighbor_rlev;
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);
  const p4est_qcoord_t nh = P4EST_QUADRANT_LEN (nlevel);
  const p4est_qcoord_t qx = q->x;
  const p4est_qcoord_t qy = q->y;
  p4est_qcoord_t      cornerx, cornery;
  p4est_quadrant_t    n;
#ifdef P4EST_DEBUG
  int                 qcid;
#endif

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (-1 <= neighbor_rlev && neighbor_rlev <= 1);
  P4EST_ASSERT (0 <= nlevel && nlevel <= P4EST_QMAXLEVEL);
  P4EST_ASSERT (node + nnum != 3);

  P4EST_QUADRANT_INIT (&n);

  switch (node) {
  case 0:
    cornerx = qx;
    cornery = qy;
    break;
  case 1:
    cornerx = qx + qh;
    cornery = qy;
    break;
  case 2:
    cornerx = qx;
    cornery = qy + qh;
    break;
  case 3:
    cornerx = qx + qh;
    cornery = qy + qh;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }

#ifdef P4EST_DEBUG
  /* Check to see if it is possible to construct the neighbor */
  qcid = p4est_quadrant_child_id (q);
  P4EST_ASSERT (neighbor_rlev >= 0 || qcid == node);
#endif

  nnode = 3 - nnum;
  n.level = (int8_t) nlevel;
  switch (nnum) {
  case 0:
    n.x = cornerx - nh;
    n.y = cornery - nh;
    break;
  case 1:
    n.x = cornerx;
    n.y = cornery - nh;
    break;
  case 2:
    n.x = cornerx - nh;
    n.y = cornery;
    break;
  case 3:
    n.x = cornerx;
    n.y = cornery;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }

  *neighbor = n;
  *neighbor_node = nnode;

  P4EST_ASSERT (p4est_quadrant_is_extended (neighbor));
}

void
p4est_order_local_vertices (p4est_t * p4est,
                            bool identify_periodic,
                            p4est_locidx_t * num_uniq_local_vertices,
                            p4est_locidx_t * quadrant_to_local_vertex)
{
  const int           rank = p4est->mpirank;
  p4est_connectivity_t *conn = p4est->connectivity;
  int                 qcid, transform;
  int                 neighbor_node;
  int                 face, corner, nnum, rlev, tree_corner;
  int                 neighbor_proc;
  bool                face_contact[4];
  bool                quad_contact[4];
  p4est_topidx_t      first_local_tree = p4est->first_local_tree;
  p4est_topidx_t      last_local_tree = p4est->last_local_tree;
  p4est_topidx_t      jt, num_trees = conn->num_trees;
  p4est_locidx_t      Ntotal = 0;
  p4est_locidx_t      il, Ncells = p4est->local_num_quadrants;
  p4est_locidx_t      vertex_num;
  p4est_locidx_t      lqid;
  p4est_locidx_t      neighbor_tree;
  size_t              ctree;
  size_t              zz, numz_quadrants;
  ssize_t             lnid;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  p4est_locidx_t     *tree_offset;
  p4est_tree_t       *tree, *ntree;
  sc_array_t         *trees = p4est->trees;
  sc_array_t         *quadrants, ctransforms;
  p4est_quadrant_t    neighbor, cneighbor;
  p4est_quadrant_t   *q;
  p4est_corner_transform_t *ct;

  P4EST_ASSERT (p4est_is_valid (p4est));

  P4EST_QUADRANT_INIT (&neighbor);
  P4EST_QUADRANT_INIT (&cneighbor);

  sc_array_init (&ctransforms, sizeof (p4est_corner_transform_t));

  /* figure out the offset of each tree into the local element id */
  tree_offset = P4EST_ALLOC_ZERO (p4est_locidx_t, num_trees);
  if (first_local_tree >= 0) {
    tree_offset[first_local_tree] = 0;
    for (jt = first_local_tree; jt < last_local_tree; ++jt) {
      tree = sc_array_index (trees, jt);
      tree_offset[jt + 1] = tree_offset[jt] + tree->quadrants.elem_count;
    }
  }
  else {
    P4EST_ASSERT (first_local_tree == -1 && last_local_tree == -2);
  }

  /* Initialize vertex list to all -1.  This way we know which values
   * get set because legitimate values are >= 0.
   */
  for (il = 0; il < 4 * Ncells; ++il) {
    quadrant_to_local_vertex[il] = -1;
  }

  /* loop over all local trees to generate the connetivity list */
  for (jt = first_local_tree, vertex_num = 0, lqid = 0;
       jt <= last_local_tree; ++jt) {
    for (face = 0; face < 4; ++face) {
      face_contact[face] = (conn->tree_to_tree[4 * jt + face] != jt ||
                            (identify_periodic &&
                             (int) conn->tree_to_face[4 * jt + face] !=
                             face));
    }
    tree = sc_array_index (p4est->trees, jt);
    quadrants = &tree->quadrants;
    numz_quadrants = quadrants->elem_count;

    /* Find the neighbors of each quadrant */
    for (zz = 0; zz < numz_quadrants; ++zz, ++lqid) {
      /* this quadrant may be on the boundary with a range of processors */
      q = sc_array_index (quadrants, zz);

      /* loop over the corners of the quadrant */
      for (corner = 0; corner < 4; ++corner) {

        /* Check to see if we have a new vertex */
        if (quadrant_to_local_vertex[lqid * 4 + corner] == -1) {
          quadrant_to_local_vertex[lqid * 4 + corner] = vertex_num;

          /* loop over the possible neighbors and set the new vertex */
          for (nnum = 0; nnum < 4; ++nnum) {
            /* Don't search for the quadrant q */
            if (3 - nnum == corner)
              continue;

            qcid = p4est_quadrant_child_id (q);

            /* loop over possible neighbor sizes */
            for (rlev = -1; rlev < 2; ++rlev) {
              /* can't check for quadrants larger than the root */
              if (q->level == 0 && rlev < 0)
                continue;
              /* can't check for quadrants larger unless child id
               * and corner line up
               */
              if (qcid != corner && rlev < 0)
                continue;

              /* get possible neighbor */
              p4est_possible_node_neighbor (q, corner, nnum, rlev,
                                            &neighbor, &neighbor_node);

              if (p4est_quadrant_is_inside_root (&neighbor)) {
                /* neighbor is in the same tree */

                neighbor_proc = p4est_comm_find_owner (p4est, jt, &neighbor,
                                                       rank);

                /* Neighbor is remote so we don't number its node */
                if (neighbor_proc != rank)
                  continue;

                lnid = sc_array_bsearch (quadrants, &neighbor,
                                         p4est_quadrant_compare);
                if (lnid != -1) {
                  lnid += tree_offset[jt];
                  /* We have found a neighbor in the same tree */
                  quadrant_to_local_vertex[lnid * 4 + neighbor_node]
                    = vertex_num;

                  /* No need to check for more quadrants for this neighbor */
                  continue;
                }
              }
              else {
                /* the neighbor is in a neighboring tree or multiple
                 * if it is a neighbor across the corner of the tree
                 */

                quad_contact[0] = (neighbor.y < 0);
                quad_contact[1] = (neighbor.x >= rh);
                quad_contact[2] = (neighbor.y >= rh);
                quad_contact[3] = (neighbor.x < 0);

                if ((quad_contact[0] || quad_contact[2]) &&
                    (quad_contact[1] || quad_contact[3])) {
                  /* Neighbor is across a corner */
                  for (tree_corner = 0; tree_corner < 4; ++tree_corner) {
                    if (quad_contact[(tree_corner + 3) % 4]
                        && quad_contact[tree_corner]) {
                      break;
                    }
                  }
                  p4est_find_corner_transform (conn, jt, tree_corner,
                                               &ctransforms);
                  for (ctree = 0; ctree < ctransforms.elem_count; ++ctree) {
                    ct = sc_array_index (&ctransforms, ctree);
                    neighbor_tree = ct->ntree;

                    /* Don't use corner identification in the same tree */
                    if (!identify_periodic && neighbor_tree == jt)
                      continue;

                    cneighbor = neighbor;
                    p4est_quadrant_transform_corner (&cneighbor,
                                                     (int) ct->ncorner, true);

                    neighbor_proc = p4est_comm_find_owner (p4est,
                                                           neighbor_tree,
                                                           &cneighbor, rank);

                    /* Neighbor is remote so we don't number its node */
                    if (neighbor_proc != rank)
                      continue;

                    ntree = sc_array_index (trees, neighbor_tree);

                    lnid = sc_array_bsearch (&ntree->quadrants, &cneighbor,
                                             p4est_quadrant_compare);
                    if (lnid != -1) {
                      lnid += tree_offset[neighbor_tree];
                      neighbor_node = (int) ct->ncorner;
                      /* We have found a corner neighbor */
                      quadrant_to_local_vertex[lnid * 4 + neighbor_node]
                        = vertex_num;
                    }
                  }
                }
                else {
                  /* Neighbor is across a face */
                  for (face = 0; face < 4; ++face) {
                    if (quad_contact[face] && face_contact[face]) {
                      neighbor_tree = conn->tree_to_tree[4 * jt + face];
                      break;
                    }
                  }
                  if (face == 4) {
                    /* this quadrant ran across a face with no neighbor */
                    continue;
                  }
                  /* transform the neighbor into the other tree's
                   * coordinates
                   */
                  transform = p4est_find_face_transform (conn, jt, face);
                  p4est_quadrant_translate_face (&neighbor, face);
                  p4est_quadrant_transform_face (&neighbor, &cneighbor,
                                                 transform);

                  neighbor_proc = p4est_comm_find_owner (p4est,
                                                         neighbor_tree,
                                                         &cneighbor, rank);
                  /* Neighbor is remote so we don't number its node */
                  if (neighbor_proc != rank)
                    continue;

                  ntree = sc_array_index (trees, neighbor_tree);

                  lnid = sc_array_bsearch (&ntree->quadrants, &cneighbor,
                                           p4est_quadrant_compare);
                  if (lnid != -1) {
                    lnid += tree_offset[neighbor_tree];
                    neighbor_node = p4est_node_transform (neighbor_node,
                                                          transform);

                    /* We have found a face neighbor */
                    quadrant_to_local_vertex[lnid * 4 + neighbor_node]
                      = vertex_num;
                  }
                }
              }
            }
          }
          ++vertex_num;
        }
      }
    }
  }

  Ntotal = vertex_num;
  P4EST_FREE (tree_offset);
  sc_array_reset (&ctransforms);

  *num_uniq_local_vertices = Ntotal;
}

#endif /* !P4_TO_P8 */

/** Determine the right tree for a node and clamp it inside the domain.
 *
 * If the node is on the boundary, assign the lowest tree to own it.
 * Clamp it just inside the tree bounds if necessary.
 *
 * \param [in] p4est    The p4est to work on.
 * \param [in] treeid   Original tree index for this node.
 * \param [in] n        The node to work on.
 * \param [out] c       The clamped node in owning tree coordinates.
 *                      Its user data will be filled with owning tree id.
 * \return              Returns the owning tree id.
 */
static              p4est_topidx_t
p4est_node_canonicalize (p4est_t * p4est, p4est_topidx_t treeid,
                         const p4est_quadrant_t * n, p4est_quadrant_t * c)
{
  bool                quad_contact[2 * P4EST_DIM];
  int                 face;
  p4est_topidx_t      ntreeid, lowest;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_quadrant_t    tmpq, o;
#ifndef P4_TO_P8
  int                 transform;
#else
  int                 ftransform[9];
  p4est_topidx_t      ntreeid2;
#endif

  P4EST_ASSERT (treeid >= 0 && treeid < conn->num_trees);
  P4EST_ASSERT (p4est_quadrant_is_node (n, false));

  P4EST_QUADRANT_INIT (&tmpq);
  P4EST_QUADRANT_INIT (&o);

  p4est_node_clamp_inside (n, c);
  c->p.which_tree = lowest = treeid;

#ifndef P4_TO_P8
  quad_contact[0] = (n->y == 0);
  quad_contact[1] = (n->x == P4EST_ROOT_LEN);
  quad_contact[2] = (n->y == P4EST_ROOT_LEN);
  quad_contact[3] = (n->x == 0);
#else
  quad_contact[0] = (n->x == 0);
  quad_contact[1] = (n->x == P4EST_ROOT_LEN);
  quad_contact[2] = (n->y == 0);
  quad_contact[3] = (n->y == P4EST_ROOT_LEN);
  quad_contact[4] = (n->z == 0);
  quad_contact[5] = (n->z == P4EST_ROOT_LEN);
#endif

  for (face = 0; face < 2 * P4EST_DIM; ++face) {
    if (!quad_contact[face]) {
      /* The node is not touching this face */
      continue;
    }
    ntreeid = conn->tree_to_tree[2 * P4EST_DIM * treeid + face];
    if (ntreeid == treeid
        && ((int) conn->tree_to_face[2 * P4EST_DIM * treeid + face] ==
            face)) {
      /* The node touches a face with no neighbor */
      continue;
    }
    if (ntreeid > lowest) {
      /* This neighbor tree is higher, so we keep the ownership */
      continue;
    }
#ifndef P4_TO_P8
    /* Transform the node into the other tree's coordinates */
    transform = p4est_find_face_transform (conn, treeid, face);
    tmpq = *n;
    p4est_quadrant_translate_face (&tmpq, face);
    p4est_quadrant_transform_face (&tmpq, &o, transform);
#else
    ntreeid2 = p8est_find_face_transform (conn, treeid, face, ftransform);
    P4EST_ASSERT (ntreeid2 == ntreeid);
    p8est_quadrant_transform_face (n, &o, ftransform);
#endif
    if (ntreeid < lowest) {
      /* we have found a new owning tree */
      p4est_node_clamp_inside (&o, c);
      c->p.which_tree = lowest = ntreeid;
    }
    else {
      P4EST_ASSERT (lowest == ntreeid);
      p4est_node_clamp_inside (&o, &tmpq);
      if (p4est_quadrant_compare (&tmpq, c) < 0) {
        /* same tree (periodic) and the new position is lower than the old */
        *c = tmpq;
        c->p.which_tree = lowest;
      }
    }
  }

  return lowest;
}

void
p4est_collect_nodes (p4est_t * p4est, sc_array_t * ghost_layer)
{
  int                 k;
  int                 qcid;
  int                 face;
  size_t              zz;
  p4est_qcoord_t      hanging[P4EST_CHILDREN];
  p4est_topidx_t      jt;
  p4est_locidx_t      num_local_nodes;
  p4est_locidx_t     *local_nodes, *quad_nodes;
  p4est_tree_t       *tree;
  p4est_quadrant_t    p, n, c;
  p4est_quadrant_t   *q, *qp[2];
  sc_array_t         *quadrants;
  sc_array_t          exist_array;
#ifdef P4_TO_P8
  int                 l, edge, corner;
#ifdef P4EST_DEBUG
  int                 check;
#endif
#endif

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_ASSERT (ghost_layer != NULL);

  sc_array_init (&exist_array, sizeof (int));
  P4EST_QUADRANT_INIT (&p);
  P4EST_QUADRANT_INIT (&n);
  P4EST_QUADRANT_INIT (&c);
  qp[0] = NULL;
  qp[1] = &p;
  num_local_nodes = P4EST_CHILDREN * p4est->local_num_quadrants;
  local_nodes = P4EST_ALLOC (p4est_locidx_t, num_local_nodes);

  quad_nodes = local_nodes;
  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    tree = p4est_array_index_topidx (p4est->trees, jt);
    quadrants = &tree->quadrants;

    for (zz = 0; zz < quadrants->elem_count;
         quad_nodes += P4EST_CHILDREN, ++zz) {
      qp[0] = q = sc_array_index (quadrants, zz);
      qcid = p4est_quadrant_child_id (q);
      if (q->level > 0) {
        p4est_quadrant_parent (q, &p);
      }
#ifdef P4EST_DEBUG
      else {
        P4EST_QUADRANT_INIT (&p);
      }
#endif

      /* establish hanging node status of all corners */
      hanging[0] = hanging[1] = hanging[2] = hanging[3] = -1;
#ifdef P4_TO_P8
      hanging[4] = hanging[5] = hanging[6] = hanging[7] = -1;
#endif
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        if (k == qcid || k == P4EST_CHILDREN - 1 - qcid || q->level == 0) {
          hanging[k] = 0;
          continue;
        }
#ifndef P4_TO_P8
        if (k == p4est_hanging_corner[qcid][0]) {
          face = p4est_hanging_face[qcid][0];
        }
        else {
          P4EST_ASSERT (k == p4est_hanging_corner[qcid][1]);
          face = p4est_hanging_face[qcid][1];
        }
#else
        face = p8est_child_corner_face[qcid][k];
        if (face == -1) {
          P4EST_ASSERT (p8est_child_corner_edge[qcid][k] >= 0);
          continue;
        }
#endif
        p4est_quadrant_face_neighbor (&p, face, &n);
        if (p4est_quadrant_exists (p4est, ghost_layer, jt, &n, NULL)) {
          hanging[k] = 1;
#ifdef P4_TO_P8
#ifdef P4EST_DEBUG
          check = 0;
#endif
          for (l = 0; l < 4; ++l) {
            corner = p8est_face_corners[face][l];
            if (corner != qcid && corner != k) {
              hanging[corner] = 1;
#ifdef P4EST_DEBUG
              ++check;
#endif
            }
          }
          P4EST_ASSERT (check == 2);
#endif
        }
        else {
          hanging[k] = 0;
        }
      }
#ifdef P4_TO_P8
#ifdef P4EST_DEBUG
      check = 0;
#endif
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        if (hanging[k] == -1) {
          edge = p8est_child_corner_edge[qcid][k];
          P4EST_ASSERT (edge >= 0 && edge < 12);
          p8est_quadrant_edge_neighbor (&p, edge, &n);
          hanging[k] = (int) p4est_quadrant_exists
            (p4est, ghost_layer, jt, &n, &exist_array);
#ifdef P4EST_DEBUG
          ++check;
#endif
        }
      }
      P4EST_ASSERT (check <= 3);
#endif

      /* assign all independent and hanging nodes to the element */
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        P4EST_ASSERT (hanging[k] == 0 || hanging[k] == 1);
        p4est_quadrant_corner_node (qp[hanging[k]], k, &n);
        p4est_node_canonicalize (p4est, jt, &n, &c);

      }
    }
  }

  P4EST_FREE (local_nodes);
  sc_array_reset (&exist_array);
}

p4est_neighborhood_t *
p4est_neighborhood_new (p4est_t * p4est)
{
  bool                success;
  p4est_topidx_t      local_num_trees, flt, nt;
  p4est_locidx_t      local_num_quadrants, lsum;
  p4est_tree_t       *tree;
  p4est_neighborhood_t *nhood;
  sc_array_t          ghost_layer;

  P4EST_ASSERT (p4est_is_valid (p4est));

  sc_array_init (&ghost_layer, sizeof (p4est_quadrant_t));
  success = p4est_build_ghost_layer (p4est, &ghost_layer);
  P4EST_ASSERT (success);

  p4est_collect_nodes (p4est, &ghost_layer);
  sc_array_reset (&ghost_layer);

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

/* EOF p4est_mesh.c */
