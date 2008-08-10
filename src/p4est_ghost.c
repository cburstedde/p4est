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
#include <p4est_ghost.h>

#endif /* !P4_TO_P8 */

int
p4est_quadrant_find_owner (p4est_t * p4est, p4est_topidx_t treeid,
                           int face, const p4est_quadrant_t * q)
{
  const int           rank = p4est->mpirank;
  p4est_connectivity_t *conn = p4est->connectivity;
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
  if (face != -1) {
    P4EST_ASSERT (face >= 0 && face < 2 * P4EST_DIM);
    P4EST_ASSERT (treeid >= 0 && treeid < conn->num_trees);
    ntreeid = conn->tree_to_tree[2 * P4EST_DIM * treeid + face];
    if (ntreeid == treeid
        && ((int) conn->tree_to_face[2 * P4EST_DIM * treeid + face] ==
            face)) {
      /* This quadrant goes across a face with no neighbor */
      return -1;
    }
  }
  else {
    /* We need to determine the face ourselves */
    const p4est_qcoord_t rh = P4EST_ROOT_LEN;
    bool                quad_contact[2 * P4EST_DIM];

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

    /* Make sure we are neither tree edge nor tree corner */
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
  }

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

#ifdef P4EST_MPI

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

bool
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
    }
    return (lnid != -1);
  }
}

p4est_locidx_t
p4est_face_quadrant_exists (p4est_t * p4est, sc_array_t * ghost_layer,
                            p4est_topidx_t treeid, int *pface,
                            const p4est_quadrant_t * q, int *owner_rank)
{
  const int           rank = p4est->mpirank;
  int                 qproc;
  int                 nface, face = *pface;
  ssize_t             lnid;
  p4est_topidx_t      tqtreeid;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_quadrant_t    tq, non_existent;
  sc_array_t         *ta;
#ifndef P4_TO_P8
  int                 transform;
  p4est_quadrant_t    tempq;
#else
  int                 ftransform[9];
  p4est_topidx_t      tqtreeid2;
#endif

  P4EST_ASSERT (treeid >= 0 && 0 <= face && face < 2 * P4EST_DIM);

  P4EST_QUADRANT_INIT (&non_existent);
  if (non_existent.level == q->level) {
    return -1;
  }
  ta = NULL;

  /* q is in the unit domain */
  if (p4est_quadrant_is_inside_root (q)) {
    *pface = p4est_face_dual[face];
    *owner_rank = qproc = p4est_comm_find_owner (p4est, treeid, q, rank);
    if (qproc == rank) {
      p4est_tree_t       *tree;

      tree = p4est_array_index_topidx (p4est->trees, treeid);
      lnid = sc_array_bsearch (&tree->quadrants, q, p4est_quadrant_compare);
      return (lnid == -1) ? (p4est_locidx_t) (-1) :
        (tree->quadrants_offset + (p4est_locidx_t) lnid);
    }
    else {
      /* off processor so search in the ghost layer */
      tq = *q;
      tq.p.piggy1.which_tree = treeid;
      lnid =
        sc_array_bsearch (ghost_layer, &tq, p4est_quadrant_compare_piggy);
      return (lnid == -1) ? (p4est_locidx_t) (-1) :
        (q = sc_array_index_ssize_t (ghost_layer, lnid),
         q->p.piggy3.local_num);
    }
  }

  /* neighbor is across a tree face */
  tqtreeid = conn->tree_to_tree[2 * P4EST_DIM * treeid + face];
  nface = (int) conn->tree_to_face[2 * P4EST_DIM * treeid + face];
  *pface = nface % 6;
  if (tqtreeid == treeid && nface == face) {
    return -2;
  }

  /* transform quadrant */
#ifndef P4_TO_P8
  tempq = *q;
  transform = p4est_find_face_transform (conn, treeid, face);
  p4est_quadrant_translate_face (&tempq, face);
  p4est_quadrant_transform_face (&tempq, &tq, transform);
#else
  tqtreeid2 = p8est_find_face_transform (conn, treeid, face, ftransform);
  P4EST_ASSERT (tqtreeid == tqtreeid2);
  p8est_quadrant_transform_face (q, &tq, ftransform);
#endif

  /* find its owner and local number */
  *owner_rank = qproc = p4est_comm_find_owner (p4est, tqtreeid, &tq, rank);
  if (qproc == rank) {
    p4est_tree_t       *tqtree;

    tqtree = p4est_array_index_topidx (p4est->trees, tqtreeid);
    lnid = sc_array_bsearch (&tqtree->quadrants, &tq, p4est_quadrant_compare);
    return (lnid == -1) ? (p4est_locidx_t) (-1) :
      (tqtree->quadrants_offset + (p4est_locidx_t) lnid);
  }
  else {
    /* off processor so search in the ghost layer */
    tq.p.piggy1.which_tree = tqtreeid;
    lnid = sc_array_bsearch (ghost_layer, &tq, p4est_quadrant_compare_piggy);
    return (lnid == -1) ? (p4est_locidx_t) (-1) :
      (q = sc_array_index_ssize_t (ghost_layer, lnid), q->p.piggy3.local_num);
  }
}

/** Checks if a quadrant's face is on the boundary of the forest.
 *
 * \param [in] p4est  The forest in which to search for \a q
 * \param [in] treeid The tree id for which \a q belongs.
 * \param [in] q      The quadrant that is in question.
 * \param [in] face   The face of the quadrant that is in question.
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
  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (conn->tree_to_tree[2 * P4EST_DIM * treeid + face] != treeid ||
      (int) conn->tree_to_face[2 * P4EST_DIM * treeid + face] != face) {
    return false;
  }

  dh = P4EST_LAST_OFFSET (q->level);
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
  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (p4est_quadrant_touches_corner (q, corner, true)) {
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
                        p4est_locidx_t number, const p4est_quadrant_t * q)
{
  p4est_quadrant_t   *qold, *qnew;

  P4EST_ASSERT (treeid >= 0 && number >= 0);

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

  /* Cram the tree id and the local number into the user_data pointer */
  qnew->p.piggy3.which_tree = treeid;
  qnew->p.piggy3.local_num = number;
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
  if (!p4est_build_ghost_layer (p4est, true, &ghost_layer, NULL)) {
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
        p4est_quadrant_all_face_neighbors (q, face, n);
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

        face = p8est_child_edge_faces[qcid][edge];
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
        face = p8est_child_corner_faces[qcid][corner];
        edge = p8est_child_corner_edges[qcid][corner];
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
p4est_build_ghost_layer (p4est_t * p4est, bool include_diagonals,
                         sc_array_t * ghost_layer, int **ghost_owner)
{
#ifdef P4EST_MPI
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  int                 face, corner, zcorner;
  int                 i, ncheck, ncount;
  int                 n0_proc, n0ur_proc, n1_proc;
  int                 num_peers, peer, peer_proc;
  int                 mpiret;
  int                *gown;
  bool                maxed, failed;
  bool                full_tree[2], tree_contact[2 * P4EST_DIM];
  bool                urg[P4EST_DIM - 1];
  size_t              pz;
  p4est_topidx_t      nt;
  p4est_topidx_t      first_local_tree = p4est->first_local_tree;
  p4est_topidx_t      last_local_tree = p4est->last_local_tree;
  p4est_locidx_t      li, local_num;
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

  P4EST_GLOBAL_PRODUCTION ("Into " P4EST_STRING "_build_ghost_layer\n");
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
  local_num = 0;
  for (nt = first_local_tree; nt <= last_local_tree; ++nt) {
    tree = p4est_array_index_topidx (p4est->trees, nt);
    quadrants = &tree->quadrants;
    num_quadrants = (p4est_locidx_t) quadrants->elem_count;
    p4est_comm_tree_info (p4est, nt, full_tree, tree_contact, NULL, NULL);

    /* Find the smaller neighboring processors of each quadrant */
    for (li = 0; li < num_quadrants; ++local_num, ++li) {
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
          p4est_quadrant_half_face_neighbors (q, face, n, nur);
          ncheck = ncount = P4EST_CHILDREN / 2;
        }

        n1_proc = -1;
        for (i = 0; i < ncount; ++i) {
          n0_proc = p4est_quadrant_find_owner (p4est, nt, face, &n[i]);
          if (i < ncheck) {
            /* Note that we will always check this
             * because it prevents deadlocks
             */
            n0ur_proc = p4est_quadrant_find_owner (p4est, nt, face, &nur[i]);
            if (n0_proc != n0ur_proc) {
              P4EST_NOTICE ("Small face owner inconsistency\n");
              failed = true;
              goto failtest;
            }
          }

          if (n0_proc != rank && n0_proc >= 0 && n0_proc != n1_proc) {
            buf = sc_array_index_int (&send_bufs, n0_proc);
            p4est_add_ghost_to_buf (buf, nt, local_num, q);
            n1_proc = n0_proc;
          }
        }
      }

      /* Optionally skip edges and corners */
      if (!include_diagonals)
        continue;

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
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
            }

            if (!maxed) {
              n1_proc = *((int *) sc_array_index (&procs[1], pz));

              if (n1_proc != n0_proc && n1_proc != rank) {
                buf = sc_array_index_int (&send_bufs, n1_proc);
                p4est_add_ghost_to_buf (buf, nt, local_num, q);
              }
            }
          }
        }
        else {
          /* We are not at a tree edge so we only have two neighbors
           * either inside the tree or across a face
           */
          n0_proc = n1_proc =
            p4est_quadrant_find_owner (p4est, nt, -1, &n[0]);
          if (!maxed) {
            n1_proc = p4est_quadrant_find_owner (p4est, nt, -1, &n[1]);
            n0ur_proc = p4est_quadrant_find_owner (p4est, nt, -1, &nur[0]);
            n1ur_proc = p4est_quadrant_find_owner (p4est, nt, -1, &nur[1]);

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
            p4est_add_ghost_to_buf (buf, nt, local_num, q);
          }

          if (n1_proc != n0_proc && n1_proc != rank && n1_proc >= 0) {
            buf = sc_array_index_int (&send_bufs, n1_proc);
            p4est_add_ghost_to_buf (buf, nt, local_num, q);
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
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
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
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
            }
          }
        }
#endif
        else {
          /* We are not at a tree edge or corner so
           * we only have one corner neighbor
           */
          n0_proc = p4est_quadrant_find_owner (p4est, nt, -1, &n[0]);
          if (!maxed) {
            n0ur_proc = p4est_quadrant_find_owner (p4est, nt, -1, &nur[0]);

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
            p4est_add_ghost_to_buf (buf, nt, local_num, q);
          }
        }
      }
    }
  }
  P4EST_ASSERT (local_num == p4est->local_num_quadrants);

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
  if (ghost_owner != NULL) {
    *ghost_owner = gown = P4EST_ALLOC (int, num_ghosts);
  }

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
      if (ghost_owner != NULL) {
        for (li = 0; li < recv_counts[peer]; ++li) {
          gown[ghost_offset + li] = peer_proc;
        }
      }
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
    P4EST_ASSERT (q->p.piggy3.local_num >= 0);
    if (q2 != NULL) {
      P4EST_ASSERT (p4est_quadrant_compare_piggy (q2, q) < 0);
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

  P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_build_ghost_layer\n");

#else /* !P4EST_MPI */
  /* If we are not running with mpi then we don't need to do anything */
  sc_array_reset (ghost_layer);
#endif /* !PEST_MPI */

  return true;
}

/* EOF p4est_ghost.c */
