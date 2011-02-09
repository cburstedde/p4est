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

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_ghost.h>
#include <p4est_search.h>
#else
/* bits and communication are included in p8est_ghost.c */
#include <p8est_ghost.h>
#include <p8est_search.h>
#endif

/* htonl is in either of these two */
#ifdef P4EST_HAVE_ARPA_NET_H
#include <arpa/inet.h>
#endif
#ifdef P4EST_HAVE_NETINET_IN_H
#include <netinet/in.h>
#endif

typedef enum
{
  P4EST_GHOST_UNBALANCED_ABORT = 0,
  P4EST_GHOST_UNBALANCED_FAIL,
  P4EST_GHOST_UNBALANCED_ALLOW
}
p4est_ghost_tolerance_t;

size_t
p4est_ghost_memory_used (p4est_ghost_t * ghost)
{
  return sizeof (p4est_ghost_t) +
    sc_array_memory_used (&ghost->ghosts, 0) +
    (ghost->mpisize + 1) * sizeof (p4est_locidx_t) +
    (ghost->num_trees + 1) * sizeof (p4est_locidx_t);
}

#ifdef P4EST_MPI

static inline sc_array_t *
p4est_ghost_array_index (sc_array_t * array, int i)
{
  return (sc_array_t *) sc_array_index_int (array, i);
}

#endif

static p4est_ghost_t *p4est_ghost_new_check (p4est_t * p4est,
                                             p4est_balance_type_t btype,
                                             p4est_ghost_tolerance_t tol);

int
p4est_quadrant_find_owner (p4est_t * p4est, p4est_topidx_t treeid,
                           int face, const p4est_quadrant_t * q)
{
  const int           rank = p4est->mpirank;
  p4est_connectivity_t *conn = p4est->connectivity;
#ifdef P4EST_DEBUG
  int                 dims;
#endif
  int                 ftransform[P4EST_FTRANSFORM];
  p4est_topidx_t      ntreeid, ntreeid2;
  p4est_quadrant_t    nq;

  if (p4est_quadrant_is_inside_root (q)) {
    return p4est_comm_find_owner (p4est, treeid, q, rank);
  }

  P4EST_QUADRANT_INIT (&nq);

  /* We are outside of the unit tree */
  if (face != -1) {
    P4EST_ASSERT (face >= 0 && face < P4EST_FACES);
    P4EST_ASSERT (treeid >= 0 && treeid < conn->num_trees);
    ntreeid = conn->tree_to_tree[P4EST_FACES * treeid + face];
    if (ntreeid == treeid
        && ((int) conn->tree_to_face[P4EST_FACES * treeid + face] == face)) {
      /* This quadrant goes across a face with no neighbor */
      return -1;
    }
  }
  else {
    /* We need to determine the face ourselves */
    const p4est_qcoord_t rh = P4EST_ROOT_LEN;
    int                 quad_contact[P4EST_FACES];

    quad_contact[0] = (q->x < 0);
    quad_contact[1] = (q->x >= rh);
    quad_contact[2] = (q->y < 0);
    quad_contact[3] = (q->y >= rh);
#ifdef P4_TO_P8
    quad_contact[4] = (q->z < 0);
    quad_contact[5] = (q->z >= rh);
#endif

    /* Make sure we are neither tree edge nor tree corner */
#ifdef P4EST_DEBUG
    dims = (((quad_contact[0] || quad_contact[1]) ? 1 : 0) +
            ((quad_contact[2] || quad_contact[3]) ? 1 : 0));
#ifdef P4_TO_P8
    dims += ((quad_contact[4] || quad_contact[5]) ? 1 : 0);
#endif
    P4EST_ASSERT (dims == 1);
#endif

    ntreeid = -1;
    for (face = 0; face < P4EST_FACES; ++face) {
      if (quad_contact[face]) {
        ntreeid = conn->tree_to_tree[P4EST_FACES * treeid + face];
        if (ntreeid == treeid
            && ((int) conn->tree_to_face[P4EST_FACES * treeid + face] ==
                face)) {
          /* This quadrant goes across a face with no neighbor */
          return -1;
        }
        break;
      }
    }
    P4EST_ASSERT (face < P4EST_FACES && ntreeid >= 0);
  }

  ntreeid2 = p4est_find_face_transform (conn, treeid, face, ftransform);
  P4EST_ASSERT (ntreeid2 == ntreeid);
  p4est_quadrant_transform_face (q, &nq, ftransform);

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
                                        sc_array_t * q_procs, int *nurgood)
{
  const int           rank = p4est->mpirank;
  p4est_connectivity_t *conn = p4est->connectivity;
  int                *proc, nurproc;
  size_t              ctree;
  p4est_quadrant_t    cq;
  p4est_corner_info_t ci;
  p4est_corner_transform_t *ct;
  sc_array_t         *cta;

  P4EST_ASSERT (p4est_quadrant_is_outside_corner (q));

  P4EST_QUADRANT_INIT (&cq);

  /* Find all corners that are not myself or from a face neighbor */
  cta = &ci.corner_transforms;
  sc_array_init (cta, sizeof (p4est_corner_transform_t));
  p4est_find_corner_transform (conn, treeid, treecorner, &ci);

  sc_array_resize (q_procs, 0);
  if (nurgood != NULL) {
    *nurgood = 1;
    if (q->level == P4EST_QMAXLEVEL)
      nurgood = NULL;
  }

  for (ctree = 0; ctree < cta->elem_count; ++ctree) {
    ct = p4est_corner_array_index (cta, ctree);

    cq = *q;
    p4est_quadrant_transform_corner (&cq, (int) ct->ncorner, 1);

    proc = (int *) sc_array_push (q_procs);
    *proc = p4est_comm_find_owner (p4est, ct->ntree, &cq, rank);

    if (nurgood != NULL) {
      p4est_quadrant_last_descendant (&cq, &cq, P4EST_QMAXLEVEL);
      nurproc = p4est_comm_find_owner (p4est, ct->ntree, &cq, *proc);
      *nurgood = *nurgood && (nurproc == *proc);
    }
  }

  sc_array_reset (cta);
}

#endif /* P4EST_MPI */

ssize_t
p4est_ghost_tree_bsearch (p4est_ghost_t * ghost, p4est_topidx_t which_tree,
                          const p4est_quadrant_t * q)
{
  size_t              start, ended;
  ssize_t             result;
  sc_array_t          ghost_view;

  start = (size_t) ghost->tree_offsets[which_tree];
  ended = (size_t) ghost->tree_offsets[which_tree + 1];

  /* create a per-tree window on the ghost layer */
  sc_array_init_view (&ghost_view, &ghost->ghosts, start, ended - start);
  result = sc_array_bsearch (&ghost_view, q, p4est_quadrant_compare);

  /* and don't forget to add the window offset */
  return (result < 0) ? (ssize_t) (-1) : result + (ssize_t) start;
}

int
p4est_quadrant_exists (p4est_t * p4est, p4est_ghost_t * ghost,
                       p4est_topidx_t treeid, const p4est_quadrant_t * q,
                       sc_array_t * exists_arr)
{
  const int           rank = p4est->mpirank;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  int                 qproc, face, edge, corner;
  int                 ftransform[P4EST_FTRANSFORM];
  int                *pexists;
  int                 face_axis[3];     /* 3 not P4EST_DIM */
  int                 quad_contact[P4EST_FACES];
  int                 quad_corner;
  int                 exists;
  size_t              ctreeidz;
  ssize_t             lnid;
  p4est_topidx_t      tqtreeid;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, treeid);
  p4est_tree_t       *tqtree;
  p4est_quadrant_t    tq, non_existent;
  sc_array_t         *quadrants = &tree->quadrants;
  sc_array_t         *ta;
#ifdef P4_TO_P8
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
#endif
  p4est_corner_info_t ci;
  p4est_corner_transform_t *ct;

  if (exists_arr != NULL) {
    P4EST_ASSERT (exists_arr->elem_size == sizeof (int));
    sc_array_resize (exists_arr, 0);
  }
  P4EST_QUADRANT_INIT (&non_existent);
  ta = NULL;

  if (non_existent.level == q->level) {
    return 0;
  }

  /* q is in the unit domain */
  if (p4est_quadrant_is_inside_root (q)) {
    qproc = p4est_comm_find_owner (p4est, treeid, q, rank);
    if (qproc == rank) {
      lnid = sc_array_bsearch (quadrants, q, p4est_quadrant_compare);
    }
    else {
      lnid = p4est_ghost_tree_bsearch (ghost, treeid, q);
    }
    return (lnid != -1);
  }

  /* q is in a neighboring tree */
  quad_contact[0] = (q->x < 0);
  quad_contact[1] = (q->x >= rh);
  face_axis[0] = quad_contact[0] || quad_contact[1];
  quad_contact[2] = (q->y < 0);
  quad_contact[3] = (q->y >= rh);
  face_axis[1] = quad_contact[2] || quad_contact[3];
#ifndef P4_TO_P8
  face_axis[2] = 0;
#else
  quad_contact[4] = (q->z < 0);
  quad_contact[5] = (q->z >= rh);
  face_axis[2] = quad_contact[4] || quad_contact[5];
#endif
  quad_corner = 0;
  face = edge = corner = -1;
  P4EST_ASSERT (face_axis[0] || face_axis[1] || face_axis[2]);
  if (!face_axis[1] && !face_axis[2]) {
    face = 0 + quad_contact[1];
  }
  else if (!face_axis[0] && !face_axis[2]) {
    face = 2 + quad_contact[3];
  }
#ifdef P4_TO_P8
  else if (!face_axis[0] && !face_axis[1]) {
    face = 4 + quad_contact[5];
  }
  else if (!face_axis[0]) {
    edge = 0 + 2 * quad_contact[5] + quad_contact[3];
    quad_corner = 1;
  }
  else if (!face_axis[1]) {
    edge = 4 + 2 * quad_contact[5] + quad_contact[1];
    quad_corner = 1;
  }
  else if (!face_axis[2]) {
    edge = 8 + 2 * quad_contact[3] + quad_contact[1];
    quad_corner = 1;
  }
#endif
  else {
    corner =
#ifdef P4_TO_P8
      4 * quad_contact[5] +
#endif
      2 * quad_contact[3] + quad_contact[1];
    quad_corner = 1;
  }
  if (quad_corner) {
    /* Neighbor is across a tree edge or corner */
    P4EST_ASSERT (exists_arr != NULL);
    P4EST_ASSERT (face == -1 &&
                  ((edge >= 0 && corner == -1) ||
                   (edge == -1 && corner >= 0)));
    if (corner >= 0) {
      ta = &ci.corner_transforms;
      sc_array_init (ta, sizeof (p4est_corner_transform_t));
      p4est_find_corner_transform (conn, treeid, corner, &ci);
    }
#ifdef P4_TO_P8
    else {
      ta = &ei.edge_transforms;
      sc_array_init (ta, sizeof (p8est_edge_transform_t));
      p8est_find_edge_transform (conn, treeid, edge, &ei);
    }
#endif
    sc_array_resize (exists_arr, ta->elem_count);

    exists = 0;
    for (ctreeidz = 0; ctreeidz < ta->elem_count; ++ctreeidz) {
      if (corner >= 0) {
        ct = p4est_corner_array_index (ta, ctreeidz);
        tqtreeid = ct->ntree;
        tq = *q;
        p4est_quadrant_transform_corner (&tq, (int) ct->ncorner, 1);
      }
#ifdef P4_TO_P8
      else {
        et = p8est_edge_array_index (ta, ctreeidz);
        tqtreeid = et->ntree;
        p8est_quadrant_transform_edge (q, &tq, &ei, et, 1);
      }
      et = NULL;
#endif
      ct = NULL;

      qproc = p4est_comm_find_owner (p4est, tqtreeid, &tq, rank);

      if (qproc == rank) {
        tqtree = p4est_tree_array_index (p4est->trees, tqtreeid);
        lnid = sc_array_bsearch (&tqtree->quadrants, &tq,
                                 p4est_quadrant_compare);
      }
      else {
        lnid = p4est_ghost_tree_bsearch (ghost, tqtreeid, &tq);
      }

      /* add the existence value */
      pexists = (int *) sc_array_index (exists_arr, ctreeidz);
      *pexists = (lnid != -1);
      exists = exists || *pexists;
    }

    sc_array_reset (ta);
    return exists;
  }
  else {
    /* Neighbor is across a tree face */
    P4EST_ASSERT (face >= 0 && edge == -1 && corner == -1);
    P4EST_ASSERT (quad_contact[face]);

    tqtreeid = p4est_find_face_transform (conn, treeid, face, ftransform);
    if (tqtreeid == -1) {
      /* there is no tree neighbor across this face */
      return 0;
    }
    p4est_quadrant_transform_face (q, &tq, ftransform);

    /* find owner of the transformed quadrant */
    qproc = p4est_comm_find_owner (p4est, tqtreeid, &tq, rank);
    if (qproc == rank) {
      tqtree = p4est_tree_array_index (p4est->trees, tqtreeid);
      lnid = sc_array_bsearch (&tqtree->quadrants, &tq,
                               p4est_quadrant_compare);
    }
    else {
      lnid = p4est_ghost_tree_bsearch (ghost, tqtreeid, &tq);
    }
    return (lnid != -1);
  }
}

p4est_locidx_t
p4est_face_quadrant_exists (p4est_t * p4est, p4est_ghost_t * ghost,
                            p4est_topidx_t treeid, const p4est_quadrant_t * q,
                            int *pface, int *phang, int *owner_rank)
{
  const int           rank = p4est->mpirank;
  int                 qproc;
  int                 nface, face, orientation;
  int                 ftransform[P4EST_FTRANSFORM];
  ssize_t             lnid;
  p4est_topidx_t      tqtreeid, tqtreeid2;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_quadrant_t    tq, non_existent;
#ifdef P4_TO_P8
  int                 face_ref, face_perm;
#endif

  face = *pface;
  P4EST_ASSERT (treeid >= 0 && 0 <= face && face < P4EST_FACES);

  P4EST_QUADRANT_INIT (&non_existent);
  if (non_existent.level == q->level) {
    return -1;
  }

  /* determine the hanging face number */
  if (phang != NULL) {
    P4EST_ASSERT (*phang >= 0 && *phang < P4EST_CHILDREN);
    *phang = p4est_corner_face_corners[*phang][face];
    P4EST_ASSERT (*phang >= 0 && *phang < P4EST_HALF);
  }

  /* q is in the unit domain */
  if (p4est_quadrant_is_inside_root (q)) {
    *pface = p4est_face_dual[face];
    *owner_rank = qproc = p4est_comm_find_owner (p4est, treeid, q, rank);
    if (qproc == rank) {
      p4est_tree_t       *tree;

      tree = p4est_tree_array_index (p4est->trees, treeid);
      lnid = sc_array_bsearch (&tree->quadrants, q, p4est_quadrant_compare);
      return (lnid == -1) ? (p4est_locidx_t) (-1) :
        (tree->quadrants_offset + (p4est_locidx_t) lnid);
    }
    else {
      lnid = p4est_ghost_tree_bsearch (ghost, treeid, q);
      return (lnid == -1) ? (p4est_locidx_t) (-1) :
        (q = p4est_quadrant_array_index (&ghost->ghosts, (size_t) lnid),
         q->p.piggy3.local_num);
    }
  }

  /* neighbor is across a tree face */
  tqtreeid = conn->tree_to_tree[P4EST_FACES * treeid + face];
  nface = (int) conn->tree_to_face[P4EST_FACES * treeid + face];
  if (tqtreeid == treeid && nface == face) {
    *owner_rank = -1;
    *pface = -1;
    if (phang) {
      *phang = -1;
    }
    return -2;
  }

  /* transform the hanging face number */
  *pface = nface;
  if (phang != NULL) {
    orientation = nface / P4EST_FACES;
#ifdef P4_TO_P8
    face_ref = p8est_face_permutation_refs[face][nface % P4EST_FACES];
    face_perm = p8est_face_permutation_sets[face_ref][orientation];
    *phang = p8est_face_permutations[face_perm][*phang];
#else
    *phang = *phang ^ orientation;
#endif
  }

  /* transform quadrant */
  tqtreeid2 = p4est_find_face_transform (conn, treeid, face, ftransform);
  P4EST_ASSERT (tqtreeid == tqtreeid2);
  p4est_quadrant_transform_face (q, &tq, ftransform);

  /* find its owner and local number */
  *owner_rank = qproc = p4est_comm_find_owner (p4est, tqtreeid, &tq, rank);
  if (qproc == rank) {
    p4est_tree_t       *tqtree;

    tqtree = p4est_tree_array_index (p4est->trees, tqtreeid);
    lnid = sc_array_bsearch (&tqtree->quadrants, &tq, p4est_quadrant_compare);
    return (lnid == -1) ? (p4est_locidx_t) (-1) :
      (tqtree->quadrants_offset + (p4est_locidx_t) lnid);
  }
  else {
    lnid = p4est_ghost_tree_bsearch (ghost, tqtreeid, &tq);
    return (lnid == -1) ? (p4est_locidx_t) (-1) :
      (q =
       p4est_quadrant_array_index (&ghost->ghosts, (size_t) lnid),
       q->p.piggy3.local_num);
  }
}

/** Checks if a quadrant's face is on the boundary of the forest.
 *
 * \param [in] p4est  The forest in which to search for \a q
 * \param [in] treeid The tree to which \a q belongs.
 * \param [in] q      The quadrant that is in question.
 * \param [in] face   The face of the quadrant that is in question.
 *
 * \return true if the quadrant's face is on the boundary of the forest and
 *         false otherwise.
 */
static int
p4est_quadrant_on_face_boundary (p4est_t * p4est, p4est_topidx_t treeid,
                                 int face, const p4est_quadrant_t * q)
{
  p4est_qcoord_t      dh, xyz;
  p4est_connectivity_t *conn = p4est->connectivity;

  P4EST_ASSERT (0 <= face && face < P4EST_FACES);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (conn->tree_to_tree[P4EST_FACES * treeid + face] != treeid ||
      (int) conn->tree_to_face[P4EST_FACES * treeid + face] != face) {
    return 0;
  }

  dh = P4EST_LAST_OFFSET (q->level);
  switch (face / 2) {
  case 0:
    xyz = q->x;
    break;
  case 1:
    xyz = q->y;
    break;
#ifdef P4_TO_P8
  case 2:
    xyz = q->z;
    break;
#endif
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }
  return xyz == ((face & 0x01) ? dh : 0);
}

/** Get the smallest corner neighbor of \a q.
 *
 * Gets the smallest corner neighbor, which is half of the size assuming the
 * 2-1 constaint.
 *
 * \param [in]  q      The quadrant whose corner neighbor will be constructed.
 * \param [in]  corner The corner across which to generate the neighbor.
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
 * \param [in]  q       The quadrant whose face neighbors will be constructed.
 * \param [in]  corner  The corner across which to generate the neighbors.
 * \param [out] n[0]    Filled with the possible half-size corner neighbor
 *                      if it exists or initialized by P4EST_QUADRANT_INIT.
 * \param [out] n[1]    Filled with the possible same-size corner neighbor.
 * \param [out] n[2]    Filled with the possible double-size corner neighbor
 *                      if it exists or initialized by P4EST_QUADRANT_INIT.
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
 * \param [in] treeid The tree to which \a q belongs.
 * \param [in] q      The quadrant that is in question.
 * \param [in] corner The corner of quadrant that is in question.
 *
 * \return true if the quadrant's corner is on the boundary of the forest and
 *         false otherwise.
 */
static int
p4est_quadrant_on_corner_boundary (p4est_t * p4est, p4est_topidx_t treeid,
                                   int corner, const p4est_quadrant_t * q)
{
  p4est_connectivity_t *conn = p4est->connectivity;
  int                 face;
  int                 on_boundary = 0;
  p4est_quadrant_t    q2;
#ifdef P4_TO_P8
  int                 edge;
  p8est_edge_info_t   ei;
  sc_array_t         *eta;
#endif
  p4est_corner_info_t ci;
  sc_array_t         *cta;

  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (p4est_quadrant_touches_corner (q, corner, 1)) {
    cta = &ci.corner_transforms;
    sc_array_init (cta, sizeof (p4est_corner_transform_t));
    p4est_find_corner_transform (conn, treeid, corner, &ci);

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
    return 0;
  }

  return
    (conn->tree_to_tree[P4EST_FACES * treeid + face] == treeid &&
     (int) conn->tree_to_face[P4EST_FACES * treeid + face] == face);
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
    qold = p4est_quadrant_array_index (buf, buf->elem_count - 1);
    if (treeid == qold->p.piggy1.which_tree &&
        p4est_quadrant_compare (q, qold) == 0) {
      return;
    }
  }

  qnew = p4est_quadrant_array_push (buf);
  *qnew = *q;

  /* Cram the tree id and the local number into the user_data pointer */
  qnew->p.piggy3.which_tree = treeid;
  qnew->p.piggy3.local_num = number;
}

#endif /* P4EST_MPI */

int
p4est_is_balanced (p4est_t * p4est, p4est_balance_type_t btype)
{
  int                 zero = 0;
  int                 face, corner;
  int                 ii, qcid;
  int                *pe0, *pe1, *pe2;
  int                 failed;
  int                 e0, e1, e0b, e1b, e2, e3;
  int                 bigger_face[P4EST_FACES];
  size_t              cez, zz;
  p4est_topidx_t      nt;
  p4est_topidx_t      first_local_tree = p4est->first_local_tree;
  p4est_topidx_t      last_local_tree = p4est->last_local_tree;
  p4est_quadrant_t   *q;
  p4est_quadrant_t    n[P4EST_HALF + 2];
  p4est_tree_t       *tree;
  sc_array_t         *quadrants;
  sc_array_t          e0_a, e1_a, e2_a;
#ifdef P4_TO_P8
  int                 edge;
  int                *pe3;
  int                 bigger_edge[12];
  sc_array_t          e3_a;
#ifdef P4EST_DEBUG
  size_t              big_count[12];
#endif
#endif
  p4est_ghost_t      *gl;

  gl = p4est_ghost_new_check (p4est, btype, P4EST_GHOST_UNBALANCED_FAIL);
  if (gl == NULL) {
    return 0;
  }

  for (ii = 0; ii < P4EST_HALF + 2; ++ii) {
    P4EST_QUADRANT_INIT (&n[ii]);
  }

  failed = 0;
  sc_array_init (&e0_a, sizeof (int));
  sc_array_init (&e1_a, sizeof (int));
  sc_array_init (&e2_a, sizeof (int));
#ifdef P4_TO_P8
  sc_array_init (&e3_a, sizeof (int));
#endif

  /* loop over all local trees */
  for (nt = first_local_tree; nt <= last_local_tree; ++nt) {
    tree = p4est_tree_array_index (p4est->trees, nt);
    quadrants = &tree->quadrants;

    /* Find the neighboring processors of each quadrant */
    for (zz = 0; zz < quadrants->elem_count; ++zz) {
      q = p4est_quadrant_array_index (quadrants, zz);
      qcid = p4est_quadrant_child_id (q);

      /* Find face neighbors */
      for (face = 0; face < P4EST_FACES; ++face) {
        bigger_face[face] = 0;

        /* If q is at a boundary then it is automatically balanced */
        if (p4est_quadrant_on_face_boundary (p4est, nt, face, q)) {
          continue;
        }

        /* Do more expensive face balance checks */
        p4est_quadrant_all_face_neighbors (q, face, n);
        e0 = p4est_quadrant_exists (p4est, gl, nt, &n[0], NULL);
        e1 = p4est_quadrant_exists (p4est, gl, nt, &n[1], NULL);
#ifndef P4_TO_P8
        e0b = e1b = e0;
#else
        e0b = p4est_quadrant_exists (p4est, gl, nt, &n[2], NULL);
        e1b = p4est_quadrant_exists (p4est, gl, nt, &n[3], NULL);
#endif
        if (e0 != e1 || e0 != e0b || e0 != e1b) {
          P4EST_NOTICE ("Contradicting small face neighbors\n");
          failed = 1;
          goto failtest;
        }
        e2 = p4est_quadrant_exists (p4est, gl, nt, &n[P4EST_HALF], NULL);
        e3 = p4est_quadrant_exists (p4est, gl, nt, &n[P4EST_HALF + 1], NULL);
        if ((int) e0 + (int) e2 + (int) e3 != 1) {
          P4EST_NOTICE ("Face balance quadrant mismatch\n");
          failed = 1;
          goto failtest;
        }
        bigger_face[face] = e3;
      }

#ifndef P4_TO_P8
      if (btype == P4EST_BALANCE_FACE)
        continue;
#else
      if (btype == P8EST_BALANCE_FACE)
        continue;

      /* Find edge neighbors */
      for (edge = 0; edge < P8EST_EDGES; ++edge) {
        bigger_edge[edge] = 0;
#ifdef P4EST_DEBUG
        big_count[edge] = 0;
#endif

        /* If q is at a boundary then it is automatically balanced */
        if (p8est_quadrant_on_edge_boundary (p4est, nt, edge, q)) {
          continue;
        }

        /* Do more expensive edge balance checks */
        p8est_quadrant_get_possible_edge_neighbors (q, edge, n);
        e0 = p4est_quadrant_exists (p4est, gl, nt, &n[0], &e0_a);
        e1 = p4est_quadrant_exists (p4est, gl, nt, &n[1], &e1_a);
        if (e0 != e1 || e0_a.elem_count != e1_a.elem_count) {
          P4EST_NOTICE ("Contradicting small edge neighbors\n");
          failed = 1;
          goto failtest;
        }
        e2 = p4est_quadrant_exists (p4est, gl, nt, &n[2], &e2_a);
        e3 = p4est_quadrant_exists (p4est, gl, nt, &n[3], &e3_a);
        P4EST_ASSERT (((e0_a.elem_count == 0 && q->level == P4EST_QMAXLEVEL)
                       || e0_a.elem_count == e2_a.elem_count)
                      && ((e3_a.elem_count == 0 && (!e3 || q->level == 0))
                          || e3_a.elem_count == e2_a.elem_count));

        face = p8est_child_edge_faces[qcid][edge];
        if (face >= 0 && bigger_face[face]) {
          P4EST_ASSERT (e2_a.elem_count == 0);
          if (e0 || e2 || e3) {
            P4EST_NOTICE ("Invalid edges across hanging face\n");
            failed = 1;
            goto failtest;
          }
        }
        else {
          if (!e0 && !e2 && !e3) {
            P4EST_NOTICE ("Edge balance missing quadrants\n");
            failed = 1;
            goto failtest;
          }
          if (e2_a.elem_count == 0 && (int) e0 + (int) e2 + (int) e3 != 1) {
            P4EST_NOTICE ("Edge balance duplicate quadrants\n");
            failed = 1;
            goto failtest;
          }
          for (cez = 0; cez < e2_a.elem_count; ++cez) {
            pe0 = (e0_a.elem_count > 0) ?
              (int *) sc_array_index (&e0_a, cez) : &zero;
            pe1 = (e1_a.elem_count > 0) ?
              (int *) sc_array_index (&e1_a, cez) : &zero;
            pe2 = (int *) sc_array_index (&e2_a, cez);
            pe3 = (e3_a.elem_count > 0) ?
              (int *) sc_array_index (&e3_a, cez) : &zero;
            if (*pe0 + *pe2 + *pe3 != 1 || *pe0 != *pe1) {
              P4EST_NOTICE ("Edge balance quadrant mismatch\n");
              failed = 1;
              goto failtest;
            }
          }
        }
        bigger_edge[edge] = e3;
#ifdef P4EST_DEBUG
        big_count[edge] = e3_a.elem_count;
#endif
      }

      if (btype == P8EST_BALANCE_EDGE)
        continue;
#endif

      /* Find corner neighbors */
      for (corner = 0; corner < P4EST_CHILDREN; ++corner) {

        /* If q is at a boundary then it is automatically balanced */
        if (p4est_quadrant_on_corner_boundary (p4est, nt, corner, q)) {
          continue;
        }

        /* Do more expensive corner balance checks */
        p4est_quadrant_get_possible_corner_neighbors (q, corner, n);
        e0 = p4est_quadrant_exists (p4est, gl, nt, &n[0], &e0_a);
        e1 = p4est_quadrant_exists (p4est, gl, nt, &n[1], &e1_a);
        e2 = p4est_quadrant_exists (p4est, gl, nt, &n[2], &e2_a);
        P4EST_ASSERT (((e0_a.elem_count == 0 && q->level == P4EST_QMAXLEVEL)
                       || e0_a.elem_count == e1_a.elem_count)
                      && ((e2_a.elem_count == 0 && (!e2 || q->level == 0))
                          || e2_a.elem_count == e1_a.elem_count));

        face = p4est_child_corner_faces[qcid][corner];
#ifdef P4_TO_P8
        edge = p8est_child_corner_edges[qcid][corner];
#endif
        if (face >= 0 && bigger_face[face]) {
          P4EST_ASSERT (e1_a.elem_count == 0);
          if (e0 || e1 || e2) {
            P4EST_NOTICE ("Invalid corners across hanging face\n");
            failed = 1;
            goto failtest;
          }
        }
#ifndef P4_TO_P8
        else if (!e0 && !e1 && !e2) {
          P4EST_NOTICE ("Corner balance missing quadrants\n");
          failed = 1;
          goto failtest;
        }
#endif
        else if (e1_a.elem_count == 0) {
#ifdef P4_TO_P8
          if (edge >= 0 && bigger_edge[edge]) {
            P4EST_ASSERT (big_count[edge] == 0);
            if (e0 || e1 || e2) {
              P4EST_NOTICE ("Invalid corners across hanging edge\n");
              failed = 1;
              goto failtest;
            }
          }
          else
#endif
          if ((int) e0 + (int) e1 + (int) e2 != 1) {
            P4EST_NOTICE ("Corner balance quadrant mismatch\n");
            failed = 1;
            goto failtest;
          }
        }
        else {
#ifdef P4_TO_P8
          e3 = 0;
          if (edge >= 0 && bigger_edge[edge]) {
            P4EST_ASSERT (big_count[edge] == e1_a.elem_count);

            /* recreate the edge neighbor information */
            p4est_quadrant_parent (q, &n[3]);
            p8est_quadrant_edge_neighbor (&n[3], edge, &n[3]);
            e3 = p4est_quadrant_exists (p4est, gl, nt, &n[3], &e3_a);
            P4EST_ASSERT (e3 && big_count[edge] == e3_a.elem_count);
          }
#endif
          for (cez = 0; cez < e1_a.elem_count; ++cez) {
            pe0 = (e0_a.elem_count > 0) ?
              (int *) sc_array_index (&e0_a, cez) : &zero;
            pe1 = (int *) sc_array_index (&e1_a, cez);
            pe2 = (e2_a.elem_count > 0) ?
              (int *) sc_array_index (&e2_a, cez) : &zero;
#ifdef P4_TO_P8
            if (e3) {
              pe3 = (int *) sc_array_index (&e3_a, cez);
              if (*pe3) {
                if (*pe0 || *pe1 || *pe2) {
                  P4EST_NOTICE ("Invalid corners across hanging edge\n");
                  failed = 1;
                  goto failtest;
                }
                continue;
              }
            }
#endif
            if (*pe0 + *pe1 + *pe2 != 1) {
              P4EST_NOTICE ("Corner balance quadrant mismatch\n");
              failed = 1;
              goto failtest;
            }
          }
        }
      }

      P4EST_ASSERT (btype == P4EST_BALANCE_FULL);
    }
  }

failtest:
  sc_array_reset (&e0_a);
  sc_array_reset (&e1_a);
  sc_array_reset (&e2_a);
#ifdef P4_TO_P8
  sc_array_reset (&e3_a);
#endif
  p4est_ghost_destroy (gl);

  return !p4est_comm_sync_flag (p4est, failed, MPI_BOR);
}

static              size_t
ghost_tree_type (sc_array_t * array, size_t zindex, void *data)
{
  p4est_quadrant_t   *q;

  P4EST_ASSERT (array->elem_size == sizeof (p4est_quadrant_t));

  q = (p4est_quadrant_t *) sc_array_index (array, zindex);
  return (size_t) q->p.which_tree;
}

#ifdef P4EST_MPI

static void
p4est_ghost_test_add (p4est_t * p4est, p4est_quadrant_t * q, p4est_topidx_t t,
                      p4est_quadrant_t * nq, p4est_topidx_t nt, int32_t touch,
                      int rank, sc_array_t * send_bufs,
                      p4est_locidx_t local_num)
{
  p4est_quadrant_t    temp;
  p4est_quadrant_t   *lq, *uq;
  int64_t             next_lid, uid;
  int                 n0_proc, n1_proc, proc;
  p4est_quadrant_t   *gfp = p4est->global_first_position;
  sc_array_t         *buf;
  int32_t             rb;

  P4EST_ASSERT (q->level == nq->level);
  n0_proc = p4est_comm_find_owner (p4est, nt, nq, rank);
  P4EST_ASSERT (n0_proc >= 0);
  if (q->level == P4EST_QMAXLEVEL) {
    if (n0_proc != rank) {
      buf = p4est_ghost_array_index (send_bufs, n0_proc);
      p4est_add_ghost_to_buf (buf, t, local_num, q);
    }
    return;
  }
  p4est_quadrant_last_descendant (nq, &temp, P4EST_QMAXLEVEL);
  n1_proc = p4est_comm_find_owner (p4est, nt, &temp, n0_proc);
  P4EST_ASSERT (n1_proc >= n0_proc);
  if (n0_proc == n1_proc) {
    if (n0_proc != rank) {
      buf = p4est_ghost_array_index (send_bufs, n0_proc);
      p4est_add_ghost_to_buf (buf, t, local_num, q);
    }
    return;
  }
  for (proc = n0_proc; proc <= n1_proc; proc++) {
    if (proc == rank) {
      continue;
    }
    lq = &(gfp[proc]);
    uq = &(gfp[proc + 1]);
    /* check for empty processor */
    if (p4est_quadrant_is_equal_piggy (lq, uq)) {
      continue;
    }
    if (proc == n0_proc) {
      lq = NULL;
    }
    else {
      P4EST_ASSERT (p4est_quadrant_is_valid (lq));
      P4EST_ASSERT (lq->p.which_tree == nt);
      P4EST_ASSERT (p4est_quadrant_is_ancestor (nq, lq) ||
                    p4est_quadrant_is_equal (nq, lq));
    }
    if (proc == n1_proc) {
      uq = NULL;
    }
    else {
      P4EST_ASSERT (p4est_quadrant_is_valid (uq));
      P4EST_ASSERT (uq->p.which_tree == nt);
      P4EST_ASSERT (p4est_quadrant_is_ancestor (nq, uq) ||
                    p4est_quadrant_is_equal (nq, uq));
      next_lid = p4est_quadrant_linear_id (uq, P4EST_QMAXLEVEL);
      P4EST_ASSERT (next_lid > 0);
      uid = next_lid - 1;
      uq = &temp;
      p4est_quadrant_set_morton (uq, P4EST_QMAXLEVEL, uid);
      P4EST_ASSERT (p4est_quadrant_is_valid (uq));
    }
#ifdef P4EST_DEBUG
    if (lq != NULL && uq != NULL) {
      P4EST_ASSERT (p4est_quadrant_compare (lq, uq) <= 0);
    }
#endif
    rb = p4est_find_range_boundaries (lq, uq, (int) q->level,
#ifdef P4_TO_P8
                                      NULL,
#endif
                                      NULL, NULL);
    if (rb & touch) {
      buf = p4est_ghost_array_index (send_bufs, proc);
      p4est_add_ghost_to_buf (buf, t, local_num, q);
    }
  }
}

#endif /* P4EST_MPI */

static p4est_ghost_t *
p4est_ghost_new_check (p4est_t * p4est, p4est_balance_type_t btype,
                       p4est_ghost_tolerance_t tol)
{
  const p4est_topidx_t num_trees = p4est->connectivity->num_trees;
  const int           num_procs = p4est->mpisize;
#ifdef P4EST_MPI
  const int           rank = p4est->mpirank;
  MPI_Comm            comm = p4est->mpicomm;
  p4est_connectivity_t *conn = p4est->connectivity;
  int                 face, corner;
  int                 nface, ncheck, ncount;
  int                 i;
  int                 n0_proc, n0ur_proc, n1_proc;
  int                 num_peers, peer, peer_proc;
  int                 mpiret;
  int                 maxed, failed;
  int                 full_tree[2], tree_contact[2 * P4EST_DIM];
  int                 urg[P4EST_DIM - 1];
  size_t              pz, zz;
  p4est_topidx_t      first_local_tree = p4est->first_local_tree;
  p4est_topidx_t      last_local_tree = p4est->last_local_tree;
#ifdef P4EST_DEBUG
  p4est_locidx_t      li;
#endif
  p4est_locidx_t      local_num;
  p4est_locidx_t      num_ghosts, ghost_offset, skipped;
  p4est_locidx_t     *send_counts, *recv_counts;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  p4est_quadrant_t    n[P4EST_HALF], nur[P4EST_HALF];
  sc_array_t          send_bufs;
  sc_array_t          procs[P4EST_DIM - 1];
  sc_array_t         *buf, *quadrants;
  MPI_Request        *recv_request, *send_request;
  MPI_Status         *recv_status, *send_status;
  MPI_Request        *recv_load_request, *send_load_request;
  MPI_Status         *recv_load_status, *send_load_status;
#ifdef P4_TO_P8
  int                 edge, nedge;
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  sc_array_t         *eta;
  size_t              etree;
  int                 o, ref, set;
  int                 c0, c1;
  int                 nc0, nc1;
  int                 oppedge;
  int                 n1ur_proc;
#endif
#ifdef P4EST_DEBUG
  p4est_quadrant_t   *q2;
#endif
  int32_t             touch;
  p4est_topidx_t      nnt;
  int                 ftransform[P4EST_FTRANSFORM];
  p4est_corner_info_t ci;
  p4est_corner_transform_t *ct;
  sc_array_t         *cta;
  size_t              ctree;
#endif
  size_t             *ppz;
  sc_array_t          split;
  sc_array_t         *ghost_layer;
  p4est_topidx_t      nt;
  p4est_ghost_t      *gl;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING "_ghost_new %s\n",
                            p4est_balance_type_string (btype));

  gl = P4EST_ALLOC (p4est_ghost_t, 1);
  gl->mpisize = num_procs;
  gl->num_trees = num_trees;
  ghost_layer = &gl->ghosts;
  sc_array_init (ghost_layer, sizeof (p4est_quadrant_t));
  gl->tree_offsets = P4EST_ALLOC (p4est_locidx_t, num_trees + 1);
  gl->proc_offsets = P4EST_ALLOC (p4est_locidx_t, num_procs + 1);
  gl->proc_offsets[0] = 0;
#ifndef P4EST_MPI
  gl->proc_offsets[1] = 0;
#else
#ifdef P4_TO_P8
  eta = &ei.edge_transforms;
#endif
  cta = &ci.corner_transforms;

  for (i = 0; i < P4EST_HALF; ++i) {
    P4EST_QUADRANT_INIT (&n[i]);
    P4EST_QUADRANT_INIT (&nur[i]);
  }

  failed = 0;
  for (i = 0; i < P4EST_DIM - 1; ++i) {
    sc_array_init (&procs[i], sizeof (int));
  }
  skipped = 0;

  /* allocate empty send buffers */
  sc_array_init (&send_bufs, sizeof (sc_array_t));
  sc_array_resize (&send_bufs, (size_t) num_procs);
  for (i = 0; i < num_procs; ++i) {
    buf = p4est_ghost_array_index (&send_bufs, i);
    sc_array_init (buf, sizeof (p4est_quadrant_t));
  }

  /* loop over all local trees */
  local_num = 0;
  for (nt = first_local_tree; nt <= last_local_tree; ++nt) {
    tree = p4est_tree_array_index (p4est->trees, nt);
    quadrants = &tree->quadrants;
    p4est_comm_tree_info (p4est, nt, full_tree, tree_contact, NULL, NULL);

    /* Find the smaller neighboring processors of each quadrant */
    for (zz = 0; zz < quadrants->elem_count; ++local_num, ++zz) {
      q = p4est_quadrant_array_index (quadrants, zz);

      if (p4est_comm_neighborhood_owned
          (p4est, nt, full_tree, tree_contact, q)) {
        /* The 3x3 neighborhood of q is owned by this processor */
        ++skipped;
        continue;
      }

      /* Find smaller face neighbors */
      for (face = 0; face < 2 * P4EST_DIM; ++face) {
        if (tol < P4EST_GHOST_UNBALANCED_ALLOW) {
          if (q->level == P4EST_QMAXLEVEL) {
            p4est_quadrant_face_neighbor (q, face, &n[0]);
            ncheck = 0;
            ncount = 1;
          }
          else {
            p4est_quadrant_half_face_neighbors (q, face, n, nur);
            ncheck = ncount = P4EST_HALF;
          }

          n1_proc = -1;
          for (i = 0; i < ncount; ++i) {
            n0_proc = p4est_quadrant_find_owner (p4est, nt, face, &n[i]);
            if (i < ncheck) {
              /* Note that we will always check this
               * because it prevents deadlocks
               */
              n0ur_proc = p4est_quadrant_find_owner (p4est, nt, face,
                                                     &nur[i]);
              if (n0_proc != n0ur_proc) {
                P4EST_NOTICE ("Small face owner inconsistency\n");
                failed = 1;
                goto failtest;
              }
            }

            if (n0_proc != rank && n0_proc >= 0 && n0_proc != n1_proc) {
              buf = p4est_ghost_array_index (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
              n1_proc = n0_proc;
            }
          }
        }
        else {
          p4est_quadrant_face_neighbor (q, face, &n[0]);
          if (p4est_quadrant_is_inside_root (&n[0])) {
            nface = face ^ 1;
            touch = ((int32_t) 1 << nface);
            p4est_ghost_test_add (p4est, q, nt, &n[0], nt, touch, rank,
                                  &send_bufs, local_num);
          }
          else {
            nnt = p4est_find_face_transform (conn, nt, face, ftransform);
            if (nnt < 0) {
              continue;
            }
            nface = (int) conn->tree_to_face[nt * P4EST_FACES + face];
            nface %= P4EST_FACES;
            touch = ((int32_t) 1 << nface);
            p4est_quadrant_transform_face (&n[0], &n[1], ftransform);
            p4est_ghost_test_add (p4est, q, nt, &n[1], nnt, touch, rank,
                                  &send_bufs, local_num);
          }
        }
      }

#ifndef P4_TO_P8
      if (btype == P4EST_BALANCE_FACE)
        continue;
#else
      if (btype == P8EST_BALANCE_FACE)
        continue;

      /* Find smaller edge neighbors */
      for (edge = 0; edge < 12; ++edge) {
        if (tol < P4EST_GHOST_UNBALANCED_ALLOW) {
          if (q->level == P4EST_QMAXLEVEL) {
            p8est_quadrant_edge_neighbor (q, edge, &n[0]);
            maxed = 1;
          }
          else {
            p8est_quadrant_get_half_edge_neighbors (q, edge, n, nur);
            maxed = 0;
          }

          /* Check to see if we are a tree edge neighbor */
          P4EST_ASSERT (!p4est_quadrant_is_outside_corner (&n[0]));
          if (p8est_quadrant_is_outside_edge (&n[0])) {
            p8est_quadrant_find_tree_edge_owners (p4est, nt, edge,
                                                  &n[0], &procs[0], &urg[0]);
            if (!maxed) {
              p8est_quadrant_find_tree_edge_owners (p4est, nt, edge,
                                                    &n[1], &procs[1],
                                                    &urg[1]);
              P4EST_ASSERT (procs[0].elem_count == procs[1].elem_count);

              if (!urg[0] || !urg[1]) {
                P4EST_NOTICE ("Tree edge owner inconsistency\n");
                failed = 1;
                goto failtest;
              }
            }

            /* Then we have to loop over multiple neighbors */
            for (pz = 0; pz < procs[0].elem_count; ++pz) {
              n0_proc = *((int *) sc_array_index (&procs[0], pz));

              if (n0_proc != rank) {
                buf = p4est_ghost_array_index (&send_bufs, n0_proc);
                p4est_add_ghost_to_buf (buf, nt, local_num, q);
              }

              if (!maxed) {
                n1_proc = *((int *) sc_array_index (&procs[1], pz));

                if (n1_proc != n0_proc && n1_proc != rank) {
                  buf = p4est_ghost_array_index (&send_bufs, n1_proc);
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
                failed = 1;
                goto failtest;
              }
            }

            if (n0_proc != rank && n0_proc >= 0) {
              buf = p4est_ghost_array_index (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
            }

            if (n1_proc != n0_proc && n1_proc != rank && n1_proc >= 0) {
              buf = p4est_ghost_array_index (&send_bufs, n1_proc);
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
            }
          }
        }
        else {
          p8est_quadrant_edge_neighbor (q, edge, &n[0]);
          if (p4est_quadrant_is_inside_root (&n[0])) {
            nedge = edge ^ 3;
            touch = ((int32_t) 1 << (6 + nedge));
            p4est_ghost_test_add (p4est, q, nt, &n[0], nt, touch, rank,
                                  &send_bufs, local_num);
          }
          else if (p4est_quadrant_is_outside_face (&n[0])) {
            P4EST_ASSERT (p4est_quadrant_is_extended (&n[0]));
            face = -1;
            if (n[0].x < 0 || n[0].x >= P4EST_ROOT_LEN) {
              face = p8est_edge_faces[edge][0];
            }
            else if (n[0].z < 0 || n[0].z >= P4EST_ROOT_LEN) {
              face = p8est_edge_faces[edge][1];
            }
            else if (n[0].y < 0) {
              face = 2;
            }
            else {
              face = 3;
            }
            nnt = p4est_find_face_transform (conn, nt, face, ftransform);
            if (nnt < 0) {
              continue;
            }
            P4EST_ASSERT (face >= 0);
            P4EST_ASSERT (p8est_edge_face_corners[edge][face][0] != -1);
            if (p8est_edge_faces[edge][0] == face) {
              oppedge = edge ^ 2;
              P4EST_ASSERT (p8est_edge_faces[oppedge][0] == face);
            }
            else {
              oppedge = edge ^ 1;
              P4EST_ASSERT (p8est_edge_faces[oppedge][1] == face);
            }
            nface = (int) conn->tree_to_face[nt * P4EST_FACES + face];
            o = nface / P4EST_FACES;
            nface %= P4EST_FACES;
            ref = p8est_face_permutation_refs[face][nface];
            set = p8est_face_permutation_sets[ref][o];
            c0 = p8est_edge_face_corners[oppedge][face][0];
            c1 = p8est_edge_face_corners[oppedge][face][1];
            nc0 = p8est_face_permutations[set][c0];
            nc1 = p8est_face_permutations[set][c1];
            nc0 = p8est_face_corners[nface][nc0];
            nc1 = p8est_face_corners[nface][nc1];
            nedge = p8est_child_corner_edges[nc0][nc1];
            touch = ((int32_t) 1 << (6 + nedge));
            p4est_quadrant_transform_face (&n[0], &n[1], ftransform);
            p4est_ghost_test_add (p4est, q, nt, &n[1], nnt, touch, rank,
                                  &send_bufs, local_num);
          }
          else {
            P4EST_ASSERT (p8est_quadrant_is_outside_edge (&n[0]));
            sc_array_init (eta, sizeof (p8est_edge_transform_t));
            p8est_find_edge_transform (conn, nt, edge, &ei);
            for (etree = 0; etree < eta->elem_count; etree++) {
              et = p8est_edge_array_index (eta, etree);
              p8est_quadrant_transform_edge (&n[0], &n[1], &ei, et, 1);
              nnt = et->ntree;
              nedge = (int) et->nedge;
              touch = ((int32_t) 1 << (6 + nedge));
              p4est_ghost_test_add (p4est, q, nt, &n[1], nnt, touch, rank,
                                    &send_bufs, local_num);
            }
            sc_array_reset (eta);
          }
        }
      }

      if (btype == P8EST_BALANCE_EDGE)
        continue;
#endif

      /* Find smaller corner neighbors */
      for (corner = 0; corner < P4EST_CHILDREN; ++corner) {
        if (tol < P4EST_GHOST_UNBALANCED_ALLOW) {
          if (q->level == P4EST_QMAXLEVEL) {
            p4est_quadrant_corner_neighbor (q, corner, &n[0]);
            maxed = 1;
          }
          else {
            p4est_quadrant_get_half_corner_neighbor (q, corner, &n[0],
                                                     &nur[0]);
            maxed = 0;
          }

          /* Check to see if we are a tree corner neighbor */
          if (p4est_quadrant_is_outside_corner (&n[0])) {
            /* Then we have to loop over multiple corner neighbors */
            p4est_quadrant_find_tree_corner_owners (p4est, nt, corner, &n[0],
                                                    &procs[0], &urg[0]);
            if (!urg[0]) {
              P4EST_NOTICE ("Tree corner owner inconsistency\n");
              failed = 1;
              goto failtest;
            }

            for (pz = 0; pz < procs[0].elem_count; ++pz) {
              n0_proc = *((int *) sc_array_index (&procs[0], pz));

              if (n0_proc != rank) {
                buf = p4est_ghost_array_index (&send_bufs, n0_proc);
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
              failed = 1;
              goto failtest;
            }

            /* Then we have to loop over multiple edge neighbors */
            for (pz = 0; pz < procs[0].elem_count; ++pz) {
              n0_proc = *((int *) sc_array_index (&procs[0], pz));

              if (n0_proc != rank) {
                buf = p4est_ghost_array_index (&send_bufs, n0_proc);
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
                failed = 1;
                goto failtest;
              }
            }

            if (n0_proc != rank && n0_proc >= 0) {
              buf = p4est_ghost_array_index (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
            }
          }
        }
        else {
          p4est_quadrant_corner_descendant (q, &n[1], corner,
                                            P4EST_QMAXLEVEL);
          p4est_quadrant_corner_neighbor (&n[1], corner, &n[0]);
          if (p4est_quadrant_is_inside_root (&n[0])) {
            n0_proc = p4est_comm_find_owner (p4est, nt, &n[0], rank);
            P4EST_ASSERT (n0_proc >= 0);
            if (n0_proc != rank) {
              buf = p4est_ghost_array_index (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
            }
          }
          else if (p4est_quadrant_is_outside_face (&n[0])) {
            if (n[0].x < 0 || n[0].x >= P4EST_ROOT_LEN) {
              face = p4est_corner_faces[corner][0];
            }
#ifdef P4_TO_P8
            else if (n[0].y < 0 || n[0].y >= P4EST_ROOT_LEN) {
              face = p4est_corner_faces[corner][1];
            }
#endif
            else {
              face = p4est_corner_faces[corner][P4EST_DIM - 1];
            }
            nnt = p4est_find_face_transform (conn, nt, face, ftransform);
            if (nnt < 0) {
              continue;
            }
            p4est_quadrant_transform_face (&n[0], &n[1], ftransform);
            n0_proc = p4est_comm_find_owner (p4est, nnt, &n[1], rank);
            if (n0_proc != rank) {
              buf = p4est_ghost_array_index (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
            }
          }
#ifdef P4_TO_P8
          else if (p8est_quadrant_is_outside_edge_extra (&n[0], &edge)) {
            sc_array_init (eta, sizeof (p8est_edge_transform_t));
            p8est_find_edge_transform (conn, nt, edge, &ei);
            for (etree = 0; etree < eta->elem_count; etree++) {
              et = p8est_edge_array_index (eta, etree);
              p8est_quadrant_transform_edge (&n[0], &n[1], &ei, et, 1);
              nnt = et->ntree;
              n0_proc = p4est_comm_find_owner (p4est, nnt, &n[1], rank);
              if (n0_proc != rank) {
                buf = p4est_ghost_array_index (&send_bufs, n0_proc);
                p4est_add_ghost_to_buf (buf, nt, local_num, q);
              }
            }
            sc_array_reset (eta);
          }
#endif
          else {
            sc_array_init (cta, sizeof (p4est_corner_transform_t));
            p4est_find_corner_transform (conn, nt, corner, &ci);
            for (ctree = 0; ctree < cta->elem_count; ++ctree) {
              ct = p4est_corner_array_index (cta, ctree);
              p4est_quadrant_transform_corner (&n[0], (int) ct->ncorner, 1);
              nnt = ct->ntree;
              n0_proc = p4est_comm_find_owner (p4est, nnt, &n[0], rank);
              if (n0_proc != rank) {
                buf = p4est_ghost_array_index (&send_bufs, n0_proc);
                p4est_add_ghost_to_buf (buf, nt, local_num, q);
              }
            }
            sc_array_reset (cta);
          }
        }
      }

      P4EST_ASSERT (btype == P4EST_BALANCE_FULL);
    }
  }
  P4EST_ASSERT (local_num == p4est->local_num_quadrants);

failtest:
  if (tol == P4EST_GHOST_UNBALANCED_FAIL) {
    if (p4est_comm_sync_flag (p4est, failed, MPI_BOR)) {
      for (i = 0; i < num_procs; ++i) {
        buf = p4est_ghost_array_index (&send_bufs, i);
        sc_array_reset (buf);
      }
      sc_array_reset (&send_bufs);
      for (i = 0; i < P4EST_DIM - 1; ++i) {
        sc_array_reset (&procs[i]);
      }
      p4est_ghost_destroy (gl);

      return NULL;
    }
  }
  else if (tol == P4EST_GHOST_UNBALANCED_ABORT) {
    SC_CHECK_ABORT (!failed, "Ghost layer");
  }

  /* Count the number of peers that I send to and receive from */
  for (i = 0, num_peers = 0; i < num_procs; ++i) {
    buf = p4est_ghost_array_index (&send_bufs, i);
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
    buf = p4est_ghost_array_index (&send_bufs, i);
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
    buf = p4est_ghost_array_index (&send_bufs, i);
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
    buf = p4est_ghost_array_index (&send_bufs, i);
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
    /* proc_offsets[0] is set at beginning of this function */
    gl->proc_offsets[i + 1] = ghost_offset;
  }
  P4EST_ASSERT (ghost_offset == num_ghosts);

  /* Send the ghosts */
  for (i = 0, peer = 0; i < num_procs; ++i) {
    buf = p4est_ghost_array_index (&send_bufs, i);
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
    q = p4est_quadrant_array_index (ghost_layer, (size_t) li);
    P4EST_ASSERT (p4est_quadrant_is_valid (q));
    P4EST_ASSERT (q->p.piggy1.which_tree >= 0 &&
                  q->p.piggy1.which_tree < num_trees);
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
    buf = p4est_ghost_array_index (&send_bufs, i);
    sc_array_reset (buf);
  }
  sc_array_reset (&send_bufs);
  for (i = 0; i < P4EST_DIM - 1; ++i) {
    sc_array_reset (&procs[i]);
  }
#endif /* P4EST_MPI */

  /* calculate tree offsets */
  sc_array_init (&split, sizeof (size_t));
  sc_array_split (ghost_layer, &split,
                  (size_t) num_trees, ghost_tree_type, NULL);
  P4EST_ASSERT (split.elem_count == (size_t) num_trees + 1);
  for (nt = 0; nt <= num_trees; ++nt) {
    ppz = (size_t *) sc_array_index (&split, (size_t) nt);
    gl->tree_offsets[nt] = *ppz;
#ifdef P4EST_DEBUG
    if (nt > 0) {
      p4est_locidx_t      lk;
      p4est_quadrant_t   *q3;

      for (lk = gl->tree_offsets[nt - 1]; lk < gl->tree_offsets[nt]; ++lk) {
        q3 = p4est_quadrant_array_index (ghost_layer, (size_t) lk);
        SC_CHECK_ABORT (q3->p.which_tree == nt - 1, "Ghost tree offset");
      }
    }
#endif
  }
  sc_array_reset (&split);
  P4EST_ASSERT (gl->tree_offsets[0] == 0);
  P4EST_ASSERT (gl->proc_offsets[0] == 0);

  P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_ghost_new\n");
  return gl;
}

p4est_ghost_t      *
p4est_ghost_new (p4est_t * p4est, p4est_balance_type_t btype)
{
  return p4est_ghost_new_check (p4est, btype, P4EST_GHOST_UNBALANCED_ALLOW);
}

void
p4est_ghost_destroy (p4est_ghost_t * ghost)
{
  sc_array_reset (&ghost->ghosts);

  P4EST_FREE (ghost->tree_offsets);
  P4EST_FREE (ghost->proc_offsets);

  P4EST_FREE (ghost);
}

unsigned
p4est_ghost_checksum (p4est_t * p4est, p4est_ghost_t * ghost)
{
  unsigned            crc;
  uint32_t           *check;
  size_t              zz, csize, qcount, offset;
  size_t              nt1, np1, local_count;
  sc_array_t         *quadrants, *checkarray;
  p4est_quadrant_t   *q;

  quadrants = &ghost->ghosts;
  qcount = quadrants->elem_count;
  nt1 = (size_t) p4est->connectivity->num_trees + 1;
  np1 = (size_t) p4est->mpisize + 1;

  P4EST_ASSERT (quadrants->elem_size == sizeof (p4est_quadrant_t));

  csize = sizeof (uint32_t);
  checkarray = sc_array_new (csize);

  local_count = qcount * (P4EST_DIM + 3) + nt1 + np1;
  sc_array_resize (checkarray, local_count);

  /* checksum ghost quadrants */
  for (zz = 0; zz < qcount; ++zz) {
    q = p4est_quadrant_array_index (quadrants, zz);
    P4EST_ASSERT (p4est_quadrant_is_valid (q));
    check = (uint32_t *) sc_array_index (checkarray, zz * (P4EST_DIM + 3));
    check[0] = htonl ((uint32_t) q->x);
    check[1] = htonl ((uint32_t) q->y);
#ifdef P4_TO_P8
    check[2] = htonl ((uint32_t) q->z);
#endif
    check[P4EST_DIM] = htonl ((uint32_t) q->level);
    check[P4EST_DIM + 1] = htonl ((uint32_t) q->p.piggy3.which_tree);
    check[P4EST_DIM + 2] = htonl ((uint32_t) q->p.piggy3.local_num);
  }

  /* checksum tree_offsets */
  offset = qcount * (P4EST_DIM + 3);
  for (zz = 0; zz < nt1; ++zz) {
    check = (uint32_t *) sc_array_index (checkarray, offset + zz);
    *check = htonl ((uint32_t) ghost->tree_offsets[zz]);
  }

  /* checksum proc_offsets */
  offset += nt1;
  for (zz = 0; zz < np1; ++zz) {
    check = (uint32_t *) sc_array_index (checkarray, offset + zz);
    *check = htonl ((uint32_t) ghost->proc_offsets[zz]);
  }
  P4EST_ASSERT (offset + zz == local_count);

  /* compute parallel checksum */
  crc = sc_array_checksum (checkarray);
  sc_array_destroy (checkarray);

  return p4est_comm_checksum (p4est, crc, csize * local_count);
}
