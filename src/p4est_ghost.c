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
#include <p4est_communication.h>
#include <p4est_ghost.h>
#include <p4est_search.h>
#include <p4est_lnodes.h>
#include <p4est_algorithms.h>
#else
/* bits and communication are included in p8est_ghost.c */
#include <p8est_ghost.h>
#include <p8est_search.h>
#include <p8est_lnodes.h>
#include <p8est_algorithms.h>
#endif
#include <sc_search.h>

/* htonl is in either of these three */
#ifdef P4EST_HAVE_ARPA_NET_H
#include <arpa/inet.h>
#endif
#ifdef P4EST_HAVE_NETINET_IN_H
#include <netinet/in.h>
#endif
#if defined P4EST_HAVE_WINSOCK2_H || defined _WIN32
#include <winsock2.h>
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

#ifdef P4EST_ENABLE_MPI

static inline sc_array_t *
p4est_ghost_array_index (sc_array_t * array, int i)
{
  return (sc_array_t *) sc_array_index_int (array, i);
}

#endif

p4est_ghost_t      *
p4est_ghost_new_local (p4est_t * p4est, p4est_connect_type_t ctype)
{
  p4est_ghost_t      *ghost;
  p4est_topidx_t      ntpo;
  int                 Ppo;

  /* assert validity of input parameters */
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_ASSERT (P4EST_CONNECT_SELF <= ctype && ctype <= P4EST_CONNECT_FULL);

  /* leave mirror_proc_mirrors and mirror_proc_front* at NULL */
  ghost = P4EST_ALLOC_ZERO (p4est_ghost_t, 1);

  /* ghost meta information */
  Ppo = (ghost->mpisize = p4est->mpisize) + 1;
  ntpo = (ghost->num_trees = p4est->connectivity->num_trees) + 1;
  ghost->btype = ctype;

  /* the ghost and mirror quadrants themselves */
  sc_array_init (&ghost->ghosts, sizeof (p4est_quadrant_t));
  sc_array_init (&ghost->mirrors, sizeof (p4est_quadrant_t));

  /* offsets into ghosts and mirrors grouped by tree */
  ghost->tree_offsets = P4EST_ALLOC_ZERO (p4est_locidx_t, ntpo);
  ghost->mirror_tree_offsets = P4EST_ALLOC_ZERO (p4est_locidx_t, ntpo);

  /* offsets into ghosts and mirrors grouped by process */
  ghost->proc_offsets = P4EST_ALLOC_ZERO (p4est_locidx_t, Ppo);
  ghost->mirror_proc_offsets = P4EST_ALLOC_ZERO (p4est_locidx_t, Ppo);

  /* this ghost layer is valid */
  P4EST_ASSERT (p4est_ghost_is_valid (p4est, ghost));
  return ghost;
}

static p4est_ghost_t *p4est_ghost_new_check (p4est_t * p4est,
                                             p4est_connect_type_t btype,
                                             p4est_ghost_tolerance_t tol);

int
p4est_quadrant_find_owner (p4est_t * p4est, p4est_topidx_t treeid,
                           int face, const p4est_quadrant_t * q)
{
  const int           rank = p4est->mpirank;
  p4est_connectivity_t *conn = p4est->connectivity;
#ifdef P4EST_ENABLE_DEBUG
  int                 dims;
#endif
  int                 ftransform[P4EST_FTRANSFORM];
  p4est_topidx_t      ntreeid;
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
#ifdef P4EST_ENABLE_DEBUG
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

  P4EST_EXECUTE_ASSERT_TOPIDX
    (p4est_find_face_transform (conn, treeid, face, ftransform), ntreeid);
  p4est_quadrant_transform_face (q, &nq, ftransform);

  return p4est_comm_find_owner (p4est, ntreeid, &nq, rank);
}

#ifdef P4EST_ENABLE_MPI

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

#endif /* P4EST_ENABLE_MPI */

/* Returns true for matching proc and tree range, false otherwise. */
static int
p4est_ghost_check_range (p4est_ghost_t * ghost,
                         int which_proc, p4est_topidx_t which_tree,
                         size_t *pstart, size_t *pended)
{
  size_t              start = 0;
  size_t              ended = ghost->ghosts.elem_count;

  if (ghost->ghosts.elem_count == 0) {
    *pstart = *pended = 0;
    return 0;
  }

  if (which_proc != -1) {
    P4EST_ASSERT (0 <= which_proc && which_proc < ghost->mpisize);
    start = SC_MAX (start, (size_t) ghost->proc_offsets[which_proc]);
    ended = SC_MIN (ended, (size_t) ghost->proc_offsets[which_proc + 1]);
  }
  if (which_tree != -1) {
    P4EST_ASSERT (0 <= which_tree && which_tree < ghost->num_trees);
    start = SC_MAX (start, (size_t) ghost->tree_offsets[which_tree]);
    ended = SC_MIN (ended, (size_t) ghost->tree_offsets[which_tree + 1]);
  }

  *pstart = start;
  *pended = ended;
  return start < ended;
}

ssize_t
p4est_ghost_bsearch (p4est_ghost_t * ghost,
                     int which_proc, p4est_topidx_t which_tree,
                     const p4est_quadrant_t * q)
{
  size_t              start, ended;

  if (p4est_ghost_check_range (ghost, which_proc, which_tree, &start, &ended)) {
    ssize_t             result;
    sc_array_t          ghost_view;

    /* create a per-tree window on the ghost layer */
    sc_array_init_view (&ghost_view, &ghost->ghosts, start, ended - start);
    result = sc_array_bsearch (&ghost_view, q, p4est_quadrant_compare);

    /* and don't forget to add the window offset */
    return (result < 0) ? (ssize_t) (-1) : result + (ssize_t) start;
  }
  else {
    P4EST_ASSERT (p4est_quadrant_is_valid (q));
    return -1;
  }
}

ssize_t
p4est_ghost_contains (p4est_ghost_t * ghost,
                      int which_proc, p4est_topidx_t which_tree,
                      const p4est_quadrant_t * q)
{
  size_t              start, ended;

  if (p4est_ghost_check_range (ghost, which_proc, which_tree, &start, &ended)) {
    size_t              nmemb = ended - start - 1;
    size_t              result;
    sc_array_t          ghost_view;
    p4est_quadrant_t   *qresult;

    /* create a per-tree window on the ghost layer */
    sc_array_init_view (&ghost_view, &ghost->ghosts, start, ended - start);
    result = sc_bsearch_range (q, ghost_view.array,
                               nmemb, sizeof (p4est_quadrant_t),
                               p4est_quadrant_compare);
    qresult = p4est_quadrant_array_index (&ghost_view, result);

    /* and don't forget to add the window offset */
    return !(p4est_quadrant_is_equal (qresult, q) ||
             p4est_quadrant_is_ancestor (qresult, q)) ?
      (ssize_t) (-1) : (ssize_t) (result + start);
  }
  else {
    P4EST_ASSERT (p4est_quadrant_is_valid (q));
    return -1;
  }
}

int
p4est_quadrant_exists (p4est_t * p4est, p4est_ghost_t * ghost,
                       p4est_topidx_t treeid, const p4est_quadrant_t * q,
                       sc_array_t * exists_arr,
                       sc_array_t * rproc_arr, sc_array_t * rquad_arr)
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
  p4est_quadrant_t    tq, non_existent, *rquad;
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
  if (rproc_arr != NULL) {
    P4EST_ASSERT (rproc_arr->elem_size == sizeof (int));
    sc_array_resize (rproc_arr, 0);
  }
  if (rquad_arr != NULL) {
    P4EST_ASSERT (rquad_arr->elem_size == sizeof (p4est_quadrant_t));
    sc_array_resize (rquad_arr, 0);
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
      lnid = p4est_ghost_bsearch (ghost, qproc, treeid, q);
      P4EST_ASSERT (lnid == -1 ||
                    (ghost->proc_offsets[qproc] <= lnid &&
                     lnid < ghost->proc_offsets[qproc + 1]));
    }
    if (rproc_arr != NULL) {
      *(int *) sc_array_push (rproc_arr) = qproc;
    }
    if (rquad_arr != NULL) {
      rquad = p4est_quadrant_array_push_copy (rquad_arr, q);
      rquad->p.piggy3.which_tree = treeid;
      rquad->p.piggy3.local_num = (p4est_locidx_t) lnid;
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
        lnid = p4est_ghost_bsearch (ghost, qproc, tqtreeid, &tq);
        P4EST_ASSERT (lnid == -1 ||
                      (ghost->proc_offsets[qproc] <= lnid &&
                       lnid < ghost->proc_offsets[qproc + 1]));
      }
      if (rproc_arr != NULL) {
        *(int *) sc_array_push (rproc_arr) = qproc;
      }
      if (rquad_arr != NULL) {
        rquad = p4est_quadrant_array_push_copy (rquad_arr, &tq);
        rquad->p.piggy3.which_tree = tqtreeid;
        rquad->p.piggy3.local_num = (p4est_locidx_t) lnid;
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
      lnid = p4est_ghost_bsearch (ghost, qproc, tqtreeid, &tq);
      P4EST_ASSERT (lnid == -1 ||
                    (ghost->proc_offsets[qproc] <= lnid &&
                     lnid < ghost->proc_offsets[qproc + 1]));
    }
    if (rproc_arr != NULL) {
      *(int *) sc_array_push (rproc_arr) = qproc;
    }
    if (rquad_arr != NULL) {
      rquad = p4est_quadrant_array_push_copy (rquad_arr, &tq);
      rquad->p.piggy3.which_tree = tqtreeid;
      rquad->p.piggy3.local_num = (p4est_locidx_t) lnid;
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
  p4est_topidx_t      tqtreeid;
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
      lnid = p4est_ghost_bsearch (ghost, qproc, treeid, q);
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
  P4EST_EXECUTE_ASSERT_TOPIDX
    (p4est_find_face_transform (conn, treeid, face, ftransform), tqtreeid);
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
    lnid = p4est_ghost_bsearch (ghost, qproc, tqtreeid, &tq);
    return (lnid == -1) ? (p4est_locidx_t) (-1) :
      (q =
       p4est_quadrant_array_index (&ghost->ghosts, (size_t) lnid),
       q->p.piggy3.local_num);
  }
}

/** Get the smallest corner neighbor of \a q.
 *
 * Gets the smallest corner neighbor, which is half of the size assuming the
 * 2-1 constraint.
 *
 * \param [in]  q      The quadrant whose corner neighbor will be constructed.
 * \param [in]  corner The corner across which to generate the neighbor.
 * \param [out] n0     Filled with the smallest corner neighbor, which is
 *                     half of the size assuming the 2-1 constraint.
 * \param [out] n0ur   If not NULL, it is filled with smallest quadrant
 *                     that fits in the upper right corner of \a n0.
 */
static void
p4est_quadrant_get_half_corner_neighbor (const p4est_quadrant_t * q,
                                         int corner,
                                         p4est_quadrant_t * n0,
                                         p4est_quadrant_t * n0ur)
{
  p4est_quadrant_half_corner_neighbor (q, corner, n0);
  if (n0ur != NULL) {
    const p4est_qcoord_t qh_2 = P4EST_QUADRANT_LEN (q->level + 1);
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

int
p4est_is_balanced (p4est_t * p4est, p4est_connect_type_t btype)
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
#ifdef P4EST_ENABLE_DEBUG
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
        e0 = p4est_quadrant_exists (p4est, gl, nt, &n[0], NULL, NULL, NULL);
        e1 = p4est_quadrant_exists (p4est, gl, nt, &n[1], NULL, NULL, NULL);
#ifndef P4_TO_P8
        e0b = e1b = e0;
#else
        e0b = p4est_quadrant_exists (p4est, gl, nt, &n[2], NULL, NULL, NULL);
        e1b = p4est_quadrant_exists (p4est, gl, nt, &n[3], NULL, NULL, NULL);
#endif
        if (e0 != e1 || e0 != e0b || e0 != e1b) {
          P4EST_NOTICE ("Contradicting small face neighbors\n");
          failed = 1;
          goto failtest;
        }
        e2 = p4est_quadrant_exists (p4est, gl, nt, &n[P4EST_HALF],
                                    NULL, NULL, NULL);
        e3 = p4est_quadrant_exists (p4est, gl, nt, &n[P4EST_HALF + 1],
                                    NULL, NULL, NULL);
        if ((int) e0 + (int) e2 + (int) e3 != 1) {
          P4EST_NOTICE ("Face balance quadrant mismatch\n");
          failed = 1;
          goto failtest;
        }
        bigger_face[face] = e3;
      }

      if (btype == P4EST_CONNECT_FACE)
        continue;

#ifdef P4_TO_P8

      /* Find edge neighbors */
      for (edge = 0; edge < P8EST_EDGES; ++edge) {
        bigger_edge[edge] = 0;
#ifdef P4EST_ENABLE_DEBUG
        big_count[edge] = 0;
#endif

        /* If q is at a boundary then it is automatically balanced */
        if (p8est_quadrant_on_edge_boundary (p4est, nt, edge, q)) {
          continue;
        }

        /* Do more expensive edge balance checks */
        p8est_quadrant_get_possible_edge_neighbors (q, edge, n);
        e0 = p4est_quadrant_exists (p4est, gl, nt, &n[0], &e0_a, NULL, NULL);
        e1 = p4est_quadrant_exists (p4est, gl, nt, &n[1], &e1_a, NULL, NULL);
        if (e0 != e1 || e0_a.elem_count != e1_a.elem_count) {
          P4EST_NOTICE ("Contradicting small edge neighbors\n");
          failed = 1;
          goto failtest;
        }
        e2 = p4est_quadrant_exists (p4est, gl, nt, &n[2], &e2_a, NULL, NULL);
        e3 = p4est_quadrant_exists (p4est, gl, nt, &n[3], &e3_a, NULL, NULL);
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
#ifdef P4EST_ENABLE_DEBUG
        big_count[edge] = e3_a.elem_count;
#endif
      }

      if (btype == P8EST_CONNECT_EDGE)
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
        e0 = p4est_quadrant_exists (p4est, gl, nt, &n[0], &e0_a, NULL, NULL);
        e1 = p4est_quadrant_exists (p4est, gl, nt, &n[1], &e1_a, NULL, NULL);
        e2 = p4est_quadrant_exists (p4est, gl, nt, &n[2], &e2_a, NULL, NULL);
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
            e3 =
              p4est_quadrant_exists (p4est, gl, nt, &n[3], &e3_a, NULL, NULL);
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

      P4EST_ASSERT (btype == P4EST_CONNECT_FULL);
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

  return !p4est_comm_sync_flag (p4est, failed, sc_MPI_BOR);
}

static              size_t
ghost_tree_type (sc_array_t * array, size_t zindex, void *data)
{
  p4est_quadrant_t   *q;

  P4EST_ASSERT (array->elem_size == sizeof (p4est_quadrant_t));

  q = (p4est_quadrant_t *) sc_array_index (array, zindex);
  return (size_t) q->p.which_tree;
}

#ifdef P4EST_ENABLE_MPI

static              size_t
ghost_proc_type (sc_array_t * array, size_t zindex, void *data)
{
  p4est_quadrant_t   *q;
  p4est_t            *p4est = (p4est_t *) data;
  int                 proc;

  P4EST_ASSERT (array->elem_size == sizeof (p4est_quadrant_t));

  q = (p4est_quadrant_t *) sc_array_index (array, zindex);
  proc = p4est_comm_find_owner (p4est, q->p.which_tree, q, 0);
  P4EST_ASSERT (proc >= 0 && proc < p4est->mpisize);
  P4EST_ASSERT (p4est_comm_is_owner (p4est, q->p.which_tree, q, proc));
  return (size_t) proc;
}

/** This adds a quadrant to the end of a buffer.
 *
 * It crams the tree id into the user_data field of the quadrant in
 * the buffer and only adds the quadrant to the end of the buffer if
 * it is unique.
 *
 * \param [in,out] buf    \a q is added to the end if it is not already there.
 * \param [in,out] q      the quadrant to be added.  The \c user_data field
 *                        is filled with \a treeid.
 * \param [in]            treeid the tree id of \a q.
 * \return                true if the ghost was added, false if duplicate.
 */
static int
p4est_add_ghost_to_buf (sc_array_t * buf, p4est_topidx_t treeid,
                        p4est_locidx_t number, const p4est_quadrant_t * q)
{
  p4est_quadrant_t   *qold, *qnew;

  P4EST_ASSERT (treeid >= 0 && number >= 0);

  /* Check to see if the quadrant already is last in the array */
  if (buf->elem_count > 0) {
    qold = p4est_quadrant_array_index (buf, buf->elem_count - 1);
    if (treeid == qold->p.piggy3.which_tree &&
        p4est_quadrant_is_equal (q, qold)) {
      return 0;
    }
  }

  qnew = p4est_quadrant_array_push_copy (buf, q);

  /* Cram the tree id and the local number into the user_data pointer */
  qnew->p.piggy3.which_tree = treeid;
  qnew->p.piggy3.local_num = number;

  return 1;
}

/** Data structure that contains temporary mirror information */
typedef struct p4est_ghost_mirror
{
  int                 mpisize, mpirank;
  int                 known;    /* was this mirror added before? */
  p4est_locidx_t      sum_all_procs;    /* sum of mirrors by processor */
  sc_array_t         *send_bufs;        /* lives in p4est_ghost_new_check */
  sc_array_t         *mirrors;  /* lives in p4est_ghost_t */
  sc_array_t         *offsets_by_proc;  /* a p4est_locidx_t array per proc */
}
p4est_ghost_mirror_t;

/** Initialize temporary mirror storage */
static void
p4est_ghost_mirror_init (p4est_ghost_t * ghost, int mpirank,
                         sc_array_t * send_bufs, p4est_ghost_mirror_t * m)
{
  int                 p;

  m->mpisize = ghost->mpisize;
  m->mpirank = mpirank;
  /* m->known is left undefined: it needs to be set to 0 for every quadrant */
  m->sum_all_procs = 0;

  m->send_bufs = send_bufs;
  P4EST_ASSERT (m->send_bufs->elem_size == sizeof (sc_array_t));
  P4EST_ASSERT (m->send_bufs->elem_count == (size_t) m->mpisize);

  m->mirrors = &ghost->mirrors;
  P4EST_ASSERT (m->mirrors->elem_size == sizeof (p4est_quadrant_t));
  P4EST_ASSERT (m->mirrors->elem_count == 0);

  m->offsets_by_proc = P4EST_ALLOC (sc_array_t, ghost->mpisize);
  for (p = 0; p < ghost->mpisize; ++p) {
    sc_array_init (m->offsets_by_proc + p, sizeof (p4est_locidx_t));
  }
}

/** Potentially record a quadrant that is to be sent as a mirror
 * \param [in] m      The temporary data structure to work on.
 * \param [in] treeid The tree number looped through by the current rank.
 * \param [in] q      The quadrant currently looked at by current rank.
 * \param [in] p      The rank that \a q should be sent to.
 */
static void
p4est_ghost_mirror_add (p4est_ghost_mirror_t * m, p4est_topidx_t treeid,
                        p4est_locidx_t number, p4est_quadrant_t * q, int p)
{
  sc_array_t         *buf;
  p4est_locidx_t     *num;
  p4est_quadrant_t   *qnew;

  P4EST_ASSERT (p != m->mpirank);
  P4EST_ASSERT (0 <= p && p < m->mpisize);

  if (!m->known) {
    /* add this quadrant to the mirror array */
    qnew = p4est_quadrant_array_push_copy (m->mirrors, q);

    /* cram the tree id and the local number into the user_data pointer */
    qnew->p.piggy3.which_tree = treeid;
    qnew->p.piggy3.local_num = number;

    m->known = 1;
  }

  buf = p4est_ghost_array_index (m->send_bufs, p);
  if (p4est_add_ghost_to_buf (buf, treeid, number, q)) {
    P4EST_ASSERT (m->mirrors->elem_count > 0);

    num = (p4est_locidx_t *) sc_array_push (m->offsets_by_proc + p);
    *num = (p4est_locidx_t) (m->mirrors->elem_count - 1);
    ++m->sum_all_procs;
  }
}

/** Populate the mirror fields in the ghost layer with final data.
 * The elements in the temporary p4est_ghost_mirror_t structure are freed. */
static void
p4est_ghost_mirror_reset (p4est_ghost_t * ghost, p4est_ghost_mirror_t * m,
                          int populate)
{
  int                 p;
  p4est_locidx_t     *mpm;
  p4est_locidx_t      pcount, sum_all_procs = 0;

  P4EST_ASSERT (ghost->mirror_proc_mirrors == NULL);

  /* if we did not run into failtest, populate the mirrors */
  if (populate) {
    mpm = ghost->mirror_proc_mirrors =
      P4EST_ALLOC (p4est_locidx_t, m->sum_all_procs);
    for (p = 0; p < ghost->mpisize; ++p) {
      pcount = (p4est_locidx_t) m->offsets_by_proc[p].elem_count;
      P4EST_ASSERT (p != m->mpirank || pcount == 0);
      memcpy (mpm + sum_all_procs, m->offsets_by_proc[p].array,
              pcount * sizeof (p4est_locidx_t));
      ghost->mirror_proc_offsets[p] = sum_all_procs;
      sum_all_procs += pcount;
    }
    P4EST_ASSERT (sum_all_procs == m->sum_all_procs);
    ghost->mirror_proc_offsets[p] = sum_all_procs;
  }

  /* clean up memory regardless */
  for (p = 0; p < ghost->mpisize; ++p) {
    sc_array_reset (m->offsets_by_proc + p);
  }
  P4EST_FREE (m->offsets_by_proc);
  memset (m, 0, sizeof (p4est_ghost_mirror_t));
}

static void
p4est_ghost_test_add (p4est_t * p4est, p4est_ghost_mirror_t * m,
                      p4est_quadrant_t * q, p4est_topidx_t t,
                      p4est_quadrant_t * nq, p4est_topidx_t nt,
                      int32_t touch, int rank,
#if 0
                      sc_array_t * send_bufs,
#endif
                      p4est_locidx_t local_num)
{
  p4est_quadrant_t    temp;
  p4est_quadrant_t   *lq, *uq;
#ifdef P4EST_ENABLE_DEBUG
  p4est_quadrant_t    debug_quad;
  p4est_lid_t         next_lid, uid, temp_lid;
#endif
  int                 n0_proc, n1_proc, proc;
  p4est_quadrant_t   *gfp = p4est->global_first_position;
#if 0
  sc_array_t         *buf;
#endif
  int32_t             rb;

  P4EST_ASSERT (q->level == nq->level);
  n0_proc = p4est_comm_find_owner (p4est, nt, nq, rank);
  P4EST_ASSERT (n0_proc >= 0);
  if (q->level == P4EST_QMAXLEVEL) {
    if (n0_proc != rank) {
#if 0
      buf = p4est_ghost_array_index (send_bufs, n0_proc);
      p4est_add_ghost_to_buf (buf, t, local_num, q);
#endif
      p4est_ghost_mirror_add (m, t, local_num, q, n0_proc);
    }
    return;
  }
  p4est_quadrant_last_descendant (nq, &temp, P4EST_QMAXLEVEL);
  n1_proc = p4est_comm_find_owner (p4est, nt, &temp, n0_proc);
  P4EST_ASSERT (n1_proc >= n0_proc);
  if (n0_proc == n1_proc) {
    if (n0_proc != rank) {
#if 0
      buf = p4est_ghost_array_index (send_bufs, n0_proc);
      p4est_add_ghost_to_buf (buf, t, local_num, q);
#endif
      p4est_ghost_mirror_add (m, t, local_num, q, n0_proc);
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

      p4est_quadrant_predecessor (uq, &temp);
      uq = &temp;
      P4EST_ASSERT (p4est_quadrant_is_valid (uq));

#ifdef P4EST_ENABLE_DEBUG
      p4est_quadrant_copy (&(gfp[proc + 1]), &debug_quad);
      p4est_quadrant_linear_id_ext128 (&debug_quad, P4EST_QMAXLEVEL,
                                       &next_lid);
      p4est_lid_set_zero (&temp_lid);
      P4EST_ASSERT (p4est_lid_compare (&next_lid, &temp_lid) > 0);

      p4est_lid_set_one (&temp_lid);
      p4est_lid_sub (&next_lid, &temp_lid, &uid);
      p4est_quadrant_set_morton_ext128 (&debug_quad, P4EST_QMAXLEVEL, &uid);
      P4EST_ASSERT (p4est_quadrant_is_valid (&debug_quad));
      P4EST_ASSERT (p4est_quadrant_is_equal (uq, &debug_quad));
#endif
    }
#ifdef P4EST_ENABLE_DEBUG
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
#if 0
      buf = p4est_ghost_array_index (send_bufs, proc);
      p4est_add_ghost_to_buf (buf, t, local_num, q);
#endif
      p4est_ghost_mirror_add (m, t, local_num, q, proc);
    }
  }
}

#endif /* P4EST_ENABLE_MPI */

static p4est_ghost_t *
p4est_ghost_new_check (p4est_t * p4est, p4est_connect_type_t btype,
                       p4est_ghost_tolerance_t tol)
{
  const p4est_topidx_t num_trees = p4est->connectivity->num_trees;
  const int           num_procs = p4est->mpisize;
#ifdef P4EST_ENABLE_MPI
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
#ifdef P4EST_ENABLE_DEBUG
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
  MPI_Request        *recv_load_request, *send_load_request;
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
#ifdef P4EST_ENABLE_DEBUG
  p4est_quadrant_t   *q2;
#endif
  int                 ftransform[P4EST_FTRANSFORM];
  int32_t             touch;
  p4est_topidx_t      nnt;
  p4est_corner_info_t ci;
  p4est_corner_transform_t *ct;
  sc_array_t         *cta;
  size_t              ctree;
  p4est_ghost_mirror_t m;
#endif
  size_t             *ppz;
  sc_array_t          split;
  sc_array_t         *ghost_layer;
  p4est_topidx_t      nt;
  p4est_ghost_t      *gl;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING "_ghost_new %s\n",
                            p4est_connect_type_string (btype));
  p4est_log_indent_push ();

  gl = P4EST_ALLOC (p4est_ghost_t, 1);
  gl->mpisize = num_procs;
  gl->num_trees = num_trees;
  gl->btype = btype;

  ghost_layer = &gl->ghosts;
  sc_array_init (ghost_layer, sizeof (p4est_quadrant_t));
  gl->tree_offsets = P4EST_ALLOC (p4est_locidx_t, num_trees + 1);
  gl->proc_offsets = P4EST_ALLOC (p4est_locidx_t, num_procs + 1);

  sc_array_init (&gl->mirrors, sizeof (p4est_quadrant_t));
  gl->mirror_tree_offsets = P4EST_ALLOC (p4est_locidx_t, num_trees + 1);
  gl->mirror_proc_mirrors = NULL;
  gl->mirror_proc_offsets = P4EST_ALLOC (p4est_locidx_t, num_procs + 1);
  gl->mirror_proc_fronts = NULL;
  gl->mirror_proc_front_offsets = NULL;

  gl->proc_offsets[0] = 0;
  gl->mirror_proc_offsets[0] = 0;
#ifndef P4EST_ENABLE_MPI
  gl->proc_offsets[1] = 0;
  gl->mirror_proc_offsets[1] = 0;
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

  /* initialize structure to keep track of mirror quadrants */
  p4est_ghost_mirror_init (gl, p4est->mpirank, &send_bufs, &m);

  /* loop over all local trees */
  local_num = 0;
  for (nt = 0; nt < first_local_tree; ++nt) {
    /* does nothing if this processor is empty */
    gl->mirror_tree_offsets[nt] = 0;
  }
  for (nt = first_local_tree; nt <= last_local_tree; ++nt) {
    /* does nothing if this processor is empty */
    tree = p4est_tree_array_index (p4est->trees, nt);
    quadrants = &tree->quadrants;
    p4est_comm_tree_info (p4est, nt, full_tree, tree_contact, NULL, NULL);
    gl->mirror_tree_offsets[nt] = (p4est_locidx_t) gl->mirrors.elem_count;

    /* Find the smaller neighboring processors of each quadrant */
    for (zz = 0; zz < quadrants->elem_count; ++local_num, ++zz) {
      q = p4est_quadrant_array_index (quadrants, zz);
      m.known = 0;

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
#if 0
              buf = p4est_ghost_array_index (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
#endif
              p4est_ghost_mirror_add (&m, nt, local_num, q, n0_proc);
              n1_proc = n0_proc;
            }
          }
        }
        else {
          p4est_quadrant_face_neighbor (q, face, &n[0]);
          if (p4est_quadrant_is_inside_root (&n[0])) {
            nface = face ^ 1;
            touch = ((int32_t) 1 << nface);
            p4est_ghost_test_add (p4est, &m, q, nt, &n[0], nt, touch, rank,
                                  local_num);
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
            p4est_ghost_test_add (p4est, &m, q, nt, &n[1], nnt, touch, rank,
                                  local_num);
          }
        }
      }

      if (btype == P4EST_CONNECT_FACE)
        continue;

#ifdef P4_TO_P8

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
#if 0
                buf = p4est_ghost_array_index (&send_bufs, n0_proc);
                p4est_add_ghost_to_buf (buf, nt, local_num, q);
#endif
                p4est_ghost_mirror_add (&m, nt, local_num, q, n0_proc);
              }

              if (!maxed) {
                n1_proc = *((int *) sc_array_index (&procs[1], pz));

                if (n1_proc != n0_proc && n1_proc != rank) {
#if 0
                  buf = p4est_ghost_array_index (&send_bufs, n1_proc);
                  p4est_add_ghost_to_buf (buf, nt, local_num, q);
#endif
                  p4est_ghost_mirror_add (&m, nt, local_num, q, n1_proc);
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
#if 0
              buf = p4est_ghost_array_index (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
#endif
              p4est_ghost_mirror_add (&m, nt, local_num, q, n0_proc);
            }

            if (n1_proc != n0_proc && n1_proc != rank && n1_proc >= 0) {
#if 0
              buf = p4est_ghost_array_index (&send_bufs, n1_proc);
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
#endif
              p4est_ghost_mirror_add (&m, nt, local_num, q, n1_proc);
            }
          }
        }
        else {
          p8est_quadrant_edge_neighbor (q, edge, &n[0]);
          if (p4est_quadrant_is_inside_root (&n[0])) {
            nedge = edge ^ 3;
            touch = ((int32_t) 1 << (6 + nedge));
            p4est_ghost_test_add (p4est, &m, q, nt, &n[0], nt, touch, rank,
                                  local_num);
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
              P4EST_ASSERT (p8est_edge_faces[edge][1] == face);
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
            p4est_ghost_test_add (p4est, &m, q, nt, &n[1], nnt, touch, rank,
                                  local_num);
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
              p4est_ghost_test_add (p4est, &m, q, nt, &n[1], nnt, touch, rank,
                                    local_num);
            }
            sc_array_reset (eta);
          }
        }
      }

      if (btype == P8EST_CONNECT_EDGE)
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
#if 0
                buf = p4est_ghost_array_index (&send_bufs, n0_proc);
                p4est_add_ghost_to_buf (buf, nt, local_num, q);
#endif
                p4est_ghost_mirror_add (&m, nt, local_num, q, n0_proc);
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
#if 0
                buf = p4est_ghost_array_index (&send_bufs, n0_proc);
                p4est_add_ghost_to_buf (buf, nt, local_num, q);
#endif
                p4est_ghost_mirror_add (&m, nt, local_num, q, n0_proc);
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
#if 0
              buf = p4est_ghost_array_index (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
#endif
              p4est_ghost_mirror_add (&m, nt, local_num, q, n0_proc);
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
#if 0
              buf = p4est_ghost_array_index (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
#endif
              p4est_ghost_mirror_add (&m, nt, local_num, q, n0_proc);
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
#if 0
              buf = p4est_ghost_array_index (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, local_num, q);
#endif
              p4est_ghost_mirror_add (&m, nt, local_num, q, n0_proc);
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
#if 0
                buf = p4est_ghost_array_index (&send_bufs, n0_proc);
                p4est_add_ghost_to_buf (buf, nt, local_num, q);
#endif
                p4est_ghost_mirror_add (&m, nt, local_num, q, n0_proc);
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
#if 0
                buf = p4est_ghost_array_index (&send_bufs, n0_proc);
                p4est_add_ghost_to_buf (buf, nt, local_num, q);
#endif
                p4est_ghost_mirror_add (&m, nt, local_num, q, n0_proc);
              }
            }
            sc_array_reset (cta);
          }
        }
      }

      P4EST_ASSERT (btype == P4EST_CONNECT_FULL);
    }
  }
  P4EST_ASSERT (local_num == p4est->local_num_quadrants);
  for (nt = SC_MAX (p4est->last_local_tree + 1, 0); nt <= num_trees; ++nt) {
    /* needs to cover all trees if this processor is empty */
    /* needs to run inclusive on num_trees */
    gl->mirror_tree_offsets[nt] = (p4est_locidx_t) gl->mirrors.elem_count;
  }

failtest:
  if (tol == P4EST_GHOST_UNBALANCED_FAIL) {
    if (p4est_comm_sync_flag (p4est, failed, MPI_BOR)) {
      p4est_ghost_mirror_reset (gl, &m, 0);

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
  send_request = P4EST_ALLOC (MPI_Request, 2 * num_peers);

  recv_counts = P4EST_ALLOC (p4est_locidx_t, 2 * num_peers);
  send_counts = recv_counts + num_peers;

  recv_load_request = recv_request + num_peers;
  send_load_request = send_request + num_peers;

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

  /* The mirrors can be assembled here since they are defined on the sender */
  p4est_ghost_mirror_reset (gl, &m, 1);

  /* Wait for the counts */
  if (num_peers > 0) {
    mpiret = sc_MPI_Waitall (num_peers, recv_request, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Waitall (num_peers, send_request, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }

#ifdef P4EST_ENABLE_DEBUG
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
    mpiret =
      sc_MPI_Waitall (num_peers, recv_load_request, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    mpiret =
      sc_MPI_Waitall (num_peers, send_load_request, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }

  /* Clean up */
  P4EST_FREE (recv_counts);

#ifdef P4EST_ENABLE_DEBUG
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
  P4EST_FREE (send_request);

  for (i = 0; i < num_procs; ++i) {
    buf = p4est_ghost_array_index (&send_bufs, i);
    sc_array_reset (buf);
  }
  sc_array_reset (&send_bufs);
  for (i = 0; i < P4EST_DIM - 1; ++i) {
    sc_array_reset (&procs[i]);
  }
#endif /* P4EST_ENABLE_MPI */

  /* calculate tree offsets */
  sc_array_init (&split, sizeof (size_t));
  sc_array_split (ghost_layer, &split,
                  (size_t) num_trees, ghost_tree_type, NULL);
  P4EST_ASSERT (split.elem_count == (size_t) num_trees + 1);
  for (nt = 0; nt <= num_trees; ++nt) {
    ppz = (size_t *) sc_array_index (&split, (size_t) nt);
    gl->tree_offsets[nt] = *ppz;
#ifdef P4EST_ENABLE_DEBUG
    if (nt > 0) {
      p4est_locidx_t      lk;
      p4est_quadrant_t   *q3;

      for (lk = gl->tree_offsets[nt - 1]; lk < gl->tree_offsets[nt]; ++lk) {
        q3 = p4est_quadrant_array_index (ghost_layer, (size_t) lk);
        SC_CHECK_ABORT (q3->p.which_tree == nt - 1, "Ghost tree offset");
      }
    }
#endif
#ifndef P4EST_ENABLE_MPI
    gl->mirror_tree_offsets[nt] = 0;
#endif
  }
  sc_array_reset (&split);
  P4EST_ASSERT (gl->tree_offsets[0] == 0);
  P4EST_ASSERT (gl->proc_offsets[0] == 0);

  gl->mirror_proc_fronts = gl->mirror_proc_mirrors;
  gl->mirror_proc_front_offsets = gl->mirror_proc_offsets;

  P4EST_ASSERT (p4est_ghost_is_valid (p4est, gl));

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_ghost_new\n");
  return gl;
}

p4est_ghost_t      *
p4est_ghost_new (p4est_t * p4est, p4est_connect_type_t btype)
{
  return p4est_ghost_new_check (p4est, btype, P4EST_GHOST_UNBALANCED_ALLOW);
}

void
p4est_ghost_destroy (p4est_ghost_t * ghost)
{
  sc_array_reset (&ghost->ghosts);
  P4EST_FREE (ghost->tree_offsets);
  P4EST_FREE (ghost->proc_offsets);

  if (ghost->mirror_proc_fronts != ghost->mirror_proc_mirrors) {
    P4EST_ASSERT (ghost->mirror_proc_front_offsets !=
                  ghost->mirror_proc_offsets);
    P4EST_FREE (ghost->mirror_proc_fronts);
    P4EST_FREE (ghost->mirror_proc_front_offsets);
  }

  sc_array_reset (&ghost->mirrors);
  P4EST_FREE (ghost->mirror_tree_offsets);
  P4EST_FREE (ghost->mirror_proc_mirrors);
  P4EST_FREE (ghost->mirror_proc_offsets);

  P4EST_FREE (ghost);
}

unsigned
p4est_ghost_checksum (p4est_t * p4est, p4est_ghost_t * ghost)
{
  unsigned            crc;
  uint32_t           *check;
#ifdef P4_TO_P8
  int                 level_difference;
#endif
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
#ifndef P4_TO_P8
    check[0] = htonl ((uint32_t) q->x);
    check[1] = htonl ((uint32_t) q->y);
#else
    if (q->level <= P4EST_OLD_QMAXLEVEL) {
      /* shift the quadrant coordinates to ensure backward compatibility */
      level_difference = P4EST_MAXLEVEL - P4EST_OLD_MAXLEVEL;
      /* *INDENT-OFF* */
      check[0] =
        htonl ((q->x < 0) ? -(((uint32_t) -q->x) >> level_difference) :
                              (((uint32_t) q->x) >> level_difference));
      check[1] =
        htonl ((q->y < 0) ? -(((uint32_t) -q->y) >> level_difference) :
                              (((uint32_t) q->y) >> level_difference));
      check[2] =
        htonl ((q->z < 0) ? -(((uint32_t) -q->z) >> level_difference) :
                              (((uint32_t) q->z) >> level_difference));
      /* *INDENT-ON* */
    }
    else {
      check[0] = htonl ((uint32_t) q->x);
      check[1] = htonl ((uint32_t) q->y);
      check[2] = htonl ((uint32_t) q->z);
    }
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

void
p4est_ghost_exchange_data (p4est_t * p4est, p4est_ghost_t * ghost,
                           void *ghost_data)
{
  p4est_ghost_exchange_data_end (p4est_ghost_exchange_data_begin
                                 (p4est, ghost, ghost_data));
}

p4est_ghost_exchange_t *
p4est_ghost_exchange_data_begin (p4est_t * p4est, p4est_ghost_t * ghost,
                                 void *ghost_data)
{
  size_t              zz;
  size_t              data_size;
#ifdef P4EST_ENABLE_DEBUG
  p4est_topidx_t      prev_tree;
#endif
  p4est_topidx_t      which_tree;
  p4est_locidx_t      which_quad;
  p4est_quadrant_t   *mirror, *q;
  p4est_tree_t       *tree;
  p4est_ghost_exchange_t *exc;
  void              **mirror_data;

  /* allocate temporary storage */
  mirror_data = P4EST_ALLOC (void *, ghost->mirrors.elem_count);

  data_size = p4est->data_size == 0 ? sizeof (void *) : p4est->data_size;
#ifdef P4EST_ENABLE_DEBUG
  prev_tree = -1;
#endif
  for (zz = 0; zz < ghost->mirrors.elem_count; ++zz) {
    mirror = p4est_quadrant_array_index (&ghost->mirrors, zz);
    which_tree = mirror->p.piggy3.which_tree;
    P4EST_ASSERT (p4est->first_local_tree <= which_tree &&
                  which_tree <= p4est->last_local_tree);
    P4EST_ASSERT (prev_tree <= which_tree);
#ifdef P4EST_ENABLE_DEBUG
    prev_tree = which_tree;
#endif
    tree = p4est_tree_array_index (p4est->trees, which_tree);
    which_quad = mirror->p.piggy3.local_num - tree->quadrants_offset;
    P4EST_ASSERT (0 <= which_quad &&
                  which_quad < (p4est_locidx_t) tree->quadrants.elem_count);
    q = p4est_quadrant_array_index (&tree->quadrants, which_quad);
    mirror_data[zz] =
      p4est->data_size == 0 ? &q->p.user_data : q->p.user_data;
  }

  /* delegate the rest of the work */
  exc = p4est_ghost_exchange_custom_begin (p4est, ghost, data_size,
                                           mirror_data, ghost_data);
  P4EST_ASSERT (exc->is_custom);
  P4EST_ASSERT (!exc->is_levels);
  exc->is_custom = 0;

  /* the mirror_data is copied before sending so it can be freed */
  P4EST_FREE (mirror_data);

  /* return message buffers */
  return exc;
}

void
p4est_ghost_exchange_data_end (p4est_ghost_exchange_t * exc)
{
  /* don't confuse this function with p4est_ghost_exchange_custom_end */
  P4EST_ASSERT (!exc->is_custom);
  P4EST_ASSERT (!exc->is_levels);

  /* delegate the rest of the work, including freeing the context */
  exc->is_custom = 1;
  p4est_ghost_exchange_custom_end (exc);
}

void
p4est_ghost_exchange_custom (p4est_t * p4est, p4est_ghost_t * ghost,
                             size_t data_size,
                             void **mirror_data, void *ghost_data)
{
  p4est_ghost_exchange_custom_end (p4est_ghost_exchange_custom_begin
                                   (p4est, ghost, data_size,
                                    mirror_data, ghost_data));
}

p4est_ghost_exchange_t *
p4est_ghost_exchange_custom_begin (p4est_t * p4est, p4est_ghost_t * ghost,
                                   size_t data_size,
                                   void **mirror_data, void *ghost_data)
{
  const int           num_procs = p4est->mpisize;
  int                 mpiret;
  int                 q;
  char               *mem, **sbuf;
  p4est_locidx_t      ng_excl, ng_incl, ng, theg;
  p4est_locidx_t      mirr;
  p4est_ghost_exchange_t *exc;
  sc_MPI_Request     *r;

  /* initialize transient storage */
  exc = P4EST_ALLOC_ZERO (p4est_ghost_exchange_t, 1);
  exc->is_custom = 1;
  exc->p4est = p4est;
  exc->ghost = ghost;
  exc->minlevel = 0;
  exc->maxlevel = P4EST_QMAXLEVEL;
  exc->data_size = data_size;
  exc->ghost_data = ghost_data;
  sc_array_init (&exc->requests, sizeof (sc_MPI_Request));
  sc_array_init (&exc->sbuffers, sizeof (char *));

  /* return early if there is nothing to do */
  if (data_size == 0) {
    return exc;
  }

  /* receive data from other processors */
  ng_excl = 0;
  for (q = 0; q < num_procs; ++q) {
    ng_incl = ghost->proc_offsets[q + 1];
    ng = ng_incl - ng_excl;
    P4EST_ASSERT (ng >= 0);
    if (ng > 0) {
      r = (sc_MPI_Request *) sc_array_push (&exc->requests);
      mpiret = sc_MPI_Irecv ((char *) ghost_data + ng_excl * data_size,
                             ng * data_size, sc_MPI_BYTE, q,
                             P4EST_COMM_GHOST_EXCHANGE, p4est->mpicomm, r);
      SC_CHECK_MPI (mpiret);
      ng_excl = ng_incl;
    }
  }
  P4EST_ASSERT (ng_excl == (p4est_locidx_t) ghost->ghosts.elem_count);

  /* send data to other processors */
  ng_excl = 0;
  for (q = 0; q < num_procs; ++q) {
    ng_incl = ghost->mirror_proc_offsets[q + 1];
    ng = ng_incl - ng_excl;
    P4EST_ASSERT (ng >= 0);
    if (ng > 0) {
      /* every peer populates its own send buffer */
      sbuf = (char **) sc_array_push (&exc->sbuffers);
      mem = *sbuf = P4EST_ALLOC (char, ng * data_size);
      for (theg = 0; theg < ng; ++theg) {
        mirr = ghost->mirror_proc_mirrors[ng_excl + theg];
        P4EST_ASSERT (0 <= mirr && (size_t) mirr < ghost->mirrors.elem_count);
        memcpy (mem, mirror_data[mirr], data_size);
        mem += data_size;
      }
      r = (sc_MPI_Request *) sc_array_push (&exc->requests);
      mpiret = sc_MPI_Isend (*sbuf, ng * data_size, sc_MPI_BYTE, q,
                             P4EST_COMM_GHOST_EXCHANGE, p4est->mpicomm, r);
      SC_CHECK_MPI (mpiret);
      ng_excl = ng_incl;
    }
  }

  /* we are done posting the messages */
  return exc;
}

void
p4est_ghost_exchange_custom_end (p4est_ghost_exchange_t * exc)
{
  int                 mpiret;
  size_t              zz;
  char              **sbuf;

  /* don't confuse this function with p4est_ghost_exchange_data_end */
  P4EST_ASSERT (exc->is_custom);

  /* don't confuse it with p4est_ghost_exchange_custom_levels_end either */
  P4EST_ASSERT (!exc->is_levels);

  /* wait for messages to complete and clean up */
  mpiret = sc_MPI_Waitall (exc->requests.elem_count, (sc_MPI_Request *)
                           exc->requests.array, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_reset (&exc->requests);
  for (zz = 0; zz < exc->sbuffers.elem_count; ++zz) {
    sbuf = (char **) sc_array_index (&exc->sbuffers, zz);
    P4EST_FREE (*sbuf);
  }
  sc_array_reset (&exc->sbuffers);

  /* free the store */
  P4EST_FREE (exc);
}

void
p4est_ghost_exchange_custom_levels (p4est_t * p4est, p4est_ghost_t * ghost,
                                    int minlevel, int maxlevel,
                                    size_t data_size,
                                    void **mirror_data, void *ghost_data)
{
  p4est_ghost_exchange_custom_levels_end
    (p4est_ghost_exchange_custom_levels_begin (p4est, ghost,
                                               minlevel, maxlevel, data_size,
                                               mirror_data, ghost_data));
}

p4est_ghost_exchange_t *
p4est_ghost_exchange_custom_levels_begin (p4est_t * p4est,
                                          p4est_ghost_t * ghost,
                                          int minlevel, int maxlevel,
                                          size_t data_size,
                                          void **mirror_data,
                                          void *ghost_data)
{
  const int           num_procs = p4est->mpisize;
  int                 mpiret;
  int                 q;
  int                *theq, *qactive, *qbuffer;
  char               *mem, **rbuf, **sbuf;
  p4est_locidx_t      ng_excl, ng_incl, ng, theg;
  p4est_locidx_t      lmatches;
  p4est_locidx_t      mirr;
  p4est_quadrant_t   *g, *m;
  p4est_ghost_exchange_t *exc;
  sc_MPI_Request     *r;

  if (minlevel <= 0 && maxlevel >= P4EST_QMAXLEVEL) {
    /* this case can be processed by a more specialized function */
    exc = p4est_ghost_exchange_custom_begin (p4est, ghost, data_size,
                                             mirror_data, ghost_data);
    P4EST_ASSERT (exc->is_custom);
    P4EST_ASSERT (!exc->is_levels);
    exc->is_levels = 1;

    /* the completion function will have to switch for this case */
    return exc;
  }

  /* initialize transient storage */
  exc = P4EST_ALLOC_ZERO (p4est_ghost_exchange_t, 1);
  exc->is_custom = 1;
  exc->is_levels = 1;
  exc->p4est = p4est;
  exc->ghost = ghost;
  exc->minlevel = minlevel;
  exc->maxlevel = maxlevel;
  exc->data_size = data_size;
  exc->ghost_data = ghost_data;
  sc_array_init (&exc->requests, sizeof (sc_MPI_Request));
  sc_array_init (&exc->rrequests, sizeof (sc_MPI_Request));
  sc_array_init (&exc->rbuffers, sizeof (char *));
  sc_array_init (&exc->sbuffers, sizeof (char *));

  /* return early if there is nothing to do */
  if (data_size == 0 || minlevel > maxlevel) {
    return exc;
  }
  qactive = exc->qactive = P4EST_ALLOC (int, num_procs);
  qbuffer = exc->qbuffer = P4EST_ALLOC (int, num_procs);

  /* receive data from other processors */
  ng_excl = 0;
  for (q = 0; q < num_procs; ++q) {
    qactive[q] = -1;
    qbuffer[q] = -1;
    ng_incl = ghost->proc_offsets[q + 1];
    ng = ng_incl - ng_excl;
    P4EST_ASSERT (ng >= 0);
    if (ng > 0) {
      P4EST_ASSERT (q != p4est->mpirank);
      /* run through ghosts to count the matching level quadrants */
      for (lmatches = 0, theg = 0; theg < ng; ++theg) {
        g = p4est_quadrant_array_index (&ghost->ghosts, ng_excl + theg);
        if (minlevel <= (int) g->level && (int) g->level <= maxlevel) {
          ++lmatches;
        }
      }
      if (lmatches > 0) {
        theq = qactive + exc->rrequests.elem_count;
        r = (sc_MPI_Request *) sc_array_push (&exc->rrequests);
        if (lmatches < ng) {
          /* every peer populates its own receive buffer */
          *theq = q;
          qbuffer[q] = (int) exc->rbuffers.elem_count;
          rbuf = (char **) sc_array_push (&exc->rbuffers);
          *rbuf = P4EST_ALLOC (char, lmatches * data_size);
          mpiret = sc_MPI_Irecv (*rbuf, lmatches * data_size, sc_MPI_BYTE, q,
                                 P4EST_COMM_GHOST_EXCHANGE, p4est->mpicomm,
                                 r);
        }
        else {
          /* use the ghost data memory as is */
          *theq = -1;
          mpiret = sc_MPI_Irecv ((char *) ghost_data + ng_excl * data_size,
                                 ng * data_size, sc_MPI_BYTE, q,
                                 P4EST_COMM_GHOST_EXCHANGE, p4est->mpicomm,
                                 r);
        }
        SC_CHECK_MPI (mpiret);
      }
      ng_excl = ng_incl;
    }
  }
  P4EST_ASSERT (ng_excl == (p4est_locidx_t) ghost->ghosts.elem_count);

  /* send data to other processors */
  ng_excl = 0;
  for (q = 0; q < num_procs; ++q) {
    ng_incl = ghost->mirror_proc_offsets[q + 1];
    ng = ng_incl - ng_excl;
    P4EST_ASSERT (ng >= 0);
    if (ng > 0) {
      P4EST_ASSERT (q != p4est->mpirank);
      /* run through mirrors to count the matching level quadrants */
      for (lmatches = 0, theg = 0; theg < ng; ++theg) {
        mirr = ghost->mirror_proc_mirrors[ng_excl + theg];
        m = p4est_quadrant_array_index (&ghost->mirrors, mirr);
        if (minlevel <= (int) m->level && (int) m->level <= maxlevel) {
          ++lmatches;
        }
      }
      if (lmatches > 0) {
        /* every peer populates its own send buffer */
        sbuf = (char **) sc_array_push (&exc->sbuffers);
        mem = *sbuf = P4EST_ALLOC (char, lmatches * data_size);
        for (theg = 0; theg < ng; ++theg) {
          mirr = ghost->mirror_proc_mirrors[ng_excl + theg];
          m = p4est_quadrant_array_index (&ghost->mirrors, mirr);
          if (minlevel <= (int) m->level && (int) m->level <= maxlevel) {
            memcpy (mem, mirror_data[mirr], data_size);
            mem += data_size;
          }
        }
        r = (sc_MPI_Request *) sc_array_push (&exc->requests);
        mpiret = sc_MPI_Isend (*sbuf, lmatches * data_size, sc_MPI_BYTE, q,
                               P4EST_COMM_GHOST_EXCHANGE, p4est->mpicomm, r);
        SC_CHECK_MPI (mpiret);
      }
      ng_excl = ng_incl;
    }
  }

  /* we are done posting messages */
  return exc;
}

void
p4est_ghost_exchange_custom_levels_end (p4est_ghost_exchange_t * exc)
{
  p4est_ghost_t      *ghost = exc->ghost;
#ifdef P4EST_ENABLE_DEBUG
  p4est_t            *p4est = exc->p4est;
  const int           num_procs = p4est->mpisize;
#endif
  const int           minlevel = exc->minlevel;
  const int           maxlevel = exc->maxlevel;
  const size_t        data_size = exc->data_size;
  int                 mpiret;
  int                 i, expected, remaining, received, *peers;
  int                 q;
  char              **rbuf, **sbuf;
  size_t              zz;
  p4est_locidx_t      ng_excl, ng_incl, ng, theg;
  p4est_locidx_t      lmatches;
  p4est_quadrant_t   *g;

  /* make sure that the begin function matches the end function */
  P4EST_ASSERT (exc->is_custom);
  P4EST_ASSERT (exc->is_levels);

  /* check whether we have used the specialized function */
  if (minlevel <= 0 && maxlevel >= P4EST_QMAXLEVEL) {
    exc->is_levels = 0;
    p4est_ghost_exchange_custom_end (exc);
    return;
  }

  /* wait for receives and copy data into the proper result array */
  peers = P4EST_ALLOC (int, exc->rrequests.elem_count);
  expected = remaining = (int) exc->rrequests.elem_count;
  while (remaining > 0) {
    mpiret =
      sc_MPI_Waitsome (expected, (sc_MPI_Request *) exc->rrequests.array,
                       &received, peers, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (received != sc_MPI_UNDEFINED);
    P4EST_ASSERT (received > 0);
    for (i = 0; i < received; ++i) {
      P4EST_ASSERT (0 <= peers[i] &&
                    peers[i] < (int) exc->rrequests.elem_count);
      q = exc->qactive[peers[i]];
      if (q >= 0) {
        P4EST_ASSERT (q != p4est->mpirank && q < num_procs);
        ng_excl = ghost->proc_offsets[q];
        ng_incl = ghost->proc_offsets[q + 1];
        ng = ng_incl - ng_excl;
        P4EST_ASSERT (ng > 0);
        /* run through ghosts to copy the matching level quadrants' data */
        rbuf = (char **) sc_array_index_int (&exc->rbuffers, exc->qbuffer[q]);
        for (lmatches = 0, theg = 0; theg < ng; ++theg) {
          g = p4est_quadrant_array_index (&ghost->ghosts, ng_excl + theg);
          if (minlevel <= (int) g->level && (int) g->level <= maxlevel) {
            memcpy ((char *) exc->ghost_data + (ng_excl + theg) * data_size,
                    *rbuf + lmatches * data_size, data_size);
            ++lmatches;
          }
        }
        P4EST_FREE (*rbuf);
        exc->qactive[peers[i]] = -1;
        exc->qbuffer[q] = -1;
      }
    }
    remaining -= received;
  }
  P4EST_FREE (peers);
  P4EST_FREE (exc->qactive);
  P4EST_FREE (exc->qbuffer);
  sc_array_reset (&exc->rrequests);
  sc_array_reset (&exc->rbuffers);

  /* wait for sends and clean up */
  mpiret = sc_MPI_Waitall (exc->requests.elem_count, (sc_MPI_Request *)
                           exc->requests.array, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_reset (&exc->requests);
  for (zz = 0; zz < exc->sbuffers.elem_count; ++zz) {
    sbuf = (char **) sc_array_index (&exc->sbuffers, zz);
    P4EST_FREE (*sbuf);
  }
  sc_array_reset (&exc->sbuffers);

  /* free temporary storage */
  P4EST_FREE (exc);
}

#ifdef P4EST_ENABLE_MPI

static void
p4est_ghost_expand_insert (p4est_quadrant_t * q, p4est_topidx_t t,
                           p4est_locidx_t idx, sc_array_t * send_bufs,
                           int target, int owner, int is_ghost)
{
  sc_array_t         *send_buf =
    (sc_array_t *) sc_array_index_int (send_bufs, target);
  p4est_quadrant_t   *qp;
  /* add to mirrors */

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  qp = p4est_quadrant_array_push (send_buf);
  qp->x = q->x;
  qp->y = q->y;
#ifdef P4_TO_P8
  qp->z = q->z;
#endif
  qp->level = q->level;
  qp->p.piggy3.which_tree = t;
  if (is_ghost) {
    sc_array_t         *from_buf =
      (sc_array_t *) sc_array_index_int (send_bufs, owner);
    p4est_quadrant_t   *qp2;

    P4EST_ASSERT (owner != target);
    P4EST_ASSERT (q->p.piggy3.which_tree == t);
    qp->p.piggy3.local_num = q->p.piggy3.local_num;

    qp2 = p4est_quadrant_array_push (from_buf);
    qp2->x = q->x;
    qp2->y = q->y;
#ifdef P4_TO_P8
    qp2->z = q->z;
#endif
    qp2->level = q->level;
    qp2->p.piggy1.which_tree = t;
    qp2->p.piggy1.owner_rank = target;  /* yes, we're putting the target in the
                                           owner field */
    P4EST_ASSERT (p4est_quadrant_is_valid (qp2));
  }
  else {
    qp->p.piggy3.local_num = idx;
  }
  P4EST_ASSERT (p4est_quadrant_is_valid (qp));
}

static void
p4est_ghost_expand_kernel (p4est_topidx_t t, p4est_quadrant_t * mq,
                           p4est_topidx_t nt, p4est_quadrant_t * nq,
                           sc_array_t * quads, int quads_is_ghost,
                           sc_array_t * mview, sc_array_t * pview,
                           p4est_connect_type_t btype, int point,
                           sc_array_t * tempquads, sc_array_t * temptrees,
                           p4est_connectivity_t * conn,
                           p4est_locidx_t pview_offset,
                           p4est_locidx_t quads_offset,
                           p4est_t * p4est, int target,
                           sc_array_t * send_bufs)
{
  ssize_t             fidx, lidx;
  size_t              zz, zy;
  int                 owner = -1;

  /* we want a list of all quadrants in quads that overlap nq */

  /* this will return a quad that overlaps nq, but not necessarily the first
   * */
  fidx = sc_array_bsearch (quads, nq, p4est_quadrant_disjoint);
  if (fidx < 0) {
    /* nothing to be done */
    return;
  }

  lidx = fidx;

  /* walk fidx back to find the first quad that overlaps nq */
  while (fidx > 0) {
    p4est_quadrant_t   *testq = p4est_quadrant_array_index (quads, (size_t)
                                                            fidx - 1);

    if (p4est_quadrant_disjoint (testq, nq)) {
      P4EST_ASSERT (p4est_quadrant_compare (testq, nq) < 0);
      break;
    }
    fidx--;
  }

  /* walk lidx forward to find the last quad that overlaps nq */
  while (lidx < (ssize_t) quads->elem_count - 1) {
    p4est_quadrant_t   *testq = p4est_quadrant_array_index (quads, (size_t)
                                                            lidx + 1);

    if (p4est_quadrant_disjoint (testq, nq)) {
      P4EST_ASSERT (p4est_quadrant_compare (testq, nq) > 0);
      break;
    }
    lidx++;
  }

  /* for every overlapping quadrant, test to see if mq overlaps the neighbor
   * */
  for (zz = (size_t) fidx; zz <= (size_t) lidx; zz++) {
    ssize_t             midx;
    p4est_topidx_t      nnt;
    p4est_quadrant_t   *p = p4est_quadrant_array_index (quads, zz);
    p4est_quadrant_t    np;

    P4EST_ASSERT (p4est_quadrant_overlaps (p, nq));
    /* if we are searching local quads, avoid trying to send mirrors back to
     * the procs for which they are already ghosts */
    if (!quads_is_ghost) {
      midx = sc_array_bsearch (mview, p, p4est_quadrant_compare);

      /* p is a mirror, but is it a mirror for this proc?  search in pview */
      if (midx >= 0) {
        ssize_t             pidx;
        p4est_locidx_t      lmidx;

        lmidx = midx + pview_offset;
        pidx = sc_array_bsearch (pview, &lmidx, p4est_locidx_compare);
        if (pidx >= 0) {
          /* this is already a mirror that this proc knows about */
          continue;
        }
      }
    }
    else {
      /* if we are searching for ghost quads, don't send back one owned by the
       * target proc */
      owner = p4est_comm_find_owner (p4est, nt, p, target);
      if (owner == target) {
        continue;
      }
    }

    /* now create the appropriate neighbor and test for overlaps */
    if (btype == P4EST_CONNECT_FACE) {
      nnt = p4est_quadrant_face_neighbor_extra (p, nt, point, &np, NULL,
                                                conn);
      P4EST_ASSERT (nnt == nt || nnt == t);
      if (nnt == t && p4est_quadrant_overlaps (mq, &np)) {
        p4est_ghost_expand_insert (p, nt, (p4est_locidx_t) zz +
                                   quads_offset, send_bufs, target, owner,
                                   quads_is_ghost);
      }
    }
#ifdef P4_TO_P8
    else if (btype == P8EST_CONNECT_EDGE) {
      p8est_quadrant_edge_neighbor_extra (p, nt, point, tempquads, temptrees,
                                          NULL, conn);

      for (zy = 0; zy < tempquads->elem_count; zy++) {
        nnt = *((p4est_topidx_t *) sc_array_index (temptrees, zy));

        if (nnt == t) {
          p4est_quadrant_t   *tempq =
            p4est_quadrant_array_index (tempquads, zy);

          if (p4est_quadrant_overlaps (mq, tempq)) {
            p4est_ghost_expand_insert (p, nt, (p4est_locidx_t) zz +
                                       quads_offset, send_bufs, target, owner,
                                       quads_is_ghost);
            break;
          }
        }
      }
      sc_array_reset (tempquads);
      sc_array_reset (temptrees);
    }
#endif
    else {
      P4EST_ASSERT (btype == P4EST_CONNECT_CORNER);

      p4est_quadrant_corner_neighbor_extra (p, nt, point, tempquads,
                                            temptrees, NULL, conn);

      for (zy = 0; zy < tempquads->elem_count; zy++) {
        nnt = *((p4est_topidx_t *) sc_array_index (temptrees, zy));

        if (nnt == t) {
          p4est_quadrant_t   *tempq =
            p4est_quadrant_array_index (tempquads, zy);

          if (p4est_quadrant_overlaps (mq, tempq)) {
            /* add to mirrors */
            p4est_ghost_expand_insert (p, nt, (p4est_locidx_t) zz +
                                       quads_offset, send_bufs, target, owner,
                                       quads_is_ghost);
            break;
          }
        }
      }
      sc_array_reset (tempquads);
      sc_array_reset (temptrees);
    }
  }
}

static void
p4est_ghost_expand_int (p4est_topidx_t t, p4est_quadrant_t * mq,
                        p4est_topidx_t nt, p4est_quadrant_t * nq,
                        sc_array_t * pview, p4est_connect_type_t btype,
                        int point, sc_array_t * tempquads,
                        sc_array_t * temptrees, int target,
                        p4est_t * p4est, p4est_ghost_t * ghost,
                        sc_array_t * send_bufs)
{
  sc_array_t          mview;
  sc_array_t          gview;

  sc_array_init_view (&mview, &ghost->mirrors, ghost->mirror_tree_offsets[nt],
                      ghost->mirror_tree_offsets[nt + 1] -
                      ghost->mirror_tree_offsets[nt]);
  sc_array_init_view (&gview, &ghost->ghosts, ghost->tree_offsets[nt],
                      ghost->tree_offsets[nt + 1] - ghost->tree_offsets[nt]);

  if (nt >= p4est->first_local_tree && nt <= p4est->last_local_tree) {
    p4est_tree_t       *ntree = p4est_tree_array_index (p4est->trees, nt);

    p4est_ghost_expand_kernel (t, mq, nt, nq, &ntree->quadrants, 0, &mview,
                               pview, btype, point, tempquads, temptrees,
                               p4est->connectivity,
                               ghost->mirror_tree_offsets[nt],
                               ntree->quadrants_offset,
                               p4est, target, send_bufs);
  }

  p4est_ghost_expand_kernel (t, mq, nt, nq, &gview, 1, &mview, pview, btype,
                             point, tempquads, temptrees, p4est->connectivity,
                             ghost->mirror_tree_offsets[nt],
                             ghost->tree_offsets[nt],
                             p4est, target, send_bufs);

  sc_array_reset (&mview);
  sc_array_reset (&gview);
}

static int
p4est_quadrant_compare_piggy_proc (const void *a, const void *b)
{
  const p4est_quadrant_t *A = (const p4est_quadrant_t *) a;
  const p4est_quadrant_t *B = (const p4est_quadrant_t *) b;
  int                 ret = p4est_quadrant_compare_piggy (a, b);

  if (ret) {
    return ret;
  }
  return (A->p.piggy1.owner_rank - B->p.piggy1.owner_rank);
}

#endif /* P4EST_ENABLE_MPI */

static void
p4est_ghost_expand_internal (p4est_t * p4est, p4est_lnodes_t * lnodes,
                             p4est_ghost_t * ghost)
{
#ifdef P4EST_ENABLE_MPI
  int                 p;
  int                 mpisize = p4est->mpisize;
  int                 mpirank = p4est->mpirank;
  MPI_Comm            comm = p4est->mpicomm;
  p4est_connect_type_t btype = ghost->btype;
  sc_array_t         *mirrors = &ghost->mirrors;
  sc_array_t         *new_mirrors;
  p4est_locidx_t     *mirror_tree_offsets = ghost->mirror_tree_offsets;
  p4est_locidx_t     *mirror_proc_mirrors = ghost->mirror_proc_mirrors;
  p4est_locidx_t     *mirror_proc_offsets = ghost->mirror_proc_offsets;
  p4est_locidx_t     *proc_offsets = ghost->proc_offsets;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_locidx_t     *mpf, *mpfo;
  p4est_locidx_t     *send_counts, *recv_counts;
  MPI_Request        *recv_request, *send_request;
  MPI_Request        *recv_load_request, *send_load_request;
  int                 num_peers;
  int                 peer;
  int                 mpiret;
  sc_array_t         *send_bufs, *buf;
  size_t              zz, *ppz;
  p4est_topidx_t      t;
  sc_array_t         *nmpma, *nmpfa;
  p4est_locidx_t      old_num_ghosts, num_new_ghosts, ghost_offset;
  p4est_locidx_t      old_num_mirrors, new_count, new_num_mirrors;
  sc_array_t         *ghost_layer = &ghost->ghosts;
  sc_array_t          split;
  sc_array_t         *tempquads;
  sc_array_t         *temptrees;
  sc_array_t         *tempquads2;
  sc_array_t         *temptrees2;
  sc_array_t         *npoints;
  p4est_locidx_t     *ntq_offset = NULL;
  p4est_locidx_t     *node_to_quad = NULL;
  p4est_topidx_t     *node_to_tree = NULL;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING "_ghost_expand %s\n",
                            p4est_connect_type_string (btype));
  p4est_log_indent_push ();

  tempquads = sc_array_new (sizeof (p4est_quadrant_t));
  temptrees = sc_array_new (sizeof (p4est_topidx_t));
  tempquads2 = sc_array_new (sizeof (p4est_quadrant_t));
  temptrees2 = sc_array_new (sizeof (p4est_topidx_t));
  npoints = sc_array_new (sizeof (int));

  /* if lnodes, build node_to_quad */
  if (lnodes) {
    p4est_gloidx_t    **mirror_data;
    p4est_locidx_t      il, num_mirrors, qid, K = lnodes->num_local_elements;
    p4est_locidx_t      G = (p4est_locidx_t) ghost->ghosts.elem_count;
    p4est_locidx_t      N = lnodes->num_local_nodes;
    int                 v, nid;
    int                 vnodes = lnodes->vnodes;
    int                *node_to_quad_count;
    int                *quad_counted;
    p4est_gloidx_t     *quad_to_node_global;
    p4est_gloidx_t     *quad_to_node_global_ghost;
    p4est_topidx_t      flt, llt, t;

    quad_to_node_global = P4EST_ALLOC (p4est_gloidx_t, K * vnodes);
    quad_to_node_global_ghost = P4EST_ALLOC (p4est_gloidx_t, G * vnodes);
#ifdef P4EST_ENABLE_DEBUG
    memset (quad_to_node_global_ghost, -1,
            vnodes * G * sizeof (p4est_gloidx_t));
#endif

    mirror_data = P4EST_ALLOC (p4est_gloidx_t *, ghost->mirrors.elem_count);

    for (qid = 0; qid < K * vnodes; qid++) {
      quad_to_node_global[qid] =
        p4est_lnodes_global_index (lnodes, lnodes->element_nodes[qid]);
    }

    num_mirrors = (p4est_locidx_t) ghost->mirrors.elem_count;
    for (il = 0; il < num_mirrors; il++) {
      p4est_quadrant_t   *q;

      q = p4est_quadrant_array_index (&ghost->mirrors, il);

      qid = q->p.piggy3.local_num;
      mirror_data[il] = &quad_to_node_global[qid * vnodes];
    }

    p4est_ghost_exchange_custom (p4est, ghost,
                                 (size_t) vnodes * sizeof (p4est_gloidx_t),
                                 (void **) mirror_data,
                                 quad_to_node_global_ghost);
    P4EST_FREE (mirror_data);
    P4EST_FREE (quad_to_node_global);

    /* convert global back to local */
    for (qid = 0; qid < G; qid++) {
      for (v = 0; v < vnodes; v++) {
        p4est_gloidx_t      gnid;

        gnid = quad_to_node_global_ghost[qid * vnodes + v];
#ifdef P4EST_ENABLE_DEBUG
        P4EST_ASSERT (gnid >= 0);
#endif
        nid = -1;
        if ((gnid >= lnodes->global_offset) &&
            (gnid < (lnodes->global_offset + lnodes->owned_count))) {
          nid = gnid - lnodes->global_offset;
        }
        else {
          sc_array_t          view;
          ssize_t             idx;

          sc_array_init_data (&view, lnodes->nonlocal_nodes,
                              sizeof (p4est_gloidx_t),
                              (size_t) (lnodes->num_local_nodes -
                                        lnodes->owned_count));

          idx = sc_array_bsearch (&view, &gnid, p4est_gloidx_compare);
          if (idx >= 0) {
            nid = idx + lnodes->owned_count;
          }
        }
        P4EST_ASSERT (nid == -1
                      || p4est_lnodes_global_index (lnodes, nid) == gnid);
        quad_to_node_global_ghost[qid * vnodes + v] = nid;
      }
    }

    node_to_quad_count = P4EST_ALLOC_ZERO (int, N);
    quad_counted = P4EST_ALLOC_ZERO (int, K);

    /* first count fronts */
    for (p = 0; p < mpisize; p++) {
      p4est_locidx_t      ilstart, ilend;

      ilstart = ghost->mirror_proc_front_offsets[p];
      ilend = ghost->mirror_proc_front_offsets[p + 1];

      for (il = ilstart; il < ilend; il++) {
        p4est_quadrant_t   *q;
        int                 v;

        qid = ghost->mirror_proc_fronts[il];
        q = p4est_quadrant_array_index (&ghost->mirrors, (size_t) qid);
        qid = q->p.piggy3.local_num;
        if (!quad_counted[qid]) {
          quad_counted[qid] = 1;
          for (v = 0; v < vnodes; v++) {
            nid = lnodes->element_nodes[qid * vnodes + v];
            node_to_quad_count[nid]++;
          }
        }
      }
    }

    /* count the rest of the quads and ghosts, only incrementing nodes that
     * have a front quad adjacent */
    for (qid = 0; qid < K; qid++) {
      if (!quad_counted[qid]) {
        quad_counted[qid] = 1;
        for (v = 0; v < vnodes; v++) {
          nid = lnodes->element_nodes[qid * vnodes + v];
          if (node_to_quad_count[nid]) {
            node_to_quad_count[nid]++;
          }
        }
      }
    }
    for (qid = 0; qid < G; qid++) {
      for (v = 0; v < vnodes; v++) {
        nid = (p4est_locidx_t) quad_to_node_global_ghost[qid * vnodes + v];
        if (nid >= 0) {
          if (node_to_quad_count[nid]) {
            node_to_quad_count[nid]++;
          }
        }
      }
    }
    P4EST_FREE (quad_counted);

    ntq_offset = P4EST_ALLOC (p4est_locidx_t, N + 1);
    ntq_offset[0] = 0;
    for (nid = 0; nid < N; nid++) {
      ntq_offset[nid + 1] = ntq_offset[nid] + node_to_quad_count[nid];
    }
    node_to_quad = P4EST_ALLOC (p4est_locidx_t, ntq_offset[N]);
    node_to_tree = P4EST_ALLOC (p4est_topidx_t, ntq_offset[N]);
#ifdef P4EST_ENABLE_DEBUG
    memset (node_to_quad, -1, ntq_offset[N] * sizeof (p4est_locidx_t));
    memset (node_to_tree, -1, ntq_offset[N] * sizeof (p4est_topidx_t));
#endif

    memset (node_to_quad_count, 0, N * sizeof (int));

    flt = p4est->first_local_tree;
    llt = p4est->last_local_tree;
    for (qid = 0, t = flt; t <= llt; t++) {
      p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, t);
      p4est_locidx_t      nquads =
        (p4est_locidx_t) tree->quadrants.elem_count;

      for (il = 0; il < nquads; il++, qid++) {
        for (v = 0; v < vnodes; v++) {
          p4est_locidx_t      ilstart;
          nid = lnodes->element_nodes[qid * vnodes + v];
          ilstart = ntq_offset[nid];
          if (ntq_offset[nid + 1] > ilstart) {
            node_to_quad[ilstart + node_to_quad_count[nid]] = qid;
            node_to_tree[ilstart + node_to_quad_count[nid]++] = t;
          }
        }
      }
    }
    for (qid = 0; qid < G; qid++) {
      p4est_quadrant_t   *q =
        p4est_quadrant_array_index (&ghost->ghosts, (size_t) qid);

      for (v = 0; v < vnodes; v++) {
        p4est_locidx_t      ilstart;
        nid = (p4est_locidx_t) quad_to_node_global_ghost[qid * vnodes + v];
        if (nid >= 0) {
          ilstart = ntq_offset[nid];
          if (ntq_offset[nid + 1] > ilstart) {
            node_to_quad[ilstart + node_to_quad_count[nid]] = -(qid + 1);
            node_to_tree[ilstart + node_to_quad_count[nid]++] =
              q->p.piggy3.which_tree;
          }
        }
      }
    }
#ifdef P4EST_ENABLE_DEBUG
    for (qid = 0; qid < ntq_offset[N]; qid++) {
      P4EST_ASSERT (node_to_tree[qid] >= 0);
    }
#endif

    P4EST_FREE (node_to_quad_count);
    P4EST_FREE (quad_to_node_global_ghost);
  }

  /* post recvs */
  for (p = 0, num_peers = 0; p < mpisize; p++) {
    if (mirror_proc_offsets[p + 1] != mirror_proc_offsets[p]) {
      /* this is an important assertion: if any proc is part of my ghost
       * layer, I am part of its ghost layer */
      P4EST_ASSERT (proc_offsets[p + 1] != proc_offsets[p]);
      num_peers++;
    }
  }
  recv_request = P4EST_ALLOC (MPI_Request, 2 * num_peers);
  send_request = P4EST_ALLOC (MPI_Request, 2 * num_peers);

  recv_counts = P4EST_ALLOC (p4est_locidx_t, 2 * num_peers);
  send_counts = recv_counts + num_peers;

  recv_load_request = recv_request + num_peers;
  send_load_request = send_request + num_peers;

  send_bufs = sc_array_new_size (sizeof (sc_array_t), mpisize);
  for (p = 0; p < mpisize; p++) {
    buf = (sc_array_t *) sc_array_index (send_bufs, p);
    sc_array_init (buf, sizeof (p4est_quadrant_t));
  }

  for (p = 0, peer = 0; p < mpisize; p++) {
    if (mirror_proc_offsets[p + 1] != mirror_proc_offsets[p]) {
      P4EST_ASSERT (p != mpirank);
      P4EST_LDEBUGF ("ghost layer expand post count receive from %d\n", p);
      mpiret = MPI_Irecv (recv_counts + peer, 1, P4EST_MPI_LOCIDX, p,
                          P4EST_COMM_GHOST_EXPAND_COUNT, comm, recv_request +
                          peer);
      SC_CHECK_MPI (mpiret);
      peer++;
    }
  }

  if (ghost->mirror_proc_fronts == ghost->mirror_proc_mirrors) {
    /* create the fronts: the last quads added to the mirrors */
    P4EST_ASSERT (ghost->mirror_proc_front_offsets ==
                  ghost->mirror_proc_offsets);

    ghost->mirror_proc_fronts = P4EST_ALLOC (p4est_locidx_t,
                                             mirror_proc_offsets[mpisize]);
    memcpy (ghost->mirror_proc_fronts, mirror_proc_mirrors,
            sizeof (p4est_locidx_t) * mirror_proc_offsets[mpisize]);

    ghost->mirror_proc_front_offsets = P4EST_ALLOC (p4est_locidx_t,
                                                    mpisize + 1);
    memcpy (ghost->mirror_proc_front_offsets, mirror_proc_offsets,
            sizeof (p4est_locidx_t) * (mpisize + 1));
  }
  mpf = ghost->mirror_proc_fronts;
  mpfo = ghost->mirror_proc_front_offsets;

  /* for every proc */
  for (p = 0; p < mpisize; p++) {
    /* get mirror_proc_offsets */
    p4est_locidx_t      first_mirror = mpfo[p];
    p4est_locidx_t      end_mirror = mpfo[p + 1];
    size_t              zm;
    sc_array_t          pview;

    if (mirror_proc_offsets[p + 1] == mirror_proc_offsets[p]) {
      continue;
    }

    sc_array_init_data (&pview, mirror_proc_mirrors + mirror_proc_offsets[p],
                        sizeof (p4est_locidx_t),
                        mirror_proc_offsets[p + 1] - mirror_proc_offsets[p]);

    /* for every mirror */
    P4EST_ASSERT (first_mirror >= 0 && end_mirror >= 0);
    for (zm = (size_t) first_mirror; zm < (size_t) end_mirror; zm++) {
      int                 f, c;
#ifdef P4_TO_P8
      int                 e;
#endif
      p4est_quadrant_t   *mq = p4est_quadrant_array_index (mirrors,
                                                           (size_t) mpf[zm]);

      t = mq->p.piggy3.which_tree;

      if (lnodes) {
        /* construct adjacency via lnodes */
        int                 v, vnodes = lnodes->vnodes;
        p4est_locidx_t      qid = mq->p.piggy3.local_num;

        for (v = 0; v < vnodes; v++) {
          p4est_locidx_t      nid = lnodes->element_nodes[qid * vnodes + v];
          p4est_locidx_t      qstart, qend, il;
          int                 owner = -1;

          qstart = ntq_offset[nid];
          qend = ntq_offset[nid + 1];
          P4EST_ASSERT (qend > qstart);

          for (il = qstart; il < qend; il++) {
            p4est_locidx_t      nqid = node_to_quad[il];
            p4est_topidx_t      nt = node_to_tree[il];
            p4est_quadrant_t    qtemp;
            p4est_quadrant_t   *q;
            int                 already_in_mirrors = 0;

            P4EST_ASSERT ((double) nqid == (double) nqid);
            P4EST_ASSERT ((double) qid == (double) qid);
            if (nqid == qid) {
              continue;
            }
            if (nqid < 0) {
              /* ghost */
              size_t              idx = (size_t) - (nqid + 1);
              q = p4est_quadrant_array_index (&ghost->ghosts, idx);
              owner = p4est_comm_find_owner (p4est, nt, q, mpirank);
            }
            else {
              /* local */
              ssize_t             idx;
              p4est_tree_t       *tree =
                p4est_tree_array_index (p4est->trees, nt);
              q =
                p4est_quadrant_array_index (&tree->quadrants,
                                            nqid - tree->quadrants_offset);

              owner = mpirank;
              qtemp = *q;

              qtemp.p.piggy3.which_tree = nt;
              qtemp.p.piggy3.local_num = qid;

              idx =
                sc_array_bsearch (mirrors, &qtemp,
                                  p4est_quadrant_compare_piggy);
              if (idx >= 0) {
                p4est_locidx_t      key = (p4est_locidx_t) idx;

                idx = sc_array_bsearch (&pview, &key, p4est_locidx_compare);

                if (idx >= 0) {
                  already_in_mirrors = 1;
                }
              }
            }
            if (!already_in_mirrors && owner != p) {
              p4est_ghost_expand_insert (q, nt, nqid < 0 ? -(nqid + 1) : nqid,
                                         send_bufs, p, owner,
                                         nqid < 0 ? 1 : 0);
            }
          }
        }
      }
      else {
        /* for every face */
        for (f = 0; f < P4EST_FACES; f++) {
          p4est_quadrant_t    nq;
          p4est_topidx_t      nt;
          int                 nf;

          nf = (f ^ 1);
          nt = t;

          /* get the neighbor */
          nt = p4est_quadrant_face_neighbor_extra (mq, t, f, &nq, &nf, conn);
          if (nt == -1) {
            continue;
          }

          nf = nf % P4EST_FACES;

          /* add any quadrant that overlaps this neighbor and touches the mirror
           * to the buffer */
          p4est_ghost_expand_int (t, mq, nt, &nq, &pview, P4EST_CONNECT_FACE,
                                  nf, tempquads, temptrees, p, p4est, ghost,
                                  send_bufs);
        }
        if (btype == P4EST_CONNECT_FACE) {
          continue;
        }

#ifdef P4_TO_P8
        /* for every edge */
        for (e = 0; e < P8EST_EDGES; e++) {
          p4est_quadrant_t   *nq;
          p4est_topidx_t      nt;
          int                 ne;

          nt = t;

          /* get the neighbors */
          p8est_quadrant_edge_neighbor_extra (mq, t, e, tempquads2,
                                              temptrees2, npoints, conn);

          /* for every neighbor */
          for (zz = 0; zz < tempquads2->elem_count; zz++) {
            nq = p4est_quadrant_array_index (tempquads2, zz);
            nt = *((p4est_locidx_t *) sc_array_index (temptrees2, zz));
            ne = *((int *) sc_array_index (npoints, zz));
            ne = ne % P8EST_EDGES;

            /* add any quadrant that overlaps this neighbor and touches the mirror
             * to the buffer */
            p4est_ghost_expand_int (t, mq, nt, nq, &pview, P8EST_CONNECT_EDGE,
                                    ne, tempquads, temptrees, p, p4est, ghost,
                                    send_bufs);

          }
          sc_array_reset (tempquads2);
          sc_array_reset (temptrees2);
          sc_array_reset (npoints);
        }
        if (btype == P8EST_CONNECT_EDGE) {
          continue;
        }
#endif
        /* for every corner */
        for (c = 0; c < P4EST_CHILDREN; c++) {
          p4est_quadrant_t   *nq;
          p4est_topidx_t      nt;
          int                 nc;

          nt = t;

          /* get the neighbors */
          p4est_quadrant_corner_neighbor_extra (mq, t, c, tempquads2,
                                                temptrees2, npoints, conn);

          /* for every neighbor */
          for (zz = 0; zz < tempquads2->elem_count; zz++) {
            nq = p4est_quadrant_array_index (tempquads2, zz);
            nt = *((p4est_locidx_t *) sc_array_index (temptrees2, zz));
            nc = *((int *) sc_array_index (npoints, zz));

            /* add any quadrant that overlaps this neighbor and touches the mirror
             * to the buffer */
            p4est_ghost_expand_int (t, mq, nt, nq, &pview,
                                    P4EST_CONNECT_CORNER, nc, tempquads,
                                    temptrees, p, p4est, ghost, send_bufs);
          }
          sc_array_reset (tempquads2);
          sc_array_reset (temptrees2);
          sc_array_reset (npoints);
        }
      }
    }

    sc_array_reset (&pview);
  }
  if (lnodes) {
    P4EST_FREE (ntq_offset);
    P4EST_FREE (node_to_quad);
    P4EST_FREE (node_to_tree);
  }
  sc_array_destroy (tempquads);
  sc_array_destroy (temptrees);
  sc_array_destroy (tempquads2);
  sc_array_destroy (temptrees2);
  sc_array_destroy (npoints);

  /* Send the counts of ghosts that are going to be sent */
  new_count = 0;
  for (p = 0, peer = 0; p < mpisize; p++) {
    buf = (sc_array_t *) sc_array_index_int (send_bufs, p);

    if (mirror_proc_offsets[p + 1] == mirror_proc_offsets[p]) {
      continue;
    }

    if (buf->elem_count) {
      sc_array_sort (buf, p4est_quadrant_compare_piggy);
      sc_array_uniq (buf, p4est_quadrant_compare_piggy_proc);
    }
    send_counts[peer] = (p4est_locidx_t) buf->elem_count;
    new_count += send_counts[peer];
    P4EST_ASSERT (p != mpirank);
    P4EST_LDEBUGF ("ghost layer expand post count send to %d\n", p);
    mpiret = MPI_Isend (send_counts + peer, 1, P4EST_MPI_LOCIDX, p,
                        P4EST_COMM_GHOST_EXPAND_COUNT, comm, send_request +
                        peer);
    SC_CHECK_MPI (mpiret);
    peer++;
  }
  P4EST_ASSERT (peer == num_peers);
  P4EST_VERBOSEF ("Total new ghosts to send %lld\n", (long long) new_count);

  /* Wait for the counts */
  if (num_peers > 0) {
    mpiret = sc_MPI_Waitall (num_peers, recv_request, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Waitall (num_peers, send_request, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }

#ifdef P4EST_ENABLE_DEBUG
  for (p = 0; p < num_peers; ++p) {
    P4EST_ASSERT (recv_request[p] == MPI_REQUEST_NULL);
  }
  for (p = 0; p < num_peers; ++p) {
    P4EST_ASSERT (send_request[p] == MPI_REQUEST_NULL);
  }
#endif

  /* Count ghosts */
  for (peer = 0, num_new_ghosts = 0; peer < num_peers; ++peer) {
    num_new_ghosts += recv_counts[peer];        /* same type */
  }
  P4EST_VERBOSEF ("Total new ghosts to receive %lld\n",
                  (long long) num_new_ghosts);

  /* Allocate space for the ghosts */
  old_num_ghosts = (p4est_locidx_t) ghost_layer->elem_count;
  sc_array_resize (ghost_layer, (size_t) (old_num_ghosts + num_new_ghosts));

  /* Post receives for the ghosts */
  for (p = 0, peer = 0, ghost_offset = old_num_ghosts; p < mpisize; p++) {

    if (mirror_proc_offsets[p + 1] == mirror_proc_offsets[p]) {
      continue;
    }

    if (recv_counts[peer]) {
      P4EST_LDEBUGF
        ("ghost layer expand post ghost receive %lld quadrants from %d\n",
         (long long) recv_counts[peer], p);
      mpiret =
        MPI_Irecv (ghost_layer->array +
                   ghost_offset * sizeof (p4est_quadrant_t),
                   (int) (recv_counts[peer] * sizeof (p4est_quadrant_t)),
                   MPI_BYTE, p, P4EST_COMM_GHOST_EXPAND_LOAD, comm,
                   recv_load_request + peer);

      SC_CHECK_MPI (mpiret);
      ghost_offset += recv_counts[peer];
    }
    else {
      recv_load_request[peer] = MPI_REQUEST_NULL;
    }
    peer++;
  }
  P4EST_ASSERT (ghost_offset == old_num_ghosts + num_new_ghosts);

  /* Send the ghosts */
  for (p = 0, peer = 0; p < mpisize; p++) {
    if (mirror_proc_offsets[p + 1] == mirror_proc_offsets[p]) {
      continue;
    }
    buf = (sc_array_t *) sc_array_index (send_bufs, p);
    if (buf->elem_count > 0) {
      P4EST_ASSERT ((p4est_locidx_t) buf->elem_count == send_counts[peer]);
      P4EST_LDEBUGF
        ("ghost layer expand post ghost send %lld quadrants to %d\n",
         (long long) send_counts[peer], p);
      mpiret =
        MPI_Isend (buf->array,
                   (int) (send_counts[peer] * sizeof (p4est_quadrant_t)),
                   MPI_BYTE, p, P4EST_COMM_GHOST_EXPAND_LOAD, comm,
                   send_load_request + peer);
      SC_CHECK_MPI (mpiret);
    }
    else {
      send_load_request[peer] = MPI_REQUEST_NULL;
    }
    peer++;
  }

  /* Wait for everything */
  if (num_peers > 0) {
    mpiret =
      sc_MPI_Waitall (num_peers, recv_load_request, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    mpiret =
      sc_MPI_Waitall (num_peers, send_load_request, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }

#ifdef P4EST_ENABLE_DEBUG
  for (p = 0; p < num_peers; p++) {
    P4EST_ASSERT (recv_load_request[p] == MPI_REQUEST_NULL);
  }
  for (p = 0; p < num_peers; p++) {
    P4EST_ASSERT (send_load_request[p] == MPI_REQUEST_NULL);
  }
#endif

  /* Clean up */
  P4EST_FREE (recv_counts);
  P4EST_FREE (recv_request);
  P4EST_FREE (send_request);

  /* sift bridges out of send buffers so that we can reuse buffers when
   * updating mirrors*/
  for (p = 0; p < mpisize; p++) {
    size_t              count;
    p4est_quadrant_t   *q1, *q2;

    buf = (sc_array_t *) sc_array_index_int (send_bufs, p);
    count = buf->elem_count;

    if (!count) {
      continue;
    }

    q1 = p4est_quadrant_array_index (buf, 0);
    for (zz = 0; zz < buf->elem_count; zz++) {
      q2 = p4est_quadrant_array_index (buf, zz);
      P4EST_ASSERT (p4est_quadrant_is_valid (q2));
      if (p4est_comm_is_owner (p4est, q2->p.which_tree, q2, mpirank)) {
#ifdef P4EST_ENABLE_DEBUG
        ssize_t             idx;

        idx = sc_array_bsearch (mirrors, q2, p4est_quadrant_compare_piggy);

        if (idx >= 0) {
          ssize_t             idx2;
          sc_array_t          pview;
          p4est_locidx_t      locidx = (p4est_locidx_t) idx;
          sc_array_init_data (&pview,
                              mirror_proc_mirrors + mirror_proc_offsets[p],
                              sizeof (p4est_locidx_t),
                              mirror_proc_offsets[p + 1] -
                              mirror_proc_offsets[p]);
          idx2 = sc_array_bsearch (&pview, &locidx, p4est_locidx_compare);
          P4EST_ASSERT (idx2 < 0);
        }
#endif
        if (q1 != q2) {
          *(q1++) = *q2;
        }
        else {
          q1++;
        }
      }
      else {
        count--;
      }
    }
    P4EST_LDEBUGF ("ghost layer expand sending %lld new non-bridges to %d\n",
                   (long long) count, p);
    sc_array_resize (buf, count);
  }

  if (num_new_ghosts) {
    p4est_quadrant_t   *q1, *q2;

    q1 = p4est_quadrant_array_index (ghost_layer, old_num_ghosts);

    /* sift out the bridges */
    for (zz = old_num_ghosts; zz < ghost_layer->elem_count; zz++) {
      q2 = p4est_quadrant_array_index (ghost_layer, zz);

      if (p4est_comm_is_owner (p4est, q2->p.which_tree, q2, mpirank)) {
        /* this is a bridge: it is already a mirror for a proc other
         * than p */
        int                 target = q2->p.piggy1.owner_rank;
        ssize_t             idx, idx2;
        p4est_locidx_t      locidx;
        sc_array_t          pview;

        P4EST_ASSERT (0 <= target && target < mpisize);
        P4EST_ASSERT (target != mpirank);

        num_new_ghosts--;
        idx = sc_array_bsearch (mirrors, q2, p4est_quadrant_compare_piggy);
        P4EST_ASSERT (idx >= 0);
        /* does the target already know about this ? */
        locidx = (p4est_locidx_t) idx;
        sc_array_init_data (&pview,
                            mirror_proc_mirrors + mirror_proc_offsets[target],
                            sizeof (p4est_locidx_t),
                            mirror_proc_offsets[target + 1] -
                            mirror_proc_offsets[target]);
        idx2 = sc_array_bsearch (&pview, &locidx, p4est_locidx_compare);
        sc_array_reset (&pview);

        if (idx2 < 0) {
          /* if the target doesn't already know about it, put it in send_bufs
           * */
          p4est_quadrant_t   *q3;

          q3 = p4est_quadrant_array_index (mirrors, (size_t) idx);
          P4EST_ASSERT (p4est_quadrant_is_equal_piggy (q2, q3));
          buf = (sc_array_t *) sc_array_index_int (send_bufs, target);
          (void) p4est_quadrant_array_push_copy (buf, q3);
        }
      }
      else {
        if (q1 != q2) {
          *(q1++) = *q2;
        }
        else {
          q1++;
        }
      }
    }

    P4EST_ASSERT (num_new_ghosts >= 0);

    for (p = 0; p < mpisize; p++) {
      buf = (sc_array_t *) sc_array_index_int (send_bufs, p);

      sc_array_sort (buf, p4est_quadrant_compare_piggy);
      sc_array_uniq (buf, p4est_quadrant_compare_piggy);
    }

    sc_array_resize (ghost_layer, (size_t) (old_num_ghosts + num_new_ghosts));
    if (num_new_ghosts) {
      /* update the ghost layer */
      sc_array_sort (ghost_layer, p4est_quadrant_compare_piggy);
      sc_array_uniq (ghost_layer, p4est_quadrant_compare_piggy);

      num_new_ghosts = ghost_layer->elem_count - old_num_ghosts;

      if (num_new_ghosts) {

        P4EST_LDEBUGF ("ghost layer expand receiving %lld new non-bridges\n",
                       (long long) num_new_ghosts);
        /* calculate tree offsets */
        sc_array_init (&split, sizeof (size_t));
        sc_array_split (ghost_layer, &split,
                        (size_t) conn->num_trees, ghost_tree_type, NULL);
        P4EST_ASSERT (split.elem_count == (size_t) conn->num_trees + 1);
        for (t = 0; t <= conn->num_trees; ++t) {
          ppz = (size_t *) sc_array_index (&split, (size_t) t);
          ghost->tree_offsets[t] = *ppz;
        }
        sc_array_reset (&split);

        /* calculate proc offsets */
        sc_array_init (&split, sizeof (size_t));
        sc_array_split (ghost_layer, &split,
                        (size_t) mpisize, ghost_proc_type, p4est);
        P4EST_ASSERT (split.elem_count == (size_t) mpisize + 1);
        for (p = 0; p <= mpisize; p++) {
          ppz = (size_t *) sc_array_index (&split, (size_t) p);
          proc_offsets[p] = (p4est_locidx_t) (*ppz);
        }
        sc_array_reset (&split);
      }
    }
  }

  /* we're going to build new mirrors structures, then destroy the old
   * structures later */
  new_mirrors = sc_array_new_size (mirrors->elem_size, mirrors->elem_count);
  sc_array_copy (new_mirrors, mirrors);
  old_num_mirrors = (p4est_locidx_t) mirrors->elem_count;
  for (p = 0; p < mpisize; p++) {
    /* add all of the potentially new mirrors */
    buf = (sc_array_t *) sc_array_index_int (send_bufs, p);

    if (!buf->elem_count) {
      continue;
    }
    else {
      size_t              oldsize = new_mirrors->elem_count;
      sc_array_resize (new_mirrors, oldsize + buf->elem_count);
      memcpy (new_mirrors->array + oldsize * new_mirrors->elem_size,
              buf->array, buf->elem_count * buf->elem_size);
    }
  }
  sc_array_sort (new_mirrors, p4est_quadrant_compare_piggy);
  sc_array_uniq (new_mirrors, p4est_quadrant_compare_piggy);
  new_num_mirrors = (p4est_locidx_t) new_mirrors->elem_count;
  P4EST_ASSERT (new_num_mirrors >= old_num_mirrors);

  if (new_num_mirrors > old_num_mirrors) {
    /* update mirror_tree_offsets */
    sc_array_init (&split, sizeof (size_t));
    sc_array_split (new_mirrors, &split,
                    (size_t) conn->num_trees, ghost_tree_type, NULL);
    P4EST_ASSERT (split.elem_count == (size_t) conn->num_trees + 1);
    for (t = 0; t <= conn->num_trees; ++t) {
      ppz = (size_t *) sc_array_index (&split, (size_t) t);
      mirror_tree_offsets[t] = *ppz;
    }
    sc_array_reset (&split);
  }

  /* update mirror_proc_fronts */
  nmpfa = sc_array_new (sizeof (p4est_locidx_t));
  for (p = 0; p < mpisize; p++) {
    size_t              offset = nmpfa->elem_count;
    buf = (sc_array_t *) sc_array_index_int (send_bufs, p);

    sc_array_resize (nmpfa, offset + buf->elem_count);
    mpfo[p] = (p4est_locidx_t) offset;

    for (zz = 0; zz < buf->elem_count; zz++) {
      ssize_t             idx;
      p4est_quadrant_t   *q1 = p4est_quadrant_array_index (buf, zz);
      p4est_locidx_t     *lp =
        (p4est_locidx_t *) sc_array_index (nmpfa, offset + zz);

      idx = sc_array_bsearch (new_mirrors, q1, p4est_quadrant_compare_piggy);
      P4EST_ASSERT (idx >= 0);
      *lp = (p4est_locidx_t) idx;
    }

    sc_array_reset (buf);
  }
  mpfo[mpisize] = nmpfa->elem_count;
  P4EST_FREE (mpf);
  ghost->mirror_proc_fronts = mpf =
    P4EST_ALLOC (p4est_locidx_t, nmpfa->elem_count);
  memcpy (mpf, nmpfa->array, nmpfa->elem_size * nmpfa->elem_count);
  sc_array_destroy (nmpfa);
  sc_array_destroy (send_bufs);

  /* update mirror_proc_mirrors */
  nmpma = sc_array_new (sizeof (p4est_locidx_t));
  for (p = 0; p < mpisize; p++) {
    size_t              frontsize = mpfo[p + 1] - mpfo[p];
    size_t              offset = nmpma->elem_count;
    p4est_locidx_t      old_offset = mirror_proc_offsets[p];
    p4est_locidx_t      old_count = mirror_proc_offsets[p + 1] - old_offset;

    P4EST_ASSERT (old_count >= 0);
    mirror_proc_offsets[p] = offset;

    P4EST_LDEBUGF
      ("ghost layer expanded with proc %d: send %lld receive %lld\n",
       p, (long long) (old_count + frontsize),
       (long long) (proc_offsets[p + 1] - proc_offsets[p]));
    sc_array_resize (nmpma, offset + old_count + frontsize);
    memcpy (nmpma->array + nmpma->elem_size * offset,
            mpf + mpfo[p], sizeof (p4est_locidx_t) * frontsize);

    if (old_count) {
      sc_array_t          pview;
      for (zz = 0; zz < (size_t) old_count; zz++) {
        ssize_t             idx;
        p4est_quadrant_t   *q1 = p4est_quadrant_array_index (mirrors,
                                                             mirror_proc_mirrors
                                                             [zz +
                                                              old_offset]);
        p4est_locidx_t     *lp =
          (p4est_locidx_t *) sc_array_index (nmpma, offset + frontsize + zz);

        idx =
          sc_array_bsearch (new_mirrors, q1, p4est_quadrant_compare_piggy);
        P4EST_ASSERT (idx >= 0);
        *lp = (p4est_locidx_t) idx;
      }

      sc_array_init_view (&pview, nmpma, offset, old_count + frontsize);
      sc_array_sort (&pview, p4est_locidx_compare);
      sc_array_reset (&pview);
    }
  }
  mirror_proc_offsets[mpisize] = nmpma->elem_count;
  P4EST_FREE (mirror_proc_mirrors);
  ghost->mirror_proc_mirrors = mirror_proc_mirrors = P4EST_ALLOC
    (p4est_locidx_t, nmpma->elem_count);
  memcpy (mirror_proc_mirrors, nmpma->array,
          nmpma->elem_size * nmpma->elem_count);
  sc_array_destroy (nmpma);

  sc_array_resize (mirrors, new_mirrors->elem_count);
  sc_array_copy (mirrors, new_mirrors);
  sc_array_destroy (new_mirrors);

#ifdef P4EST_ENABLE_DEBUG
  for (p = 0; p < mpisize; p++) {
    int                 ghost_count =
      ghost->proc_offsets[p + 1] - ghost->proc_offsets[p];
    int                 mirror_count =
      ghost->mirror_proc_offsets[p + 1] - ghost->mirror_proc_offsets[p];

    P4EST_ASSERT ((mirror_count && ghost_count)
                  || (!mirror_count && !ghost_count));

  }
#endif
  P4EST_ASSERT (p4est_ghost_is_valid (p4est, ghost));

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_ghost_expand\n");
#endif
}

void
p4est_ghost_expand (p4est_t * p4est, p4est_ghost_t * ghost)
{
  p4est_ghost_expand_internal (p4est, NULL, ghost);
}

void
p4est_ghost_expand_by_lnodes (p4est_t * p4est, p4est_lnodes_t * lnodes,
                              p4est_ghost_t * ghost)
{
  p4est_ghost_expand_internal (p4est, lnodes, ghost);
}

int
p4est_ghost_is_valid (p4est_t * p4est, p4est_ghost_t * ghost)
{
  const p4est_topidx_t num_trees = ghost->num_trees;
  const int           mpisize = ghost->mpisize;
  int                 i, mpiret, retval;
  size_t              view_length, proc_length;
  p4est_locidx_t      proc_offset;
  sc_array_t          array, *workspace, *requests;
  uint64_t           *checksums_recv, *checksums_send;

  /* check if the last entries of the offset arrays are the element count
   * of ghosts/mirrors array. */
  if ((size_t) ghost->tree_offsets[num_trees] != ghost->ghosts.elem_count
      || (size_t) ghost->proc_offsets[mpisize] != ghost->ghosts.elem_count
      || (size_t) ghost->mirror_tree_offsets[num_trees] !=
      ghost->mirrors.elem_count) {
    return 0;
  }

  /* check if quadrants in ghost and mirror layer are
   * in p4est_quadrant_compare_piggy order.
   * Also check if tree_offsets, proc_offsets, mirror_tree_offsets
   * and mirror_proc_offsets are sorted.
   */
  if (!sc_array_is_sorted (&ghost->ghosts, p4est_quadrant_compare_piggy) ||
      !sc_array_is_sorted (&ghost->mirrors, p4est_quadrant_compare_piggy)
      || !sc_array_is_sorted (&ghost->mirrors,
                              p4est_quadrant_compare_local_num)) {
    return 0;
  }
  sc_array_init_data (&array, ghost->tree_offsets, sizeof (p4est_locidx_t),
                      num_trees + 1);
  if (!sc_array_is_sorted (&array, p4est_locidx_compare))
    return 0;
  sc_array_init_data (&array, ghost->proc_offsets, sizeof (p4est_locidx_t),
                      mpisize + 1);
  if (!sc_array_is_sorted (&array, p4est_locidx_compare))
    return 0;
  sc_array_init_data (&array, ghost->mirror_tree_offsets,
                      sizeof (p4est_locidx_t), num_trees + 1);
  if (!sc_array_is_sorted (&array, p4est_locidx_compare))
    return 0;
  sc_array_init_data (&array, ghost->mirror_proc_offsets,
                      sizeof (p4est_locidx_t), mpisize + 1);
  if (!sc_array_is_sorted (&array, p4est_locidx_compare))
    return 0;

  /* check if local number in piggy3 data member of the quadrants in ghost is
   * ascending within each rank.
   */
  for (i = 0; i < mpisize; i++) {
    proc_offset = ghost->proc_offsets[i];
    view_length = (size_t) (ghost->proc_offsets[i + 1] - proc_offset);
    sc_array_init_view (&array, &ghost->ghosts, (size_t) proc_offset,
                        view_length);
    if (!sc_array_is_sorted (&array, p4est_quadrant_compare_local_num)) {
      return 0;
    }
  }

  /* check if mirror_proc_offsets is ascending within each rank
   */
  for (i = 0; i < mpisize; i++) {
    proc_offset = ghost->mirror_proc_offsets[i];
    proc_length = (size_t) (ghost->mirror_proc_offsets[i + 1] - proc_offset);
    sc_array_init_data (&array, ghost->mirror_proc_mirrors + proc_offset,
                        sizeof (p4est_locidx_t), proc_length);
    if (!sc_array_is_sorted (&array, p4est_locidx_compare)) {
      return 0;
    }
  }

  /* compare checksums of ghosts with checksums of mirrors */
  checksums_recv = P4EST_ALLOC (uint64_t, mpisize);
  checksums_send = P4EST_ALLOC (uint64_t, mpisize);
  requests = sc_array_new (sizeof (sc_MPI_Request));
  workspace = sc_array_new (sizeof (p4est_quadrant_t));
  for (i = 0; i < mpisize; i++) {
    p4est_locidx_t      count;
    sc_MPI_Request     *req;

    proc_offset = ghost->proc_offsets[i];
    count = ghost->proc_offsets[i + 1] - proc_offset;

    if (count) {
      req = (sc_MPI_Request *) sc_array_push (requests);
      mpiret = sc_MPI_Irecv (&checksums_recv[i], 1, sc_MPI_LONG_LONG_INT, i,
                             P4EST_COMM_GHOST_CHECKSUM, p4est->mpicomm, req);
      SC_CHECK_MPI (mpiret);
    }

    proc_offset = ghost->mirror_proc_offsets[i];
    count = ghost->mirror_proc_offsets[i + 1] - proc_offset;

    if (count) {
      p4est_locidx_t      jl;

      sc_array_truncate (workspace);

      for (jl = proc_offset; jl < proc_offset + count; jl++) {
        p4est_locidx_t      idx;
        p4est_quadrant_t   *q1;

        idx = ghost->mirror_proc_mirrors[jl];

        q1 = p4est_quadrant_array_index (&ghost->mirrors, (size_t) idx);
        (void) p4est_quadrant_array_push_copy (workspace, q1);
      }

      checksums_send[i] =
        (uint64_t) p4est_quadrant_checksum (workspace, NULL, 0);

      req = (sc_MPI_Request *) sc_array_push (requests);
      mpiret = sc_MPI_Isend (&checksums_send[i], 1, sc_MPI_LONG_LONG_INT, i,
                             P4EST_COMM_GHOST_CHECKSUM, p4est->mpicomm, req);
      SC_CHECK_MPI (mpiret);
    }
  }

  mpiret = sc_MPI_Waitall (requests->elem_count, (sc_MPI_Request *)
                           requests->array, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy (workspace);
  sc_array_destroy (requests);
  P4EST_FREE (checksums_send);

  retval = 1;
  for (i = 0; i < mpisize; i++) {
    p4est_locidx_t      count;

    proc_offset = ghost->proc_offsets[i];
    count = ghost->proc_offsets[i + 1] - proc_offset;

    if (count) {
      sc_array_t          view;
      uint64_t            thiscrc;

      sc_array_init_view (&view, &ghost->ghosts, (size_t) proc_offset,
                          (size_t) count);

      thiscrc = (uint64_t) p4est_quadrant_checksum (&view, NULL, 0);
      if (thiscrc != checksums_recv[i]) {
        P4EST_LERRORF ("Ghost layer checksum mismatch: "
                       "proc %d, my checksum %llu, their checksum %llu\n",
                       i, (long long unsigned) thiscrc,
                       (long long unsigned) checksums_recv[i]);
        retval = 0;
      }
    }
  }
  P4EST_FREE (checksums_recv);
  return retval;
}
