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

#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_mesh.h>

/** Gets the procid of the owner of \a q.
 *
 * \param [in] p4est  The forest in which to search for \a q
 * \param [in] treeid The tree id for which \a q belongs.
 * \param [in] q      The quadrant that is being searched for.
 *
 * \return Procid of the owner of \a q or -1 if the qudrant does not exist in
 *         the mesh.
 *
 * \warning Does not work for tree corner neighbors.
 */
static int
p4est_quadrant_find_owner (p4est_t * p4est, p4est_locidx_t treeid,
                           p4est_quadrant_t * q)
{
  const int           rank = p4est->mpirank;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  int                 owner = -1;
  p4est_connectivity_t *conn = p4est->connectivity;
  int                 quad_contact[4];
  int                 face, transform;
  p4est_locidx_t      ntreeid;
  p4est_quadrant_t    tmpq = *q, nq;

  P4EST_QUADRANT_INIT (&nq);

  if (p4est_quadrant_is_inside_root (q)) {
    owner = p4est_comm_find_owner (p4est, treeid, q, rank);
  }
  else {
    quad_contact[0] = (q->y < 0);
    quad_contact[1] = (q->x >= rh);
    quad_contact[2] = (q->y >= rh);
    quad_contact[3] = (q->x < 0);

    /* Make sure we are not a tree corner */
    P4EST_ASSERT (!((quad_contact[0] || quad_contact[2]) &&
                    (quad_contact[1] || quad_contact[3])));

    for (face = 0; face < 4; ++face) {
      if (quad_contact[face]
          && (conn->tree_to_tree[4 * treeid + face] != treeid
              || ((int) conn->tree_to_face[4 * treeid + face] != face))) {
        ntreeid = conn->tree_to_tree[4 * treeid + face];
        break;
      }
    }
    if (face == 4) {
      /* This quadrant does not exist in the mesh */
      owner = -1;
    }
    else {
      transform = p4est_find_face_transform (conn, treeid, face);
      p4est_quadrant_translate_face (&tmpq, face);
      p4est_quadrant_transform_face (&tmpq, &nq, transform);

      owner = p4est_comm_find_owner (p4est, ntreeid, &nq, rank);
    }
  }

  return owner;
}

/** Gets the procids of the owners of \a q.
 *
 * For a quadrant across the corner of a tree has possibly multiple
 * trees in which it lives, and thus multiple procs.
 *
 * \param [in]     p4est      The forest in which to search for \a q.
 * \param [in]     treeid     The tree id for which \a q belongs.
 * \param [in]     treecorner The corner of the tree \a q is across from.
 * \param [in]     q          The quadrant that is being searched for.
 * \param [in,out] qprocs     Starts as an initialize array and ends with
 *                            the list of processors that \a q belongs too.
 */
static void
p4est_quadrant_find_tree_corner_owners (p4est_t * p4est,
                                        p4est_topidx_t treeid,
                                        int treecorner,
                                        p4est_quadrant_t * q,
                                        sc_array_t * q_procs)
{
  const int           rank = p4est->mpirank;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  int                 cproc, *proc;
  size_t              ctree;
  p4est_topidx_t      ctreeid;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_quadrant_t    cq;
  sc_array_t          ctransforms, *cta = &ctransforms;
  p4est_corner_transform_t *ct;

  P4EST_QUADRANT_INIT (&cq);

  P4EST_ASSERT (((q->y < 0) || (q->y >= rh)) && ((q->x >= rh) || (q->x < 0)));

  sc_array_init (cta, sizeof (p4est_corner_transform_t));
  p4est_find_corner_transform (conn, treeid, treecorner, cta);

  sc_array_resize (q_procs, 0);

  for (ctree = 0; ctree < cta->elem_count; ++ctree) {
    ct = sc_array_index (cta, ctree);
    ctreeid = ct->ntree;

    if (ctreeid == treeid)
      continue;

    cq = *q;
    p4est_quadrant_transform_corner (&cq, (int) ct->ncorner, true);

    cproc = p4est_comm_find_owner (p4est, ctreeid, &cq, rank);

    proc = sc_array_push (q_procs);
    *proc = cproc;
  }

  sc_array_reset (cta);
}

/** Get the possible face neighbors of \a q.
 *
 * Gets the all face neighbors, possible assuming the 2-1 constant.
 * If the larger quadrant doesn't exisit than it returned
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
 * \param [out] n0     Filled with the first possible face neighbor, which is
 *                     half of the size if it exists or initialized to
 *                     P4EST_QUADRANT_INIT.
 * \param [out] n1     Filled with the second possible face neighbor, which is
 *                     half of the size if it exists or initialized to
 *                     P4EST_QUADRANT_INIT.
 * \param [out] n2     Filled with the face neighbor, which is the same size.
 * \param [out] n3     Filled with the face neighbor, which is twice the size
 *                     if it exists or initialized to P4EST_QUADRANT_INIT.
 *
 */
static void
p4est_quadrant_get_possible_face_neighbors (p4est_quadrant_t * q,
                                            int face,
                                            p4est_quadrant_t * n0,
                                            p4est_quadrant_t * n1,
                                            p4est_quadrant_t * n2,
                                            p4est_quadrant_t * n3)
{
  p4est_qcoord_t      qh = P4EST_QUADRANT_LEN (q->level);
  p4est_qcoord_t      ph = P4EST_QUADRANT_LEN (q->level - 1);
  p4est_qcoord_t      qh_2 = P4EST_QUADRANT_LEN (q->level + 1);
  int                 qcid = p4est_quadrant_child_id (q);
  int                 rqcid = p4est_corner_to_zorder[qcid];

  if (q->level >= P4EST_MAXLEVEL) {
    P4EST_QUADRANT_INIT (n0);
    P4EST_QUADRANT_INIT (n2);
  }
  else {
    *n0 = *q;
    *n1 = *q;

    n0->level += 1;
    n1->level += 1;

    switch (face) {
    case 0:
      n0->y -= qh_2;

      n1->x += qh_2;
      n1->y -= qh_2;
      break;
    case 1:
      n0->x += qh;

      n1->x += qh;
      n1->y += qh_2;
      break;
    case 2:
      n0->y += qh;

      n1->x += qh_2;
      n1->y += qh;
      break;
    case 3:
      n0->x -= qh_2;

      n1->x -= qh_2;
      n1->y += qh_2;
      break;
    default:
      SC_CHECK_NOT_REACHED ();
      break;
    }
  }

  *n2 = *q;
  switch (face) {
  case 0:
    n2->y -= qh;
    break;
  case 1:
    n2->x += qh;
    break;
  case 2:
    n2->y += qh;
    break;
  case 3:
    n2->x -= qh;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }

  /* Check to see if the larger element exists */
  if (((face != rqcid) && (face != ((rqcid + 3) % 4))) || (q->level == 0)) {
    P4EST_QUADRANT_INIT (n3);
  }
  else {
    p4est_quadrant_parent (q, n3);
    switch (face) {
    case 0:
      n3->y -= ph;
      break;
    case 1:
      n3->x += ph;
      break;
    case 2:
      n3->y += ph;
      break;
    case 3:
      n3->x -= ph;
      break;
    default:
      SC_CHECK_NOT_REACHED ();
      break;
    }
  }
}

/** Get the possible corner neighbors of \a q.
 *
 * Gets the corner face neighbors, possible assuming the 2-1 constant.
 * If the larger quadrant doesn't exisit than it returned as initialized by
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
 * \param [out] n0     Filled with the possible corner neighbor, which is
 *                     half of the size.
 * \param [out] n1     Filled with the face neighbor, which is the same size.
 * \param [out] n2     Filled with the face neighbor, which is twice the size
 *                     if it exists or initialized to P4EST_QUADRANT_INIT.
 *
 */
static void
p4est_quadrant_get_possible_corner_neighbors (p4est_quadrant_t * q,
                                              int corner,
                                              p4est_quadrant_t * n0,
                                              p4est_quadrant_t * n1,
                                              p4est_quadrant_t * n2)
{
  p4est_qcoord_t      qh = P4EST_QUADRANT_LEN (q->level);
  p4est_qcoord_t      ph = P4EST_QUADRANT_LEN (q->level - 1);
  p4est_qcoord_t      qh_2 = P4EST_QUADRANT_LEN (q->level + 1);
  int                 qcid = p4est_quadrant_child_id (q);

  P4EST_ASSERT (q->level != P4EST_MAXLEVEL);

  if (q->level >= P4EST_MAXLEVEL) {
    P4EST_QUADRANT_INIT (n0);
  }
  else {
    *n0 = *q;
    n0->level += 1;

    switch (corner) {
    case 0:
      n0->x -= qh_2;
      n0->y -= qh_2;
      break;
    case 1:
      n0->x += qh_2;
      n0->y -= qh_2;
      break;
    case 2:
      n0->x -= qh_2;
      n0->y += qh_2;
      break;
    case 3:
      n0->x += qh_2;
      n0->y += qh_2;
      break;
    default:
      SC_CHECK_NOT_REACHED ();
      break;
    }
  }

  *n1 = *q;
  switch (corner) {
  case 0:
    n1->x -= qh;
    n1->y -= qh;
    break;
  case 1:
    n1->x += qh;
    n1->y -= qh;
    break;
  case 2:
    n1->x -= qh;
    n1->y += qh;
    break;
  case 3:
    n1->x += qh;
    n1->y += qh;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }

  /* Check to see if the larger element exists */
  if ((corner != qcid) || (q->level == 0)) {
    P4EST_QUADRANT_INIT (n2);
  }
  else {
    p4est_quadrant_parent (q, n2);
    switch (corner) {
    case 0:
      n2->x -= ph;
      n2->y -= ph;
      break;
    case 1:
      n2->x += ph;
      n2->y -= ph;
      break;
    case 2:
      n2->x -= ph;
      n2->y += ph;
      break;
    case 3:
      n2->x += ph;
      n2->y += ph;
      break;
    default:
      SC_CHECK_NOT_REACHED ();
      break;
    }
  }
}

/** Checks if quadrant exists in the local forest or the ghost layer.
 *
 * For quadrants across tree corners it checks that the quadrant exists in
 * all of the corner neighbors
 *
 * \param [in]  p4est        The forest in which to search for \a q
 * \param [in]  ghost_layer  The ghost layer in which to search for \a q
 * \param [in]  treeid       The tree id for which \a q belongs.
 * \param [in]  q            The quadrant that is being searched for.
 * \param [out] exists_arr   Filled for tree corner cases if the qudrant
 *                           exists in the multiple tree corner neighbors.
 *                           An entry for each corner neighbor is set to 1 if
 *                           it exists in the local forest or ghost_layer.
 *
 * \return 1 if the quadrant exists in the local forest or in the ghost_layer,
 *         and 0 if doesn't exist in either
 */
static int
p4est_quadrant_exists (p4est_t * p4est,
                       sc_array_t * ghost_layer,
                       p4est_topidx_t treeid, p4est_quadrant_t * q,
                       sc_array_t * exists_arr)
{
  const int           rank = p4est->mpirank;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  int                 exists = 0, *pexists;
  int                 qproc, face, transform, tree_corner;
  int                 quad_contact[4];
  size_t              ctreeidz, num_ctrees = 0;
  ssize_t             lnid;
  p4est_topidx_t      tqtreeid = -1;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_tree_t       *tree = p4est_array_index_topidx (p4est->trees, treeid);
  p4est_tree_t       *tqtree;
  p4est_quadrant_t    tempq, tq, non_existant;
  sc_array_t         *quadrants = &tree->quadrants;
  sc_array_t          ctransforms, *cta = &ctransforms;
  p4est_corner_transform_t *ct;

  P4EST_QUADRANT_INIT (&non_existant);

  if (non_existant.level == q->level) {
    exists = 0;
  }
  else if (p4est_quadrant_is_inside_root (q)) {
    qproc = p4est_comm_find_owner (p4est, treeid, q, rank);
    if (qproc == rank) {
      lnid = sc_array_bsearch (quadrants, q, p4est_quadrant_compare);
      exists = (lnid != -1) ? 1 : 0;
    }
    else {
      /* off processor so search in the ghost layer */
      tq = *q;
      P4EST_ASSERT (sizeof (long) >= sizeof (p4est_locidx_t));
      tq.p.piggy1.which_tree = treeid;
      lnid = sc_array_bsearch (ghost_layer, &tq,
                               p4est_quadrant_compare_piggy);
      exists = (lnid != -1) ? 2 : 0;
    }
  }
  else {
    /* q is in a neighboring tree */
    quad_contact[0] = (q->y < 0);
    quad_contact[1] = (q->x >= rh);
    quad_contact[2] = (q->y >= rh);
    quad_contact[3] = (q->x < 0);

    if ((quad_contact[0] || quad_contact[2]) &&
        (quad_contact[1] || quad_contact[3])) {
      /* Neighbor is across a tree corner */
      for (tree_corner = 0; tree_corner < 4; ++tree_corner) {
        if (quad_contact[(tree_corner + 3) % 4]
            && quad_contact[tree_corner]) {
          break;
        }
      }
      sc_array_init (cta, sizeof (p4est_corner_transform_t));
      p4est_find_corner_transform (conn, treeid, tree_corner, cta);

      sc_array_resize (exists_arr, cta->elem_count);

      num_ctrees = 0;
      for (ctreeidz = 0; ctreeidz < cta->elem_count; ++ctreeidz) {
        ct = sc_array_index (cta, ctreeidz);
        tqtreeid = ct->ntree;

        /* Don't use corner identification in the same tree */
        if (tqtreeid == treeid)
          continue;

        tq = *q;
        p4est_quadrant_transform_corner (&tq, (int) ct->ncorner, true);

        qproc = p4est_comm_find_owner (p4est, tqtreeid, &tq, rank);

        if (qproc == rank) {
          tqtree = p4est_array_index_topidx (p4est->trees, tqtreeid);
          lnid = sc_array_bsearch (&tqtree->quadrants, &tq,
                                   p4est_quadrant_compare);
          exists = (lnid != -1) ? 1 : 0;
        }
        else {
          P4EST_ASSERT (sizeof (long) >= sizeof (p4est_locidx_t));
          tq.p.piggy1.which_tree = tqtreeid;
          lnid = sc_array_bsearch (ghost_layer, &tq,
                                   p4est_quadrant_compare_piggy);
          exists = (lnid != -1) ? 1 : 0;
        }

        /* add the existance value */
        pexists = sc_array_index (exists_arr, num_ctrees);
        *pexists = exists;
        ++num_ctrees;
      }
      sc_array_resize (exists_arr, num_ctrees);
      sc_array_reset (cta);
      exists = -1;
    }
    else {
      /* Neighbor is across a tree face */
      for (face = 0; face < 4; ++face) {
        if (quad_contact[face] &&
            ((conn->tree_to_tree[4 * treeid + face] != treeid)
             || ((int) conn->tree_to_face[4 * treeid + face] != face))) {
          tqtreeid = conn->tree_to_tree[4 * treeid + face];
          break;
        }
      }
      if (face == 4) {
        /* this quadrant ran across a face with no neighbor */
        exists = 0;
      }
      else {
        /* transform the neighbor into the other tree's
         * coordinates
         */
        tempq = *q;
        transform = p4est_find_face_transform (conn, treeid, face);
        p4est_quadrant_translate_face (&tempq, face);
        p4est_quadrant_transform_face (&tempq, &tq, transform);

        qproc = p4est_comm_find_owner (p4est, tqtreeid, &tq, rank);

        if (qproc == rank) {
          tqtree = p4est_array_index_topidx (p4est->trees, tqtreeid);

          lnid = sc_array_bsearch (&tqtree->quadrants, &tq,
                                   p4est_quadrant_compare);
          exists = (lnid != -1) ? 1 : 0;
        }
        else {
          /* off processor so search in the ghost layer */

          P4EST_ASSERT (sizeof (long) >= sizeof (p4est_locidx_t));
          tq.p.piggy1.which_tree = tqtreeid;
          lnid = sc_array_bsearch (ghost_layer, &tq,
                                   p4est_quadrant_compare_piggy);
          exists = (lnid != -1) ? 1 : 0;
        }
      }
    }
  }

  return exists;
}

/** Checks if a quadrant's face is on the boundary of the forest.
 *
 * \param [in] p4est  The forest in which to search for \a q
 * \param [in] treeid The tree id for which \a q belongs.
 * \param [in] q      The quadrant that is in questioned.
 * \param [in] face   The face of quadrant that is in question.
 *
 * \return 1 if the quadrant's face is on the boundary of the forest and
 *         0 otherwise.
 */
static              bool
p4est_quadrant_on_face_boundary (p4est_t * p4est, p4est_topidx_t treeid,
                                 p4est_quadrant_t * q, int face)
{
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);
  p4est_connectivity_t *conn = p4est->connectivity;
  bool                on_boundary;
  bool                tree_boundary =
    ((conn->tree_to_tree[4 * treeid + face] == treeid) &&
     ((int) conn->tree_to_face[4 * treeid + face] == face));

  P4EST_ASSERT (p4est_quadrant_is_inside_root (q));

  switch (face) {
  case 0:
    on_boundary = (tree_boundary && (q->y == 0));
    break;
  case 1:
    on_boundary = (tree_boundary && (q->x == rh - qh));
    break;
  case 2:
    on_boundary = (tree_boundary && (q->y == rh - qh));
    break;
  case 3:
    on_boundary = (tree_boundary && (q->x == 0));
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }

  return on_boundary;
}

/** Checks if a quadrant's corner is on the boundary of the forest.
 *
 * This means that the quadrant's tree doesn't have any non face neighbors.
 *
 * \param [in] p4est  The forest in which to search for \a q
 * \param [in] treeid The tree id for which \a q belongs.
 * \param [in] q      The quadrant that is in questioned.
 * \param [in] corner The corner of quadrant that is in question.
 *
 * \return 1 if the quadrant's corner is on the boundary of the forest and
 *         0 otherwise.
 */
static              bool
p4est_quadrant_on_corner_boundary (p4est_t * p4est, p4est_topidx_t treeid,
                                   p4est_quadrant_t * q, int corner)
{
  int                 rcorner = p4est_corner_to_zorder[corner];
  bool                on_boundary = false;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  p4est_topidx_t      corner_trees;
  p4est_topidx_t      vertex, ntree1, ntree2, ctree, ntree;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_qcoord_t      cx = q->x, cy = q->y;
  p4est_qcoord_t      qh = P4EST_QUADRANT_LEN (q->level);

  P4EST_ASSERT (p4est_quadrant_is_inside_root (q));

  switch (corner) {
  case 0:
    cx -= qh;
    cy -= qh;
    break;
  case 1:
    cx += qh;
    cy -= qh;
    break;
  case 2:
    cx -= qh;
    cy += qh;
    break;
  case 3:
    cx += qh;
    cy += qh;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }

  if (((cy < 0) || (cy >= rh)) && ((cx >= rh) || (cx < 0))) {

    vertex = conn->tree_to_vertex[4 * treeid + rcorner];
    corner_trees = conn->vtt_offset[vertex + 1] - conn->vtt_offset[vertex];

    ntree1 = conn->tree_to_tree[4 * treeid + (rcorner + 3) % 4];
    ntree2 = conn->tree_to_tree[4 * treeid + rcorner];

    on_boundary = true;
    for (ctree = 0; ctree < corner_trees; ++ctree) {
      ntree = conn->vertex_to_tree[conn->vtt_offset[vertex] + ctree];
      if (!(ntree == treeid || ntree == ntree1 || ntree == ntree2)) {
        on_boundary = false;
        break;
      }
    }
  }
  else {
    on_boundary = (cy < 0) || (cy >= rh) || (cx >= rh) || (cx < 0);
  }

  return on_boundary;
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
 * \param [out] n0ur   Filled with smallest quadrant that fits in the
 *                     upper right corner of \a n0.
 */
static void
p4est_quadrant_get_half_corner_neighbors (p4est_quadrant_t * q, int corner,
                                          p4est_quadrant_t * n0,
                                          p4est_quadrant_t * n0ur)
{
  p4est_qcoord_t      th = P4EST_QUADRANT_LEN (P4EST_MAXLEVEL);
  p4est_qcoord_t      qh = P4EST_QUADRANT_LEN (q->level);
  p4est_qcoord_t      qh_2 = P4EST_QUADRANT_LEN (q->level + 1);

  *n0 = *q;

  P4EST_ASSERT (n0->level != P4EST_MAXLEVEL);
  n0->level += 1;

  switch (corner) {
  case 0:
    n0->x -= qh_2;
    n0->y -= qh_2;
    break;
  case 1:
    n0->x += qh;
    n0->y -= qh_2;
    break;
  case 2:
    n0->x -= qh_2;
    n0->y += qh;
    break;
  case 3:
    n0->x += qh;
    n0->y += qh;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }

  n0ur->x = n0->x + qh_2 - th;
  n0ur->y = n0->y + qh_2 - th;
  n0ur->level = P4EST_MAXLEVEL;
}

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
 * \param [out] n0     Filled with the first possible face neighbor, which is
 *                     half of the size assuming the 2-1 constaint.
 * \param [out] n0ur   Filled with smallest quadrant that fits in the
 *                     upper right corner of \a n0.
 * \param [out] n1     Filled with the second possible face neighbor, which is
 *                     half of the size assuming the 2-1 constaint.
 * \param [out] n1ur   Filled with smallest quadrant that fits in the
 *                     upper right corner of \a n1.
 *
 */
static void
p4est_quadrant_get_half_face_neighbors (p4est_quadrant_t * q, int face,
                                        p4est_quadrant_t * n0,
                                        p4est_quadrant_t * n0ur,
                                        p4est_quadrant_t * n1,
                                        p4est_quadrant_t * n1ur)
{
  p4est_qcoord_t      th = P4EST_QUADRANT_LEN (P4EST_MAXLEVEL);
  p4est_qcoord_t      qh = P4EST_QUADRANT_LEN (q->level);
  p4est_qcoord_t      qh_2 = P4EST_QUADRANT_LEN (q->level + 1);

  *n0 = *q;
  *n1 = *q;

  P4EST_ASSERT (n0->level != P4EST_MAXLEVEL);
  P4EST_ASSERT (n1->level != P4EST_MAXLEVEL);
  n0->level += 1;
  n1->level += 1;

  switch (face) {
  case 0:
    n0->y -= qh_2;

    n1->x += qh_2;
    n1->y -= qh_2;
    break;
  case 1:
    n0->x += qh;

    n1->x += qh;
    n1->y += qh_2;
    break;
  case 2:
    n0->y += qh;

    n1->x += qh_2;
    n1->y += qh;
    break;
  case 3:
    n0->x -= qh_2;

    n1->x -= qh_2;
    n1->y += qh_2;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }

  n0ur->x = n0->x + qh_2 - th;
  n0ur->y = n0->y + qh_2 - th;
  n0ur->level = P4EST_MAXLEVEL;

  n1ur->x = n1->x + qh_2 - th;
  n1ur->y = n1->y + qh_2 - th;
  n1ur->level = P4EST_MAXLEVEL;
}

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
                        p4est_quadrant_t * q)
{
  bool                add_to_proc = true;
  p4est_quadrant_t   *qold, *qnew;

  /* Check to see if the quadrant already exists in the array */
  if (buf->elem_count > 0) {
    qold = sc_array_index (buf, buf->elem_count - 1);
    if (p4est_quadrant_compare_piggy (q, qold) == 0) {
      add_to_proc = false;
    }
  }

  if (add_to_proc) {
    qnew = sc_array_push (buf);
    *qnew = *q;

    /* Cram the tree id into the user_data pointer */
    P4EST_ASSERT (sizeof (long) >= sizeof (p4est_locidx_t));
    qnew->p.piggy1.which_tree = treeid;
  }
}

bool
p4est_is_balanced (p4est_t * p4est)
{
  int                 face, corner;
  int                 qcid;
  int                 e0, e1, e2, e3;
  int                *pe0, *pe1, *pe2;
  bool                is_balanced = true;
  bool                is_face_balanced = false;
  bool                is_corner_balanced = false;
  size_t              cez;
  p4est_locidx_t      li, lj;
  p4est_locidx_t      num_quadrants;
  p4est_locidx_t      first_local_tree = p4est->first_local_tree;
  p4est_locidx_t      last_local_tree = p4est->last_local_tree;
  p4est_quadrant_t   *q;
  p4est_quadrant_t    n0, n1, n2, n3;
  p4est_tree_t       *tree;
  sc_array_t          ghost_layer;
  sc_array_t         *quadrants;
  sc_array_t          e0_a, e1_a, e2_a;

  P4EST_QUADRANT_INIT (&n0);
  P4EST_QUADRANT_INIT (&n1);
  P4EST_QUADRANT_INIT (&n2);
  P4EST_QUADRANT_INIT (&n3);

  sc_array_init (&ghost_layer, sizeof (p4est_quadrant_t));
  sc_array_init (&e0_a, sizeof (int));
  sc_array_init (&e1_a, sizeof (int));
  sc_array_init (&e2_a, sizeof (int));

  p4est_build_ghost_layer (p4est, &ghost_layer);

  /* loop over all local trees */
  for (lj = first_local_tree; lj <= last_local_tree; ++lj) {
    tree = sc_array_index (p4est->trees, lj);
    quadrants = &tree->quadrants;
    num_quadrants = quadrants->elem_count;

    /* Find the neighboring processors of each quadrant */
    for (li = 0; li < num_quadrants; ++li) {
      q = sc_array_index (quadrants, li);
      qcid = p4est_quadrant_child_id (q);

      /* Find Face Neighbors */
      for (face = 0; face < 4; ++face) {
        is_face_balanced = false;

        /* If q is at a boundary then it is automatically balanced */
        if (p4est_quadrant_on_face_boundary (p4est, lj, q, face)) {
          is_face_balanced = !is_face_balanced;
        }
        else {
          p4est_quadrant_get_possible_face_neighbors (q, face,
                                                      &n0, &n1, &n2, &n3);

          e0 = p4est_quadrant_exists (p4est, &ghost_layer, lj, &n0, &e0_a);
          e1 = p4est_quadrant_exists (p4est, &ghost_layer, lj, &n1, &e0_a);
          e2 = p4est_quadrant_exists (p4est, &ghost_layer, lj, &n2, &e0_a);
          e3 = p4est_quadrant_exists (p4est, &ghost_layer, lj, &n3, &e0_a);

          P4EST_ASSERT (e0 != -1 && e1 != -1 && e2 != -1 && e3 != -1);

          /* Check existance of the half size quadrants */
          if (e0 && e1)
            is_face_balanced = !is_face_balanced;

          /* Check existance of the same size quadrant */
          if (e2)
            is_face_balanced = !is_face_balanced;

          /* Check existance of the double size quadrant */
          if (e3)
            is_face_balanced = !is_face_balanced;
        }

        is_balanced = (is_balanced && is_face_balanced);
      }

      /* Find Corner Neighbors */
      for (corner = 0; corner < 4; ++corner) {
        is_corner_balanced = false;

        /* If q is at a boundary then it is automatically balanced */
        if (corner != qcid) {
          is_corner_balanced = !is_corner_balanced;
        }
        else if (p4est_quadrant_on_corner_boundary (p4est, lj, q, corner)) {
          is_corner_balanced = !is_corner_balanced;
        }
        else {
          p4est_quadrant_get_possible_corner_neighbors (q, corner,
                                                        &n0, &n1, &n2);

          e0 = p4est_quadrant_exists (p4est, &ghost_layer, lj, &n0, &e0_a);
          e1 = p4est_quadrant_exists (p4est, &ghost_layer, lj, &n1, &e1_a);
          e2 = p4est_quadrant_exists (p4est, &ghost_layer, lj, &n2, &e2_a);

          if (e0 == -1 && e1 == -1 && e2 == -1) {
            /* We are dealing with a tree neighbor */
            P4EST_ASSERT (e0_a.elem_count == e1_a.elem_count
                          && e0_a.elem_count == e2_a.elem_count);
            is_corner_balanced = true;
            for (cez = 0; cez < e0_a.elem_count; ++cez) {
              pe0 = sc_array_index (&e0_a, cez);
              pe1 = sc_array_index (&e1_a, cez);
              pe2 = sc_array_index (&e2_a, cez);

              /* check to see if one and only one of the neighboring quadrants
               * exists
               */
              e3 = *pe0 ^ *pe1 ^ *pe2;

              /* Corner is only balanced if all of the cornering quadrants are
               * the right size
               */
              is_corner_balanced = is_corner_balanced & e3;
            }
          }
          else {
            P4EST_ASSERT (e0 != -1 && e1 != -1 && e2 != -1);

            /* Check existance of the half size quadrants */
            if (e0)
              is_corner_balanced = !is_corner_balanced;

            /* Check existance of the same size quadrant */
            if (e1)
              is_corner_balanced = !is_corner_balanced;

            /* Check existance of the double size quadrant */
            if (e2)
              is_corner_balanced = !is_corner_balanced;
          }
        }
        is_balanced = (is_balanced && is_corner_balanced);
      }
    }
  }

  sc_array_reset (&ghost_layer);
  sc_array_reset (&e0_a);
  sc_array_reset (&e1_a);
  sc_array_reset (&e2_a);

  return is_balanced;
}

void
p4est_build_ghost_layer (p4est_t * p4est, sc_array_t * ghost_layer)
{
#ifdef P4EST_MPI
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  int                 face, corner, treecorner;
  int                 i;
  int                 n0_proc, n0ur_proc, n1_proc, n1ur_proc;
  int                 num_peers, peer, peer_proc;
  int                 mpiret;
  size_t              pz;
  p4est_topidx_t      nt;
  p4est_topidx_t      first_local_tree = p4est->first_local_tree;
  p4est_topidx_t      last_local_tree = p4est->last_local_tree;
  p4est_locidx_t      li;
  p4est_locidx_t      num_quadrants;
  p4est_locidx_t      num_ghosts;
  p4est_locidx_t     *peer_counts;
  p4est_locidx_t      count;
  p4est_locidx_t      ghost_offset;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  p4est_quadrant_t    n0, n0ur, n1, n1ur;
  sc_array_t          send_bufs;
  sc_array_t          procs, urprocs;
  sc_array_t         *buf;
  sc_array_t         *quadrants;
  MPI_Comm            comm = p4est->mpicomm;
  MPI_Request        *recv_request, *send_request;
  MPI_Status         *recv_status, *send_status;
  MPI_Request        *recv_load_request, *send_load_request;
  MPI_Status         *recv_load_status, *send_load_status;

  P4EST_QUADRANT_INIT (&n0);
  P4EST_QUADRANT_INIT (&n0ur);
  P4EST_QUADRANT_INIT (&n1);
  P4EST_QUADRANT_INIT (&n1ur);

  sc_array_init (&procs, sizeof (int));
  sc_array_init (&urprocs, sizeof (int));

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
    num_quadrants = quadrants->elem_count;

    /* Find the neighboring processors of each quadrant */
    for (li = 0; li < num_quadrants; ++li) {
      q = sc_array_index (quadrants, li);

      /* Find Face Neighbors */
      for (face = 0; face < 4; ++face) {
        p4est_quadrant_get_half_face_neighbors (q, face, &n0, &n0ur,
                                                &n1, &n1ur);

        n0_proc = p4est_quadrant_find_owner (p4est, nt, &n0);
        n1_proc = p4est_quadrant_find_owner (p4est, nt, &n1);

        n0ur_proc = p4est_quadrant_find_owner (p4est, nt, &n0ur);
        n1ur_proc = p4est_quadrant_find_owner (p4est, nt, &n1ur);

        /* Note that we will always check this because it is cheap
         * and prevents deadlocks
         */
        SC_CHECK_ABORT (n0_proc == n0ur_proc, "Non reciprocal communication");
        SC_CHECK_ABORT (n1_proc == n1ur_proc, "Non reciprocal communication");

        if (n0_proc != rank && n0_proc >= 0) {
          buf = sc_array_index_int (&send_bufs, n0_proc);
          p4est_add_ghost_to_buf (buf, nt, q);
        }

        if (n1_proc != rank && n1_proc >= 0 && n0_proc != n1_proc) {
          buf = sc_array_index_int (&send_bufs, n1_proc);
          p4est_add_ghost_to_buf (buf, nt, q);
        }
      }

      /* Find Corner Neighbors */
      for (corner = 0; corner < 4; ++corner) {
        p4est_quadrant_get_half_corner_neighbors (q, corner, &n0, &n0ur);

        /* Check to see if we are a tree corner neighbor */
        if (((n0.y < 0) || (n0.y >= rh)) && ((n0.x >= rh) || (n0.x < 0))) {
          /* We have to loop over multiple neighbors if we are at
           * a tree corner
           */
          treecorner = p4est_corner_to_zorder[corner];
          p4est_quadrant_find_tree_corner_owners (p4est, nt, treecorner,
                                                  &n0, &procs);
          p4est_quadrant_find_tree_corner_owners (p4est, nt, treecorner,
                                                  &n0ur, &urprocs);
          for (pz = 0; pz < procs.elem_count; ++pz) {
            n0_proc = *((int *) sc_array_index (&procs, pz));
            n0ur_proc = *((int *) sc_array_index (&urprocs, pz));

            /* Note that we will always check this because it is cheap
             * and prevents deadlocks
             */
            SC_CHECK_ABORT (n0_proc == n0ur_proc,
                            "Non reciprocal communication");

            if (n0_proc != rank) {
              buf = sc_array_index_int (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, nt, q);
            }
          }
        }
        else {
          /* We are not at a tree corner so we only have one
           * corner neighbor
           */
          n0_proc = p4est_quadrant_find_owner (p4est, nt, &n0);
          n0ur_proc = p4est_quadrant_find_owner (p4est, nt, &n0ur);
          /* Note that we will always check this because it is cheap
           * and prevents deadlocks
           */
          SC_CHECK_ABORT (n0_proc == n0ur_proc,
                          "Non reciprocal communication");

          if (n0_proc != rank && n0_proc >= 0) {
            buf = sc_array_index_int (&send_bufs, n0_proc);
            p4est_add_ghost_to_buf (buf, nt, q);
          }
        }
      }
    }
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

  peer_counts = P4EST_ALLOC (p4est_locidx_t, num_peers);

  recv_load_request = recv_request + num_peers;
  recv_load_status = recv_status + num_peers;

  send_load_request = send_request + num_peers;
  send_load_status = send_status + num_peers;

  /* Post receives for the counts of ghosts to be received */
  for (i = 0, peer = 0; i < num_procs; ++i) {
    buf = sc_array_index_int (&send_bufs, i);
    if (buf->elem_count > 0) {
      peer_proc = i;
      P4EST_LDEBUGF ("ghost layer post count receive from %d\n", peer_proc);
      mpiret = MPI_Irecv (peer_counts + peer, 1, P4EST_MPI_LOCIDX,
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
      P4EST_LDEBUGF ("ghost layer post count send to %d\n", peer_proc);
      count = (p4est_locidx_t) buf->elem_count;
      mpiret = MPI_Isend (&count, 1, P4EST_MPI_LOCIDX, peer_proc,
                          P4EST_COMM_GHOST_COUNT, comm, send_request + peer);
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

  /* Count Ghosts */
  for (peer = 0, num_ghosts = 0; peer < num_peers; ++peer) {
    P4EST_ASSERT (peer_counts[peer] > 0);
    num_ghosts += peer_counts[peer];
  }

  /* Allocate space for the ghosts */
  sc_array_resize (ghost_layer, num_ghosts);

  /* Post receives for the ghosts */
  for (i = 0, peer = 0, ghost_offset = 0; i < num_procs; ++i) {
    buf = sc_array_index_int (&send_bufs, i);
    if (buf->elem_count > 0) {
      peer_proc = i;
      P4EST_LDEBUGF
        ("ghost layer post ghost receive %lld quadrants from %d\n",
         (long long) peer_counts[peer], peer_proc);
      mpiret =
        MPI_Irecv (ghost_layer->array +
                   ghost_offset * sizeof (p4est_quadrant_t),
                   (int) (peer_counts[peer] * sizeof (p4est_quadrant_t)),
                   MPI_CHAR, peer_proc, P4EST_COMM_GHOST_LOAD, comm,
                   recv_load_request + peer);
      SC_CHECK_MPI (mpiret);
      ghost_offset += peer_counts[peer];
      ++peer;
    }
  }

  /* Send the ghosts */
  for (i = 0, peer = 0; i < num_procs; ++i) {
    buf = sc_array_index_int (&send_bufs, i);
    if (buf->elem_count > 0) {
      peer_proc = i;
      count = (p4est_locidx_t) buf->elem_count;
      P4EST_LDEBUGF ("ghost layer post ghost send %lld quadrants to %d\n",
                     (long long) count, peer_proc);
      mpiret =
        MPI_Isend (&buf->array, (int) (count * sizeof (p4est_quadrant_t)),
                   MPI_CHAR, peer_proc, P4EST_COMM_GHOST_LOAD, comm,
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
  P4EST_FREE (peer_counts);

#ifdef P4EST_DEBUG
  for (i = 0; i < num_peers; ++i) {
    P4EST_ASSERT (recv_load_request[i] == MPI_REQUEST_NULL);
  }
  for (i = 0; i < num_peers; ++i) {
    P4EST_ASSERT (send_load_request[i] == MPI_REQUEST_NULL);
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
  sc_array_reset (&procs);
  sc_array_reset (&urprocs);

#else
  /* If we are not running with mpi then we don't need to do anything */
  sc_array_reset (ghost_layer);
#endif
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

void
p4est_possible_node_neighbor (const p4est_quadrant_t * q, int node,
                              int nnum, int neighbor_rlev,
                              p4est_quadrant_t * neighbor, int *neighbor_node)
{
  int                 nnode;
  const int           nlevel = (int) q->level + neighbor_rlev;
  const p4est_qcoord_t qh =
    (p4est_qcoord_t) (1 << (P4EST_MAXLEVEL - q->level));
  const p4est_qcoord_t nh = (p4est_qcoord_t) (1 << (P4EST_MAXLEVEL - nlevel));
  const p4est_qcoord_t qx = q->x;
  const p4est_qcoord_t qy = q->y;
  p4est_qcoord_t      cornerx, cornery;
  p4est_quadrant_t    n;
#ifdef P4EST_DEBUG
  int                 qcid;
#endif

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (-1 <= neighbor_rlev && neighbor_rlev <= 1);
  P4EST_ASSERT (0 <= nlevel && nlevel <= P4EST_MAXLEVEL);

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

/* EOF p4est_mesh.c */
