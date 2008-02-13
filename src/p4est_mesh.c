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
  int                 owner = -1;
  p4est_connectivity_t *conn = p4est->connectivity;
  int                 quad_contact[4];
  int                 face_contact[4];
  int                 face, transform;
  p4est_locidx_t      ntreeid;
  p4est_quadrant_t    tmpq = *q, nq;

  P4EST_QUADRANT_INIT (&nq);

  if (p4est_quadrant_is_inside (q)) {
    owner = p4est_comm_find_owner (p4est, treeid, q, rank);
  }
  else {
    quad_contact[0] = (q->y < 0);
    quad_contact[1] = (q->x >= P4EST_ROOT_LEN);
    quad_contact[2] = (q->y >= P4EST_ROOT_LEN);
    quad_contact[3] = (q->x < 0);

    /* Make sure we are not a tree corner */
    P4EST_ASSERT (!((quad_contact[0] || quad_contact[2]) &&
                    (quad_contact[1] || quad_contact[3])));

    for (face = 0; face < 4; ++face) {
      if (quad_contact[face]
          && (conn->tree_to_tree[4 * treeid + face] != treeid
              || (conn->tree_to_face[4 * treeid + face] != face))) {
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
      p4est_quadrant_translate (&tmpq, face);
      p4est_quadrant_transform (&tmpq, &nq, transform);

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
                                        p4est_locidx_t treeid,
                                        int treecorner,
                                        p4est_quadrant_t * q,
                                        p4est_array_t * q_procs)
{
  const int           rank = p4est->mpirank;
  int                 zcorner, cproc, *proc;
  p4est_locidx_t      ctree;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_array_t       corner_info;
  p4est_array_init (&corner_info, sizeof (p4est_corner_info_t));
  p4est_find_corner_info (conn, treeid, treecorner, &corner_info);
  p4est_corner_info_t *ci;
  p4est_locidx_t      ctreeid;
  p4est_quadrant_t    cq;

  P4EST_QUADRANT_INIT (&cq);

  P4EST_ASSERT (((q->y < 0) || (q->y >= P4EST_ROOT_LEN)) &&
                ((q->x >= P4EST_ROOT_LEN) || (q->x < 0)));

  p4est_array_resize (q_procs, 0);

  for (ctree = 0; ctree < corner_info.elem_count; ++ctree) {
    ci = p4est_array_index (&corner_info, ctree);
    ctreeid = ci->ntree;

    if (ctreeid == treeid)
      continue;

    zcorner = p4est_corner_to_zorder[ci->ncorner];
    cq = *q;
    p4est_quadrant_corner (&cq, zcorner, 1);

    cproc = p4est_comm_find_owner (p4est, ctreeid, &cq, rank);

    p4est_array_resize (q_procs, ctree + 1);
    proc = p4est_array_index (q_procs, ctree);
    *proc = cproc;
  }
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
    P4EST_ASSERT_NOT_REACHED ();
    break;
  }

  n0ur->x = n0->x - th;
  n0ur->y = n0->y - th;
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
    P4EST_ASSERT_NOT_REACHED ();
    break;
  }

  n0ur->x = n0->x - th;
  n0ur->y = n0->y - th;
  n0ur->level = P4EST_MAXLEVEL;

  n1ur->x = n1->x - th;
  n1ur->y = n1->y - th;
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
p4est_add_ghost_to_buf (p4est_array_t * buf, p4est_locidx_t treeid,
                        p4est_quadrant_t * q)
{
  int                 add_to_proc = 1;
  p4est_quadrant_t   *qold, *qnew;

  /* Cram the tree id into the user_data pointer */
  P4EST_ASSERT (sizeof (long) >= sizeof (p4est_locidx_t));
  q->user_data = (void *) (long) treeid;

  /* Check to see if the quadrant already exists in the array */
  if (buf->elem_count > 0) {
    qold = p4est_array_index (buf, buf->elem_count - 1);
    if (p4est_quadrant_compare_piggy (q, qold) == 0) {
      add_to_proc = 0;
    }
  }

  if (add_to_proc) {
    p4est_array_resize (buf, buf->elem_count + 1);
    qnew = p4est_array_index (buf, buf->elem_count - 1);
    *qnew = *q;
  }
}

int
p4est_is_balanced (p4est_t * p4est)
{
  int                 is_balanced = 1;

  return is_balanced;
}

void
p4est_build_ghost_layer (p4est_t * p4est, p4est_array_t * ghost_layer)
{
#ifdef HAVE_MPI
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  int                 face, corner, rlev, tree_corner;
  int                 i, p;
  p4est_locidx_t      li, lj;
  p4est_array_t      *trees = p4est->trees;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_array_t      *quadrants, corner_info;
  p4est_locidx_t      num_trees = conn->num_trees;
  p4est_locidx_t      Ncells = p4est->local_num_quadrants;
  p4est_locidx_t      first_local_tree = p4est->first_local_tree;
  p4est_locidx_t      last_local_tree = p4est->last_local_tree;
  p4est_locidx_t      num_quadrants;
  p4est_locidx_t      num_ghosts;
  p4est_locidx_t     *peer_counts;
  p4est_locidx_t     *peer_offsets;
  p4est_array_t       send_bufs;
  p4est_array_t       procs, urprocs;
  p4est_array_t      *buf;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  p4est_quadrant_t    n0, n0ur, n1, n1ur;
  int                 n0_proc, n0ur_proc, n1_proc, n1ur_proc;
  int                 num_peers, peer, peer_proc;
  int                 mpiret;
  MPI_Comm            comm = p4est->mpicomm;
  MPI_Request        *recv_request, *send_request;
  MPI_Status         *recv_status, *send_status;

  P4EST_QUADRANT_INIT (&n0);
  P4EST_QUADRANT_INIT (&n0ur);
  P4EST_QUADRANT_INIT (&n1);
  P4EST_QUADRANT_INIT (&n1ur);

  p4est_array_init (&procs, sizeof (int));
  p4est_array_init (&urprocs, sizeof (int));

  /* allocate empty send buffers */
  p4est_array_init (&send_bufs, sizeof (p4est_array_t));
  p4est_array_resize (&send_bufs, num_procs);
  for (i = 0; i < num_procs; ++i) {
    buf = p4est_array_index (&send_bufs, i);
    p4est_array_init (buf, sizeof (p4est_quadrant_t));
  }

  /* loop over all local trees */
  for (lj = first_local_tree; lj <= last_local_tree; ++lj) {
    tree = p4est_array_index (p4est->trees, lj);
    quadrants = &tree->quadrants;
    num_quadrants = quadrants->elem_count;

    /* Find the neighboring processors of each quadrant */
    for (li = 0; li < num_quadrants; ++li) {
      q = p4est_array_index (quadrants, li);

      /* Find Face Neighbors */
      for (face = 0; face < 4; ++face) {
        p4est_quadrant_get_half_face_neighbors (q, face, &n0, &n0ur,
                                                &n1, &n1ur);

        n0_proc = p4est_quadrant_find_owner (p4est, lj, &n0);
        n1_proc = p4est_quadrant_find_owner (p4est, lj, &n1);

        n0ur_proc = p4est_quadrant_find_owner (p4est, lj, &n0ur);
        n1ur_proc = p4est_quadrant_find_owner (p4est, lj, &n1ur);

        /* Note that we will always check this because it is cheap
         * and prevents deadlocks
         */
        P4EST_CHECK_ABORT (n0_proc == n0ur_proc,
                           "Non reciprocal communication");
        P4EST_CHECK_ABORT (n1_proc == n1ur_proc,
                           "Non reciprocal communication");

        if (n0_proc != rank && n0_proc >= 0) {
          buf = p4est_array_index (&send_bufs, n0_proc);
          p4est_add_ghost_to_buf (buf, lj, q);
        }

        if (n1_proc != rank && n1_proc >= 0 && n0_proc != n1_proc) {
          buf = p4est_array_index (&send_bufs, n1_proc);
          p4est_add_ghost_to_buf (buf, lj, q);
        }
      }

      /* Find Corner Neighbors */
      for (corner = 0; corner < 4; ++corner) {
        p4est_quadrant_get_half_corner_neighbors (q, corner, &n0, &n0ur);

        /* Check to see if we are a tree corner neighbor */
        if (((n0.y < 0) || (n0.y >= P4EST_ROOT_LEN)) &&
            ((n0.x >= P4EST_ROOT_LEN) || (n0.x < 0))) {
          /* We have to loop over multiple neighbors if we are at
           * a tree corner
           */
          p4est_quadrant_find_tree_corner_owners (p4est, lj, corner,
                                                  &n0, &procs);
          p4est_quadrant_find_tree_corner_owners (p4est, lj, corner,
                                                  &n0ur, &urprocs);
          for (p = 0; p < procs.elem_count; ++p) {
            n0_proc = *((int *) p4est_array_index (&procs, p));
            n0ur_proc = *((int *) p4est_array_index (&urprocs, p));

            /* Note that we will always check this because it is cheap
             * and prevents deadlocks
             */
            P4EST_CHECK_ABORT (n0_proc == n0ur_proc,
                               "Non reciprocal communication");

            if (n0_proc != rank) {
              buf = p4est_array_index (&send_bufs, n0_proc);
              p4est_add_ghost_to_buf (buf, lj, q);
            }
          }
        }
        else {
          /* We are not at a tree corner so we only have one
           * corner neighbor
           */
          n0_proc = p4est_quadrant_find_owner (p4est, lj, &n0);
          n0ur_proc = p4est_quadrant_find_owner (p4est, lj, &n0ur);
          /* Note that we will always check this because it is cheap
           * and prevents deadlocks
           */
          P4EST_CHECK_ABORT (n0_proc == n0ur_proc,
                             "Non reciprocal communication");

          if (n0_proc != rank && n0_proc >= 0) {
            buf = p4est_array_index (&send_bufs, n0_proc);
            p4est_add_ghost_to_buf (buf, lj, q);
          }
        }
      }
    }
  }

  /* Count the number of peers that I send to and receive from */
  for (i = 0, num_peers = 0; i < num_procs; ++i) {
    buf = p4est_array_index (&send_bufs, i);
    if (buf->elem_count > 0)
      ++num_peers;
  }

  recv_request = P4EST_ALLOC (MPI_Request, num_peers);
  P4EST_CHECK_ALLOC (recv_request);
  recv_status = P4EST_ALLOC (MPI_Status, num_peers);
  P4EST_CHECK_ALLOC (recv_status);

  send_request = P4EST_ALLOC (MPI_Request, num_peers);
  P4EST_CHECK_ALLOC (send_request);
  send_status = P4EST_ALLOC (MPI_Status, num_peers);
  P4EST_CHECK_ALLOC (send_status);

  peer_counts = P4EST_ALLOC (p4est_locidx_t, num_peers);
  P4EST_CHECK_ALLOC (peer_counts);

  /* Post receives for the counts of ghosts to be received */
  for (i = 0, peer = 0; i < num_procs; ++i) {
    buf = p4est_array_index (&send_bufs, i);
    if (buf->elem_count > 0) {
      peer_proc = i;
      P4EST_DEBUGF ("ghost layer post count receive from %d\n", peer_proc);
      mpiret = MPI_Irecv (peer_counts + peer, 1, P4EST_MPI_LOCIDX,
                          peer_proc, P4EST_COMM_GHOST_COUNT, comm,
                          recv_request + peer);
      P4EST_CHECK_MPI (mpiret);
      ++peer;
    }
  }

  /* Send the counts of ghosts that are going to be sent */
  for (i = 0, peer = 0; i < num_procs; ++i) {
    buf = p4est_array_index (&send_bufs, i);
    if (buf->elem_count > 0) {
      peer_proc = i;
      P4EST_DEBUGF ("ghost layer post count sent to %d\n", peer_proc);
      mpiret = MPI_Isend (&buf->elem_count, 1, P4EST_MPI_LOCIDX, peer_proc,
                          P4EST_COMM_GHOST_COUNT, comm, send_request + peer);
      P4EST_CHECK_MPI (mpiret);
      ++peer;
    }
  }

  /* Wait for the counts */
  mpiret = MPI_Waitall (peer_proc, recv_request, recv_status);
  P4EST_CHECK_MPI (mpiret);

  /* Allocate space for the ghosts */

  /* Post receives for the ghosts */

  /* Send the ghosts */

  /* Clean up */
  P4EST_FREE (peer_counts);

#ifdef P4EST_HAVE_DEBUG
  for (i = 0; i < num_peers; ++i) {
    P4EST_ASSERT (recv_request[i] == MPI_REQUEST_NULL);
  }
  for (i = 0; i < num_peers; ++i) {
    P4EST_ASSERT (send_request[i] == MPI_REQUEST_NULL);
  }
#endif
  P4EST_FREE (recv_request);
  P4EST_FREE (recv_status);
  P4EST_FREE (send_request);
  P4EST_FREE (send_status);

  for (i = 0; i < num_procs; ++i) {
    buf = p4est_array_index (&send_bufs, i);
    p4est_array_reset (buf);
  }
  p4est_array_reset (&send_bufs);
  p4est_array_reset (&procs);
  p4est_array_reset (&urprocs);

#else
  /* If we are not running with mpi then we don't need to do anything */
  p4est_array_reset (ghost_layer);
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
#ifdef P4EST_HAVE_DEBUG
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

#ifdef P4EST_HAVE_DEBUG
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
