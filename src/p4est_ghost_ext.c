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
#include <p4est_ghost_ext.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#else
#include <p8est_ghost_ext.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#endif

#ifdef P4EST_ENABLE_MPI

void
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

void
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
    qnew = p4est_quadrant_array_push (m->mirrors);
    *qnew = *q;

    /* cram the tree id and the local number into the user_data pointer */
    qnew->p.piggy3.which_tree = treeid;
    qnew->p.piggy3.local_num = number;

    m->known = 1;
  }

  buf = p4est_ghost_array_index_int (m->send_bufs, p);
  if (p4est_ghost_add_to_buf (buf, treeid, number, q)) {
    P4EST_ASSERT (m->mirrors->elem_count > 0);

    num = (p4est_locidx_t *) sc_array_push (m->offsets_by_proc + p);
    *num = (p4est_locidx_t) (m->mirrors->elem_count - 1);
    ++m->sum_all_procs;
  }
}

sc_array_t         *
p4est_ghost_array_index_int (sc_array_t * array, int i)
{
  return (sc_array_t *) sc_array_index_int (array, i);
}

int
p4est_ghost_add_to_buf (sc_array_t * buf, p4est_topidx_t treeid,
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

  qnew = p4est_quadrant_array_push (buf);
  *qnew = *q;

  /* Cram the tree id and the local number into the user_data pointer */
  qnew->p.piggy3.which_tree = treeid;
  qnew->p.piggy3.local_num = number;

  return 1;
}

/** Populate the mirror fields in the ghost layer with final data.
 * The elements in the temporary p4est_ghost_mirror_t structure are freed. */
void
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

void
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
      buf = p4est_ghost_array_index_int (send_bufs, n0_proc);
      p4est_ghost_add_to_buf (buf, t, local_num, q);
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
      buf = p4est_ghost_array_index_int (send_bufs, n0_proc);
      p4est_ghost_add_to_buf (buf, t, local_num, q);
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
      buf = p4est_ghost_array_index_int (send_bufs, proc);
      p4est_ghost_add_to_buf (buf, t, local_num, q);
#endif
      p4est_ghost_mirror_add (m, t, local_num, q, proc);
    }
  }
}

#endif /* P4EST_ENABLE_MPI */

size_t
p4est_ghost_tree_type (sc_array_t * array, size_t zindex, void *data)
{
  p4est_quadrant_t   *q;

  P4EST_ASSERT (array->elem_size == sizeof (p4est_quadrant_t));

  q = (p4est_quadrant_t *) sc_array_index (array, zindex);
  return (size_t) q->p.which_tree;
}
