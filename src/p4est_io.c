/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Copyright (C) 2012 Carsten Burstedde
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
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_io.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_io.h>
#endif
#include <sc_search.h>

sc_array_t         *
p4est_deflate_quadrants (p4est_t * p4est, sc_array_t ** data)
{
  const size_t        qsize = sizeof (p4est_qcoord_t);
  const size_t        dsize = p4est->data_size;
  size_t              qtreez, qz;
  sc_array_t         *qarr, *darr;
  p4est_topidx_t      tt;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  p4est_qcoord_t     *qap;
  char               *dap;

  qarr = sc_array_new_size (qsize,
                            (P4EST_DIM + 1) * p4est->local_num_quadrants);
  qap = (p4est_qcoord_t *) qarr->array;
  darr = NULL;
  dap = NULL;
  if (data != NULL) {
    P4EST_ASSERT (dsize > 0);
    darr = sc_array_new_size (dsize, p4est->local_num_quadrants);
    dap = darr->array;
  }
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    qtreez = tree->quadrants.elem_count;
    for (qz = 0; qz < qtreez; ++qz) {
      q = p4est_quadrant_array_index (&tree->quadrants, qz);
      *qap++ = q->x;
      *qap++ = q->y;
#ifdef P4_TO_P8
      *qap++ = q->z;
#endif
      *qap++ = (p4est_qcoord_t) q->level;
      if (data != NULL) {
        memcpy (dap, q->p.user_data, dsize);
        dap += dsize;
      }
    }
  }
  P4EST_ASSERT ((void *) qap ==
                qarr->array + qarr->elem_size * qarr->elem_count);
  if (data != NULL) {
    P4EST_ASSERT (dap == darr->array + darr->elem_size * darr->elem_count);
    *data = darr;
  }
  return qarr;
}

p4est_t            *
p4est_inflate (sc_MPI_Comm mpicomm, p4est_connectivity_t * connectivity,
               const p4est_gloidx_t * global_first_quadrant,
               const p4est_gloidx_t * pertree,
               sc_array_t * quadrants, sc_array_t * data, void *user_pointer)
{
  const p4est_gloidx_t *gfq;
  int                 i;
  int                 mpiret;
  int                 num_procs, rank;
  p4est_topidx_t      num_trees, jt;
  p4est_gloidx_t      gkey, gtreeskip, gtreeremain, gquadremain;
  p4est_t            *p4est;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
#ifdef P4EST_ENABLE_DEBUG
  int                 p;
#endif
  int8_t              ql, tml;
  size_t              dsize;
  size_t              gk1, gk2;
  size_t              qz, zqoffset, zqthistree;
  p4est_qcoord_t     *qap;
  char               *dap;

  P4EST_GLOBAL_PRODUCTION ("Into " P4EST_STRING "_inflate\n");
  p4est_log_indent_push ();

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));
  P4EST_ASSERT (global_first_quadrant != NULL);
  P4EST_ASSERT (pertree != NULL);
  P4EST_ASSERT (quadrants != NULL);
  P4EST_ASSERT (quadrants->elem_size == sizeof (p4est_qcoord_t));
  /* data may be NULL, in this case p4est->data_size will be 0 */
  /* user_pointer may be anything, we don't look at it */

  /* retrieve MPI information */
  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* assign some data members */
  p4est = P4EST_ALLOC_ZERO (p4est_t, 1);
  p4est->mpicomm = mpicomm;
  p4est->mpisize = num_procs;
  p4est->mpirank = rank;
  dsize = p4est->data_size = data == NULL ? 0 : data->elem_size;
  dap = (char *) (data == NULL ? NULL : data->array);
  qap = (p4est_locidx_t *) quadrants->array;
  p4est->user_pointer = user_pointer;
  p4est->connectivity = connectivity;
  num_trees = connectivity->num_trees;

  /* create global first quadrant offsets */
  gfq = p4est->global_first_quadrant =
    P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  memcpy (p4est->global_first_quadrant, global_first_quadrant,
          (num_procs + 1) * sizeof (p4est_gloidx_t));
#ifdef P4EST_ENABLE_DEBUG
  P4EST_ASSERT (gfq[0] == 0);
  for (p = 0; p < num_procs; ++p) {
    P4EST_ASSERT (gfq[p] <= gfq[p + 1]);
  }
  P4EST_ASSERT (pertree[0] == 0);
  for (jt = 0; jt < num_trees; ++jt) {
    P4EST_ASSERT (pertree[jt] <= pertree[jt + 1]);
  }
  P4EST_ASSERT (gfq[num_procs] == pertree[num_trees]);
#endif
  gquadremain = gfq[rank + 1] - gfq[rank];
  p4est->local_num_quadrants = (p4est_locidx_t) gquadremain;
  p4est->global_num_quadrants = gfq[num_procs];
  P4EST_ASSERT (quadrants->elem_count ==
                (P4EST_DIM + 1) * (size_t) p4est->local_num_quadrants);
  P4EST_ASSERT (data == NULL || data->elem_count ==
                (size_t) p4est->local_num_quadrants);

  /* allocate memory pools */
  if (dsize > 0) {
    p4est->user_data_pool = sc_mempool_new (dsize);
  }
  else {
    p4est->user_data_pool = NULL;
  }
  p4est->quadrant_pool = sc_mempool_new (sizeof (p4est_quadrant_t));

  /* find the first and last tree on this processor */
  if (p4est->local_num_quadrants > 0) {
    gkey = gfq[rank];
    gk1 = sc_bsearch_range (&gkey, pertree, num_trees,
                            sizeof (p4est_gloidx_t), p4est_gloidx_compare);
    P4EST_ASSERT (gk1 < (size_t) num_trees);
    gtreeskip = gkey - pertree[gk1];
    gkey = gfq[rank + 1] - 1;
    gk2 = sc_bsearch_range (&gkey, pertree, num_trees,
                            sizeof (p4est_gloidx_t), p4est_gloidx_compare);
    P4EST_ASSERT (gk1 <= gk2 && gk2 < (size_t) num_trees);
    p4est->first_local_tree = (p4est_topidx_t) gk1;
    p4est->last_local_tree = (p4est_topidx_t) gk2;
  }
  else {
    gtreeskip = 0;
    p4est->first_local_tree = -1;
    p4est->last_local_tree = -2;
  }

  /* populate trees */
  zqoffset = 0;
  gquadremain = p4est->local_num_quadrants;
  p4est->trees = sc_array_new_size (sizeof (p4est_tree_t), num_trees);
  for (jt = 0; jt < num_trees; ++jt) {
    /* all trees need at least some basic setup */
    tree = p4est_tree_array_index (p4est->trees, jt);
    sc_array_init (&tree->quadrants, sizeof (p4est_quadrant_t));
    P4EST_QUADRANT_INIT (&tree->first_desc);
    P4EST_QUADRANT_INIT (&tree->last_desc);
    tree->quadrants_offset = (p4est_locidx_t) zqoffset;
    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }
    for (; i <= P4EST_MAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = -1;
    }
    q = NULL;
    tree->maxlevel = 0;
    if (jt >= p4est->first_local_tree && jt <= p4est->last_local_tree) {
      /* this tree has local quadrants */
      gtreeremain = pertree[jt + 1] - pertree[jt] - gtreeskip;
      P4EST_ASSERT (gtreeremain > 0 && gquadremain > 0);
      zqthistree = (size_t) SC_MIN (gtreeremain, gquadremain);
      P4EST_ASSERT (zqthistree > 0);
      sc_array_resize (&tree->quadrants, zqthistree);
      tml = 0;
      for (qz = 0; qz < zqthistree; ++qz) {
        q = p4est_quadrant_array_index (&tree->quadrants, qz);
        P4EST_QUADRANT_INIT (q);
        q->x = *qap++;
        q->y = *qap++;
#ifdef P4_TO_P8
        q->z = *qap++;
#endif
/* *INDENT-OFF* HORRIBLE indent bug */
        q->level = ql = (int8_t) *qap++;
/* *INDENT-ON* */
        P4EST_ASSERT (ql >= 0 && ql <= P4EST_QMAXLEVEL);
        ++tree->quadrants_per_level[ql];
        tml = SC_MAX (tml, ql);
        p4est_quadrant_init_data (p4est, jt, q, NULL);
        if (data != NULL) {
          memcpy (q->p.user_data, dap, dsize);
          dap += dsize;
        }
        if (qz == 0) {
          p4est_quadrant_first_descendant (q, &tree->first_desc,
                                           P4EST_QMAXLEVEL);
        }
      }
      p4est_quadrant_last_descendant (q, &tree->last_desc, P4EST_QMAXLEVEL);
      tree->maxlevel = tml;
      zqoffset += zqthistree;
      gquadremain -= (p4est_gloidx_t) zqthistree;
      gtreeskip = 0;
    }
  }
  P4EST_ASSERT (zqoffset == (size_t) p4est->local_num_quadrants);
  P4EST_ASSERT (gquadremain == 0);

  /* communicate partition information */
  p4est->global_first_position =
    P4EST_ALLOC (p4est_quadrant_t, num_procs + 1);
  p4est_comm_global_partition (p4est, NULL);

  /* print more statistics */
  P4EST_VERBOSEF ("total local quadrants %lld\n",
                  (long long) p4est->local_num_quadrants);

  P4EST_ASSERT (p4est_is_valid (p4est));
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_inflate\n");

  return p4est;
}
