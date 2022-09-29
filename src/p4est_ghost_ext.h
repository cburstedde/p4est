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

/** \file p4est_ghost_ext.h
 *
 * passing quadrants and data to neighboring processes
 *
 * \ingroup p4est
 */

#ifndef P4EST_GHOST_EXT_H
#define P4EST_GHOST_EXT_H

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

SC_EXTERN_C_BEGIN;

typedef enum
{
  P4EST_GHOST_UNBALANCED_ABORT = 0,
  P4EST_GHOST_UNBALANCED_FAIL,
  P4EST_GHOST_UNBALANCED_ALLOW
}
p4est_ghost_tolerance_t;

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

/** Builds the ghost layer using recursive top down approach.
 *
 * This will gather the quadrants from each neighboring proc to build
 * one layer of face and corner based ghost elements around the ones they own.
 *
 * \param [in] p4est            The forest for which the ghost layer will be
 *                              generated.
 * \param [in] btype            Which ghosts to include (across face, corner
 *                              or full).
 * \return                      A fully initialized ghost layer.
 */
p4est_ghost_t      *p4est_ghost_new_ext (p4est_t * p4est,
                                         p4est_connect_type_t btype);

#ifdef P4EST_ENABLE_MPI

sc_array_t         *p4est_ghost_array_index (sc_array_t * array, int i);

void                p4est_quadrant_find_tree_corner_owners (p4est_t * p4est,
                                                            p4est_topidx_t
                                                            treeid,
                                                            int treecorner,
                                                            const
                                                            p4est_quadrant_t *
                                                            q,
                                                            sc_array_t *
                                                            q_procs,
                                                            int *nurgood);

size_t              ghost_proc_type (sc_array_t * array, size_t zindex,
                                     void *data);

int                 p4est_add_ghost_to_buf (sc_array_t * buf,
                                            p4est_topidx_t treeid,
                                            p4est_locidx_t number,
                                            const p4est_quadrant_t * q);

void                p4est_ghost_mirror_init (p4est_ghost_t * ghost,
                                             int mpirank,
                                             sc_array_t * send_bufs,
                                             p4est_ghost_mirror_t * m);

void                p4est_ghost_mirror_add (p4est_ghost_mirror_t * m,
                                            p4est_topidx_t treeid,
                                            p4est_locidx_t number,
                                            p4est_quadrant_t * q, int p);

void                p4est_ghost_mirror_reset (p4est_ghost_t * ghost,
                                              p4est_ghost_mirror_t * m,
                                              int populate);

void                p4est_ghost_test_add (p4est_t * p4est,
                                          p4est_ghost_mirror_t * m,
                                          p4est_quadrant_t * q,
                                          p4est_topidx_t t,
                                          p4est_quadrant_t * nq,
                                          p4est_topidx_t nt, int32_t touch,
                                          int rank,
#if 0
                                          sc_array_t * send_bufs,
#endif
                                          p4est_locidx_t local_num);

#endif

SC_EXTERN_C_END;

#endif
