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

/*
 * This is an internal header file.
 * It is included from a C file that may be preceded by p4est_to_p8est.h.
 */

#ifndef P4EST_SPHERES_GLOBAL_H
#define P4EST_SPHERES_GLOBAL_H

#include <sc_random.h>
#include <sc_statistics.h>
#ifndef P4_TO_P8
#include "p4est_spheres.h"
#else
#include "p8est_spheres.h"
#endif /* P4_TO_P8 */

SC_EXTERN_C_BEGIN;

typedef struct qu_data
{
  int                 set_refine;
}
qu_data_t;

typedef struct sph_item
{
  p4est_sphere_t      sph;
  p4est_gloidx_t      gid;
}
sph_item_t;

typedef struct sr_buf
{
  int                 rank;
  sc_array_t         *items;
}
sr_buf_t;

typedef struct spheres_global
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize, mpirank;

  int                 minlevel;
  int                 maxlevel;
  int                 scaling;

  int                 repetitions;

  int                 write_vtk;
  int                 mpiwrap;

  double              rmax;
  double              thickness;
  double              lfraction;
  double              spherelems;

  const char         *prefix;

  p4est_connectivity_t *conn;
  p4est_t            *p4est;

  sc_rand_state_t     rstate;

  /* data of generated spheres */
  sc_array_t         *sphr;
  p4est_locidx_t      lsph;
  p4est_gloidx_t      gsoff;
  sc_array_t         *lcounts;
  sc_array_t         *goffsets;

  /* metadata for partition search */
  int                 last_to_rank;
  p4est_sphere_t      box;
  sc_array_t         *to_procs;
  sc_array_t         *sphere_procs;
  sr_buf_t           *last_to_proc;

  /* metadata for pattern reversal */
  int                 ntop, nint, nbot;
  int                 notify_alltoall;
  int                *last_payload;
  sc_array_t         *notify;
  sc_array_t         *payload;

  /* send and receive messages */
  int                 num_to_procs, num_from_procs;
  sc_array_t         *to_requests, *from_requests;
  sc_array_t         *from_procs;

  /* refine and partition */
  p4est_locidx_t      lqindex, lqindex_refined;
  p4est_locidx_t      lsph_offset;
  sc_array_t         *lcounts_refined;

  /* performance statistics */
  int                 num_stats;
  sc_array_t         *stats;
}
spheres_global_t;

SC_EXTERN_C_END;

#endif /* !P4EST_SPHERES_GLOBAL_H */
