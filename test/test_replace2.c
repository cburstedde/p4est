/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2012 The University of Texas System
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
#include <p4est_extended.h>
#include <p4est_communication.h>
#include <p4est_vtk.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_communication.h>
#include <p8est_vtk.h>
#endif

#ifndef P4_TO_P8
static int          refine_level = 5;
#else
static int          refine_level = 3;
#endif

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 cid;

  if (which_tree == 2 || which_tree == 3) {
    return 0;
  }

  cid = p4est_quadrant_child_id (quadrant);

  if (cid == P4EST_CHILDREN - 1 ||
      (quadrant->x >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2) &&
       quadrant->y >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#ifdef P4_TO_P8
       && quadrant->z >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#endif
      )) {
    return 1;
  }
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && cid == 2) {
    return 1;
  }
  if (quadrant->x == P4EST_QUADRANT_LEN (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->y >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

static int
coarsen_fn (p4est_t * p4est, p4est_topidx_t which_tree,
            p4est_quadrant_t * q[])
{
  SC_CHECK_ABORT (p4est_quadrant_is_familypv (q), "Coarsen invocation");

  return q[0]->y < P4EST_ROOT_LEN / 2;
}

static void
replace_fn (p4est_t * p4est, p4est_topidx_t which_tree,
            int num_outgoing, p4est_quadrant_t * outgoing[],
            int num_incoming, p4est_quadrant_t * incoming[])
{
  p4est_quadrant_t  **fam;
  p4est_quadrant_t   *p;

  P4EST_ASSERT (num_outgoing + num_incoming == P4EST_CHILDREN + 1);
  P4EST_ASSERT (num_outgoing == 1 || num_incoming == 1);

  if (num_outgoing == 1) {
    p = outgoing[0];
    fam = incoming;
  }
  else {
    p = incoming[0];
    fam = outgoing;
  }

  SC_CHECK_ABORT (p4est_quadrant_is_familypv (fam),
                  P4EST_STRING "_replace_t family is not a family");
  SC_CHECK_ABORT (p4est_quadrant_is_parent (p, fam[0]),
                  P4EST_STRING
                  "_replace_t incoming and outgoing don't align");
}

int
main (int argc, char **argv)
{
  int                 mpirank, mpisize;
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* create connectivity and forest structures */
#ifdef P4_TO_P8
  connectivity = p8est_connectivity_new_rotcubes ();
#else
  connectivity = p4est_connectivity_new_star ();
#endif
  p4est = p4est_new_ext (mpicomm, connectivity, 15, 0, 0, 1, NULL, NULL);
  p4est_refine_ext (p4est, 1, P4EST_QMAXLEVEL, refine_fn, NULL, replace_fn);
  p4est_coarsen_ext (p4est, 1, 0, coarsen_fn, NULL, replace_fn);
  p4est_balance_ext (p4est, P4EST_CONNECT_FULL, NULL, replace_fn);

  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
