/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2009 Carsten Burstedde, Lucas Wilcox,
                     Toby Isaac.

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
#include <p4est_bits.h>
#include <p4est_algorithms.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#else
#include <p8est_bits.h>
#include <p8est_algorithms.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
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

int
main (int argc, char **argv)
{
  MPI_Comm            mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  sc_array_t          ghost_layer;
  int                *ghost_owner;
  bool                success;
  int                 ntests;
  int                 i, j;
  p4est_lnodes_t     *lnodes;

#ifndef P4_TO_P8
  ntests = 3;
#else
  ntests = 4;
#endif

  /* initialize MPI */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, true, true, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  for (i = 0; i < ntests; i++) {
    /* create connectivity and forest structures */
    switch (i) {
#ifndef P4_TO_P8
    case 0:
      connectivity = p4est_connectivity_new_moebius ();
      break;
    case 1:
      connectivity = p4est_connectivity_new_star ();
      break;
    default:
      connectivity = p4est_connectivity_new_periodic ();
      break;
#else
    case 0:
      connectivity = p8est_connectivity_new_periodic ();
      break;
    case 1:
      connectivity = p8est_connectivity_new_rotwrap ();
      break;
    case 2:
      connectivity = p8est_connectivity_new_rotcubes ();
      break;
    default:
      connectivity = p8est_connectivity_new_shell ();
      break;
#endif
    }
    p4est = p4est_new (mpicomm, connectivity, 15, 0, NULL, NULL);

    /* refine to make the number of elements interesting */
    p4est_refine (p4est, true, refine_fn, NULL);

    /* balance the forest */
#ifndef P4_TO_P8
    p4est_balance (p4est, P4EST_BALANCE_FULL, NULL);
#else
    p4est_balance (p4est, P8EST_BALANCE_FULL, NULL);
#endif

    /* do a uniform partition */
    p4est_partition (p4est, NULL);

    sc_array_init (&ghost_layer, sizeof (p4est_quadrant_t));
    success = p4est_build_ghost_layer (p4est, P4EST_BALANCE_FULL,
                                       &ghost_layer, &ghost_owner);
    P4EST_ASSERT (success);

    for (j = 1; j <= 4; j++) {
      P4EST_GLOBAL_PRODUCTIONF ("Begin lnodes test %d:%d\n", i, j);
      lnodes = p4est_lnodes_new (p4est, &ghost_layer, j);
      p4est_lnodes_destroy (lnodes);
      P4EST_GLOBAL_PRODUCTIONF ("End lnodes test %d:%d\n", i, j);
    }

    /* clean up */
    sc_array_reset (&ghost_layer);
    P4EST_FREE (ghost_owner);

    p4est_destroy (p4est);
    p4est_connectivity_destroy (connectivity);
  }

  /* exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_lnodes2.c */
