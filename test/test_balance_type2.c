/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2008 Carsten Burstedde, Lucas Wilcox.

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
#include <p4est.h>
#include <p4est_bits.h>
#include <p4est_ghost.h>
#else
#include <p8est.h>
#include <p8est_bits.h>
#include <p8est_ghost.h>
#endif

#ifndef P4_TO_P8
static const int    refine_level = 5;
#else
static const int    refine_level = 3;
#endif
static const int    ueber_level = 10;

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
    return quadrant->level < ueber_level;
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
  int                 size, rank;
  unsigned            crcF, crcC;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  p4est_t            *p4estF, *p4estC;
#ifdef P4_TO_P8
  unsigned            crcE;
  p4est_t            *p4estE;
#endif

  /* initialize */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpicomm, &size);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (rank, sc_generic_abort, &mpicomm, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* create forest and refine */
#ifndef P4_TO_P8
  connectivity = p4est_connectivity_new_star ();
#else
  connectivity = p8est_connectivity_new_rotcubes ();
#endif
  p4est = p4est_new (mpicomm, connectivity, 0, 0, NULL, NULL);
  p4est_refine (p4est, true, refine_fn, NULL);

  /* test face balance */
  p4estF = p4est_copy (p4est, false);
#ifndef P4_TO_P8
  p4est_balance (p4estF, P4EST_BALANCE_FACE, NULL);
#else
  p4est_balance (p4estF, P8EST_BALANCE_FACE, NULL);
#endif
  crcF = p4est_checksum (p4estF);
  P4EST_GLOBAL_INFOF ("Face balance with %lld quadrants and crc 0x%x\n",
                      (long long) p4estF->global_num_quadrants, crcF);

#ifdef P4_TO_P8
  /* test edge balance */
  p4estE = p4est_copy (p4est, true);
  p4est_balance (p4estF, P8EST_BALANCE_EDGE, NULL);
  p4est_balance (p4estE, P8EST_BALANCE_EDGE, NULL);
  crcE = p4est_checksum (p4estE);
  SC_CHECK_ABORT (crcE == p4est_checksum (p4estF), "mismatch A");
  P4EST_GLOBAL_INFOF ("Edge balance with %lld quadrants and crc 0x%x\n",
                      (long long) p4estE->global_num_quadrants, crcE);
#endif

  /* test corner balance */
  p4estC = p4est_copy (p4est, true);
#ifndef P4_TO_P8
  p4est_balance (p4estF, P4EST_BALANCE_CORNER, NULL);
  p4est_balance (p4estC, P4EST_BALANCE_CORNER, NULL);
#else
  p4est_balance (p4estF, P8EST_BALANCE_CORNER, NULL);
  p4est_balance (p4estC, P8EST_BALANCE_CORNER, NULL);
#endif
  crcC = p4est_checksum (p4estC);
  SC_CHECK_ABORT (crcC == p4est_checksum (p4estF), "mismatch B");
  P4EST_GLOBAL_INFOF ("Corner balance with %lld quadrants and crc 0x%x\n",
                      (long long) p4estC->global_num_quadrants, crcC);

  /* destroy forests and connectivity */
  p4est_destroy (p4est);
  p4est_destroy (p4estF);
#ifdef P4_TO_P8
  p4est_destroy (p4estE);
#endif
  p4est_destroy (p4estC);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_balance_type2.c */
