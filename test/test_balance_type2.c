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
#include <p4est_extended.h>
#include <p4est_ghost.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
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
    return (int) quadrant->level < ueber_level;
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

static unsigned
test_checksum (p4est_t * p4est, int have_zlib)
{
  return have_zlib ? p4est_checksum (p4est) : 0;
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 have_zlib;
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
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &size);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* establish parallel logging */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* check for ZLIB usability */
  if (!(have_zlib = p4est_have_zlib ())) {
    P4EST_GLOBAL_LERROR
      ("Not found a working ZLIB installation: ignoring CRCs\n");
  }

  /* create forest and refine */
#ifndef P4_TO_P8
  connectivity = p4est_connectivity_new_star ();
#else
  connectivity = p8est_connectivity_new_rotcubes ();
#endif
  p4est = p4est_new_ext (mpicomm, connectivity, 0, 0, 0, 0, NULL, NULL);
  p4est_refine (p4est, 1, refine_fn, NULL);

  /* test face balance */
  p4estF = p4est_copy (p4est, 0);
#ifndef P4_TO_P8
  p4est_balance (p4estF, P4EST_CONNECT_FACE, NULL);
#else
  p4est_balance (p4estF, P8EST_CONNECT_FACE, NULL);
#endif
  crcF = test_checksum (p4estF, have_zlib);
  P4EST_GLOBAL_INFOF ("Face balance with %lld quadrants and crc 0x%08x\n",
                      (long long) p4estF->global_num_quadrants, crcF);

#ifdef P4_TO_P8
  /* test edge balance */
  p4estE = p4est_copy (p4est, 1);
  p4est_balance (p4estF, P8EST_CONNECT_EDGE, NULL);
  p4est_balance (p4estE, P8EST_CONNECT_EDGE, NULL);
  crcE = test_checksum (p4estE, have_zlib);
  SC_CHECK_ABORT (crcE == test_checksum (p4estF, have_zlib), "mismatch A");
  P4EST_GLOBAL_INFOF ("Edge balance with %lld quadrants and crc 0x%08x\n",
                      (long long) p4estE->global_num_quadrants, crcE);
#endif

  /* test corner balance */
  p4estC = p4est_copy (p4est, 1);
#ifndef P4_TO_P8
  p4est_balance (p4estF, P4EST_CONNECT_CORNER, NULL);
  p4est_balance (p4estC, P4EST_CONNECT_CORNER, NULL);
#else
  p4est_balance (p4estF, P8EST_CONNECT_CORNER, NULL);
  p4est_balance (p4estC, P8EST_CONNECT_CORNER, NULL);
#endif
  crcC = test_checksum (p4estC, have_zlib);
  SC_CHECK_ABORT (crcC == test_checksum (p4estF, have_zlib), "mismatch B");
  P4EST_GLOBAL_INFOF ("Corner balance with %lld quadrants and crc 0x%08x\n",
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

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
