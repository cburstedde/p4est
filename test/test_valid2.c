/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
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

#ifndef P4EST_BACKWARD_DEALII
#define P4EST_BACKWARD_DEALII
#endif

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_nodes.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_nodes.h>
#include <p8est_vtk.h>
#endif

#ifndef P4_TO_P8
static const int    refine_level = 5;
#else
static const int    refine_level = 4;
#endif

static p4est_balance_type_t
check_backward_compatibility (void)
{
  p4est_balance_type_t b;

  b = P4EST_CONNECT_FULL;
  return b;
}

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 cid;
  int                 endlevel = refine_level + (int) (which_tree % 3);

  if ((int) quadrant->level >= endlevel)
    return 0;

  cid = p4est_quadrant_child_id (quadrant);
  if (cid == 0 || cid == 3
#ifdef P4_TO_P8
      || cid == 6
#endif
    ) {
    return 1;
  }

  return 0;
}

static int
coarsen_fn (p4est_t * p4est, p4est_topidx_t which_tree,
            p4est_quadrant_t * q[])
{
  int                 pid;
  p4est_quadrant_t    p;

  SC_CHECK_ABORT (p4est_quadrant_is_familypv (q), "Coarsen invocation");

  if (q[0]->level <= 2)
    return 0;

  p4est_quadrant_parent (q[0], &p);
  pid = p4est_quadrant_child_id (&p);

  return pid == 3;
}

static void
check_all (sc_MPI_Comm mpicomm, p4est_connectivity_t * conn,
           const char *vtkname, unsigned crc_expected, unsigned gcrc_expected)
{
  int                 mpiret;
  unsigned            crc_computed, gcrc_computed;
  long long           lsize[3], gsize[3];
  size_t              size_conn, size_p4est, size_ghost;
  p4est_t            *p4est;
  p4est_nodes_t      *nodes;
  p4est_ghost_t      *ghost;

  P4EST_GLOBAL_STATISTICSF ("Testing configuration %s\n", vtkname);

  p4est = p4est_new_ext (mpicomm, conn, 0, 0, 0, 0, NULL, NULL);
  p4est_refine (p4est, 1, refine_fn, NULL);
  p4est_coarsen (p4est, 1, coarsen_fn, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  p4est_partition (p4est, 0, NULL);
  p4est_vtk_write_file (p4est, NULL, vtkname);

  crc_computed = p4est_checksum (p4est);
  P4EST_GLOBAL_STATISTICSF ("Forest checksum 0x%08x\n", crc_computed);
  if (p4est->mpisize == 2 && p4est->mpirank == 0) {
    SC_CHECK_ABORT (crc_computed == crc_expected, "Forest checksum mismatch");
  }

  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

  /* compute total size of forest storage */
  size_conn = p4est_connectivity_memory_used (conn);
  size_p4est = p4est_memory_used (p4est);
  size_ghost = p4est_ghost_memory_used (ghost);
  lsize[0] = (long long) size_conn;
  lsize[1] = (long long) size_p4est;
  lsize[2] = (long long) size_ghost;
  mpiret = sc_MPI_Reduce (lsize, gsize, 3, sc_MPI_LONG_LONG_INT, sc_MPI_SUM,
                          0, mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_INFOF ("Global byte sizes: %lld %lld %lld\n",
                      gsize[0], gsize[1], gsize[2]);

  gcrc_computed = p4est_ghost_checksum (p4est, ghost);
  P4EST_GLOBAL_STATISTICSF ("Ghost checksum 0x%08x\n", gcrc_computed);
  if (p4est->mpisize == 2 && p4est->mpirank == 0) {
    SC_CHECK_ABORT (gcrc_computed == gcrc_expected,
                    "Ghost checksum mismatch");
  }

  nodes = p4est_nodes_new (p4est, ghost);
  p4est_nodes_destroy (nodes);
  p4est_ghost_destroy (ghost);

  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 size, rank;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &size);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  (void) check_backward_compatibility ();

#ifndef P4_TO_P8
  check_all (mpicomm, p4est_connectivity_new_unitsquare (),
             "test_unitsquare", 0xef45243bU, 0xbc5d0907U);
  check_all (mpicomm, p4est_connectivity_new_rotwrap (),
             "test_rotwrap2", 0x266d2739U, 0x29a31248U);
  check_all (mpicomm, p4est_connectivity_new_corner (),
             "test_corner", 0x9dad92ccU, 0x937b27afU);
  check_all (mpicomm, p4est_connectivity_new_moebius (),
             "test_moebius", 0xbbc10f7fU, 0x09b6319eU);
  check_all (mpicomm, p4est_connectivity_new_star (),
             "test_star", 0xfb28233fU, 0x8e8a32b3);
#else
  check_all (mpicomm, p8est_connectivity_new_unitcube (),
             "test_unitcube", 0x2574801fU, 0x312559a7U);
  check_all (mpicomm, p8est_connectivity_new_periodic (),
             "test_periodic3", 0xdc7e8a93U, 0x0787ca2dU);
  check_all (mpicomm, p8est_connectivity_new_rotwrap (),
             "test_rotwrap", 0xa675888dU, 0x626cbe90U);
  check_all (mpicomm, p8est_connectivity_new_twocubes (),
             "test_twocubes", 0x7188978aU, 0x4124bcabU);
  check_all (mpicomm, p8est_connectivity_new_twowrap (),
             "test_twowrap", 0x8e3f994cU, 0x9dd49e94);
  check_all (mpicomm, p8est_connectivity_new_rotcubes (),
             "test_rotcubes", 0xc0e1b235U, 0x974af07a);
  check_all (mpicomm, p8est_connectivity_new_shell (),
             "test_shell", 0x558723a2U, 0x4dedf35eU);
  check_all (mpicomm,
             p8est_connectivity_new_brick (2, 3, 4, 0, 0, 1),
             "test_brick", 0x82174e14U, 0x211da6c5);
#endif

  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
