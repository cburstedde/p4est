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

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#endif
#include <sc_options.h>

/* #define LOADCONN_VTK */

static int          refine_level;
static int          level_shift;

static int
refine_fractal (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  int                 qid;

  if ((int) q->level >= refine_level) {
    return 0;
  }
  if ((int) q->level < refine_level - level_shift) {
    return 1;
  }

  qid = ((int) q->level == 0 ?
         (which_tree % P4EST_CHILDREN) : p4est_quadrant_child_id (q));

  return (qid == 0 || qid == 3
#ifdef P4_TO_P8
          || qid == 5 || qid == 6
#endif
    );
}

static void
run_load (sc_MPI_Comm mpicomm, p4est_connectivity_t * conn, int level)
{
  int                 mpiret;
  double              elapsed_create, elapsed_partition, elapsed_balance;
#ifdef LOADCONN_VTK
  char                filename[BUFSIZ];
#endif
  p4est_t            *p4est;

  P4EST_GLOBAL_PRODUCTIONF ("Run load on level %d\n", level);

  /* create and refine the forest */

  mpiret = sc_MPI_Barrier (mpicomm);
  SC_CHECK_MPI (mpiret);
  elapsed_create = -sc_MPI_Wtime ();

  p4est = p4est_new_ext (mpicomm, conn, 0, level, 1, 0, NULL, NULL);

  level_shift = 4;
  refine_level = level + level_shift;
  p4est_refine (p4est, 1, refine_fractal, NULL);

  elapsed_create += sc_MPI_Wtime ();

#ifdef LOADCONN_VTK
  snprintf (filename, BUFSIZ, "loadconn%d_%02d_C", P4EST_DIM, level);
  p4est_vtk_write_file (p4est, NULL, filename);
#endif

  /* partition the forest */

  mpiret = sc_MPI_Barrier (mpicomm);
  SC_CHECK_MPI (mpiret);
  elapsed_partition = -sc_MPI_Wtime ();

  p4est_partition (p4est, 0, NULL);

  elapsed_partition += sc_MPI_Wtime ();

  /* balance the forest */

  mpiret = sc_MPI_Barrier (mpicomm);
  SC_CHECK_MPI (mpiret);
  elapsed_balance = -sc_MPI_Wtime ();

  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  elapsed_balance += sc_MPI_Wtime ();

#ifdef LOADCONN_VTK
  snprintf (filename, BUFSIZ, "loadconn%d_%02d_B", P4EST_DIM, level);
  p4est_vtk_write_file (p4est, NULL, filename);
#endif

  /* report timings */

  P4EST_GLOBAL_PRODUCTIONF ("Timings %d: %g %g %g\n", level, elapsed_create,
                            elapsed_partition, elapsed_balance);

  p4est_destroy (p4est);
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret, retval;
  int                 level;
  const char         *filename;
  p4est_connectivity_t *conn;
  sc_options_t       *opt;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "Upfront refinement level");
  retval = sc_options_parse (p4est_package_id, SC_LP_ERROR, opt, argc, argv);
  if (retval == -1 || retval + 1 != argc) {
    sc_options_print_usage (p4est_package_id, SC_LP_PRODUCTION, opt, NULL);
    sc_abort_collective ("Usage error");
  }
  filename = argv[retval];
  P4EST_LDEBUGF ("Loading %s\n", filename);
  conn = p4est_connectivity_load (filename, NULL);

  run_load (mpicomm, conn, level);

  p4est_connectivity_destroy (conn);
  sc_options_destroy (opt);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
