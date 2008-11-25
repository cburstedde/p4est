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

#ifdef P4_TO_P8
#include <p8est_bits.h>
#include <p8est_points.h>
#include <p8est_vtk.h>
#else
#include <p4est_bits.h>
#include <p4est_points.h>
#include <p4est_vtk.h>
#endif /* !P4_TO_P8 */
#include <sc_io.h>

/*
 * Usage: p4est_points <configuration> <level> <prefix>
 *        possible configurations:
 *        o unit      Refinement on the unit square.
 *        o three     Refinement on a forest with three trees.
 *        o moebius   Refinement on a 5-tree Moebius band.
 *        o star      Refinement on a 6-tree star shaped domain.
 *        o periodic  Refinement on the unit square with periodic b.c.
 */

static p4est_quadrant_t *
read_points (const char *filename, p4est_topidx_t num_trees,
             p4est_locidx_t * num_points)
{
  int                 retval;
  int                 qshift;
  unsigned            ucount, u, ignored;
  double              x, y, z;
  double              qlen;
  double             *point_buffer;
  p4est_quadrant_t   *points, *q;
  FILE               *file;

  file = fopen (filename, "rb");
  SC_CHECK_ABORTF (file != NULL, "Open file %s", filename);

  sc_fread (&ucount, sizeof (unsigned int), 1, file, "Read point count");

  point_buffer = P4EST_ALLOC (double, 3 * ucount);
  sc_fread (point_buffer, sizeof (double), (size_t) (3 * ucount), file,
            "Read points");

  retval = fclose (file);
  SC_CHECK_ABORTF (retval == 0, "Close file %s", filename);

  q = points = P4EST_ALLOC_ZERO (p4est_quadrant_t, ucount);
  qlen = (double) (1 << P4EST_QMAXLEVEL);
  qshift = P4EST_MAXLEVEL - P4EST_QMAXLEVEL;
  for (ignored = u = 0; u < ucount; ++u) {
    x = point_buffer[3 * u + 0];
    y = point_buffer[3 * u + 2];
    z = point_buffer[3 * u + 1];
    if (x < 0. || x > 1. || y < 0. || y > 1. || z < 0. || z > 1.) {
      ++ignored;
      continue;
    }

    q->x = (p4est_qcoord_t) (x * qlen) << qshift;
    q->x = SC_MIN (q->x, P4EST_ROOT_LEN - 1);
    q->y = (p4est_qcoord_t) (y * qlen) << qshift;
    q->y = SC_MIN (q->y, P4EST_ROOT_LEN - 1);
#ifdef P4_TO_P8
    q->z = (p4est_qcoord_t) (z * qlen) << qshift;
    q->z = SC_MIN (q->z, P4EST_ROOT_LEN - 1);
#endif
    q->level = P4EST_MAXLEVEL;
    q->p.which_tree =
      (p4est_topidx_t) ((double) num_trees * rand () / (RAND_MAX + 1.0));
    P4EST_ASSERT (p4est_quadrant_is_node (q, true));

    ++q;
  }
  P4EST_FREE (point_buffer);

  if (num_points != NULL) {
    *num_points = (p4est_locidx_t) (ucount - ignored);
  }

  return points;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 num_procs, rank;
  int                 maxlevel;
  bool                wrongusage;
  char                buffer[BUFSIZ];
  p4est_locidx_t      num_points, max_points;
  p4est_connectivity_t *conn;
  p4est_quadrant_t   *points;
  p4est_t            *p4est;
  MPI_Comm            mpicomm;
  const char         *usage;

  /* initialize MPI and p4est internals */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (rank, sc_generic_abort, &mpicomm, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  usage =
    "Arguments: <configuration> <level> <maxpoints> <prefix>\n"
    "   Configuration can be any of\n"
#ifndef P4_TO_P8
    "      unit|three|moebius|star|periodic\n"
#else
    "      unit|periodic|rotwrap|twocubes|rotcubes\n"
#endif
    "   Level controls the maximum depth of refinement\n"
    "   Maxpoints is the maximum number of points per quadrant\n"
    "      which applies to all quadrants above maxlevel\n"
    "      A value of 0 refines recursively to maxlevel\n"
    "      A value of -1 does no refinement at all\n"
    "   Prefix is for loading a point data file";
  wrongusage = false;
  if (!wrongusage && argc != 5) {
    wrongusage = true;
  }
  if (!wrongusage) {
#ifndef P4_TO_P8
    if (!strcmp (argv[1], "unit")) {
      conn = p4est_connectivity_new_unitsquare ();
    }
    else if (!strcmp (argv[1], "three")) {
      conn = p4est_connectivity_new_corner ();
    }
    else if (!strcmp (argv[1], "moebius")) {
      conn = p4est_connectivity_new_moebius ();
    }
    else if (!strcmp (argv[1], "star")) {
      conn = p4est_connectivity_new_star ();
    }
    else if (!strcmp (argv[1], "periodic")) {
      conn = p4est_connectivity_new_periodic ();
    }
#else
    if (!strcmp (argv[1], "unit")) {
      conn = p8est_connectivity_new_unitcube ();
    }
    else if (!strcmp (argv[1], "periodic")) {
      conn = p8est_connectivity_new_periodic ();
    }
    else if (!strcmp (argv[1], "rotwrap")) {
      conn = p8est_connectivity_new_rotwrap ();
    }
    else if (!strcmp (argv[1], "twocubes")) {
      conn = p8est_connectivity_new_twocubes ();
    }
    else if (!strcmp (argv[1], "rotcubes")) {
      conn = p8est_connectivity_new_rotcubes ();
    }
#endif
    else {
      wrongusage = true;
    }
  }
  if (!wrongusage) {
    maxlevel = atoi (argv[2]);
    if (maxlevel < 0 || maxlevel > P4EST_QMAXLEVEL) {
      wrongusage = true;
    }
  }
  if (!wrongusage) {
    max_points = (p4est_locidx_t) atoi (argv[3]);
    if (max_points < -1) {
      wrongusage = true;
    }
  }
  if (wrongusage) {
    SC_CHECK_ABORT (rank != 0, usage);
    mpiret = MPI_Barrier (mpicomm);
    SC_CHECK_MPI (mpiret);
  }

  snprintf (buffer, BUFSIZ, "%s%d_%d.pts", argv[4], rank, num_procs);
  points = read_points (buffer, conn->num_trees, &num_points);
  SC_LDEBUGF ("Read %lld points\n", (long long) num_points);

  p4est = p4est_new_points (mpicomm, conn, maxlevel, points,
                            num_points, max_points, 5, NULL, NULL);
  P4EST_FREE (points);
  p4est_vtk_write_file (p4est, P4EST_STRING "_points_created");

  p4est_partition (p4est, NULL);
  p4est_vtk_write_file (p4est, P4EST_STRING "_points_partition");

  p4est_balance (p4est, P4EST_BALANCE_FULL, NULL);
  p4est_vtk_write_file (p4est, P4EST_STRING "_points_balance");

  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF points2.c */
