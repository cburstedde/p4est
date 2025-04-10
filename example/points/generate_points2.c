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

/* Please see doc/example_points.dox for a documentation of this program. */

static void
generate_points (const char *filename,
                 p4est_connectivity_t * conn,
                 p4est_gloidx_t global_num_points, sc_MPI_Comm mpicomm)
{
  int                 mpiret;
  int                 num_procs, rank;
  int                 count;
  p4est_gloidx_t      offset_mine, offset_next;
  p4est_locidx_t      local_num_points;
  p4est_locidx_t      u;
  double             *point_buffer;
  double              theta;
#ifdef P4_TO_P8
  double              phi;
#else
  double              dtheta;
#endif
  sc_MPI_File         file_handle;
  sc_MPI_Offset       mpi_offset;

  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* open a file (create if the file does not exist) */
  mpiret = sc_io_open (mpicomm, filename, SC_IO_WRITE_CREATE,
                       sc_MPI_INFO_NULL, &file_handle);
  SC_CHECK_MPI (mpiret);

  /* local (MPI) number of points */
  local_num_points = global_num_points / num_procs;

  /* offset to first point of current MPI process */
  offset_mine = p4est_partition_cut_gloidx (global_num_points,
                                            rank, num_procs);

  /* offset to first point of successor MPI process */
  offset_next = p4est_partition_cut_gloidx (global_num_points,
                                            rank + 1, num_procs);
  local_num_points = (p4est_locidx_t) (offset_next - offset_mine);

  /* allocate buffer for point's coordinates */
  point_buffer = P4EST_ALLOC (double, 3 * local_num_points);

  /* set file offset (in bytes) for this calling process */
  /* *INDENT-OFF* HORRIBLE indent bug */
  mpi_offset = (sc_MPI_Offset) offset_mine * 3 * sizeof (double);
  /* *INDENT-ON* */

#ifndef P4_TO_P8
  /* 2D */
  dtheta = (2. * M_PI) / global_num_points;

  for (u = 0; u < local_num_points; ++u) {
    theta = (offset_mine + u) * dtheta;

    point_buffer[3 * u + 0] = 0.5 + 0.25 * cos (theta);
    point_buffer[3 * u + 1] = 0.5 + 0.25 * sin (2 * theta);
    point_buffer[3 * u + 2] = 0.0;
  }

#else

  /* 3D */
  for (u = 0; u < local_num_points; ++u) {
    theta = 2. * M_PI * rand () / (RAND_MAX + 1.0);
    phi = M_PI * rand () / (RAND_MAX + 1.0);

    point_buffer[3 * u + 0] = 0.5 + 0.25 * cos (theta) * sin (phi);
    point_buffer[3 * u + 1] = 0.5 + 0.25 * sin (theta) * sin (phi);
    point_buffer[3 * u + 2] = 0.5 + 0.25 * cos (phi);
  }
#endif

  if (rank == 0) {
    /* write the global number number of points */
    mpiret = sc_io_write_at (file_handle, 0, &global_num_points,
                             sizeof (p4est_gloidx_t), sc_MPI_BYTE, &count);
    SC_CHECK_MPI (mpiret);
    SC_CHECK_ABORT (count == (int) sizeof (p4est_gloidx_t),
                    "Write number of global points: count mismatch");
  }

  /* each MPI process writes its data for its own offset */
  mpiret =
    sc_io_write_at_all (file_handle, mpi_offset + sizeof (p4est_gloidx_t),
                        &point_buffer[0],
                        3 * local_num_points * sizeof (double), sc_MPI_BYTE,
                        &count);
  SC_CHECK_MPI (mpiret);
  SC_CHECK_ABORT (count == (int) (3 * local_num_points * sizeof (double)),
                  "Write points: count mismatch");

  P4EST_FREE (point_buffer);

  /* close the file collectively */
  mpiret = sc_io_close (&file_handle);
  SC_CHECK_MPI (mpiret);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 num_procs, rank;
  int                 wrongusage;
  char                buffer[BUFSIZ];
  p4est_gloidx_t      global_num_points;
  p4est_connectivity_t *conn;
  sc_MPI_Comm         mpicomm;
  const char         *usage;

  /* initialize MPI and p4est internals */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  usage =
    "Arguments: <configuration> <globalnumpoints> <prefix>\n"
    "   Configuration can be any of\n"
#ifndef P4_TO_P8
    "      unit|brick|three|moebius|star|periodic\n"
#else
    "      unit|brick|periodic|rotwrap|twocubes|rotcubes\n"
#endif
    "   Globalnumpoints is the total number of points generated\n"
    "      over all MPI process, >= 0\n"
    "   Prefix is for writing a point data file (using MPI-IO)\n";
  wrongusage = 0;
  if (!wrongusage && argc != 4) {
    P4EST_GLOBAL_LERROR ("Invalid command line argument count\n");
    wrongusage = 1;
  }
  conn = NULL;
  if (!wrongusage) {
#ifndef P4_TO_P8
    if (!strcmp (argv[1], "unit")) {
      conn = p4est_connectivity_new_unitsquare ();
    }
    else if (!strcmp (argv[1], "brick")) {
      conn = p4est_connectivity_new_brick (2, 3, 0, 0);
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
    else if (!strcmp (argv[1], "brick")) {
      conn = p8est_connectivity_new_brick (2, 3, 4, 0, 0, 0);
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
      P4EST_GLOBAL_LERROR ("Invalid connectivity configuration\n");
      wrongusage = 1;
    }
  }
  if (!wrongusage) {
    global_num_points = (p4est_gloidx_t) atol (argv[2]);
    if (global_num_points <= -1) {
      P4EST_GLOBAL_LERROR ("Invalid global number of points\n");
      wrongusage = 1;
    }
  }
  if (wrongusage) {
    P4EST_GLOBAL_LERROR (usage);
    sc_abort_collective ("Usage error");
  }

  SC_GLOBAL_PRODUCTIONF ("Write %lld total points\n",
                         (long long) global_num_points);
  snprintf (buffer, BUFSIZ, "%s.pts", argv[3]);
  generate_points (buffer, conn, global_num_points, mpicomm);

  /* in the present version of this program the connectivity is not used */
  p4est_connectivity_destroy (conn);

  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
