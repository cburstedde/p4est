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

#include <p4est_io.h>
#include <p4est_extended.h>

static void
write_header (size_t data_size, char *buffer, void *user)
{
  p4est_t            *p4est = (p4est_t *) user;
  int                 dummy = 42;
  memcpy (buffer, &dummy, sizeof (int));
#if 0
  sprintf (buffer, "%d\n", p4est->mpirank);
#else

#endif
}

static void
write_quad_data (size_t data_size, char *buffer, void *user)
{
  callback_context_t *callback_ct = (callback_context_t *) user;

#if 0
  /* TODO: Dimension idenpendent version. */
  snprintf (buffer, 23, "(%d,%d,%d)", callback_ct->quad->x,
            callback_ct->quad->y, callback_ct->quad->level);
  printf ("strlen: %ld\n", strlen (buffer));
#else
  int                 dummy = 42;
  memcpy (buffer, &dummy, sizeof (int));
#endif
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 rank, size;
  int                 level = 3;
  const size_t        header_size = 4;
  const size_t        quad_data_write_size = 23;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  p4est_file_context_t *fc;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &size);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  connectivity = p4est_connectivity_new_unitsquare ();

  p4est = p4est_new_ext (mpicomm, connectivity, 0, level, 1, 0, NULL, NULL);

  fc =
    p4est_file_open_create (p4est, "test_io.txt", header_size, write_header,
                            p4est);
  //p4est_file_write (fc, quad_data_write_size, write_quad_data, NULL);

  p4est_file_close (fc);

  /*clean up */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
