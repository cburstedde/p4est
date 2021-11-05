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
  sprintf (buffer, "%d\n", p4est->mpirank);     /* better snprintf */
#else

#endif
}

static void
read_header (size_t data_size, char *buffer, void *user)
{
  return;
}

static void
write_quad_data (p4est_t * p4est, sc_array_t * quad_data)
{
  p4est_locidx_t      i;
  int                *current;

  for (i = 0; i < p4est->local_num_quadrants; ++i) {
    current = (int *) sc_array_index (quad_data, i);
    *current = p4est->mpirank;
    //printf ("write: %i (%i), ", *current, p4est->mpirank);
  }
  //printf ("\n");
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 rank, size;
  int                 level = 3;
  int                *current;
  const size_t        header_size = 0;  //4;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  p4est_file_context_t *fc;
  p4est_locidx_t      i;
  sc_array_t          quad_data;
  sc_array_t          read_data;

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

  /* intialize quadrant data array */
  sc_array_init (&quad_data, sizeof (int));
  sc_array_resize (&quad_data, p4est->local_num_quadrants);

  fc =
    p4est_file_open_create (p4est, "test_io.out", header_size, write_header,
                            p4est);
  write_quad_data (p4est, &quad_data);
  p4est_file_write (fc, &quad_data);
#if 0
  MPI_Offset          offset;
  MPI_File_get_position (fc->file, &offset);
  printf ("main: [%i] offset = %lld\n", fc->p4est->mpirank, offset);
#endif

  p4est_file_close (fc);

  /* intialize read quadrant data array */
  sc_array_init (&read_data, sizeof (int));
  sc_array_resize (&read_data, p4est->local_num_quadrants);

  fc =
    p4est_file_open_read (p4est, "test_io.out", header_size, read_header,
                          NULL);
  p4est_file_read (fc, &read_data);

  p4est_file_close (fc);

  /* print read data */
  for (i = 0; i < p4est->local_num_quadrants; ++i) {
    current = (int *) sc_array_index (&read_data, i);
    printf ("%i (%i), ", *current, p4est->mpirank);
  }
  printf ("\n");

  /*clean up */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);
  sc_array_reset (&quad_data);
  sc_array_reset (&read_data);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
