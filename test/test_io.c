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

#ifdef P4EST_ENABLE_MPIIO
static void
write_header (int *header)
{
  header[0] = 1;
  header[1] = 4;
}

#if 0
static void
read_header (size_t data_size, char *buffer, void *user)
{
  return;
}
#endif

static void
write_quad_data (p4est_t * p4est, sc_array_t * quad_data)
{
  p4est_locidx_t      i;
  int                *current;

  for (i = 0; i < p4est->local_num_quadrants; ++i) {
    current = (int *) sc_array_index (quad_data, i);
    *current = p4est->mpirank;
  }
}
#endif /* !ENABLE_MPIIO */

int
main (int argc, char **argv)
{
#ifdef P4EST_ENABLE_MPIIO
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 rank, size;
  int                 level = 3;
  int                *current;
  const size_t        header_size = 8;
  size_t              si;
  long               *current_elem_size;
  int                 header[2], read_header[2];
  int                 file_io_rev, magic_num;
  char                p4est_version[16];
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  p4est_file_context_t *fc;
  p4est_locidx_t      i;
  p4est_gloidx_t      global_quad_num;
  sc_array_t          quad_data;
  sc_array_t          read_data;
  sc_array_t          elem_size;

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

  /* intialize the header */
  write_header (header);

  /* intialize quadrant data array */
  sc_array_init (&quad_data, sizeof (int));
  sc_array_resize (&quad_data, p4est->local_num_quadrants);

  fc = p4est_file_open_create (p4est, "test_io.out", header_size, header);
  write_quad_data (p4est, &quad_data);
  p4est_file_write (fc, &quad_data);

  p4est_file_close (fc);

  /* intialize read quadrant data array */
  sc_array_init (&read_data, sizeof (int));
  sc_array_resize (&read_data, p4est->local_num_quadrants);

  fc = p4est_file_open_read (p4est, "test_io.out", header_size, read_header);

  /* print read header TODO: Do not use printf */
  printf ("number of arrays = %i\nnumber of bytes per element = %i\n",
          read_header[0], read_header[1]);

  p4est_file_info (fc, &global_quad_num, p4est_version, &file_io_rev,
                   &magic_num, &elem_size);
  printf
    ("file info: number of global quadrants = %ld\np4est version = %s\nfile io revision number = %d\nmagic_num = %d, number of arrays = %ld\n",
     global_quad_num, p4est_version, file_io_rev, magic_num,
     elem_size.elem_count);
  for (si = 0; si < elem_size.elem_count; ++si) {
    current_elem_size = sc_array_index (&elem_size, si);
    printf ("Array %ld: element size %ld\n", si, *current_elem_size);
  }

  p4est_file_read (fc, &read_data);

  p4est_file_close (fc);

  /* print read data */
  for (i = 0; i < p4est->local_num_quadrants; ++i) {
    current = (int *) sc_array_index (&read_data, i);
    printf ("%i (%i), ", *current, p4est->mpirank);
  }
  printf ("\n");

  /* clean up */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);
  sc_array_reset (&quad_data);
  sc_array_reset (&read_data);
  sc_array_reset (&elem_size);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
#endif /* !ENABLE_MPIIO */

  return 0;
}
