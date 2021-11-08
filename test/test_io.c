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
#include <p4est_bits.h>

#define HEADER_INT1 42
#define HEADER_INT2 84 

#ifdef P4EST_ENABLE_MPIIO
static void
write_header (int *header)
{
  header[0] = HEADER_INT1;
  header[1] = HEADER_INT2;
}

#if 0
static void
read_header (size_t data_size, char *buffer, void *user)
{
  return;
}
#endif

static void
write_rank (p4est_t * p4est, sc_array_t * quad_data)
{
  p4est_locidx_t      i;
  int                *current;

  for (i = 0; i < p4est->local_num_quadrants; ++i) {
    current = (int *) sc_array_index (quad_data, i);
    *current = p4est->mpirank;
  }
}

#if 0
static void
write_quads (p4est_t * p4est, sc_array_t * quad_data)
{
  p4est_locidx_t      i
    p4est_tree_t * tree = p4est_tree_index (p4est->trees, 0);
  p4est_quadrant_t   *quad, *dest_quad;

  sc_array_init (quad_data, sizeof (p4est_quadrant_t));
  sc_array_resize (quad_data, tree->quadrants.elem_count);

  for (i = 0; i < tree->quadrants.elem_count; ++i)
    quad = p4est_quadrant_array_index (&tree->quadrants, i);
  dest_quad = p4est_quadrant_array_index (quad_data, i);
}
#endif

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
  p4est_tree_t       *tree;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  p4est_file_context_t *fc;
  p4est_locidx_t      i;
  p4est_gloidx_t      global_quad_num;
  sc_array_t          quad_data;
  sc_array_t          read_data;
  sc_array_t          elem_size;
  sc_array_t          quads;

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
  write_rank (p4est, &quad_data);
  p4est_file_write (fc, &quad_data);

  tree = p4est_tree_array_index (p4est->trees, 0);
  p4est_file_write (fc, &tree->quadrants);

  p4est_file_close (fc);

  /* intialize read quadrant data array */
  sc_array_init (&read_data, sizeof (int));
  sc_array_resize (&read_data, p4est->local_num_quadrants);

  sc_array_init (&quads, sizeof (p4est_quadrant_t));
  sc_array_resize (&quads, p4est->local_num_quadrants);

  fc = p4est_file_open_read (p4est, "test_io.out", header_size, read_header);

  /* check read header */
  SC_CHECK_ABORT (read_header[0] == HEADER_INT1 && read_header[1] == HEADER_INT2, "Read user-defined header");

  /* read the first data array */
  p4est_file_read (fc, &read_data);

  /* read the second data array */
  p4est_file_read (fc, &quads);

  /* check the read data */
  for (i = 0; i < p4est->local_num_quadrants; ++i) {
    p4est_quadrant_array_index (&quads, i);
    SC_CHECK_ABORT (p4est_quadrant_is_equal
                    (p4est_quadrant_array_index (&quads, i),
                     p4est_quadrant_array_index (&tree->quadrants, i)),
                    "Quadrant read");
  }

  p4est_file_close (fc);

  /* check read data of the first array */
  for (i = 0; i < p4est->local_num_quadrants; ++i) {
    current = (int *) sc_array_index (&read_data, i);
    SC_CHECK_ABORT (*current == p4est->mpirank, "Rank read");
  }

  /* append data to the existing file */
  fc = p4est_file_open_append (p4est, "test_io.out", header_size);

  for (i = 0; i < p4est->local_num_quadrants; ++i) {
    current = sc_array_index (&quad_data, i);
    *current = 42;
  }

  p4est_file_write (fc, &quad_data);

  p4est_file_close (fc);

  fc = p4est_file_open_read (p4est, "test_io.out", header_size, read_header);

  p4est_file_info (fc, &global_quad_num, p4est_version, &file_io_rev,
                   &magic_num, &elem_size);
  P4EST_GLOBAL_PRODUCTIONF
    ("file info: number of global quadrants = %ld\np4est version = %s\nfile io revision number = %d\nmagic_num = %d, number of arrays = %ld\n",
     global_quad_num, p4est_version, file_io_rev, magic_num,
     elem_size.elem_count);
  for (si = 0; si < elem_size.elem_count; ++si) {
    current_elem_size = sc_array_index (&elem_size, si);
    P4EST_GLOBAL_PRODUCTIONF ("Array %ld: element size %ld\n", si, *current_elem_size);
  }

  p4est_file_close (fc);

  /* clean up */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);
  sc_array_reset (&quad_data);
  sc_array_reset (&read_data);
  sc_array_reset (&elem_size);
  sc_array_reset (&quads);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
#endif /* !ENABLE_MPIIO */

  return 0;
}
