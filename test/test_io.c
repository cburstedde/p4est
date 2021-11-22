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

//#ifdef P4EST_ENABLE_MPIIO
static void
write_header (int *header)
{
  header[0] = HEADER_INT1;
  header[1] = HEADER_INT2;
}

static int
refine (p4est_t * p4est, p4est_topidx_t which_tree,
        p4est_quadrant_t * quadrant)
{
  return quadrant->x == 0 && quadrant->y == 0;
}

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

//#endif /* !ENABLE_MPIIO */

int
main (int argc, char **argv)
{
//#ifdef P4EST_ENABLE_MPIIO
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 rank, size;
  int                 level = 3;
  int                *current;
  char               *current_char;
  const size_t        header_size = 8;
  size_t              si, read_header_size;
  long               *current_elem_size;
  int                 header[2], read_header[2];
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
  sc_array_t          unaligned;

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

  /* Test the data array padding by provoking a number of qudrants that
   * is not divisible by 16
   */
  p4est_refine (p4est, 1, refine, NULL);

  /* intialize the header */
  write_header (header);

  /* intialize quadrant data array */
  sc_array_init (&quad_data, sizeof (int));
  sc_array_resize (&quad_data, p4est->local_num_quadrants);

  fc = p4est_file_open_create (p4est, "test_io.out", header_size, header);
  SC_CHECK_ABORT (fc != NULL, "Open create");
  write_rank (p4est, &quad_data);
  SC_CHECK_ABORT (p4est_file_write (fc, &quad_data) != NULL, "Write ranks");

  tree = p4est_tree_array_index (p4est->trees, 0);
  SC_CHECK_ABORT (p4est_file_write (fc, &tree->quadrants) != NULL,
                  "Write quadrants");

  p4est_file_close (fc);

  /* intialize read quadrant data array */
  sc_array_init (&read_data, sizeof (int));
  sc_array_resize (&read_data, p4est->local_num_quadrants);

  sc_array_init (&quads, sizeof (p4est_quadrant_t));
  sc_array_resize (&quads, p4est->local_num_quadrants);

  fc = p4est_file_open_read (p4est, "test_io.out", header_size, read_header);
  SC_CHECK_ABORT (fc != NULL, "Open read");

  /* check read header */
  SC_CHECK_ABORT (read_header[0] == HEADER_INT1
                  && read_header[1] == HEADER_INT2,
                  "Read user-defined header");

  /* read the first data array */
  SC_CHECK_ABORT (p4est_file_read (fc, &read_data) != NULL, "Read ranks");

  /* read the second data array */
  SC_CHECK_ABORT (p4est_file_read (fc, &quads) != NULL, "Read quadrants");

  /* check the read data */
  for (i = 0; i < p4est->local_num_quadrants; ++i) {
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
  SC_CHECK_ABORT (fc != NULL, "Open append");

  sc_array_init (&unaligned, 3 * sizeof (char));
  sc_array_resize (&unaligned, p4est->local_num_quadrants);

  for (i = 0; i < p4est->local_num_quadrants; ++i) {
    current_char = sc_array_index (&unaligned, i);
    current_char[0] = 'a';
    current_char[1] = 'b';
    current_char[2] = 'c';
  }

  SC_CHECK_ABORT (p4est_file_write (fc, &unaligned) != NULL,
                  "Write unaligned");

  p4est_file_close (fc);

  sc_array_init (&elem_size, sizeof (size_t));
  p4est_file_info (p4est, "test_io.out", &global_quad_num, &read_header_size,
                   &elem_size);
  P4EST_GLOBAL_PRODUCTIONF
    ("file info: number of global quadrants = %ld\n, number of arrays = %ld,\nheader_size = %ld\n",
     global_quad_num, elem_size.elem_count, read_header_size);
  for (si = 0; si < elem_size.elem_count; ++si) {
    current_elem_size = sc_array_index (&elem_size, si);
    P4EST_GLOBAL_PRODUCTIONF ("Array %ld: element size %ld\n", si,
                              *current_elem_size);
  }

  fc = p4est_file_open_read (p4est, "test_io.out", header_size, read_header);

  /* skip two data arrays */
  SC_CHECK_ABORT (p4est_file_read (fc, NULL) == NULL, "Read skip 1");
  SC_CHECK_ABORT (p4est_file_read (fc, NULL) == NULL, "Read skip 2");
  SC_CHECK_ABORT (p4est_file_read (fc, &unaligned) != NULL, "Read unaligned");

  for (i = 0; i < p4est->local_num_quadrants; ++i) {
    current_char = sc_array_index (&unaligned, i);
    SC_CHECK_ABORT (current_char[0] == 'a' &&
                    current_char[1] == 'b' &&
                    current_char[2] == 'c', "Read after array padding");
  }

  p4est_file_close (fc);

  /* clean up */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);
  sc_array_reset (&quad_data);
  sc_array_reset (&read_data);
  sc_array_reset (&elem_size);
  sc_array_reset (&quads);
  sc_array_reset (&unaligned);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
//#endif /* !ENABLE_MPIIO */

  return 0;
}
