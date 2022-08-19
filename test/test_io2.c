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
#include <p4est_io.h>
#include <p4est_extended.h>
#include <p4est_bits.h>
#else
#include <p8est_io.h>
#include <p8est_extended.h>
#include <p8est_bits.h>
#endif
#include <sc_options.h>

#define HEADER_INT1 42
#define HEADER_INT2 84

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
  return quadrant->x == 0 && quadrant->y == 0
#ifdef P4_TO_P8
    && quadrant->z == 0
#endif
    ;
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

static void
parse_file_metadata (p4est_t * p4est, char *filename)
{
  int                 mpiret, ecode, eclass, msglen;
  sc_array_t          data_sizes;
  char                msg[sc_MPI_MAX_ERROR_STRING];
  size_t              read_header_size;

  P4EST_GLOBAL_PRODUCTIONF ("Parse %s\n", filename);

  sc_array_init (&data_sizes, sizeof (size_t));
  p4est_file_info (p4est, filename, &read_header_size, &data_sizes, &ecode);
  mpiret = p4est_file_error_class (ecode, &eclass);
  SC_CHECK_MPI (mpiret);
  mpiret = p4est_file_error_string (eclass, msg, &msglen);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_LERRORF ("file_info of %s at %s:%d: %s\n",
                        filename, __FILE__, __LINE__, msg);
  P4EST_GLOBAL_PRODUCTIONF
    ("file info (%s): number of global quadrants = %ld, number of arrays = %lld, header_size = %ld\n",
     filename, p4est->global_num_quadrants,
     (unsigned long long) data_sizes.elem_count, read_header_size);
  sc_array_reset (&data_sizes);
}

/** Write some invalid files in serial to the disk to check
 * the error handling of the p4est_file_* functions.
 */
static void
write_invalid_files (p4est_t * p4est)
{
  if (p4est->mpirank == 0) {
    char                string0[P4EST_NUM_METADATA_BYTES + 1];
    FILE               *file;
    int                 ret;

    /* invalid0 */
    snprintf (string0, P4EST_NUM_METADATA_BYTES + 1,
              "%.7s\n%-23s\n%.15ld\n%.15ld\n", "p4data1",
              p4est_version (), p4est->global_num_quadrants, 8L);
    string0[P4EST_NUM_METADATA_BYTES] = '\0';

    file = fopen ("invaild0." P4EST_DATA_FILE_EXT, "w");
    ret = fprintf (file, "%s", string0);
    if ((size_t) ret != strlen (string0)) {
      P4EST_LERROR ("Could not write invaild0." P4EST_DATA_FILE_EXT);
    }
    fclose (file);

    /* invalid1 */
    snprintf (string0, P4EST_NUM_METADATA_BYTES + 1,
              "%.7s\n%-22s\n%.15ld\n%.15ld\n", P4EST_MAGIC_NUMBER,
              "A wrong version string", p4est->global_num_quadrants, 8L);
    string0[P4EST_NUM_METADATA_BYTES] = '\0';

    file = fopen ("invaild1." P4EST_DATA_FILE_EXT, "w");
    ret = fprintf (file, "%s", string0);
    if ((size_t) ret != strlen (string0)) {
      P4EST_LERROR ("Could not write invaild1." P4EST_DATA_FILE_EXT);
    }
    fclose (file);

    /* invalid2 */
    snprintf (string0, P4EST_NUM_METADATA_BYTES + 1,
              "%.7s\n%-23s\n%.15ld\n%.15ld\n", P4EST_MAGIC_NUMBER,
              p4est_version (), 42L, 8L);
    string0[P4EST_NUM_METADATA_BYTES] = '\0';

    file = fopen ("invaild2." P4EST_DATA_FILE_EXT, "w");
    ret = fprintf (file, "%s", string0);
    if ((size_t) ret != strlen (string0)) {
      P4EST_LERROR ("Could not write invaild2." P4EST_DATA_FILE_EXT);
    }
    fclose (file);

    /* invalid3 */
    snprintf (string0, P4EST_NUM_METADATA_BYTES + 1,
              "%.7s\n%-23s\n%.15ld\n%.14ld\n", P4EST_MAGIC_NUMBER,
              p4est_version (), p4est->global_num_quadrants, -8L);
    string0[P4EST_NUM_METADATA_BYTES] = '\0';

    file = fopen ("invaild3." P4EST_DATA_FILE_EXT, "w");
    ret = fprintf (file, "%s", string0);
    if ((size_t) ret != strlen (string0)) {
      P4EST_LERROR ("Could not write invaild3." P4EST_DATA_FILE_EXT);
    }
    fclose (file);
  }
  parse_file_metadata (p4est, "invaild0." P4EST_DATA_FILE_EXT);
  /* leave the version string length missmatch out due to read byte count abort */
#ifndef P4EST_ENABLE_DEBUG
  parse_file_metadata (p4est, "invaild1." P4EST_DATA_FILE_EXT);
#endif
  parse_file_metadata (p4est, "invaild2." P4EST_DATA_FILE_EXT);
  parse_file_metadata (p4est, "invaild3." P4EST_DATA_FILE_EXT);
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret, errcode, errclass;
  int                 rank, size;
  int                 level = 3;
  int                 empty_header, read_only, header_only;
  int                *current;
  char               *current_char;
  int                 header_size = 8;
  size_t              si, read_header_size;
  size_t              current_elem_size;
  int                 header[2], read_header[2];
  char                msg[sc_MPI_MAX_ERROR_STRING];
  int                 msglen;
  unsigned            checksum;
  p4est_tree_t       *tree;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  p4est_file_context_t *fc, *fc1;
  p4est_locidx_t      i;
  sc_array_t          quad_data;
  sc_array_t          read_data;
  sc_array_t          elem_size;
  sc_array_t          quads;
  sc_array_t          unaligned;
  sc_options_t       *opt;
  char                current_user_string[P4EST_NUM_USER_STRING_BYTES];

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

  /* parse command line parameters */
  opt = sc_options_new (argv[0]);
  sc_options_add_bool (opt, 'E', "empty header", &empty_header, 0,
                       "Write no user header");
  sc_options_add_bool (opt, 'R', "read only", &read_only, 0, "Only reading");
  sc_options_add_bool (opt, 'H', "header only", &header_only, 0,
                       "Write only the header");
  sc_options_parse (p4est_package_id, SC_LP_DEFAULT, opt, argc, argv);

  if (empty_header) {
    header_size = 0;
  }

#ifndef P4_TO_P8
  connectivity = p4est_connectivity_new_unitsquare ();
#else
  connectivity = p8est_connectivity_new_unitcube ();
#endif

  p4est = p4est_new_ext (mpicomm, connectivity, 0, level, 1, 0, NULL, NULL);

  /* Test the data array padding by provoking a number of quadrants that
   * is not divisible by \ref P4EST_BYTE_DIV
   */
  p4est_refine (p4est, 1, refine, NULL);

  write_invalid_files (p4est);

  /* initialize the header */
  write_header (header);

  if (!header_only) {
    /* intialize quadrant data array */
    sc_array_init (&quad_data, sizeof (int));
    sc_array_resize (&quad_data, p4est->local_num_quadrants);
  }

  tree = p4est_tree_array_index (p4est->trees, 0);
  sc_array_init (&unaligned, 3 * sizeof (char));
  sc_array_resize (&unaligned, p4est->local_num_quadrants);
  if (!read_only) {
    fc =
      p4est_file_open_create (p4est, "test_io." P4EST_DATA_FILE_EXT,
                              header_size, header, &errcode);
    SC_CHECK_ABORT (fc != NULL, "Open create");

    if (!header_only) {
      write_rank (p4est, &quad_data);
      SC_CHECK_ABORT (p4est_file_write_field
                      (fc, &quad_data, "Quadrant-wise rank data",
                       &errcode) != NULL, "Write ranks");

      SC_CHECK_ABORT (p4est_file_write_field
                      (fc, &tree->quadrants, "Quadrant data", &errcode)
                      != NULL, "Write quadrants");

      for (i = 0; i < p4est->local_num_quadrants; ++i) {
        current_char = (char *) sc_array_index (&unaligned, i);
        current_char[0] = 'a';
        current_char[1] = 'b';
        current_char[2] = 'c';
      }

      SC_CHECK_ABORT (p4est_file_write_field
                      (fc, &unaligned, "Data that needs to be padded",
                       &errcode) != NULL, "Write unaligned");

      SC_CHECK_ABORT (p4est_file_write_header
                      (fc, (size_t) header_size, header, "Header as a block",
                       &errcode), "Write header");

      checksum = p4est_checksum (p4est);

/* *INDENT-OFF* */
      SC_CHECK_ABORT (p4est_file_write_header (fc, sizeof (unsigned),
                                               &checksum, "p4est checksum",
                                               &errcode) != NULL,
                                               "Write forest checksum");
/* *INDENT-ON* */

    }

    SC_CHECK_ABORT (p4est_file_close (fc, &errcode) == 0,
                    "Close file context 1");
  }

  if (!header_only) {
    /* intialize read quadrant data array */
    sc_array_init (&read_data, sizeof (int));
    sc_array_resize (&read_data, p4est->local_num_quadrants);

    sc_array_init (&quads, sizeof (p4est_quadrant_t));
    sc_array_resize (&quads, p4est->local_num_quadrants);

    fc =
      p4est_file_open_read (p4est, "test_io." P4EST_DATA_FILE_EXT,
                            header_size, read_header, &errcode);
    SC_CHECK_ABORT (fc != NULL, "Open read 1");

    /* Try to open an already opened file to test the error code/class. */
    fc1 =
      p4est_file_open_read (p4est, "test_iot." P4EST_DATA_FILE_EXT,
                            header_size, read_header, &errcode);
    mpiret = p4est_file_error_class (errcode, &errclass);
    SC_CHECK_MPI (mpiret);
    mpiret = p4est_file_error_string (errclass, msg, &msglen);
    SC_CHECK_MPI (mpiret);
    P4EST_GLOBAL_LERRORF ("Intended error by opening a non-existing"
                          " file (but we can not gurantee non-existence)"
                          " at %s:%d: %s\n", __FILE__, __LINE__, msg);
    if (fc1 != NULL) {
      /* the file seems to be existent by accident */
      SC_CHECK_ABORT (p4est_file_close (fc1, &errcode) == 0,
                      "Close accidently opened file");
    }

    if (!empty_header) {
      /* check read header */
      SC_CHECK_ABORT (read_header[0] == HEADER_INT1
                      && read_header[1] == HEADER_INT2,
                      "Read user-defined header");
    }

    /* read the first data array */
    SC_CHECK_ABORT (p4est_file_read_field
                    (fc, &read_data, current_user_string, &errcode) != NULL,
                    "Read ranks");
    P4EST_GLOBAL_PRODUCTIONF ("Read data with user string: %s\n",
                              current_user_string);

    /* read the second data array */
    SC_CHECK_ABORT (p4est_file_read_field
                    (fc, &quads, current_user_string, &errcode) != NULL,
                    "Read quadrants");
    P4EST_GLOBAL_PRODUCTIONF ("Read data with user string: %s\n",
                              current_user_string);

    SC_CHECK_ABORT (p4est_file_close (fc, &errcode) == 0,
                    "Close file context 2");

    /* check the read data */
    for (i = 0; i < p4est->local_num_quadrants; ++i) {
      SC_CHECK_ABORT (p4est_quadrant_is_equal
                      (p4est_quadrant_array_index (&quads, i),
                       p4est_quadrant_array_index (&tree->quadrants, i)),
                      "Quadrant read");
    }

    /* check read data of the first array */
    for (i = 0; i < p4est->local_num_quadrants; ++i) {
      current = (int *) sc_array_index (&read_data, i);
      SC_CHECK_ABORT (*current == p4est->mpirank, "Rank read");
    }

  }

  sc_array_init (&elem_size, sizeof (size_t));
  SC_CHECK_ABORT (p4est_file_info
                  (p4est, "test_io." P4EST_DATA_FILE_EXT, &read_header_size,
                   &elem_size, &errcode) == sc_MPI_SUCCESS, "Get file info");
  P4EST_GLOBAL_PRODUCTIONF
    ("file info: number of global quadrants = %ld, number of arrays = %lld, header_size = %ld\n",
     p4est->global_num_quadrants, (unsigned long long) elem_size.elem_count,
     read_header_size);
  for (si = 0; si < elem_size.elem_count; ++si) {
    current_elem_size = *(size_t *) sc_array_index (&elem_size, si);
    P4EST_GLOBAL_PRODUCTIONF ("Array %ld: element size %ld\n", si,
                              current_elem_size);
  }

  if (!header_only) {
    /* zero unaligned array */
    for (i = 0; i < p4est->local_num_quadrants; ++i) {
      current_char = (char *) sc_array_index (&unaligned, i);
      current_char[0] = '\0';
      current_char[1] = '\0';
      current_char[2] = '\0';
    }
  }

  if (!header_only) {
    fc =
      p4est_file_open_read (p4est, "test_io." P4EST_DATA_FILE_EXT,
                            header_size, read_header, &errcode);
    SC_CHECK_ABORT (fc != NULL, "Open read 2");

    /* skip two data arrays */
    SC_CHECK_ABORT (p4est_file_read_field (fc, NULL, NULL, &errcode) == NULL
                    && !errcode, "Read skip 1");
    SC_CHECK_ABORT (p4est_file_read_field (fc, NULL, NULL, &errcode) == NULL
                    && !errcode, "Read skip 2");
    SC_CHECK_ABORT (p4est_file_read_field
                    (fc, &unaligned, current_user_string, &errcode) != NULL
                    && !errcode, "Read unaligned");
    P4EST_GLOBAL_PRODUCTIONF ("Read data with user string: %s\n",
                              current_user_string);

    for (i = 0; i < p4est->local_num_quadrants; ++i) {
      current_char = (char *) sc_array_index (&unaligned, i);
      SC_CHECK_ABORT (current_char[0] == 'a' &&
                      current_char[1] == 'b' &&
                      current_char[2] == 'c', "Read after array padding");
    }

    read_header[0] = -1;
    read_header[1] = -1;
    SC_CHECK_ABORT (p4est_file_read_header
                    (fc, header_size, read_header, current_user_string,
                     &errcode)
                    != NULL, "Read header block");
    P4EST_GLOBAL_PRODUCTIONF ("Read header with user string: %s\n",
                              current_user_string);
    /* check read content of the header block */
    SC_CHECK_ABORT (read_header[0] == 42
                    && read_header[1] == 84, "Read header block");

    SC_CHECK_ABORT (p4est_file_close (fc, &errcode) == 0,
                    "Close file context 3");

    /* read and check the forest checksum */
    fc =
      p4est_file_open_read (p4est, "test_io." P4EST_DATA_FILE_EXT,
                            header_size, read_header, &errcode);
    SC_CHECK_ABORT (fc != NULL, "Open read 3");

    /* skip three data fields and one header block */
    SC_CHECK_ABORT (p4est_file_read_field (fc, NULL, NULL, &errcode) == NULL
                    && !errcode, "Read skip 1");
    SC_CHECK_ABORT (p4est_file_read_field (fc, NULL, NULL, &errcode) == NULL
                    && !errcode, "Read skip 2");
    SC_CHECK_ABORT (p4est_file_read_field (fc, NULL, NULL, &errcode) == NULL
                    && !errcode, "Read skip 3");
    SC_CHECK_ABORT (p4est_file_read_header (fc, 0, NULL, NULL, &errcode) ==
                    NULL && !errcode, "Read skip header 4");

    /* read the header containing the forest checksum */
    checksum = 1;
    SC_CHECK_ABORT (p4est_file_read_header
                    (fc, sizeof (unsigned), &checksum, current_user_string,
                     &errcode) != NULL, "Read checksum");
    P4EST_GLOBAL_PRODUCTIONF ("Read header data with user string: %s\n",
                              current_user_string);

    /* check the checksum */
    SC_CHECK_ABORT (p4est_checksum (p4est) ==
                    ((p4est->mpirank == 0) ? checksum : 0),
                    "Forest checksum equality");

    SC_CHECK_ABORT (p4est_file_close (fc, &errcode) == 0,
                    "Close file context 4");
  }

  /* clean up */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  if (!header_only) {
    sc_array_reset (&quad_data);
    sc_array_reset (&read_data);
    sc_array_reset (&quads);
    sc_array_reset (&unaligned);
  }
  sc_array_reset (&elem_size);
  sc_options_destroy (opt);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
