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
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_io.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_io.h>
#endif
#include <sc_search.h>
#include <sc.h>
#ifndef P4EST_ENABLE_MPII0
#include <errno.h>
#endif

sc_array_t         *
p4est_deflate_quadrants (p4est_t * p4est, sc_array_t ** data)
{
  const size_t        qsize = sizeof (p4est_qcoord_t);
  const size_t        dsize = p4est->data_size;
  size_t              qtreez, qz;
  sc_array_t         *qarr, *darr;
  p4est_topidx_t      tt;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  p4est_qcoord_t     *qap;
  char               *dap;

  qarr = sc_array_new_size (qsize,
                            (P4EST_DIM + 1) * p4est->local_num_quadrants);
  qap = (p4est_qcoord_t *) qarr->array;
  darr = NULL;
  dap = NULL;
  if (data != NULL) {
    P4EST_ASSERT (dsize > 0);
    darr = sc_array_new_size (dsize, p4est->local_num_quadrants);
    dap = darr->array;
  }
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    qtreez = tree->quadrants.elem_count;
    for (qz = 0; qz < qtreez; ++qz) {
      q = p4est_quadrant_array_index (&tree->quadrants, qz);
      *qap++ = q->x;
      *qap++ = q->y;
#ifdef P4_TO_P8
      *qap++ = q->z;
#endif
      *qap++ = (p4est_qcoord_t) q->level;
      if (data != NULL) {
        memcpy (dap, q->p.user_data, dsize);
        dap += dsize;
      }
    }
  }
  P4EST_ASSERT ((void *) qap ==
                qarr->array + qarr->elem_size * qarr->elem_count);
  if (data != NULL) {
    P4EST_ASSERT (dap == darr->array + darr->elem_size * darr->elem_count);
    *data = darr;
  }
  return qarr;
}

p4est_t            *
p4est_inflate (sc_MPI_Comm mpicomm, p4est_connectivity_t * connectivity,
               const p4est_gloidx_t * global_first_quadrant,
               const p4est_gloidx_t * pertree,
               sc_array_t * quadrants, sc_array_t * data, void *user_pointer)
{
  const p4est_gloidx_t *gfq;
  int                 i;
  int                 num_procs, rank;
  p4est_topidx_t      num_trees, jt;
  p4est_gloidx_t      gkey, gtreeskip, gtreeremain, gquadremain;
  p4est_t            *p4est;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
#ifdef P4EST_ENABLE_DEBUG
  int                 p;
#endif
  int8_t              ql, tml;
  size_t              dsize;
  size_t              gk1, gk2;
  size_t              qz, zqoffset, zqthistree;
  p4est_qcoord_t     *qap;
  char               *dap;

  P4EST_GLOBAL_PRODUCTION ("Into " P4EST_STRING "_inflate\n");
  p4est_log_indent_push ();

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));
  P4EST_ASSERT (global_first_quadrant != NULL);
  P4EST_ASSERT (pertree != NULL);
  P4EST_ASSERT (quadrants != NULL);
  P4EST_ASSERT (quadrants->elem_size == sizeof (p4est_qcoord_t));
  /* data may be NULL, in this case p4est->data_size will be 0 */
  /* user_pointer may be anything, we don't look at it */

  /* create p4est object and assign some data members */
  p4est = P4EST_ALLOC_ZERO (p4est_t, 1);
  dsize = p4est->data_size = (data == NULL ? 0 : data->elem_size);
  dap = (char *) (data == NULL ? NULL : data->array);
  qap = (p4est_locidx_t *) quadrants->array;
  p4est->user_pointer = user_pointer;
  p4est->connectivity = connectivity;
  num_trees = connectivity->num_trees;

  /* set parallel environment */
  p4est_comm_parallel_env_assign (p4est, mpicomm);
  num_procs = p4est->mpisize;
  rank = p4est->mpirank;

  /* create global first quadrant offsets */
  gfq = p4est->global_first_quadrant =
    P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  memcpy (p4est->global_first_quadrant, global_first_quadrant,
          (num_procs + 1) * sizeof (p4est_gloidx_t));
#ifdef P4EST_ENABLE_DEBUG
  P4EST_ASSERT (gfq[0] == 0);
  for (p = 0; p < num_procs; ++p) {
    P4EST_ASSERT (gfq[p] <= gfq[p + 1]);
  }
  P4EST_ASSERT (pertree[0] == 0);
  for (jt = 0; jt < num_trees; ++jt) {
    P4EST_ASSERT (pertree[jt] <= pertree[jt + 1]);
  }
  P4EST_ASSERT (gfq[num_procs] == pertree[num_trees]);
#endif
  gquadremain = gfq[rank + 1] - gfq[rank];
  p4est->local_num_quadrants = (p4est_locidx_t) gquadremain;
  p4est->global_num_quadrants = gfq[num_procs];
  P4EST_ASSERT (quadrants->elem_count ==
                (P4EST_DIM + 1) * (size_t) p4est->local_num_quadrants);
  P4EST_ASSERT (data == NULL || data->elem_count ==
                (size_t) p4est->local_num_quadrants);

  /* allocate memory pools */
  if (dsize > 0) {
    p4est->user_data_pool = sc_mempool_new (dsize);
  }
  else {
    p4est->user_data_pool = NULL;
  }
  p4est->quadrant_pool = sc_mempool_new (sizeof (p4est_quadrant_t));

  /* find the first and last tree on this processor */
  if (p4est->local_num_quadrants > 0) {
    gkey = gfq[rank];
    gk1 = sc_bsearch_range (&gkey, pertree, num_trees,
                            sizeof (p4est_gloidx_t), p4est_gloidx_compare);
    P4EST_ASSERT (gk1 < (size_t) num_trees);
    gtreeskip = gkey - pertree[gk1];
    gkey = gfq[rank + 1] - 1;
    gk2 = sc_bsearch_range (&gkey, pertree, num_trees,
                            sizeof (p4est_gloidx_t), p4est_gloidx_compare);
    P4EST_ASSERT (gk1 <= gk2 && gk2 < (size_t) num_trees);
    p4est->first_local_tree = (p4est_topidx_t) gk1;
    p4est->last_local_tree = (p4est_topidx_t) gk2;
  }
  else {
    gtreeskip = 0;
    p4est->first_local_tree = -1;
    p4est->last_local_tree = -2;
  }

  /* populate trees */
  zqoffset = 0;
  gquadremain = p4est->local_num_quadrants;
  p4est->trees = sc_array_new_size (sizeof (p4est_tree_t), num_trees);
  for (jt = 0; jt < num_trees; ++jt) {
    /* all trees need at least some basic setup */
    tree = p4est_tree_array_index (p4est->trees, jt);
    sc_array_init (&tree->quadrants, sizeof (p4est_quadrant_t));
    P4EST_QUADRANT_INIT (&tree->first_desc);
    P4EST_QUADRANT_INIT (&tree->last_desc);
    tree->quadrants_offset = (p4est_locidx_t) zqoffset;
    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }
    for (; i <= P4EST_MAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = -1;
    }
    q = NULL;
    tree->maxlevel = 0;
    if (jt >= p4est->first_local_tree && jt <= p4est->last_local_tree) {
      /* this tree has local quadrants */
      gtreeremain = pertree[jt + 1] - pertree[jt] - gtreeskip;
      P4EST_ASSERT (gtreeremain > 0 && gquadremain > 0);
      zqthistree = (size_t) SC_MIN (gtreeremain, gquadremain);
      P4EST_ASSERT (zqthistree > 0);
      sc_array_resize (&tree->quadrants, zqthistree);
      tml = 0;
      for (qz = 0; qz < zqthistree; ++qz) {
        q = p4est_quadrant_array_index (&tree->quadrants, qz);
        P4EST_QUADRANT_INIT (q);
        q->x = *qap++;
        q->y = *qap++;
#ifdef P4_TO_P8
        q->z = *qap++;
#endif
/* *INDENT-OFF* HORRIBLE indent bug */
        q->level = ql = (int8_t) *qap++;
/* *INDENT-ON* */
        P4EST_ASSERT (ql >= 0 && ql <= P4EST_QMAXLEVEL);
        ++tree->quadrants_per_level[ql];
        tml = SC_MAX (tml, ql);
        p4est_quadrant_init_data (p4est, jt, q, NULL);
        if (data != NULL) {
          memcpy (q->p.user_data, dap, dsize);
          dap += dsize;
        }
        if (qz == 0) {
          p4est_quadrant_first_descendant (q, &tree->first_desc,
                                           P4EST_QMAXLEVEL);
        }
      }
      p4est_quadrant_last_descendant (q, &tree->last_desc, P4EST_QMAXLEVEL);
      tree->maxlevel = tml;
      zqoffset += zqthistree;
      gquadremain -= (p4est_gloidx_t) zqthistree;
      gtreeskip = 0;
    }
  }
  P4EST_ASSERT (zqoffset == (size_t) p4est->local_num_quadrants);
  P4EST_ASSERT (gquadremain == 0);

  /* communicate partition information */
  p4est->global_first_position =
    P4EST_ALLOC (p4est_quadrant_t, num_procs + 1);
  p4est_comm_global_partition (p4est, NULL);

  /* print more statistics */
  P4EST_VERBOSEF ("total local quadrants %lld\n",
                  (long long) p4est->local_num_quadrants);

  P4EST_ASSERT (p4est->revision == 0);
  P4EST_ASSERT (p4est_is_valid (p4est));
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_inflate\n");

  return p4est;
}

/** The opaque file context for for p4est data files. */
struct p4est_file_context
{
  p4est_t            *p4est;            /**< corresponding p4est */
  size_t              num_calls;        /**< redundant but for convience;
                                            counts the number of calls of
                                            write and read, respectively */
#if 0
#ifndef P4EST_ENABLE_MPIIO
  const char         *filename; /**< We need to store the path for
                                   the successive opening strategy
                                   and to use the append mode to write
                                   data without MPI. */
  FILE               *file;     /**< The file pointer for serial IO */
#else
  sc_MPI_File         file;     /**< File variable with MPI IO. */
#endif
#endif
  sc_MPI_File         file;
  sc_MPI_Offset       accessed_bytes;   /**< count only array data bytes and
                                           array metadata bytes */
};

/** This function calculates a padding string consisting of spaces.
 * We require an already allocated array pad or NULL.
 * The number of bytes in pad must be at least divisor!
 * For NULL the function only calculates the number of padding bytes.
 */
static void
get_padding_string (size_t num_bytes, size_t divisor, char *pad,
                    size_t * num_pad_bytes)
{
  P4EST_ASSERT (divisor != 0 && num_pad_bytes != NULL);

  *num_pad_bytes = (divisor - (num_bytes % divisor)) % divisor;
  if (*num_pad_bytes == 0 || *num_pad_bytes == 1) {
    /* In these cases there is no space to add new line characters
     * but this is necessary to ensure a consistent layout in a text editor
     */
    *num_pad_bytes += divisor;
  }

  P4EST_ASSERT (*num_pad_bytes > 1);
  if (pad != NULL) {
    snprintf (pad, *num_pad_bytes + 1, "\n%-*s\n", (int) *num_pad_bytes - 2,
              "");
  }
}

static int
check_file_metadata (p4est_t * p4est, const char *filename, char *user_string,
                     char *metadata)
{
  long                read_global_num_quads;
  int                 error_flag;

  P4EST_ASSERT (metadata != NULL);

  error_flag = 0;

  /* check magic number */
  if (metadata[7] != '\n') {
    if (p4est->mpirank == 0) {
      P4EST_LERROR (P4EST_STRING
                    "_io: Error reading. Wrong file header format.\n");
    }
    return P4EST_ERR_IO;
  }

  metadata[7] = '\0';
  if (strcmp (metadata, P4EST_MAGIC_NUMBER)) {
    /* TODO: check for wrong endianness */
    P4EST_LERRORF (P4EST_STRING
                   "_io: Error reading <%s>. Wrong magic number (in file = %s, magic number = %s).\n",
                   filename, metadata, P4EST_MAGIC_NUMBER);
    error_flag = 1;
  }

  /* check format of version string line */
  if (metadata[31] != '\n') {
    if (p4est->mpirank == 0) {
      P4EST_LERROR (P4EST_STRING
                    "_io: Error reading. Wrong file header format.\n");
    }
    return P4EST_ERR_IO;
  }

  metadata[31] = '\0';
  if (strlen (&metadata[8]) != 23) {
    if (p4est->mpirank == 0) {
      P4EST_LERROR (P4EST_STRING
                    "_io: Error reading. Wrong file header format.\n");
    }
    return P4EST_ERR_IO;
  }

  /* check the format of the user string */
  if (metadata[47] != '\n') {
    if (p4est->mpirank == 0) {
      P4EST_LERROR (P4EST_STRING
                    "_io: Error reading. Wrong file header format.\n");
    }
    return P4EST_ERR_IO;
  }
  /* the content of the user string is not checked */
  strncpy (user_string, &metadata[32], 15);
  user_string[15] = '\0';

  /* check number of global quadrants */
  /* there is no \n at the end of this line */
  metadata[P4EST_NUM_METADATA_BYTES] = '\0';

  if (strlen (&metadata[48]) != 16) {
    if (p4est->mpirank == 0) {
      P4EST_LERROR (P4EST_STRING
                    "_io: Error reading. Wrong file header format.\n");
    }
    return P4EST_ERR_IO;
  }

  read_global_num_quads = sc_atol (&metadata[48]);
  if ((p4est_gloidx_t) read_global_num_quads != p4est->global_num_quadrants) {
    P4EST_LERRORF (P4EST_STRING
                   "_io: Error reading <%s>. Wrong global number of quadrants (in file = %ld, in given p4est = %ld).\n",
                   filename, read_global_num_quads,
                   p4est->global_num_quadrants);
    error_flag = 1;
  }

  return (error_flag) ? P4EST_ERR_IO : sc_MPI_SUCCESS;
}

int
p4est_file_error_cleanup (sc_MPI_File * file)
{
  /* no error checking since we are called under an error condition */
  P4EST_ASSERT (file != NULL);
#ifdef P4EST_ENABLE_MPIIO
  if (*file != sc_MPI_FILE_NULL) {
#else
  if ((*file)->file != sc_MPI_FILE_NULL) {
#endif
    /* We do not use here the libsc closing function since we do not perform
     * error checking in this functiont that is only called if we had already
     * an error.
     */
#ifdef P4EST_ENABLE_MPIIO
    MPI_File_close (file);
#else
    {
#ifdef P4EST_ENABLE_MPI
      int                 rank, mpiret;
#endif

#ifdef P4EST_ENABLE_MPI
      mpiret = sc_MPI_Comm_rank ((*file)->mpicomm, &rank);
      SC_CHECK_MPI (mpiret);

      if (rank == 0) {
#endif
        fclose ((*file)->file);
        (*file)->file = NULL;
#ifdef P4EST_ENABLE_MPI
      }
#endif
    }
#endif
  }
#ifndef P4EST_ENABLE_MPIIO
  SC_FREE (*file);
#endif
  return -1;
}

p4est_file_context_t *
p4est_file_open_create (p4est_t * p4est, const char *filename,
                        char user_string[15], int *errcode)
{
  int                 mpiret, count, count_error;
  /* We enforce the padding of the file header. */
  char                metadata[P4EST_NUM_METADATA_BYTES + P4EST_BYTE_DIV + 1];
  p4est_file_context_t *file_context = P4EST_ALLOC (p4est_file_context_t, 1);

  /* Open the file and create a new file if necessary */
  mpiret =
    sc_io_open (p4est->mpicomm, filename,
                SC_WRITE_CREATE, sc_MPI_INFO_NULL, &file_context->file);
  P4EST_FILE_CHECK_OPEN (mpiret, file_context, "File open create", errcode);

  if (p4est->mpirank == 0) {
    /* write padded p4est-defined header */
    snprintf (metadata, P4EST_NUM_METADATA_BYTES + P4EST_BYTE_DIV + 1,
              "%.7s\n%-23s\n%-15s\n%.16ld\n%-14s\n", P4EST_MAGIC_NUMBER,
              p4est_version (), user_string, p4est->global_num_quadrants, "");
    mpiret =
      sc_io_write_at (file_context->file, 0, metadata,
                      P4EST_NUM_METADATA_BYTES + P4EST_BYTE_DIV, sc_MPI_BYTE,
                      &count);
    P4EST_FILE_CHECK_MPI (mpiret, "Writing the file header");
    count_error = (P4EST_NUM_METADATA_BYTES + P4EST_BYTE_DIV != count);
    P4EST_FILE_CHECK_COUNT_SERIAL (P4EST_NUM_METADATA_BYTES + P4EST_BYTE_DIV,
                                   count);
  }

  P4EST_HANDLE_MPI_ERROR (mpiret, file_context, p4est->mpicomm, errcode);

  file_context->p4est = p4est;

  P4EST_HANDLE_MPI_COUNT_ERROR (count_error, file_context, errcode);

  file_context->accessed_bytes = 0;
  file_context->num_calls = 0;

  return file_context;
}

p4est_file_context_t *
p4est_file_open_read (p4est_t * p4est, const char *filename,
                      char *user_string, int *errcode)
{
  int                 mpiret;
  int                 count, count_error;
  char                metadata[P4EST_NUM_METADATA_BYTES + 1];
  p4est_file_context_t *file_context = P4EST_ALLOC (p4est_file_context_t, 1);

  /* Open the file in the reading mode */
  mpiret =
    sc_io_open (p4est->mpicomm, filename, SC_READ,
                sc_MPI_INFO_NULL, &file_context->file);
  P4EST_FILE_CHECK_OPEN (mpiret, file_context, "File open read", errcode);

  file_context->p4est = p4est;
  file_context->accessed_bytes = 0;
  file_context->num_calls = 0;

  /* read metadata and deallocate in case of error */
  if (file_context->p4est->mpirank == 0) {
    /* read metadata on rank 0 */
    mpiret =
      sc_io_read_at (file_context->file, 0, metadata,
                     P4EST_NUM_METADATA_BYTES, sc_MPI_BYTE, &count);

    P4EST_FILE_CHECK_MPI (mpiret, "Reading metadata");
    count_error = (P4EST_NUM_METADATA_BYTES != count);
    P4EST_FILE_CHECK_COUNT_SERIAL (P4EST_NUM_METADATA_BYTES, count);

    metadata[P4EST_NUM_METADATA_BYTES] = '\0';
    /* parse metadata; we do not use file_info because we do not want a Bcast */
    mpiret = check_file_metadata (p4est, filename, user_string, metadata);
    P4EST_FILE_CHECK_MPI (mpiret, "Check file header");
  }

  /* error checking */
  P4EST_HANDLE_MPI_ERROR (mpiret, file_context, p4est->mpicomm, errcode);
  P4EST_HANDLE_MPI_COUNT_ERROR (count_error, file_context, errcode);

  /* broadcast the user string of the file */
  mpiret = sc_MPI_Bcast (user_string, 16, sc_MPI_BYTE, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  return file_context;
}

p4est_file_context_t *
p4est_file_write_header (p4est_file_context_t * fc, size_t header_size,
                         const void *header_data, char user_string[47],
                         int *errcode)
{
  size_t              num_pad_bytes;
  char                header_metadata[P4EST_NUM_FIELD_HEADER_BYTES + 1],
    pad[P4EST_MAX_NUM_PAD_BYTES];
  int                 mpiret, count, count_error;

  P4EST_ASSERT (fc != NULL);
  P4EST_ASSERT (header_data != NULL);

  if (header_size == 0) {
    /* nothing to write */
    *errcode = sc_MPI_SUCCESS;
    return NULL;
  }

#ifdef P4EST_ENABLE_MPIIO
  /* set the file size */
  mpiret = MPI_File_set_size (fc->file,
                              P4EST_NUM_METADATA_BYTES + P4EST_BYTE_DIV +
                              header_size +
                              P4EST_NUM_FIELD_HEADER_BYTES +
                              fc->accessed_bytes);
  P4EST_FILE_CHECK_NULL (mpiret, fc, "Set file size", errcode);
#else
  /* We do not perform this optimization without MPI I/O */
#endif

  num_pad_bytes = 0;
  if (fc->p4est->mpirank == 0) {
    /* header-dependent metadata */
    snprintf (header_metadata,
              P4EST_NUM_FIELD_HEADER_BYTES +
              1, "H %.13ld\n%-47s\n", header_size, user_string);

    /* write header-dependent metadata */
    mpiret =
      sc_io_write_at (fc->file,
                      fc->accessed_bytes + P4EST_NUM_METADATA_BYTES +
                      P4EST_BYTE_DIV, header_metadata,
                      P4EST_NUM_FIELD_HEADER_BYTES, sc_MPI_BYTE, &count);

    P4EST_FILE_CHECK_MPI (mpiret, "Writing header metadata");
    count_error = (P4EST_NUM_FIELD_HEADER_BYTES != count);
    P4EST_FILE_CHECK_COUNT_SERIAL (P4EST_NUM_FIELD_HEADER_BYTES, count);
  }
  P4EST_HANDLE_MPI_ERROR (mpiret, fc, fc->p4est->mpicomm, errcode);
  P4EST_HANDLE_MPI_COUNT_ERROR (count_error, fc, errcode);

  /*write header data */
  if (fc->p4est->mpirank == 0) {
    mpiret =
      sc_io_write_at (fc->file,
                      fc->accessed_bytes + P4EST_NUM_METADATA_BYTES +
                      P4EST_BYTE_DIV + P4EST_NUM_FIELD_HEADER_BYTES,
                      header_data, header_size, sc_MPI_BYTE, &count);

    P4EST_FILE_CHECK_MPI (mpiret, "Writing header data");
    count_error = ((int) header_size != count);
    P4EST_FILE_CHECK_COUNT_SERIAL (header_size, count);

    /* write padding bytes */
    get_padding_string (header_size, P4EST_BYTE_DIV, pad, &num_pad_bytes);
    mpiret =
      sc_io_write_at (fc->file,
                      fc->accessed_bytes + P4EST_NUM_METADATA_BYTES +
                      P4EST_BYTE_DIV + P4EST_NUM_FIELD_HEADER_BYTES +
                      header_size, pad, num_pad_bytes, sc_MPI_BYTE, &count);
    P4EST_FILE_CHECK_MPI (mpiret, "Writing padding bytes for header data");
    count_error = ((int) num_pad_bytes != count);
    P4EST_FILE_CHECK_COUNT_SERIAL (num_pad_bytes, count);
  }
  else {
    get_padding_string (header_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
  }

  /* This is *not* the processor local value */
  fc->accessed_bytes +=
    header_size + P4EST_NUM_FIELD_HEADER_BYTES + num_pad_bytes;
  ++fc->num_calls;

  return fc;
}

/** Collectivly read and check block metadata. */
static p4est_file_context_t *
read_block_metadata (p4est_file_context_t * fc, size_t * read_data_size,
                     size_t data_size, char block_type, char *user_string,
                     int *errcode)
{
  int                 mpiret, count, count_error;
  int                 bytes_to_read;
  char                block_metadata[P4EST_NUM_FIELD_HEADER_BYTES];

  P4EST_ASSERT (read_data_size != NULL);

  bytes_to_read =
    (user_string !=
     NULL) ? P4EST_NUM_FIELD_HEADER_BYTES : (P4EST_NUM_ARRAY_METADATA_BYTES +
                                             2);
  if (fc->p4est->mpirank == 0) {
    mpiret = sc_io_read_at (fc->file,
                            fc->accessed_bytes +
                            P4EST_NUM_METADATA_BYTES +
                            P4EST_BYTE_DIV, block_metadata,
                            bytes_to_read, sc_MPI_BYTE, &count);
    P4EST_FILE_CHECK_MPI (mpiret, "Reading block-wise metadata");
    count_error = (bytes_to_read != count);
    P4EST_FILE_CHECK_COUNT_SERIAL (bytes_to_read, count);
  }
  /* In the case of error the return value is still NULL */
  P4EST_HANDLE_MPI_ERROR (mpiret, fc, fc->p4est->mpicomm, errcode);
  P4EST_HANDLE_MPI_COUNT_ERROR (count_error, fc, errcode);

  /* broadcast block metadata to calculate correct internals on each rank */
  mpiret =
    sc_MPI_Bcast (block_metadata, bytes_to_read,
                  sc_MPI_BYTE, 0, fc->p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  /* check for given block specifying character */
  if (block_metadata[0] != block_type) {
    if (fc->p4est->mpirank == 0) {
      P4EST_LERROR (P4EST_STRING "_io: Error reading. Wrong block type.\n");
    }
    p4est_file_error_cleanup (&fc->file);
    P4EST_FREE (fc);
    *errcode = P4EST_ERR_IO;
    return NULL;
  }

  /* check '\n' to check the format */
  if (block_metadata[P4EST_NUM_ARRAY_METADATA_BYTES + 1] != '\n') {
    if (fc->p4est->mpirank == 0) {
      P4EST_LERROR (P4EST_STRING
                    "_io: Error reading. Wrong block header format.\n");
    }
    p4est_file_error_cleanup (&fc->file);
    P4EST_FREE (fc);
    *errcode = P4EST_ERR_IO;
    return NULL;
  }

  /* process the block metadata */
  block_metadata[P4EST_NUM_ARRAY_METADATA_BYTES + 1] = '\0';
  /* we cut off the block type specifier */
  *read_data_size = sc_atol (&block_metadata[2]);

  if (user_string != NULL && *read_data_size != data_size) {
    if (fc->p4est->mpirank == 0) {
      P4EST_LERRORF (P4EST_STRING
                     "_io: Error reading. Wrong array element size (in file = %ld, by parameter = %ld).\n",
                     *read_data_size, data_size);
    }
    p4est_file_error_cleanup (&fc->file);
    P4EST_FREE (fc);
    *errcode = P4EST_ERR_IO;
    return NULL;
  }

  if (user_string != NULL) {
    /* check '\n' to check the format */
    if (block_metadata[P4EST_NUM_FIELD_HEADER_BYTES - 1] != '\n') {
      if (fc->p4est->mpirank == 0) {
        P4EST_LERROR (P4EST_STRING
                      "_io: Error reading. Wrong block header format.\n");
      }
      p4est_file_error_cleanup (&fc->file);
      P4EST_FREE (fc);
      *errcode = P4EST_ERR_IO;
      return NULL;
    }

    /* null-terminate the user string of the current block */
    block_metadata[P4EST_NUM_FIELD_HEADER_BYTES - 1] = '\0';

    /* copy the user string, '\0' was already set above */
    strncpy (user_string, &block_metadata[P4EST_NUM_ARRAY_METADATA_BYTES + 2],
             P4EST_NUM_USER_STRING_BYTES);
    P4EST_ASSERT (user_string[P4EST_NUM_USER_STRING_BYTES - 1] == '\0');
  }

  return fc;
}

p4est_file_context_t *
p4est_file_read_header (p4est_file_context_t * fc,
                        size_t header_size, void *header_data,
                        char *user_string, int *errcode)
{
  int                 mpiret, count, count_error;
  size_t              num_pad_bytes, read_data_size;
#ifdef P4EST_ENABLE_MPIIO
  sc_MPI_Offset       size;
#endif

  num_pad_bytes = 0;
  if (header_data == NULL || header_size == 0) {
    /* Nothing to read but we shift our own file pointer */
    if (read_block_metadata (fc, &read_data_size, 0, 'H', NULL, errcode) ==
        NULL) {
      return NULL;
    }

    /* calculate the padding bytes for this header data */
    get_padding_string (header_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
    fc->accessed_bytes +=
      header_size + P4EST_NUM_FIELD_HEADER_BYTES + num_pad_bytes;
    ++fc->num_calls;
    *errcode = sc_MPI_SUCCESS;
    return NULL;
  }

#ifdef P4EST_ENABLE_MPIIO
  /* check file size; no sync required because the file size does not
   * change in the reading mode.
   */
  mpiret = MPI_File_get_size (fc->file, &size);
  P4EST_FILE_CHECK_NULL (mpiret, fc, "Get file size for read", errcode);
  if (((size_t) size) - P4EST_NUM_METADATA_BYTES - P4EST_BYTE_DIV -
      P4EST_NUM_FIELD_HEADER_BYTES < header_size) {
    /* report wrong file size, collectively close the file and deallocate fc */
    if (fc->p4est->mpirank == 0) {
      P4EST_LERROR (P4EST_STRING
                    "_io: Error reading. File has less bytes than the user wants to read.\n");
    }
    mpiret = p4est_file_close (fc, &mpiret);
    P4EST_FILE_CHECK_NULL (mpiret, fc,
                           P4EST_STRING "_file_read_data: close file",
                           errcode);
    return NULL;
  }
#else
  /* There is no C-standard functionality to get the file size */
#endif

  /* check the header metadata */
  if (read_block_metadata
      (fc, &read_data_size, header_size, 'H', user_string, errcode) == NULL) {
    return NULL;
  }

  if (fc->p4est->mpirank == 0) {
    mpiret = sc_io_read_at (fc->file, fc->accessed_bytes +
                            P4EST_NUM_METADATA_BYTES +
                            P4EST_NUM_FIELD_HEADER_BYTES +
                            P4EST_BYTE_DIV, header_data, (int) header_size,
                            sc_MPI_BYTE, &count);
    P4EST_FILE_CHECK_MPI (mpiret, "Reading header data");
    count_error = ((int) header_size != count);
    P4EST_FILE_CHECK_COUNT_SERIAL (header_size, count);
  }
  P4EST_HANDLE_MPI_ERROR (mpiret, fc, fc->p4est->mpicomm, errcode);
  P4EST_HANDLE_MPI_COUNT_ERROR (count_error, fc, errcode);

  /* broadcast the header data */
  mpiret =
    sc_MPI_Bcast (header_data, header_size, sc_MPI_BYTE, 0,
                  fc->p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  fc->accessed_bytes +=
    header_size + P4EST_NUM_FIELD_HEADER_BYTES + num_pad_bytes;
  ++fc->num_calls;

  return fc;
}

p4est_file_context_t *
p4est_file_write_field (p4est_file_context_t * fc, sc_array_t * quadrant_data,
                        char user_string[47], int *errcode)
{
  size_t              bytes_to_write, num_pad_bytes, array_size;
  char                array_metadata[P4EST_NUM_FIELD_HEADER_BYTES + 1],
    pad[P4EST_MAX_NUM_PAD_BYTES];
  sc_MPI_Offset       write_offset;
  int                 mpiret, count, count_error;

  P4EST_ASSERT (fc != NULL);
  P4EST_ASSERT (quadrant_data != NULL
                && quadrant_data->elem_count ==
                (size_t) fc->p4est->local_num_quadrants);

  if (quadrant_data->elem_size == 0) {
    /* nothing to write */
    *errcode = sc_MPI_SUCCESS;
    return NULL;
  }

  /* Check how many bytes we write to the disk */
  bytes_to_write = quadrant_data->elem_count * quadrant_data->elem_size;

  /* rank-dependent byte offset */
  write_offset = P4EST_NUM_METADATA_BYTES + P4EST_BYTE_DIV +
    fc->p4est->global_first_quadrant[fc->p4est->mpirank] *
    quadrant_data->elem_size;

#ifdef P4EST_ENABLE_MPIIO
  /* set the file size */
  mpiret = MPI_File_set_size (fc->file,
                              P4EST_NUM_METADATA_BYTES + P4EST_BYTE_DIV +
                              fc->p4est->global_num_quadrants *
                              quadrant_data->elem_size +
                              P4EST_NUM_FIELD_HEADER_BYTES +
                              fc->accessed_bytes);
  P4EST_FILE_CHECK_NULL (mpiret, fc, "Set file size", errcode);
#else
  /* We do not perform this optimization without MPI I/O */
#endif

  num_pad_bytes = 0;
  if (fc->p4est->mpirank == 0) {
    /* array-dependent metadata */
    snprintf (array_metadata,
              P4EST_NUM_FIELD_HEADER_BYTES +
              1, "F %.13ld\n%-47s\n", quadrant_data->elem_size, user_string);

    /* write array-dependent metadata */
    mpiret =
      sc_io_write_at (fc->file, fc->accessed_bytes + write_offset,
                      array_metadata,
                      P4EST_NUM_FIELD_HEADER_BYTES, sc_MPI_BYTE, &count);

    P4EST_FILE_CHECK_MPI (mpiret, "Writing array metadata");
    count_error = (P4EST_NUM_FIELD_HEADER_BYTES != count);
    P4EST_FILE_CHECK_COUNT_SERIAL (P4EST_NUM_FIELD_HEADER_BYTES, count);
  }
  P4EST_HANDLE_MPI_ERROR (mpiret, fc, fc->p4est->mpicomm, errcode);
  P4EST_HANDLE_MPI_COUNT_ERROR (count_error, fc, errcode);

  /* write array data */
  mpiret =
    sc_io_write_at_all (fc->file,
                        fc->accessed_bytes + write_offset +
                        P4EST_NUM_FIELD_HEADER_BYTES, quadrant_data->array,
                        bytes_to_write, sc_MPI_BYTE, &count);
  P4EST_FILE_CHECK_NULL (mpiret, fc, "Writing quadrant-wise", errcode);
  P4EST_FILE_CHECK_COUNT (bytes_to_write, count, fc, errcode);

  /** We place the padding bytes write here because for the sequential
   * IO operations the order of fwrite calls plays a role.
   */
  /* write padding bytes */
  if (fc->p4est->mpirank == 0) {
    /* Caculate and write padding bytes for array data */
    array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
    get_padding_string (array_size, P4EST_BYTE_DIV, pad, &num_pad_bytes);

    mpiret =
      sc_io_write_at (fc->file,
                      fc->accessed_bytes + P4EST_NUM_METADATA_BYTES +
                      P4EST_BYTE_DIV + array_size +
                      P4EST_NUM_FIELD_HEADER_BYTES, pad, num_pad_bytes,
                      sc_MPI_BYTE, &count);
    P4EST_FILE_CHECK_MPI (mpiret, "Writing padding bytes for a data array");
    /* We do not need to call P4EST_FILE_HANDLE_MPI_ERROR in the next
     * collective line of code since P4EST_FILE_HANDLE_MPI_ERROR was already
     * called in this scope.
     */
    count_error = ((int) num_pad_bytes != count);
    P4EST_FILE_CHECK_COUNT_SERIAL (num_pad_bytes, count);
  }
  else {
    array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
    get_padding_string (array_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
  }

  /* This is *not* the processor local value */
  fc->accessed_bytes +=
    quadrant_data->elem_size * fc->p4est->global_num_quadrants +
    P4EST_NUM_FIELD_HEADER_BYTES + num_pad_bytes;
  ++fc->num_calls;

  return fc;
}

p4est_file_context_t *
p4est_file_read_field (p4est_file_context_t * fc, sc_array_t * quadrant_data,
                       char *user_string, int *errcode)
{
  int                 count;
  size_t              bytes_to_read, num_pad_bytes, array_size,
    read_data_size;
#ifdef P4EST_ENABLE_MPIIO
  sc_MPI_Offset       size;
#endif
  int                 mpiret;

  P4EST_ASSERT (fc != NULL);

  if (quadrant_data == NULL || quadrant_data->elem_size == 0) {
    /* Nothing to read but we shift our own file pointer */
    if (read_block_metadata (fc, &read_data_size, 0, 'F', NULL, errcode) ==
        NULL) {
      return NULL;
    }

    /* calculate the padding bytes for this data array */
    array_size = fc->p4est->global_num_quadrants * read_data_size;
    get_padding_string (array_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
    fc->accessed_bytes +=
      read_data_size * fc->p4est->global_num_quadrants +
      P4EST_NUM_FIELD_HEADER_BYTES + num_pad_bytes;
    ++fc->num_calls;
    *errcode = sc_MPI_SUCCESS;
    return NULL;
  }

  P4EST_ASSERT (quadrant_data->elem_count ==
                (size_t) fc->p4est->local_num_quadrants);

  /* check how many bytes we read from the disk */
  bytes_to_read = quadrant_data->elem_count * quadrant_data->elem_size;

#ifdef P4EST_ENABLE_MPIIO
  /* check file size; no sync required because the file size does not
   * change in the reading mode.
   */
  mpiret = MPI_File_get_size (fc->file, &size);
  P4EST_FILE_CHECK_NULL (mpiret, fc, "Get file size for read", errcode);
  if (((size_t) size) - P4EST_NUM_METADATA_BYTES - P4EST_BYTE_DIV -
      P4EST_NUM_FIELD_HEADER_BYTES < bytes_to_read) {
    /* report wrong file size, collectively close the file and deallocate fc */
    if (fc->p4est->mpirank == 0) {
      P4EST_LERROR (P4EST_STRING
                    "_io: Error reading. File has less bytes than the user wants to read.\n");
    }
    mpiret = p4est_file_close (fc, &mpiret);
    P4EST_FILE_CHECK_NULL (mpiret, fc,
                           P4EST_STRING "_file_read_data: close file",
                           errcode);
    return NULL;
  }
#else
  /* There is no C-standard functionality to get the file size */
#endif

  /* check the array metadata */
  if (read_block_metadata
      (fc, &read_data_size, quadrant_data->elem_size, 'F', user_string,
       errcode) == NULL) {
    return NULL;
  }

  /* calculate the padding bytes for this data array */
  array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
  get_padding_string (array_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);

  mpiret = sc_io_read_at_all (fc->file,
                              fc->accessed_bytes +
                              P4EST_NUM_METADATA_BYTES +
                              P4EST_NUM_FIELD_HEADER_BYTES +
                              P4EST_BYTE_DIV +
                              fc->p4est->global_first_quadrant[fc->p4est->
                                                               mpirank]
                              * quadrant_data->elem_size,
                              quadrant_data->array, bytes_to_read,
                              sc_MPI_BYTE, &count);

  P4EST_FILE_CHECK_NULL (mpiret, fc, "Reading quadrant-wise", errcode);
  P4EST_FILE_CHECK_COUNT (bytes_to_read, count, fc, errcode);

  fc->accessed_bytes +=
    quadrant_data->elem_size * fc->p4est->global_num_quadrants +
    P4EST_NUM_FIELD_HEADER_BYTES + num_pad_bytes;
  ++fc->num_calls;

  return fc;
}

int
p4est_file_info (p4est_t * p4est, const char *filename,
                 char *user_string, sc_array_t * elem_size, int *errcode)
{
  int                 mpiret, eclass;
  int                 retval;
  int                 count, count_error;
  long                long_header;
  size_t              current_size, num_pad_bytes;
  char                metadata[P4EST_NUM_METADATA_BYTES + 1];
  char                block_metadata[P4EST_NUM_FIELD_HEADER_BYTES + 1];
  p4est_file_block_metadata_t *current_member;
  sc_MPI_Offset       current_position;
  sc_MPI_File         file;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (filename != NULL);
  P4EST_ASSERT (user_string != NULL);
  P4EST_ASSERT (elem_size != NULL);
  P4EST_ASSERT (elem_size->elem_size == sizeof (p4est_file_block_metadata_t));

  /* set default output values */
  sc_array_reset (elem_size);

  /* open the file in reading mode */
  *errcode = eclass = sc_MPI_SUCCESS;   /* MPI defines MPI_SUCCESS to equal 0. */
  file = sc_MPI_FILE_NULL;

  if ((retval =
       sc_io_open (p4est->mpicomm, filename, SC_READ,
                   sc_MPI_INFO_NULL, &file)) != sc_MPI_SUCCESS) {
    mpiret = sc_io_error_class (retval, &eclass);
    SC_CHECK_MPI (mpiret);
  }

  if (eclass) {
    *errcode = eclass;
    SC_FREE (file);
    return -1;
  }

  /* read file metadata on root rank */
  P4EST_ASSERT (!eclass);
  if (p4est->mpirank == 0) {
    if ((retval = sc_io_read_at (file, 0, metadata,
                                 P4EST_NUM_METADATA_BYTES, sc_MPI_BYTE,
                                 &count))
        != sc_MPI_SUCCESS) {
      mpiret = sc_io_error_class (retval, &eclass);
      SC_CHECK_MPI (mpiret);
      *errcode = eclass;
      /* There is no count error for a non-successful read. */
      count_error = 0;
    }
    else {
      count_error = (P4EST_NUM_METADATA_BYTES != count);
    }
  }
  mpiret = sc_MPI_Bcast (&eclass, 1, sc_MPI_INT, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  if (eclass) {
    p4est_file_error_cleanup (&file);
    return -1;
  }
  mpiret = sc_MPI_Bcast (&count_error, 1, sc_MPI_INT, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  if (count_error) {
    if (p4est->mpirank == 0) {
      P4EST_LERROR (P4EST_STRING
                    "_file_info: read count error for file metadata reading");
    }
    *errcode = P4EST_FILE_COUNT_ERROR;
    p4est_file_error_cleanup (&file);
    return -1;
  }

  /* broadcast file metadata to all ranks and null-terminate it */
  mpiret = sc_MPI_Bcast (metadata, P4EST_NUM_METADATA_BYTES, sc_MPI_BYTE, 0,
                         p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  metadata[P4EST_NUM_METADATA_BYTES] = '\0';

  if ((mpiret =
       check_file_metadata (p4est, filename, user_string,
                            metadata)) != sc_MPI_SUCCESS) {
    *errcode = P4EST_ERR_IO;
    return p4est_file_error_cleanup (&file);
  }

  current_position =
    (sc_MPI_Offset) (P4EST_NUM_METADATA_BYTES + P4EST_BYTE_DIV);

  /* read all data headers that we find and skip the data itself */
  if (p4est->mpirank == 0) {
    for (;;) {
      /* read block metadata for current record */
      mpiret = sc_io_read_at (file, current_position, block_metadata,
                              P4EST_NUM_FIELD_HEADER_BYTES, sc_MPI_BYTE,
                              &count);
      mpiret = sc_io_error_class (mpiret, &eclass);
      SC_CHECK_MPI (mpiret);
      *errcode = eclass;
      if (eclass) {
        return p4est_file_error_cleanup (&file);
      }
      if (P4EST_NUM_FIELD_HEADER_BYTES != count) {
        /* we did not read the correct number of bytes */
        break;
      }

      /* parse and store the element size, the block type and the user string */
      current_member =
        (p4est_file_block_metadata_t *) sc_array_push (elem_size);
      if (block_metadata[0] == 'H' || block_metadata[0] == 'F') {
        /* we want to read the block type */
        current_member->block_type = block_metadata[0];
      }
      else {
        /* the last entry is incomplete and is therefore removed */
        sc_array_rewind (elem_size, elem_size->elem_count - 1);
        /* current_member is freed if the whole array is freed */
        break;
      }

      /* check format */
      if (block_metadata[P4EST_NUM_ARRAY_METADATA_BYTES + 1] != '\n') {
        /* the last entry is incomplete and is therefore removed */
        sc_array_rewind (elem_size, elem_size->elem_count - 1);
        break;
      }

      /* process block metadata */
      block_metadata[P4EST_NUM_ARRAY_METADATA_BYTES + 1] = '\0';
      /* we cut off the block type specifier to read the data size */
      current_member->data_size = (size_t) sc_atol (&block_metadata[2]);

      /* read the user string */
      /* check '\n' to check the format */
      if (block_metadata[P4EST_NUM_FIELD_HEADER_BYTES - 1] != '\n') {
        /* the last entry is incomplete and is therefore removed */
        sc_array_rewind (elem_size, elem_size->elem_count - 1);
        break;
      }

      /* null-terminate the user string of the current block */
      block_metadata[P4EST_NUM_FIELD_HEADER_BYTES - 1] = '\0';

      /* copy the user string, '\0' was already set above */
      strncpy (current_member->user_string,
               &block_metadata[P4EST_NUM_ARRAY_METADATA_BYTES + 2],
               P4EST_NUM_USER_STRING_BYTES);
      P4EST_ASSERT (current_member->user_string
                    [P4EST_NUM_USER_STRING_BYTES - 1] == '\0');

      /* get padding bytes of the current block */
      if (current_member->block_type == 'F') {
        current_size =
          (size_t) (p4est->global_num_quadrants * current_member->data_size);
      }
      else if (current_member->block_type == 'H') {
        current_size = current_member->data_size;
      }
      else {
        /* \ref read_block_metadata checks for valid block type */
        SC_ABORT_NOT_REACHED ();
      }
      get_padding_string (current_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
      current_position +=
        P4EST_NUM_FIELD_HEADER_BYTES + current_size + num_pad_bytes;
    }
  }

  /* replicate block metadata in parallel */
  long_header = (long) elem_size->elem_count;   /* 0 on non-root */
  mpiret = sc_MPI_Bcast (&long_header, 1, sc_MPI_LONG, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  if (p4est->mpirank != 0) {
    sc_array_resize (elem_size, (size_t) long_header);
  }
  mpiret = sc_MPI_Bcast (elem_size->array,
                         elem_size->elem_count * elem_size->elem_size,
                         sc_MPI_BYTE, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  /* close the file with error checking */
  P4EST_ASSERT (!eclass);
#ifdef P4EST_ENABLE_MPIIO
  if ((retval = sc_MPI_File_close (&file)) != sc_MPI_SUCCESS) {
    mpiret = sc_io_error_class (retval, &eclass);
    SC_CHECK_MPI (mpiret);
  }
#else
  if (p4est->mpirank == 0) {
    errno = 0;
    if (fclose (file->file)) {
      mpiret = sc_io_error_class (errno, &eclass);
      SC_CHECK_MPI (mpiret);
    }
  }
  else {
    P4EST_ASSERT (file->file == NULL);
  }
  SC_FREE (file);
  mpiret = sc_MPI_Bcast (&eclass, 1, sc_MPI_INT, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
#endif
  return 0;
}

int
p4est_file_error_class (int errcode, int *errclass)
{
  /* count error exists only on p4est-level */
  if (errcode == P4EST_FILE_COUNT_ERROR) {
    *errclass = errcode;
    return sc_MPI_SUCCESS;
  }
  else {
    /* all other errcodes can be handeled by libsc */
    return sc_io_error_class (errcode, errclass);

  }
}

int
p4est_file_error_string (int errclass, char *string, int *resultlen)
{
  int                 retval;

  if (string == NULL || resultlen == NULL) {
    return sc_MPI_ERR_ARG;
  }

  /* count error exists only on p4est-level */
  if (errclass == P4EST_FILE_COUNT_ERROR) {
    if ((retval =
         snprintf (string, sc_MPI_MAX_ERROR_STRING, "%s",
                   "Read or write cout error")) < 0) {
      /* unless something goes against the current standard of snprintf */
      return sc_MPI_ERR_NO_MEM;
    }
    if (retval >= sc_MPI_MAX_ERROR_STRING) {
      retval = sc_MPI_MAX_ERROR_STRING - 1;
    }
    *resultlen = retval;
    return sc_MPI_SUCCESS;
  }
  else {
    /* all other errocodes can be handled by libsc */
    return sc_MPI_Error_string (errclass, string, resultlen);
  }
}

int
p4est_file_close (p4est_file_context_t * fc, int *errcode)
{
  P4EST_ASSERT (fc != NULL);

  int                 mpiret;

  mpiret = sc_io_close (&fc->file);
  P4EST_FILE_CHECK_INT (mpiret, "Close file", errcode);

  P4EST_FREE (fc);

  return 0;
}
