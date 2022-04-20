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
  size_t              header_size;      /**< the user-defined header */
  size_t              num_calls;        /**< redundant but for convience;
                                            counts the number of calls of
                                            write and read, respectivly */
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

  if (*num_pad_bytes > 0 && pad != NULL) {
    snprintf (pad, *num_pad_bytes + 1, "%-*s", (int) *num_pad_bytes, "");
  }
}

static int
check_file_metadata (p4est_t * p4est, size_t header_size,
                     const char *filename, char *metadata)
{
  long                read_global_num_quads, read_header_size;
  int                 count, error_flag;
  char               *parsing_arg;

  P4EST_ASSERT (metadata != NULL);

  count = 0;
  error_flag = 0;
  parsing_arg = strtok (metadata, "\n");
  if (parsing_arg == NULL) {
    P4EST_LERRORF (P4EST_STRING
                   "_io: Error reading <%s>. Could not parse the file.\n",
                   filename);
    error_flag = 1;
  }
  while (parsing_arg != NULL && count < 4) {
    if (count == 0) {
      if (strcmp (parsing_arg, P4EST_MAGIC_NUMBER)) {
        /* TODO: check for wrong endianess */
        P4EST_LERRORF (P4EST_STRING
                       "_io: Error reading <%s>. Wrong magic number (in file = %s, magic number = %s).\n",
                       filename, parsing_arg, P4EST_MAGIC_NUMBER);
        error_flag = 1;
      }
    }
    else if (count == 2) {
      read_global_num_quads = sc_atol (parsing_arg);
      if (read_global_num_quads != p4est->global_num_quadrants) {
        P4EST_LERRORF (P4EST_STRING
                       "_io: Error reading <%s>. Wrong global number of quadrants (in file = %ld, in given p4est = %ld).\n",
                       filename, read_global_num_quads,
                       p4est->global_num_quadrants);
        error_flag = 1;
      }
    }
    else if (count == 3) {
      read_header_size = sc_atol (parsing_arg);
      if (read_header_size < 0 || (size_t) read_header_size != header_size) {
        P4EST_LERRORF (P4EST_STRING
                       "_io: Error reading <%s>. Wrong header_size (in file = %ld, as parameter = %ld).\n",
                       filename, read_header_size, header_size);
        error_flag = 1;
      }
    }
    parsing_arg = strtok (NULL, "\n");
    ++count;
  }

  return error_flag;
}

p4est_file_context_t *
p4est_file_open_create (p4est_t * p4est, const char *filename,
                        size_t header_size, const void *header_data)
{
  int                 mpiret;
  char                metadata[P4EST_NUM_METADATA_BYTES + 1];
  char                pad[P4EST_BYTE_DIV];
  size_t              num_pad_bytes;
  p4est_file_context_t *file_context = P4EST_ALLOC (p4est_file_context_t, 1);

  /* Open the file and create a new file if necessary */
  mpiret =
    sc_mpi_file_open (p4est->mpicomm, filename,
                      sc_MPI_MODE_WRONLY_CREATE,
                      sc_MPI_INFO_NULL, &file_context->file);
  P4EST_FILE_CHECK_OPEN (mpiret, file_context, "File open create");

  num_pad_bytes = 0;
  if (p4est->mpirank == 0) {
    /* header_size > 0 must imply header_data != NULL */
    P4EST_ASSERT (header_size <= 0 || header_data != NULL);

    /* write p4est-defined header */
    snprintf (metadata, P4EST_NUM_METADATA_BYTES + 1,
              "%.7s\n%-23s\n%.15ld\n%.15ld\n", P4EST_MAGIC_NUMBER,
              p4est_version (), p4est->global_num_quadrants, header_size);
    mpiret =
      sc_mpi_file_write_at (file_context->file, 0, metadata,
                            P4EST_NUM_METADATA_BYTES, sc_MPI_BYTE);
    P4EST_FILE_CHECK_MPI (mpiret, "Writing the metadata");

    if (header_size != 0) {
      P4EST_ASSERT (header_data != NULL);
      /* Write the user-defined header */
      /* non-collective and blocking */
      mpiret =
        sc_mpi_file_write_at (file_context->file, P4EST_NUM_METADATA_BYTES,
                              header_data, header_size, sc_MPI_BYTE);
      P4EST_FILE_CHECK_MPI (mpiret, "Writing the header");

      /* Write padding bytes for the user-defined header */
      get_padding_string (header_size, P4EST_BYTE_DIV, pad, &num_pad_bytes);
      mpiret =
        sc_mpi_file_write_at (file_context->file,
                              P4EST_NUM_METADATA_BYTES + header_size, pad,
                              num_pad_bytes, sc_MPI_BYTE);
      P4EST_FILE_CHECK_MPI (mpiret, "Writing padding bytes for header");
    }
    else {
      /* There is no header padding */
      num_pad_bytes = 0;
    }
  }
  else {
    get_padding_string (header_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
  }

  P4EST_HANDLE_MPI_ERROR (mpiret, file_context, p4est->mpicomm);

  file_context->p4est = p4est;
  file_context->header_size = header_size + num_pad_bytes;
  file_context->accessed_bytes = 0;
  file_context->num_calls = 0;

  return file_context;
}

p4est_file_context_t *
p4est_file_open_append (p4est_t * p4est, const char *filename,
                        size_t header_size)
{
  int                 mpiret;
  size_t              num_pad_bytes;
  p4est_file_context_t *file_context = P4EST_ALLOC (p4est_file_context_t, 1);
#ifdef P4EST_ENABLE_MPIIO
  sc_MPI_Offset       file_size;
#endif

  /* We do not need the mpi append mode for MPI IO since we use our own byte counter */
  mpiret = sc_mpi_file_open (p4est->mpicomm, filename,
                             sc_MPI_MODE_WRONLY_APPEND, sc_MPI_INFO_NULL,
                             &file_context->file);
  P4EST_FILE_CHECK_OPEN (mpiret, file_context, "File open append");

  get_padding_string (header_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
  file_context->p4est = p4est;
  file_context->header_size = header_size + num_pad_bytes;
  file_context->num_calls = 0;

  /* We can not caculate the file size collectively since this result in
   * different results for different ranks due to the situation that some
   * ranks will have already written some bytes to the file during other ranks
   * calculate the number of accessed bytes.
   */
  if (file_context->p4est->mpirank == 0) {
#ifdef P4EST_ENABLE_MPIIO
    mpiret = MPI_File_get_size (file_context->file, &file_size);
    P4EST_FILE_CHECK_MPI (mpiret, "Get file size for open append");
#else
    /* There is no C-Standard functionality to get the file size */
#endif
  }
#ifdef P4EST_ENABLE_MPIIO
  P4EST_HANDLE_MPI_ERROR (mpiret, file_context, file_context->p4est->mpicomm);
  sc_MPI_Bcast (&file_size, sizeof (sc_MPI_Offset), sc_MPI_BYTE, 0,
                p4est->mpicomm);

  /* Calculate already written quadrant data array bytes */
  file_context->accessed_bytes =
    file_size - P4EST_NUM_METADATA_BYTES - file_context->header_size;
#endif

  return file_context;
}

p4est_file_context_t *
p4est_file_open_read (p4est_t * p4est, const char *filename,
                      size_t header_size, void *header_data)
{
  int                 mpiret, mpiret_sec;
  int                 error_flag;
  size_t              num_pad_bytes;
  char                metadata[P4EST_NUM_METADATA_BYTES + 1];
#ifdef P4EST_ENABLE_MPIIO
  sc_MPI_Offset       file_size;
#endif
  p4est_file_context_t *file_context = P4EST_ALLOC (p4est_file_context_t, 1);

  /* Open the file in the reading mode */
  mpiret =
    sc_mpi_file_open (p4est->mpicomm, filename, sc_MPI_MODE_RDONLY,
                      sc_MPI_INFO_NULL, &file_context->file);
  P4EST_FILE_CHECK_OPEN (mpiret, file_context, "File open read");

  get_padding_string (header_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
  file_context->p4est = p4est;
  file_context->header_size = header_size + num_pad_bytes;
  file_context->accessed_bytes = 0;
  file_context->num_calls = 0;

  /* read metadata and deallocate in case of error */
  if (file_context->p4est->mpirank == 0) {
    error_flag = 0;

#ifdef P4EST_ENABLE_MPIIO
    /* check size of the file */
    mpiret = MPI_File_get_size (file_context->file, &file_size);
    P4EST_FILE_CHECK_MPI (mpiret, "Get file size for open read");
    if (header_size > (size_t) file_size) {
      P4EST_GLOBAL_LERRORF (P4EST_STRING
                            "_io: Error reading <%s>. Header_size is bigger than the file size.\n",
                            filename);
      error_flag = 1;
    }
#else
    /* There is no C-standard functionality to get the file size */
#endif

    /* read metadata on rank 0 */
    mpiret_sec =
      sc_mpi_file_read_at (file_context->file, 0, metadata,
                           P4EST_NUM_METADATA_BYTES, sc_MPI_BYTE);
    /* combine error flags to save one broadcast */
    error_flag |= mpiret | mpiret_sec;
    P4EST_FILE_CHECK_MPI (mpiret_sec, "Reading metadata");

    metadata[P4EST_NUM_METADATA_BYTES] = '\0';
    /* parse metadata; we do not use file_info because we do not want a Bcast */
    error_flag |=
      check_file_metadata (p4est, header_size, filename, metadata);

    if (!error_flag) {
      /* read header on rank 0 and skip the metadata */
      mpiret =
        sc_mpi_file_read_at (file_context->file, P4EST_NUM_METADATA_BYTES,
                             header_data, header_size, sc_MPI_BYTE);
      /* combine error flags to save one broadcast */
      error_flag |= mpiret;
      P4EST_FILE_CHECK_MPI (mpiret, "Reading header");
    }
  }
  /* error checking */
  P4EST_HANDLE_MPI_ERROR (error_flag, file_context, p4est->mpicomm);

  /* broadcast header to all ranks */
  sc_MPI_Bcast (header_data, header_size, sc_MPI_BYTE, 0, p4est->mpicomm);

  return file_context;
}

p4est_file_context_t *
p4est_file_write (p4est_file_context_t * fc, sc_array_t * quadrant_data)
{
  size_t              bytes_to_write, num_pad_bytes, array_size;
  char                array_metadata[P4EST_NUM_ARRAY_METADATA_BYTES + 1],
    pad[P4EST_BYTE_DIV];
  sc_MPI_Offset       write_offset;
#ifdef P4EST_ENABLE_MPI
  int                 mpiret;
#endif

  P4EST_ASSERT (quadrant_data != NULL
                && quadrant_data->elem_count ==
                (size_t) fc->p4est->local_num_quadrants);

  if (quadrant_data->elem_size == 0) {
    /* nothing to write */
    return NULL;
  }

  /* Check how many bytes we write to the disk */
  bytes_to_write = quadrant_data->elem_count * quadrant_data->elem_size;

  /* rank-dependent byte offset */
  write_offset = P4EST_NUM_METADATA_BYTES + fc->header_size +
    fc->p4est->global_first_quadrant[fc->p4est->mpirank] *
    quadrant_data->elem_size;

#ifdef P4EST_ENABLE_MPIIO
  /* set the file size */
  mpiret = MPI_File_set_size (fc->file,
                              P4EST_NUM_METADATA_BYTES + fc->header_size +
                              fc->p4est->global_num_quadrants *
                              quadrant_data->elem_size +
                              P4EST_NUM_ARRAY_METADATA_BYTES +
                              fc->accessed_bytes);
  P4EST_FILE_CHECK_NULL (mpiret, "Set file size");
#else
  /* We do not perform this optimization without MPI I/O */
#endif

  num_pad_bytes = 0;
  if (fc->p4est->mpirank == 0) {
    /* array-dependent metadata */
    snprintf (array_metadata, P4EST_NUM_ARRAY_METADATA_BYTES + 1,
              "\n%.14ld\n", quadrant_data->elem_size);

    /* write array-dependent metadata */
    mpiret =
      sc_mpi_file_write_at (fc->file, fc->accessed_bytes + write_offset,
                            array_metadata, P4EST_NUM_ARRAY_METADATA_BYTES,
                            sc_MPI_BYTE);
    P4EST_FILE_CHECK_MPI (mpiret, "Writing array metadata");

#if 0
    /* The sequence of of the fwrite calls plays a role and fseek has no
     * effect in the append mode.
     */

    SC_CHECK_FOPEN_NULL (fc->file, fopen (fc->filename, "ab"));
    /* write array metadata */
    sc_fwrite (array_metadata, 1, P4EST_NUM_ARRAY_METADATA_BYTES, fc->file,
               "Writing array metadata");
    fflush (fc->file);

    /* write array data */
    sc_fwrite (quadrant_data->array, 1, bytes_to_write, fc->file,
               "Writing array data");
    fflush (fc->file);

    /* write padding bytes */
    array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
    get_padding_string (array_size, P4EST_BYTE_DIV, pad, &num_pad_bytes);
    sc_fwrite (pad, 1, num_pad_bytes, fc->file, "Writitng padding bytes");
    fflush (fc->file);
#endif
  }
  P4EST_HANDLE_MPI_ERROR (mpiret, fc, fc->p4est->mpicomm);

  /* write array data */
  mpiret =
    sc_mpi_file_write_at_all (&fc->file,
                              fc->accessed_bytes + write_offset +
                              P4EST_NUM_ARRAY_METADATA_BYTES,
                              quadrant_data->array, bytes_to_write,
                              sc_MPI_BYTE);
  P4EST_FILE_CHECK_NULL (mpiret, "Writing quadrant-wise");

  /** We place the padding bytes write here because for the sequential
   * IO operations the order of fwrite calls plays a role.
   */
  /* write padding bytes */
  if (fc->p4est->mpirank == 0) {
    /* Caculate and write padding bytes for array data */
    array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
    get_padding_string (array_size, P4EST_BYTE_DIV, pad, &num_pad_bytes);

    mpiret =
      sc_mpi_file_write_at (fc->file,
                            fc->accessed_bytes + P4EST_NUM_METADATA_BYTES +
                            fc->header_size + array_size +
                            P4EST_NUM_ARRAY_METADATA_BYTES, pad,
                            num_pad_bytes, sc_MPI_BYTE);
    P4EST_FILE_CHECK_MPI (mpiret, "Writing padding bytes for a data array");
  }
  else {
    array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
    get_padding_string (array_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
  }

  /* This is *not* the processor local value */
  fc->accessed_bytes +=
    quadrant_data->elem_size * fc->p4est->global_num_quadrants +
    P4EST_NUM_ARRAY_METADATA_BYTES + num_pad_bytes;
  ++fc->num_calls;

  return fc;
}

#ifdef P4EST_ENABLE_MPIIO
p4est_file_context_t *
p4est_file_read (p4est_file_context_t * fc, sc_array_t * quadrant_data)
{
#if !defined (P4EST_ENABLE_MPIIO) && defined (P4EST_ENABLE_MPI)
  sc_MPI_Status       status;
  int                 count;
  int                 active = (fc->p4est->mpirank == 0) ? 1 : 0;
#endif
  int                 error_flag;
  size_t              bytes_to_read, num_pad_bytes, array_size,
    read_data_size;
  char                array_metadata[P4EST_NUM_ARRAY_METADATA_BYTES + 1];
#ifdef P4EST_ENABLE_MPIIO
  sc_MPI_Offset       size;
#else
  int                 no_data_flag = 0;
#endif
#ifdef P4EST_ENABLE_MPI
  int                 mpiret;
#endif

  P4EST_ASSERT (fc != NULL);

  if (quadrant_data == NULL || quadrant_data->elem_size == 0) {
    /* Nothing to read but we shift our own file pointer */
#ifdef P4EST_ENABLE_MPIIO
    mpiret = sc_mpi_file_read_at (fc->file,
                                  fc->accessed_bytes +
                                  P4EST_NUM_METADATA_BYTES + fc->header_size,
                                  array_metadata,
                                  P4EST_NUM_ARRAY_METADATA_BYTES,
                                  sc_MPI_BYTE);
    /* In the case of error the return value is still NULL */
    P4EST_FILE_CHECK_NULL (mpiret, "Reading quadrant-wise metadata");
    array_metadata[P4EST_NUM_ARRAY_METADATA_BYTES] = '\0';
    read_data_size = sc_atol (array_metadata);
#else
    if (fc->p4est->mpirank == 0) {
      fseek (fc->file,
             fc->accessed_bytes + P4EST_NUM_METADATA_BYTES + fc->header_size,
             SEEK_SET);
      if (getc (fc->file) != EOF) {
        sc_fread (array_metadata, 1, P4EST_NUM_ARRAY_METADATA_BYTES, fc->file,
                  "Reading array metadata");
      }
      else {
        /* There is no data at this position  */
        P4EST_LERRORF (P4EST_STRING
                       "_io: The end of the file %s was reached.\n",
                       fc->filename);
        no_data_flag = 1;
        goto no_data;
      }

      array_metadata[P4EST_NUM_ARRAY_METADATA_BYTES] = '\0';
      read_data_size = sc_atol (array_metadata);
    }
    sc_MPI_Bcast (&read_data_size, sizeof (size_t), sc_MPI_BYTE, 0,
                  fc->p4est->mpicomm);
#endif

    /* calculate the padding bytes for this data array */
    array_size = fc->p4est->global_num_quadrants * read_data_size;
    get_padding_string (array_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
    fc->accessed_bytes +=
      read_data_size * fc->p4est->global_num_quadrants +
      P4EST_NUM_ARRAY_METADATA_BYTES + num_pad_bytes;
    ++fc->num_calls;
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
  P4EST_FILE_CHECK_NULL (mpiret, "Get file size for read");
  if (size - P4EST_NUM_METADATA_BYTES - fc->header_size < bytes_to_read) {
    /* report wrong file size, collectively close the file and deallocate fc */
    if (fc->p4est->mpirank == 0) {
      P4EST_LERROR (P4EST_STRING
                    "_io: Error reading. File has less bytes than the user wants to read.\n");
    }
    p4est_file_close (fc);
    return NULL;
  }
#else
  /* There is no C-standard functionality to get the file size */
#endif

  /* check the array metadata */
  error_flag = 0;
  if (fc->p4est->mpirank == 0) {
#ifdef P4EST_ENABLE_MPIIO
    mpiret = sc_mpi_file_read_at (fc->file,
                                  fc->accessed_bytes +
                                  P4EST_NUM_METADATA_BYTES + fc->header_size,
                                  array_metadata,
                                  P4EST_NUM_ARRAY_METADATA_BYTES,
                                  sc_MPI_BYTE);
    P4EST_FILE_CHECK_MPI (mpiret, "Reading quadrant-wise metadata");
#else
    fseek (fc->file,
           fc->accessed_bytes + P4EST_NUM_METADATA_BYTES + fc->header_size,
           SEEK_SET);
    if (getc (fc->file) != EOF) {
      sc_fread (array_metadata, 1, P4EST_NUM_ARRAY_METADATA_BYTES, fc->file,
                "Reading array metadata");
    }
    else {
      /* There is no data after the header in the file */
      P4EST_LERRORF (P4EST_STRING
                     "_io: The end of the file %s was reached.\n",
                     fc->filename);
      no_data_flag = 1;
      goto no_data;
    }
#endif
    array_metadata[P4EST_NUM_ARRAY_METADATA_BYTES] = '\0';
    read_data_size = sc_atol (array_metadata);
    if (read_data_size != quadrant_data->elem_size) {
      P4EST_LERRORF (P4EST_STRING
                     "_io: Error reading. Wrong array element size (in file = %ld, by parameter = %ld).\n",
                     read_data_size, quadrant_data->elem_size);
      error_flag = 1;
    }
  }
#ifdef P4EST_ENABLE_MPIIO
  P4EST_HANDLE_MPI_ERROR (mpiret, fc, fc->p4est->mpicomm);
#endif
  /* broadcast the error flag to decicde if we continue */
  sc_MPI_Bcast (&error_flag, 1, sc_MPI_INT, 0, fc->p4est->mpicomm);
  if (error_flag) {
    p4est_file_close (fc);
    return NULL;
  }

  /* calculate the padding bytes for this data array */
  array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
  get_padding_string (array_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);

#ifdef P4EST_ENABLE_MPIIO
  mpiret = sc_mpi_file_read_at_all (fc->file,
                                    fc->accessed_bytes +
                                    P4EST_NUM_METADATA_BYTES +
                                    P4EST_NUM_ARRAY_METADATA_BYTES +
                                    fc->header_size +
                                    fc->p4est->
                                    global_first_quadrant[fc->p4est->mpirank]
                                    * quadrant_data->elem_size,
                                    quadrant_data->array, bytes_to_read,
                                    sc_MPI_BYTE);
  P4EST_FILE_CHECK_NULL (mpiret, "Reading quadrant-wise");
#elif defined (P4EST_ENABLE_MPI)
  if (fc->p4est->mpirank != 0) {
    /* wait until the preceding process finished the I/O operation */
    /* receive */
    mpiret = sc_MPI_Recv (&active, 1, sc_MPI_INT,
                          fc->p4est->mpirank - 1, sc_MPI_ANY_TAG,
                          fc->p4est->mpicomm, &status);
    SC_CHECK_MPI (mpiret);
    mpiret = MPI_Get_count (&status, sc_MPI_INT, &count);
    SC_CHECK_MPI (mpiret);
    SC_CHECK_ABORT (count == 1, "MPI receive");
  }

  if (active) {
    /* process 0 must not wait */
    if (fc->p4est->mpirank != 0) {
      /* open */
      fc->file =
        sc_fopen (fc->filename, "rb", "Open for sequential MPI read");
    }

    /* read array data */
    fseek (fc->file, fc->accessed_bytes + P4EST_NUM_METADATA_BYTES +
           P4EST_NUM_ARRAY_METADATA_BYTES + fc->header_size +
           fc->p4est->global_first_quadrant[fc->p4est->mpirank] *
           quadrant_data->elem_size, SEEK_SET);
    sc_fread (quadrant_data->array, 1, bytes_to_read, fc->file,
              "Reading data array");

    /* close */
    fclose (fc->file);
    /* only update active process if there are processes left */
    if (fc->p4est->mpirank < fc->p4est->mpisize - 1) {
      /* the current process finished its I/O operations */
      P4EST_ASSERT (active == 1);
      /* send */
      mpiret = sc_MPI_Send (&active, 1, sc_MPI_INT,
                            fc->p4est->mpirank + 1, 1, fc->p4est->mpicomm);
      SC_CHECK_MPI (mpiret);
    }
  }
#else
  fseek (fc->file, fc->accessed_bytes + P4EST_NUM_METADATA_BYTES +
         P4EST_NUM_ARRAY_METADATA_BYTES + fc->header_size +
         fc->p4est->global_first_quadrant[fc->p4est->mpirank] *
         quadrant_data->elem_size, SEEK_SET);
  sc_fread (quadrant_data->array, 1, bytes_to_read, fc->file,
            "Reading data array");
#endif

  fc->accessed_bytes +=
    quadrant_data->elem_size * fc->p4est->global_num_quadrants +
    P4EST_NUM_ARRAY_METADATA_BYTES + num_pad_bytes;
  ++fc->num_calls;
#if !defined (P4EST_ENABLE_MPIIO) && defined (P4EST_ENABLE_MPI)
/* The processes have to wait here because they are not allowed to start
     * other I/O operations.
     */
  sc_MPI_Barrier (fc->p4est->mpicomm);
  /* reset the processor activity values */
  if (fc->p4est->mpirank == 0) {
    fc->file =
      sc_fopen (fc->filename, "rb", "Open after reading of one chunk");
  }
  else {
    fc->file = NULL;
  }
#endif

#ifndef P4EST_ENABLE_MPIIO
  /* Without MPI IO we can not check the file size and therefore we must
   * read the header and then jump out of reading process if we detect EOF
   */
no_data:
  if (no_data_flag) {
    return NULL;
  }
#endif

  return fc;
}
#endif

static int
p4est_file_error_cleanup (sc_MPI_File * file)
{
  /* no error checking since we are called under an error condition */
  P4EST_ASSERT (file != NULL);
#ifdef P4EST_ENABLE_MPIIO
  if (*file != sc_MPI_FILE_NULL) {
#else
  /** The MPI IO file object is a pointer itself
   * but for the sake of simplicity of the IO
   * functions in libsc we use a struct as file
   * object.
   */
  if (file != sc_MPI_FILE_NULL) {
#endif
#ifdef P4EST_ENABLE_MPIIO
    MPI_File_close (file);
#else
    /* TODO: Use here a libsc closing function */
    fclose (file->file);
    //*file = NULL;
#endif
  }
  return sc_MPI_ERR_IO;
}

#ifdef P4EST_ENABLE_MPIIO
int
p4est_file_info (p4est_t * p4est, const char *filename,
                 size_t * header_size, sc_array_t * elem_size)
{
  int                 mpiret, eclass;
#ifdef P4EST_ENABLE_MPIIO
  int                 retval, icount;
  sc_MPI_Status       mpistatus;
#endif
  long                long_header;
  size_t              current_member, num_pad_bytes, padded_header;
  char                metadata[P4EST_NUM_METADATA_BYTES + 1];
  char                array_metadata[P4EST_NUM_ARRAY_METADATA_BYTES + 1];
  char               *parsing_arg;
  sc_MPI_Offset       current_position;
  sc_MPI_File         file;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (filename != NULL);
  P4EST_ASSERT (header_size != NULL);
  P4EST_ASSERT (elem_size != NULL);
  P4EST_ASSERT (elem_size->elem_size == sizeof (size_t));

  /* set default output values */
  *header_size = 0;
  sc_array_reset (elem_size);

  /* open the file in reading mode */
  eclass = sc_MPI_SUCCESS;      /* MPI defines MPI_SUCCESS to equal 0. */
  file = sc_MPI_FILE_NULL;
#ifdef P4EST_ENABLE_MPIIO
  if ((retval =
       sc_MPI_File_open (p4est->mpicomm, filename, sc_MPI_MODE_RDONLY,
                         sc_MPI_INFO_NULL, &file)) != sc_MPI_SUCCESS) {
    mpiret = sc_mpi_file_error_class (retval, &eclass);
    SC_CHECK_MPI (mpiret);
  }
#else
  if (p4est->mpirank == 0) {
    if ((file = fopen (filename, "rb")) == NULL) {
      mpiret = sc_mpi_file_error_class (errno, &eclass);
      SC_CHECK_MPI (mpiret);
    }
  }
  mpiret = sc_MPI_Bcast (&eclass, 1, sc_MPI_INT, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
#endif
  if (eclass) {
    return eclass;
  }

  /* read file metadata on root rank */
  P4EST_ASSERT (!eclass);
  if (p4est->mpirank == 0) {
#ifdef P4EST_ENABLE_MPIIO
    if ((retval = sc_mpi_file_read_at (file, 0, metadata,
                                       P4EST_NUM_METADATA_BYTES, sc_MPI_BYTE))
        != sc_MPI_SUCCESS) {
      mpiret = sc_mpi_file_error_class (retval, &eclass);
      SC_CHECK_MPI (mpiret);
    }
#else
    if (fread (metadata, P4EST_NUM_METADATA_BYTES, 1, file) != 1) {
      eclass = sc_MPI_ERR_IO;
    }
#endif
  }
  mpiret = sc_MPI_Bcast (&eclass, 1, sc_MPI_INT, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  if (eclass) {
    return eclass;
  }

  /* broadcast file metadata to all ranks and null-terminate it */
  mpiret = sc_MPI_Bcast (metadata, P4EST_NUM_METADATA_BYTES, sc_MPI_BYTE, 0,
                         p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  metadata[P4EST_NUM_METADATA_BYTES] = '\0';

  /* split the input string for the magic and skip version string */
  parsing_arg = strtok (metadata, "\n");
  if (parsing_arg == NULL || strcmp (parsing_arg, P4EST_MAGIC_NUMBER)) {
    P4EST_GLOBAL_LERROR ("p4est_file_info: magic string mismatch\n");
    return p4est_file_error_cleanup (&file);
  }
  parsing_arg = strtok (NULL, "\n");
  if (parsing_arg == NULL || strlen (parsing_arg) != 23) {
    P4EST_GLOBAL_LERROR ("p4est_file_info: version string length mismatch\n");
    return p4est_file_error_cleanup (&file);
  }

  /* split the metadata: global number of quadrants */
  parsing_arg = strtok (NULL, "\n");
  if (parsing_arg == NULL || strlen (parsing_arg) != 15) {
    P4EST_GLOBAL_LERROR ("p4est_file_info: global count length mismatch\n");
    return p4est_file_error_cleanup (&file);
  }
  if (((p4est_gloidx_t) sc_atol (parsing_arg)) != p4est->global_num_quadrants) {
    P4EST_GLOBAL_LERROR ("p4est_file_info: global quadrant count mismatch\n");
    return p4est_file_error_cleanup (&file);
  }

  /* split the metadata: number of header bytes */
  parsing_arg = strtok (NULL, "\n");
  if (parsing_arg == NULL || strlen (parsing_arg) != 15) {
    P4EST_GLOBAL_LERROR ("p4est_file_info: header string length mismatch\n");
    return p4est_file_error_cleanup (&file);
  }
  if ((long_header = sc_atol (parsing_arg)) < 0) {
    P4EST_GLOBAL_LERROR ("p4est_file_info: header length negative\n");
    return p4est_file_error_cleanup (&file);
  }
  *header_size = (size_t) long_header;

  /* calculate the padding bytes for the user-defined header */
  get_padding_string (*header_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
  padded_header = *header_size + num_pad_bytes;
  current_position =
    (sc_MPI_Offset) (P4EST_NUM_METADATA_BYTES + padded_header);

  /* read all data headers that we find and skip the data itself */
  if (p4est->mpirank == 0) {
    for (;;) {
      /* read array metadata for current record */
#ifdef P4EST_ENABLE_MPIIO
      mpiret = MPI_File_read_at (file, current_position, array_metadata,
                                 P4EST_NUM_ARRAY_METADATA_BYTES, sc_MPI_BYTE,
                                 &mpistatus);
      P4EST_FILE_CHECK_INT (mpiret, "MPI_File_read_at on proc 0");
      mpiret = sc_MPI_Get_count (&mpistatus, sc_MPI_BYTE, &icount);
      SC_CHECK_MPI (mpiret);
      if (icount != P4EST_NUM_ARRAY_METADATA_BYTES) {
        break;
      }
#else
      if (fseek (file, current_position, SEEK_SET) ||
          fread (array_metadata, P4EST_NUM_ARRAY_METADATA_BYTES, 1,
                 file) != 1) {
        break;
      }
#endif
      array_metadata[P4EST_NUM_ARRAY_METADATA_BYTES] = '\0';

      /* parse and store the element size of the array */
      parsing_arg = strtok (array_metadata, "\n");
      if (parsing_arg == NULL
          || strlen (parsing_arg) != P4EST_NUM_ARRAY_METADATA_CHARS) {
        break;
      }
      if ((long_header = sc_atol (parsing_arg)) < 0) {
        break;
      }
      *(long *) sc_array_push (elem_size) = long_header;

      /* get padding bytes of the current array */
      current_member = (size_t) (p4est->global_num_quadrants * long_header);
      get_padding_string (current_member, P4EST_BYTE_DIV, NULL,
                          &num_pad_bytes);
      current_position +=
        P4EST_NUM_ARRAY_METADATA_BYTES + current_member + num_pad_bytes;
    }
  }

  /* replicate array metadata in parallel */
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
    mpiret = sc_mpi_file_error_class (retval, &eclass);
    SC_CHECK_MPI (mpiret);
  }
#else
  if (p4est->mpirank == 0) {
    if (fclose (file)) {
      mpiret = sc_mpi_file_error_class (errno, &eclass);
      SC_CHECK_MPI (mpiret);
    }
  }
  else {
    P4EST_ASSERT (file == NULL);
  }
  mpiret = sc_MPI_Bcast (&eclass, 1, sc_MPI_INT, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
#endif
  return eclass;
}
#endif

void
p4est_file_close (p4est_file_context_t * fc)
{
  P4EST_ASSERT (fc != NULL);
#ifdef P4EST_ENABLE_MPIIO
  /* TODO: consider cases in libsc and not here */
  sc_MPI_File_close (&fc->file);
#else
  if (fc->file.file != NULL) {
    fclose (fc->file.file);
  }
#endif
  P4EST_FREE (fc);
}
