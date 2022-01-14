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
#if !defined (P4EST_ENABLE_MPII0) && !defined (P4EST_ENABLE_MPI)
/* We report errors without MPI I/O and MPI but without MPI we abort */
#include <errno.h>
#endif

#if !defined (P4EST_ENABLE_MPII0) && !defined (P4EST_ENABLE_MPI)
extern int          errno;
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

#ifndef P4_TO_P8

struct p4est_file_context
{
  p4est_t            *p4est;
  size_t              header_size;      /* only the user-defined header */
  size_t              num_calls;        /* redundant but for convience;
                                           counts the number of calls of
                                           write and read repectivly */
#ifndef P4EST_ENABLE_MPIIO
  const char         *filename; /* we need to store the path for
                                   the successive opening strategy
                                   and to use the append mode to write
                                   data without MPI */
  FILE               *file;
#else
  sc_MPI_File         file;
  sc_MPI_Offset       accessed_bytes;   /* count only array data bytes and
                                           array metadata bytes */
#endif
};

#if 0
/* *INDENT-OFF* */
static int
 p4est_file_info_extra (p4est_file_context_t * fc,
                        p4est_gloidx_t * global_num_quads,
                        char p4est_version[24],
                        char magic_num[8],size_t *header_size,
                        sc_array_t * elem_size, long max_num_arrays);
/* *INDENT-ON* */
#endif

/** This function calculates a padding string consisting of spaces.
 * We require an already allocated array pad or NULL.
 * For NULL the function only calculates the number of padding bytes.
 */
static void
get_padding_string (size_t num_bytes, size_t divisor, char *pad,
                    size_t *num_pad_bytes)
{
  P4EST_ASSERT (divisor != 0 && num_pad_bytes != NULL);

  *num_pad_bytes = (divisor - (num_bytes % divisor)) % divisor;

  if (*num_pad_bytes > 0 && pad != NULL) {
    snprintf (pad, *num_pad_bytes + 1, "%-*s", *((int *) num_pad_bytes), "");
  }
}

#ifdef P4EST_MPI_WRITE_NOT_READY

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
  P4EST_ASSERT (parsing_arg != NULL);
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
#ifdef P4EST_ENABLE_MPIIO
  int                 mpiret;
#endif
  char                metadata[P4EST_NUM_METADATA_BYTES + 1];
  char                pad[P4EST_BYTE_DIV];
  size_t              num_pad_bytes;
  p4est_file_context_t *file_context = P4EST_ALLOC (p4est_file_context_t, 1);

  /* Open the file and create a new file if necessary */
#ifdef P4EST_ENABLE_MPIIO
  mpiret =
    sc_mpi_open (p4est->mpicomm, filename,
                 sc_MPI_MODE_WRONLY | sc_MPI_MODE_CREATE, sc_MPI_INFO_NULL,
                 &file_context->file, "File open create");
  P4EST_FILE_CHECK_OPEN (mpiret, file_context, "File open create");
#elif defined (P4EST_ENABLE_MPI)
  /* serialize the I/O operations */
  /* active flag is set later in \ref p4est_file_write */
  file_context->filename = filename;
  if (p4est->mpirank == 0) {
    file_context->file =
      sc_fopen (file_context->filename, "wb", "File open create");
  }
  else {
    file_context->file = NULL;
  }
#else
  /* no MPI */
  file_context->filename = filename;
  SC_CHECK_FOPEN_NULL (file_context->file, fopen (filename, "wb"));
#endif

  if (p4est->mpirank == 0) {
    /* header_size > 0 must imply header_data != NULL */
    P4EST_ASSERT (header_size <= 0 || header_data != NULL);

    /* write p4est-defined header */
    snprintf (metadata, P4EST_NUM_METADATA_BYTES + 1,
              "%.7s\n%-23s\n%.15ld\n%.15ld\n", P4EST_MAGIC_NUMBER,
              p4est_version (), p4est->global_num_quadrants, header_size);
#ifdef P4EST_ENABLE_MPIIO
    mpiret =
      sc_mpi_write (file_context->file, metadata, P4EST_NUM_METADATA_BYTES,
                    sc_MPI_BYTE, "Writing the metadata");
    P4EST_FILE_CHECK_MPI (mpiret, "Writing the metadata");
#else
    /* this works with and without MPI */
    sc_fwrite (metadata, 1, P4EST_NUM_METADATA_BYTES, file_context->file,
               "Writing the metadata");
    fflush (file_context->file);
#endif

    if (header_size != 0) {
      P4EST_ASSERT (header_data != NULL);
      /* Write the user-defined header */
      /* non-collective and blocking */
#ifdef P4EST_ENABLE_MPIIO
      mpiret = sc_mpi_write (file_context->file, header_data,
                             header_size, sc_MPI_BYTE, "Writing the header");
      P4EST_FILE_CHECK_MPI (mpiret, "Writing the header");
#else
      /* this works with and without MPI */
      sc_fwrite (header_data, 1, header_size, file_context->file,
                 "Writing the header");
      fflush (file_context->file);
#endif

      /* Write padding bytes for the user-defined header */
      get_padding_string (header_size, P4EST_BYTE_DIV, pad, &num_pad_bytes);
#ifdef P4EST_ENABLE_MPIIO
      mpiret = sc_mpi_write (file_context->file, pad,
                             num_pad_bytes, sc_MPI_BYTE,
                             "Writing padding bytes for header");
      P4EST_FILE_CHECK_MPI (mpiret, "Writing padding bytes for header");
#else
      /* this works with and without MPI */
      sc_fwrite (pad, 1, num_pad_bytes, file_context->file,
                 "Writing padding bytes for header");
      fflush (file_context->file);
      /* file is closed to use the append mode to write data */
      fclose (file_context->file);
#endif
    }
    else {
      /* There is no header padding */
      num_pad_bytes = 0;
    }
  }
  else {
    get_padding_string (header_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
  }

#ifdef P4EST_ENABLE_MPIIO
  P4EST_HANDLE_MPI_ERROR (mpiret, file_context, p4est->mpicomm);
#endif

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
#ifdef P4EST_ENABLE_MPIIO
  int                 mpiret;
#endif
  size_t              num_pad_bytes;
  p4est_file_context_t *file_context = P4EST_ALLOC (p4est_file_context_t, 1);
#ifdef P4EST_ENABLE_MPIIO
  sc_MPI_Offset       file_size;
#endif

#ifdef P4EST_ENABLE_MPIIO
  /* We do not need the mpi append mode since we use our own byte counter */
  mpiret = sc_mpi_open (p4est->mpicomm, filename,
                        sc_MPI_MODE_WRONLY, sc_MPI_INFO_NULL,
                        &file_context->file, "File open append");
  P4EST_FILE_CHECK_OPEN (mpiret, file_context, "File open append");
#elif defined (P4EST_ENABLE_MPI)
  /* the file is opened in rank-order in \ref p4est_file_write */
  file_context->filename = filename;
  file_context->file = NULL;
#else
  /* no MPI */
  file_context->filename = filename;
  file_context->file = NULL;
#endif

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
    mpiret = sc_mpi_get_file_size (file_context->file, &file_size,
                                   "Get file size for open append");
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
#ifdef P4EST_ENABLE_MPIIO
  int                 mpiret, mpiret_sec;
#endif
  int                 error_flag;
  size_t              num_pad_bytes;
  char                metadata[P4EST_NUM_METADATA_BYTES + 1];
#ifdef P4EST_ENABLE_MPIIO
  sc_MPI_Offset       file_size;
#endif
  p4est_file_context_t *file_context = P4EST_ALLOC (p4est_file_context_t, 1);

#ifdef P4EST_ENABLE_MPIIO
  /* Open the file in the reading mode */
  mpiret =
    sc_mpi_open (p4est->mpicomm, filename, sc_MPI_MODE_RDONLY,
                 sc_MPI_INFO_NULL, &file_context->file, "File open read");
  P4EST_FILE_CHECK_OPEN (mpiret, file_context, "File open read");
#elif defined (P4EST_ENABLE_MPI)
  file_context->filename = filename;
  if (p4est->mpirank == 0) {
    file_context->file =
      sc_fopen (file_context->filename, "rb", "File open read");
  }
  else {
    file_context->file = NULL;
  }
#else
  SC_CHECK_FOPEN_NULL (file_context->file, fopen (filename, "rb"));
  file_context->filename = filename;
#endif

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
    mpiret = sc_mpi_get_file_size (file_context->file, &file_size,
                                   "Get file size for open read");
    P4EST_FILE_CHECK_MPI (mpiret, "Get file size for open read");
    if (header_size > (size_t) file_size) {
      P4EST_LERRORF (P4EST_STRING
                     "_io: Error reading <%s>. Header_size is bigger than the file size.\n",
                     filename);
      error_flag = 1;
    }
#else
    /* There is no C-standard functionality to get the file size */
#endif

#ifdef P4EST_ENABLE_MPIIO
    /* read metadata on rank 0 */
    mpiret_sec =
      sc_mpi_read_at (file_context->file, 0, metadata,
                      P4EST_NUM_METADATA_BYTES, sc_MPI_BYTE,
                      "Reading metadata");
    /* combine error flags to save one broadcast */
    error_flag |= mpiret | mpiret_sec;
    P4EST_FILE_CHECK_MPI (mpiret_sec, "Reading metadata");
#else
    /* file pointer is already set to 0 above */
    sc_fread (metadata, 1, P4EST_NUM_METADATA_BYTES, file_context->file,
              "Reading file metadata");
#endif
    metadata[P4EST_NUM_METADATA_BYTES] = '\0';
    /* parse metadata; we do not use file_info because we do not want a Bcast */
    error_flag |=
      check_file_metadata (p4est, header_size, filename, metadata);

    if (!error_flag) {
#ifdef P4EST_ENABLE_MPIIO
      /* read header on rank 0 and skip the metadata */
      mpiret =
        sc_mpi_read_at (file_context->file, P4EST_NUM_METADATA_BYTES,
                        header_data, header_size, sc_MPI_BYTE,
                        "Reading header");
      /* combine error flags to save one broadcast */
      error_flag |= mpiret;
      P4EST_FILE_CHECK_MPI (mpiret, "Reading header");
#else
      fseek (file_context->file, P4EST_NUM_METADATA_BYTES, SEEK_SET);
      sc_fread (header_data, 1, header_size, file_context->file,
                "Reading header");
#endif
    }
  }
  /* error checking */
#ifdef P4EST_ENABLE_MPIIO
  P4EST_HANDLE_MPI_ERROR (error_flag, file_context, p4est->mpicomm);
#endif

  /* broadcast header to all ranks */
  sc_MPI_Bcast (header_data, header_size, sc_MPI_BYTE, 0, p4est->mpicomm);

  return file_context;
}

p4est_file_context_t *
p4est_file_write (p4est_file_context_t * fc, sc_array_t * quadrant_data)
{
#if !defined (P4EST_ENABLE_MPIIO) && defined (P4EST_ENABLE_MPI)
  sc_MPI_Status       status;
  int                 count;
  int                 active = (fc->p4est->mpirank == 0) ? 1 : 0;
#endif
  size_t              bytes_to_write, num_pad_bytes, array_size;
  char                array_metadata[P4EST_NUM_ARRAY_METADATA_BYTES + 1],
    pad[P4EST_BYTE_DIV];
#ifdef P4EST_ENABLE_MPIIO
  sc_MPI_Offset       write_offset;
#endif
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

#ifdef P4EST_ENABLE_MPIIO
  /* rank-dependent byte offset */
  write_offset = P4EST_NUM_METADATA_BYTES + fc->header_size +
    fc->p4est->global_first_quadrant[fc->p4est->mpirank] *
    quadrant_data->elem_size;
#endif

#ifdef P4EST_ENABLE_MPIIO
  /* set the file size */
  mpiret = sc_mpi_set_file_size (fc->file,
                                 P4EST_NUM_METADATA_BYTES + fc->header_size +
                                 fc->p4est->global_num_quadrants *
                                 quadrant_data->elem_size +
                                 P4EST_NUM_ARRAY_METADATA_BYTES +
                                 fc->accessed_bytes, "Set file size");
  P4EST_FILE_CHECK_NULL (mpiret, "Set file size");
#else
  /* We do not perform this optimization without MPI I/O */
#endif

  if (fc->p4est->mpirank == 0) {
    /* array-dependent metadata */
    snprintf (array_metadata, P4EST_NUM_ARRAY_METADATA_BYTES + 1,
              "\n%.14ld\n", quadrant_data->elem_size);
#ifdef P4EST_ENABLE_MPIIO
    mpiret = sc_mpi_write_at (fc->file, fc->accessed_bytes + write_offset,
                              array_metadata, P4EST_NUM_ARRAY_METADATA_BYTES,
                              sc_MPI_BYTE, "Writing array metadata");
    P4EST_FILE_CHECK_MPI (mpiret, "Writing array metadata");
#elif defined (P4EST_ENABLE_MPI)
    fc->file = sc_fopen (fc->filename, "ab", "Open before writing metadata");
    /* write array metadata */
    sc_fwrite (array_metadata, 1, P4EST_NUM_ARRAY_METADATA_BYTES, fc->file,
               "Writing array metadata");
    fflush (fc->file);
#else
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

#ifdef P4EST_ENABLE_MPIIO
    /* Caculate and write padding bytes for array data */
    array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
    get_padding_string (array_size, P4EST_BYTE_DIV, pad, &num_pad_bytes);

    mpiret =
      sc_mpi_write_at (fc->file,
                       fc->accessed_bytes + P4EST_NUM_METADATA_BYTES +
                       fc->header_size + array_size +
                       P4EST_NUM_ARRAY_METADATA_BYTES, pad, num_pad_bytes,
                       sc_MPI_BYTE, "Writing padding bytes for a data array");
    P4EST_FILE_CHECK_MPI (mpiret, "Writing padding bytes for a data array");
#endif
  }
  else {
    array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
    get_padding_string (array_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
  }
#ifdef P4EST_ENABLE_MPIIO
  P4EST_HANDLE_MPI_ERROR (mpiret, fc, fc->p4est->mpicomm);
#endif

#ifdef P4EST_ENABLE_MPIIO
  mpiret =
    sc_mpi_write_at_all (fc->file,
                         fc->accessed_bytes + write_offset +
                         P4EST_NUM_ARRAY_METADATA_BYTES, quadrant_data->array,
                         bytes_to_write, sc_MPI_BYTE,
                         "Writing quadrant-wise");
  P4EST_FILE_CHECK_NULL (mpiret, "Writing quadrant-wise");
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
        sc_fopen (fc->filename, "ab", "Open for seqeuential MPI data write");
    }

    /* write array data */
    sc_fwrite (quadrant_data->array, 1, bytes_to_write, fc->file,
               "Writitng array data");
    fflush (fc->file);

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
  /* The case without MPI was already considered above */
#endif

  /* This is *not* the processor local value */
  fc->accessed_bytes +=
    quadrant_data->elem_size * fc->p4est->global_num_quadrants +
    P4EST_NUM_ARRAY_METADATA_BYTES + num_pad_bytes;
  ++fc->num_calls;
#if !defined (P4EST_ENABLE_MPIIO) && defined (P4EST_ENABLE_MPI)
  /* The processes have to wait here because they are not allowed to start
   * other I/O operations.
   */
  sc_MPI_Barrier (fc->p4est->mpicomm);
  if (fc->p4est->mpirank == 0) {
    fc->file = sc_fopen (fc->filename, "ab", "Open to write padding bytes");
    /* write padding bytes */
    array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
    get_padding_string (array_size, P4EST_BYTE_DIV, pad, &num_pad_bytes);
    sc_fwrite (pad, 1, num_pad_bytes, fc->file, "Writing padding bytes");
    fflush (fc->file);
    /* close */
    fclose (fc->file);
  }

  if (fc->p4est->mpirank == 0) {
    fc->file = sc_fopen (fc->filename, "ab", "Open after write one chunk");
  }
  else {
    fc->file = NULL;
  }
#endif

  return fc;
}

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
    mpiret = sc_mpi_read_at (fc->file,
                             fc->accessed_bytes + P4EST_NUM_METADATA_BYTES +
                             fc->header_size, array_metadata,
                             P4EST_NUM_ARRAY_METADATA_BYTES, sc_MPI_BYTE,
                             "Reading quadrant-wise metadata");
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
  mpiret = sc_mpi_get_file_size (fc->file, &size, "Get file size for read");
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
    mpiret = sc_mpi_read_at (fc->file,
                             fc->accessed_bytes + P4EST_NUM_METADATA_BYTES +
                             fc->header_size, array_metadata,
                             P4EST_NUM_ARRAY_METADATA_BYTES, sc_MPI_BYTE,
                             "Reading quadrant-wise metadata");
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
  sc_MPI_Bcast (&error_flag, sizeof (int), sc_MPI_BYTE, 0,
                fc->p4est->mpicomm);
  if (error_flag) {
    p4est_file_close (fc);
    return NULL;
  }

  /* calculate the padding bytes for this data array */
  array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
  get_padding_string (array_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);

#ifdef P4EST_ENABLE_MPIIO
  mpiret = sc_mpi_read_at_all (fc->file,
                               fc->accessed_bytes + P4EST_NUM_METADATA_BYTES +
                               P4EST_NUM_ARRAY_METADATA_BYTES +
                               fc->header_size +
                               fc->p4est->global_first_quadrant[fc->
                                                                p4est->mpirank]
                               * quadrant_data->elem_size,
                               quadrant_data->array, bytes_to_read,
                               sc_MPI_BYTE, "Reading quadrant-wise");
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

#endif /* P4EST_MPI_WRITE_NOT_READY */

/* always executed on process 0 exclusively */
static int
fill_elem_size (p4est_t * p4est,
#ifdef SC_ENABLE_MPIIO
   sc_MPI_File file,
#else
  FILE *file,
#endif
   size_t header_size, sc_array_t * elem_size, long max_num_arrays)
{
  char                array_metadata[P4EST_NUM_ARRAY_METADATA_BYTES + 1];
  char               *parsing_arg;
  size_t              num_pad_bytes, array_size;
  long               *new_elem;
#ifndef P4EST_ENABLE_MPIIO
  long                current_position;
#else
  int                 mpiret;
  sc_MPI_Offset       size;
  sc_MPI_Offset       current_position;
#endif

  /* read the array metadata */
  P4EST_ASSERT (elem_size->elem_size == sizeof (size_t));
  P4EST_ASSERT (elem_size->elem_count == 0);

#ifdef P4EST_ENABLE_MPIIO
  /* get the file size */
  mpiret = sc_mpi_get_file_size (file, &size, "Get file size for info_extra");
  SC_CHECK_MPI_VERBOSE (mpiret, "Get file size for info_extra");
  if (mpiret != sc_MPI_SUCCESS) {
    return -1;
  }
#else
  /* There is no C-standard functionality to get the file size */
#endif
  /* We read the metadata and we skip the user-defined header */
  current_position = P4EST_NUM_METADATA_BYTES + header_size;
#ifdef P4EST_ENABLE_MPIIO
  while (current_position < size) {
    mpiret = sc_mpi_read_at (file, current_position, array_metadata,
                             P4EST_NUM_ARRAY_METADATA_BYTES, sc_MPI_BYTE);
    SC_CHECK_MPI_VERBOSE (mpiret, "Reading array metadata");
    if (mpiret != sc_MPI_SUCCESS) {
      return mpiret;
    }
#else
  while (!fseek (file, current_position, SEEK_SET)
         && fread (array_metadata, 1, P4EST_NUM_ARRAY_METADATA_BYTES,
                   file) == P4EST_NUM_ARRAY_METADATA_BYTES
         && ((max_num_arrays < 0) ? 1
             : (elem_size->elem_count < (size_t) max_num_arrays))) {
#endif
    array_metadata[P4EST_NUM_ARRAY_METADATA_BYTES] = '\0';
    /* parse and store the element size of the array */
    parsing_arg = strtok (array_metadata, "\n");
    P4EST_ASSERT (parsing_arg != NULL);
    new_elem = (long *) sc_array_push (elem_size);
    *new_elem = sc_atol (parsing_arg);
    P4EST_ASSERT (*new_elem > 0);

    /* get padding bytes of the current array */
    array_size = p4est->global_num_quadrants * *new_elem;
    get_padding_string (array_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
    current_position +=
      array_size + P4EST_NUM_ARRAY_METADATA_BYTES + num_pad_bytes;
  }
  return sc_MPI_SUCCESS;
}

/* unused */
#if 0
static int
p4est_file_info_extra (p4est_file_context_t * fc,
                       p4est_gloidx_t * global_num_quads,
                       char p4est_version[24],
                       char magic_num[8], size_t *header_size,
                       sc_array_t * elem_size, long max_num_arrays)
{
#ifdef P4EST_ENABLE_MPIIO
  int                 mpiret, fillret;
#endif
  int                 read_file_metadata, count;
  char                metadata[P4EST_NUM_METADATA_BYTES + 1];
  char               *parsing_arg;
  size_t              current_member;

  P4EST_ASSERT (fc != NULL);

  read_file_metadata = global_num_quads != NULL || p4est_version != NULL
    || magic_num != NULL || header_size != NULL;
  if (fc->p4est->mpirank == 0) {
    if (read_file_metadata) {
#ifdef P4EST_ENABLE_MPIIO
      /* read metadata on rank 0 */
      mpiret =
        sc_mpi_read_at (fc->file, 0, metadata, P4EST_NUM_METADATA_BYTES,
                        sc_MPI_BYTE, "Reading metadata");
      SC_CHECK_MPI_VERBOSE (mpiret, "Reading metadata");
#else
      sc_fread (metadata, 1, P4EST_NUM_METADATA_BYTES, fc->file,
                "Reading file metadata");
#endif
    }
    metadata[P4EST_NUM_METADATA_BYTES] = '\0';
    if (elem_size != NULL) {
#ifdef P4EST_ENABLE_MPIIO
      fillret =
        fill_elem_size (fc->p4est, fc->file, fc->header_size, elem_size,
                        max_num_arrays);
#else
      fill_elem_size (fc->p4est, fc->file, fc->header_size, elem_size,
                      max_num_arrays);
#endif
    }
  }
#ifdef P4EST_ENABLE_MPIIO
  if (read_file_metadata) {
    if (fc->p4est->mpirank == 0) {
      mpiret |= fillret;
    }
    sc_MPI_Bcast (&mpiret, sizeof (int), sc_MPI_BYTE, 0, fc->p4est->mpicomm);
    if (mpiret != sc_MPI_SUCCESS) {
      return mpiret;
    }
  }
#endif

  if (read_file_metadata) {
    /* broadcast to all ranks */
    /* file metadata */
    sc_MPI_Bcast (metadata, P4EST_NUM_METADATA_BYTES, sc_MPI_BYTE, 0,
                  fc->p4est->mpicomm);
  }

  if (elem_size != NULL) {
    /* array metadata */
    current_member = elem_size->elem_size;
    sc_MPI_Bcast (&current_member, sizeof (size_t), sc_MPI_BYTE, 0,
                  fc->p4est->mpicomm);
    if (fc->p4est->mpirank != 0) {
      sc_array_init (elem_size, current_member);
    }
    current_member = elem_size->elem_count;
    sc_MPI_Bcast (&current_member, sizeof (size_t), sc_MPI_BYTE, 0,
                  fc->p4est->mpicomm);
    if (fc->p4est->mpirank != 0) {
      sc_array_resize (elem_size, current_member);
    }
    sc_MPI_Bcast (elem_size->array,
                  elem_size->elem_count * elem_size->elem_size, sc_MPI_BYTE,
                  0, fc->p4est->mpicomm);
  }

  if (read_file_metadata) {
    /* split the input string */
    count = 0;
    parsing_arg = strtok (metadata, "\n");
    P4EST_ASSERT (parsing_arg != NULL);
    while (parsing_arg != NULL && count < 4) {
      if (magic_num != NULL && count == 0) {
        *magic_num = sc_atoi (parsing_arg);
      }
      else if (p4est_version != NULL && count == 1) {
        strcpy (p4est_version, parsing_arg);
      }
      else if (global_num_quads != NULL && count == 2) {
        *global_num_quads = sc_atol (parsing_arg);
      }
      else if (header_size != NULL && count == 3) {
        *header_size = sc_atol (parsing_arg);
      }
      parsing_arg = strtok (NULL, "\n");
      ++count;
    }
  }
  return sc_MPI_SUCCESS;
}
#endif

int
p4est_file_info (p4est_t * p4est, const char *filename,
                 p4est_gloidx_t * global_num_quadrants,
                 size_t *header_size, sc_array_t * elem_size)
{
  int                 retval;
#ifndef P4EST_ENABLE_MPIIO
  FILE               *file;
#else
  sc_MPI_File         file;
#endif
#ifdef P4EST_ENABLE_MPI
  int                 mpiret;
#endif
  int                 count;
  char                metadata[P4EST_NUM_METADATA_BYTES];
  char               *parsing_arg;
  size_t              current_member, num_pad_bytes, padded_header;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (filename != NULL);
  P4EST_ASSERT (global_num_quadrants != NULL);
  P4EST_ASSERT (header_size != NULL);
  P4EST_ASSERT (elem_size != NULL);
  P4EST_ASSERT (elem_size->elem_size == sizeof (size_t));

  /* set default output values */
  *global_num_quadrants = 0;
  *header_size = 0;
  sc_array_reset (elem_size);

#ifdef P4EST_ENABLE_MPIIO
  /* open the file in reading mode */
  if ((mpiret = sc_MPI_File_open (p4est->mpicomm, filename, sc_MPI_MODE_RDONLY,
                                 sc_MPI_INFO_NULL, &file)) != sc_MPI_SUCCESS) {
    /* TODO: rather do this on root rank only? */
    /* TODO: rather do not print anything and copy string to caller */
    SC_CHECK_MPI_VERBOSE (mpiret, "p4est_file_info open");
    return -1;
  }
#else
  retval = 0;
  file = NULL;
  if (p4est->mpirank == 0) {
    if ((file = fopen (filename, "rb")) == NULL) {
      /* TODO: print/copy an error string */
      retval = -1;
    }
  }
  mpiret = sc_MPI_Bcast (&retval, 1, sc_MPI_INT, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  if (retval) {
    return -1;
  }
#endif

  /* read file metadata on root rank */
  retval = 0;
  if (p4est->mpirank == 0) {
#ifdef P4EST_ENABLE_MPIIO
    if ((mpiret = sc_mpi_read_at (file, 0, metadata,
                                  P4EST_NUM_METADATA_BYTES, sc_MPI_BYTE))
        != sc_MPI_SUCCESS) {
      SC_CHECK_MPI_VERBOSE (mpiret, "Reading metadata for file_info");
      /* TODO record error message: read failed on rank 0 */
      retval = -1;
    }
#else
    P4EST_ASSERT (file != NULL);
    if (fread (metadata, P4EST_NUM_METADATA_BYTES, 1, file) != 1) {
      retval = -1;
    }
#endif
  }
  mpiret = sc_MPI_Bcast (&retval, 1, sc_MPI_INT, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  if (retval) {
    /* Close file on error.  No checking of further errors. */
#ifdef P4EST_ENABLE_MPIIO
    MPI_File_close (&file);
#else
    if (p4est->mpirank == 0) {
      P4EST_ASSERT (file != NULL);
      fclose (file);
    }
    else {
      P4EST_ASSERT (file == NULL);
    }
#endif
    return -1;
  }

  /* broadcast file metadata to all ranks */
  mpiret = sc_MPI_Bcast (metadata, P4EST_NUM_METADATA_BYTES, sc_MPI_BYTE, 0,
                         p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  /* split the input string */
  count = 0;
  parsing_arg = strtok (metadata, "\n");
  P4EST_ASSERT (parsing_arg != NULL);
  while (parsing_arg != NULL && count < 4) {
    if (count == 2) {
      *global_num_quadrants = sc_atol (parsing_arg);
    }
    else if (count == 3) {
      *header_size = sc_atol (parsing_arg);
    }
    parsing_arg = strtok (NULL, "\n");
    ++count;
  }

  /* calculate the padding bytes for the user-defined header */
  get_padding_string (*header_size, P4EST_BYTE_DIV, NULL, &num_pad_bytes);
  padded_header = *header_size + num_pad_bytes;

  if (p4est->mpirank == 0) {
#ifdef P4EST_ENABLE_MPIIO
    mpiret = fill_elem_size (p4est, file, padded_header, elem_size, -1);
#else
    fill_elem_size (p4est, file, padded_header, elem_size, -1);
#endif
  }
#ifdef P4EST_ENABLE_MPIIO
  sc_MPI_Bcast (&mpiret, sizeof (int), sc_MPI_BYTE, 0, p4est->mpicomm);
  if (mpiret != sc_MPI_SUCCESS) {
    return mpiret;
  }
#endif

  /* array metadata */
  current_member = elem_size->elem_size;
  sc_MPI_Bcast (&current_member, sizeof (size_t), sc_MPI_BYTE, 0,
                p4est->mpicomm);
  if (p4est->mpirank != 0) {
    sc_array_init (elem_size, current_member);
  }
  current_member = elem_size->elem_count;
  sc_MPI_Bcast (&current_member, sizeof (size_t), sc_MPI_BYTE, 0,
                p4est->mpicomm);
  if (p4est->mpirank != 0) {
    sc_array_resize (elem_size, current_member);
  }
  sc_MPI_Bcast (elem_size->array,
                elem_size->elem_count * elem_size->elem_size, sc_MPI_BYTE, 0,
                p4est->mpicomm);

#ifdef P4EST_ENABLE_MPIIO
  sc_MPI_File_close (&file);
#else
  if (p4est->mpirank == 0) {
    fclose (file);
  }
#endif
  return sc_MPI_SUCCESS;
}

void
p4est_file_close (p4est_file_context_t * fc)
{
  P4EST_ASSERT (fc != NULL);
#ifdef P4EST_ENABLE_MPIIO
  sc_MPI_File_close (&fc->file);
#else
  if (fc->file != NULL) {
    fclose (fc->file);
  }
#endif
  P4EST_FREE (fc);
}

//#endif /* !ENABLE_MPIIO */
#endif /* P4_TO_P8 */
