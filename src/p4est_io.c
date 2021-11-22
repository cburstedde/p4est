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

#define MAGIC_NUMBER 0x123456   /* TODO: compare to other p4est magic num */
#define MAGIC_NUMBER_HTONL 1446253056
#define NUM_METADATA_BYTES 64
#define NUM_ARRAY_METADATA_BYTES 16
#define BYTE_DIV 16
#define FILE_IO_REV 0

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

#if !defined (P4EST_ENABLE_MPIIO) || !defined (P4EST_ENABLE_MPI)
#undef sc_MPI_Offset
#undef sc_MPI_File
#define sc_MPI_Offset long long
#define sc_MPI_File FILE*
#endif

struct p4est_file_context
{
  p4est_t            *p4est;
  sc_MPI_File         file;
  sc_MPI_Offset       accessed_bytes;   /* count only array data bytes and
                                           array metadata bytes */
  size_t              header_size;      /* only the user-defined header */
  size_t              num_calls;        /* redundant but for convience;
                                           counts the number of calls of
                                           write and read repectivly */
#if !defined (P4EST_ENABLE_MPIIO) && defined (P4EST_ENABLE_MPI)
  int                 active;   /* an int that indicates which process
                                   currently opened the file */
  const char         *filename; /* we need to store the path for
                                   the successive opening strategy */
#endif
};

static void
 
 
 
 
p4est_file_info_extra (p4est_file_context_t * fc,
                       p4est_gloidx_t * global_num_quads,
                       char p4est_version[16], int *file_io_rev,
                       int *magic_num, sc_array_t * elem_size);

#if !defined (P4EST_ENABLE_MPI) || !defined (P4EST_ENABLE_MPIIO)
static long long
get_file_size (FILE * f)
{
  int                 ret_val;
  long long           file_size;

  ret_val = fseek (f, 0L, SEEK_END);    /* not okay with the C standard */
  if (ret_val != 0)
    perror ("Error occurred while file seek.\n");
  SC_CHECK_ABORT (ret_val == 0, "seek end of file");
  file_size = ftell (f);
  rewind (f);
  return file_size;
}
#endif

/** This function calculates a padding string consisting of spaces.
 * We require an already allocated array pad or NULL.
 * For NULL the function only calculates the number of padding bytes.
 */
static void
get_padding_string (size_t num_bytes, size_t divisor, char *pad,
                    size_t * num_pad_bytes)
{
  P4EST_ASSERT (divisor != 0 && num_pad_bytes != NULL);

  *num_pad_bytes = (divisor - (num_bytes % divisor)) % divisor;

  if (*num_pad_bytes > 0 && pad != NULL) {
    snprintf (pad, *num_pad_bytes + 1, "%-*s", *((int *) num_pad_bytes), "");
  }
}

static int
check_file_metadata (p4est_t * p4est, size_t header_size,
                     const char *filename, char *metadata)
{
  int                 read_magic_num, read_file_io_rev;
  long                read_global_num_quads, read_header_size;
  int                 count, error_flag;
  char               *parsing_arg;

  P4EST_ASSERT (metadata != NULL);

  count = 0;
  error_flag = 0;
  parsing_arg = strtok (metadata, "\n");
  P4EST_ASSERT (parsing_arg != NULL);
  while (parsing_arg != NULL) {
    if (count == 0) {
      read_magic_num = sc_atoi (parsing_arg);
      if (read_magic_num != MAGIC_NUMBER) {
        /* check for wrong endianess */
        if (read_magic_num == MAGIC_NUMBER_HTONL) {
          P4EST_LERRORF (P4EST_STRING
                         "_io: Error reading <%s>. Maybe wrong endianness because read magic number == htonl (actual magic number).\n",
                         filename);
        }
        P4EST_LERRORF (P4EST_STRING
                       "_io: Error reading <%s>. Wrong magic number (in file = %d, magic number = %d).\n",
                       filename, read_magic_num, MAGIC_NUMBER);
        error_flag = 1;
      }
    }
    else if (count == 2) {
      read_file_io_rev = sc_atoi (parsing_arg);
      if (read_file_io_rev != FILE_IO_REV) {
        P4EST_LERRORF (P4EST_STRING
                       "_io: Error reading <%s>. Wrong file io revision (in file = %d, used file io rev. = %d).\n",
                       filename, read_file_io_rev, FILE_IO_REV);
        error_flag = 1;
      }
    }
    else if (count == 3) {
      read_global_num_quads = sc_atol (parsing_arg);
      if (read_global_num_quads != p4est->global_num_quadrants) {
        P4EST_LERRORF (P4EST_STRING
                       "_io: Error reading <%s>. Wrong global number of quadrants (in file = %ld, in given p4est = %ld).\n",
                       filename, read_global_num_quads,
                       p4est->global_num_quadrants);
        error_flag = 1;
      }
    }
    else if (count == 4) {
      read_header_size = sc_atol (parsing_arg);;
      if (read_header_size != header_size) {
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
  char                metadata[NUM_METADATA_BYTES + 1];
  char                pad[BYTE_DIV];
  size_t              num_pad_bytes;
  p4est_file_context_t *file_context = P4EST_ALLOC (p4est_file_context_t, 1);

  /* Open the file and create a new file if necessary */
#ifdef P4EST_ENABLE_MPIIO
  sc_mpi_open (p4est->mpicomm, filename,
               sc_MPI_MODE_WRONLY | sc_MPI_MODE_CREATE, sc_MPI_INFO_NULL,
               &file_context->file, "File open create");
#elif defined (P4EST_ENABLE_MPI)
  /* serialize the I/O operations */
  /* set active flag */
  file_context->filename = filename;
  if (p4est->mpirank == 0) {
    file_context->file = fopen (filename, "wb");
    file_context->active = 1;
  }
  else {
    file_context->active = 0;
  }
#else
  /* no MPI */
  file_context->file = fopen (filename, "wb");
#endif

  if (p4est->mpirank == 0) {
    /* header_size > 0 must imply header_data != NULL */
    P4EST_ASSERT (header_size <= 0 || header_data != NULL);

    /* write application-defined header */
    snprintf (metadata, NUM_METADATA_BYTES + 1,
              "%d\n%.15s\n%.7d\n%.15ld\n%.15ld\n", MAGIC_NUMBER,
              p4est_version (), FILE_IO_REV, p4est->global_num_quadrants,
              header_size);
#ifdef P4EST_ENABLE_MPIIO
    sc_mpi_write (file_context->file, metadata, NUM_METADATA_BYTES,
                  sc_MPI_BYTE, "Writing the metadata");
#else
    /* this works with and without MPI */
    fwrite (metadata, 1, NUM_METADATA_BYTES, file_context->file);
    fflush (file_context->file);
#endif

    if (header_size != 0) {
      P4EST_ASSERT (header_data != NULL);
      /* Write the user-defined header */
      /* non-collective and blocking */
#ifdef P4EST_ENABLE_MPIIO
      sc_mpi_write (file_context->file, header_data,
                    header_size, sc_MPI_BYTE, "Writing the header");
#else
      /* this works with and without MPI */
      fwrite (header_data, 1, header_size, file_context->file);
      fflush (file_context->file);
#endif

      /* Write padding bytes for the user-defined header */
      get_padding_string (header_size, BYTE_DIV, pad, &num_pad_bytes);
#ifdef P4EST_ENABLE_MPIIO
      sc_mpi_write (file_context->file, pad,
                    num_pad_bytes, sc_MPI_BYTE,
                    "Writing padding bytes for header");
#else
      /* this works with and without MPI */
      fwrite (pad, 1, num_pad_bytes, file_context->file);
      fflush (file_context->file);
#endif
    }
  }
  else {
    get_padding_string (header_size, BYTE_DIV, NULL, &num_pad_bytes);
  }
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
  size_t              num_pad_bytes;
  p4est_file_context_t *file_context = P4EST_ALLOC (p4est_file_context_t, 1);
  sc_MPI_Offset       file_size;

#ifdef P4EST_ENABLE_MPIIO
  /* We do not need the mpi append mode since we use our own byte counter */
  sc_mpi_open (p4est->mpicomm, filename,
               sc_MPI_MODE_WRONLY, sc_MPI_INFO_NULL,
               &file_context->file, "File open append");
#elif defined (P4EST_ENABLE_MPI)
  /* the file is opened in rank-order in \ref p4est_file_write */
  file_context->filename = filename;
  if (p4est->mpirank == 0) {
    file_context->file = fopen (filename, "ab");
    file_context->active = 1;
  }
  else {
    file_context->active = 0;
  }
#else
  /* no MPI */
  file_context->file = fopen (filename, "ab");
#endif

  get_padding_string (header_size, BYTE_DIV, NULL, &num_pad_bytes);
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
    sc_mpi_get_file_size (file_context->file, &file_size,
                          "Get file size for open append");
#else
    /* this works with and without MPI */
    file_size = get_file_size (file_context->file);
#endif
  }
  sc_MPI_Bcast (&file_size, sizeof (sc_MPI_Offset), sc_MPI_BYTE, 0,
                p4est->mpicomm);

  /* Calculate already written quadrant data array bytes */
  file_context->accessed_bytes =
    file_size - NUM_METADATA_BYTES - file_context->header_size;

  return file_context;
}

p4est_file_context_t *
p4est_file_open_read (p4est_t * p4est, const char *filename,
                      size_t header_size, void *header_data)
{
  int                 error_flag;
  size_t              num_pad_bytes;
  char                metadata[NUM_METADATA_BYTES + 1];
  sc_MPI_Offset       file_size;
  p4est_file_context_t *file_context = P4EST_ALLOC (p4est_file_context_t, 1);

#ifdef P4EST_ENABLE_MPIIO
  /* Open the file in the reading mode */
  sc_mpi_open (p4est->mpicomm, filename, sc_MPI_MODE_RDONLY, sc_MPI_INFO_NULL,
               &file_context->file, "File open read");
#elif defined (P4EST_ENABLE_MPI)
  file_context->filename = filename;
  if (p4est->mpirank == 0) {
    file_context->file = fopen (filename, "rb");
    file_context->active = 1;
  }
  else {
    file_context->active = 0;
  }
#else
  file_context->file = fopen (filename, "rb");
#endif

  get_padding_string (header_size, BYTE_DIV, NULL, &num_pad_bytes);
  file_context->p4est = p4est;
  file_context->header_size = header_size + num_pad_bytes;
  file_context->accessed_bytes = 0;
  file_context->num_calls = 0;

  /* read metadata and deallocate in case of error */
  if (file_context->p4est->mpirank == 0) {
    error_flag = 0;

#ifdef P4EST_ENABLE_MPIIO
    /* check size of the file */
    sc_mpi_get_file_size (file_context->file, &file_size,
                          "Get file size for open read");
#else
    file_size = get_file_size (file_context->file);
#endif
    if (header_size > file_size) {
      P4EST_LERRORF (P4EST_STRING
                     "_io: Error reading <%s>. Header_size is bigger than the file size.\n",
                     filename);
      error_flag = 1;
    }

#ifdef P4EST_ENABLE_MPIIO
    /* read metadata on rank 0 */
    sc_mpi_read_at (file_context->file, 0, metadata, NUM_METADATA_BYTES,
                    sc_MPI_BYTE, "Reading metadata");
#else
    /* file pointer is already set to 0 above */
    fread (metadata, 1, NUM_METADATA_BYTES, file_context->file);
#endif
    metadata[NUM_METADATA_BYTES] = '\0';
    /* parse metadata; we do not use file_info because we do not want a Bcast */
    error_flag |=
      check_file_metadata (p4est, header_size, filename, metadata);

    if (!error_flag) {
#ifdef P4EST_ENABLE_MPIIO
      /* read header on rank 0 and skip the metadata */
      sc_mpi_read_at (file_context->file, NUM_METADATA_BYTES, header_data,
                      header_size, sc_MPI_BYTE, "Reading header");
#else
      fseek (file_context->file, NUM_METADATA_BYTES, SEEK_SET);
      fread (header_data, 1, header_size, file_context->file);
#endif
    }
  }

  sc_MPI_Bcast (&error_flag, sizeof (int), sc_MPI_BYTE, 0, p4est->mpicomm);

  if (!error_flag) {
    /* broadcast to all ranks */
    sc_MPI_Bcast (header_data, header_size, sc_MPI_BYTE, 0, p4est->mpicomm);
  }
  else {
    /* error case */
    p4est_file_close (file_context);
    file_context = NULL;
  }

  return file_context;
}

p4est_file_context_t *
p4est_file_write (p4est_file_context_t * fc, sc_array_t * quadrant_data)
{
#if !defined (P4EST_ENABLE_MPIIO) && defined (P4EST_ENABLE_MPI)
  sc_MPI_Status       status;
#endif
  size_t              bytes_to_write, num_pad_bytes, array_size;
  char                array_metadata[NUM_ARRAY_METADATA_BYTES + 1],
    pad[BYTE_DIV];
  sc_MPI_Offset       write_offset;

  P4EST_ASSERT (quadrant_data != NULL
                && quadrant_data->elem_count ==
                fc->p4est->local_num_quadrants);

  if (quadrant_data->elem_size == 0) {
    /* nothing to write */
    return NULL;
  }

  /* Check how many bytes we write to the disk */
  bytes_to_write = quadrant_data->elem_count * quadrant_data->elem_size;

  /* rank-dependent byte offset */
  write_offset = NUM_METADATA_BYTES + fc->header_size +
    fc->p4est->global_first_quadrant[fc->p4est->mpirank] *
    quadrant_data->elem_size;

#ifdef P4EST_ENABLE_MPIIO
  /* set the file size */
  sc_mpi_set_file_size (fc->file,
                        NUM_METADATA_BYTES + fc->header_size +
                        fc->p4est->global_num_quadrants *
                        quadrant_data->elem_size + NUM_ARRAY_METADATA_BYTES +
                        fc->accessed_bytes, "Set file size");
#else
  /* We do not perform this optimization without MPI I/O */
#endif

  if (fc->p4est->mpirank == 0) {
    /* array-dependent metadata */
    snprintf (array_metadata, NUM_ARRAY_METADATA_BYTES + 1, "\n%.14ld\n",
              quadrant_data->elem_size);
#ifdef P4EST_ENABLE_MPIIO
    sc_mpi_write_at (fc->file, fc->accessed_bytes + write_offset,
                     array_metadata, NUM_ARRAY_METADATA_BYTES, sc_MPI_BYTE,
                     "Writing array metadata");
#elif defined (P4EST_ENABLE_MPI)
    /* file already opened on process 0 */
    /* write array metadata */
    fseek (fc->file, fc->accessed_bytes + write_offset, SEEK_SET);
    fwrite (array_metadata, 1, NUM_ARRAY_METADATA_BYTES, fc->file);
    fflush (fc->file);
#else
    /* If the user fills the quadrant data array with chars the sequence of
     * of the fwrite calls plays a role. */

    /* write array metadata */
    fseek (fc->file, fc->accessed_bytes + write_offset, SEEK_SET);
    fwrite (array_metadata, 1, NUM_ARRAY_METADATA_BYTES, fc->file);
    fflush (fc->file);

    /* write array data */
    fseek (fc->file, fc->accessed_bytes + write_offset +
           NUM_ARRAY_METADATA_BYTES, SEEK_SET);
    fwrite (quadrant_data->array, 1, bytes_to_write, fc->file);
    fflush (fc->file);

    /* write padding bytes */
    array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
    get_padding_string (array_size, BYTE_DIV, pad, &num_pad_bytes);
    fseek (fc->file, fc->accessed_bytes + NUM_METADATA_BYTES +
           fc->header_size + array_size + NUM_ARRAY_METADATA_BYTES, SEEK_SET);
    fwrite (pad, 1, num_pad_bytes, fc->file);
    fflush (fc->file);
#endif

#ifdef P4EST_ENABLE_MPIIO
    /* Caculate and write padding bytes for array data */
    array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
    get_padding_string (array_size, BYTE_DIV, pad, &num_pad_bytes);

    sc_mpi_write_at (fc->file,
                     fc->accessed_bytes + NUM_METADATA_BYTES +
                     fc->header_size + array_size + NUM_ARRAY_METADATA_BYTES,
                     pad, num_pad_bytes, sc_MPI_BYTE,
                     "Writing padding bytes for a data array");
#endif
  }
  else {
    array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
    get_padding_string (array_size, BYTE_DIV, NULL, &num_pad_bytes);
  }

#ifdef P4EST_ENABLE_MPIIO
  sc_mpi_write_at_all (fc->file,
                       fc->accessed_bytes + write_offset +
                       NUM_ARRAY_METADATA_BYTES, quadrant_data->array,
                       bytes_to_write, sc_MPI_BYTE, "Writing quadrant-wise");
#elif defined (P4EST_ENABLE_MPI)
  if (fc->p4est->mpirank != 0) {
    /* wait until the preceding process finished the I/O operation */
    /* TODO: check the return value */
    /* receive */
    sc_MPI_Recv (&fc->active, 1, sc_MPI_INT,
                 fc->p4est->mpirank - 1, sc_MPI_ANY_TAG, fc->p4est->mpicomm,
                 &status);
  }

  if (fc->active) {
    /* process 0 must not wait */
    if (fc->p4est->mpirank != 0) {
      /* open */
      fc->file = fopen (fc->filename, "ab");
    }

    /* write array data */
    fseek (fc->file, fc->accessed_bytes + write_offset +
           NUM_ARRAY_METADATA_BYTES, SEEK_SET);
    fwrite (quadrant_data->array, 1, bytes_to_write, fc->file);
    fflush (fc->file);

    /* close */
    fclose (fc->file);

    /* only update active process if there are processes left */
    if (fc->p4est->mpirank < fc->p4est->mpisize - 1) {
      /* the current process finished its I/O operations */
      P4EST_ASSERT (fc->active == 1);
      /* send */
      sc_MPI_Send (&fc->active, 1, sc_MPI_INT,
                   fc->p4est->mpirank + 1, 1, fc->p4est->mpicomm);
    }
  }
#endif

  /* This is *not* the processor local value */
  fc->accessed_bytes +=
    quadrant_data->elem_size * fc->p4est->global_num_quadrants +
    NUM_ARRAY_METADATA_BYTES + num_pad_bytes;
  ++fc->num_calls;
#if !defined (P4EST_ENABLE_MPIIO) && defined (P4EST_ENABLE_MPI)
  /* The processes have to wait here because they are not allowed to start
   * other I/O operations.
   */
  sc_MPI_Barrier (fc->p4est->mpicomm);
  if (fc->p4est->mpirank == 0) {
    fc->file = fopen (fc->filename, "ab");
    /* write padding bytes */
    array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
    get_padding_string (array_size, BYTE_DIV, pad, &num_pad_bytes);
    fseek (fc->file, fc->accessed_bytes + NUM_METADATA_BYTES +
           fc->header_size + array_size + NUM_ARRAY_METADATA_BYTES, SEEK_SET);
    fwrite (pad, 1, num_pad_bytes, fc->file);
    fflush (fc->file);
    /* close */
    fclose (fc->file);
  }

  /* reset the processor activity values */
  if (fc->p4est->mpirank == 0) {
    fc->active = 1;
  }
  else {
    fc->active = 0;
  }
  fc->file = (fc->p4est->mpirank == 0) ? fopen (fc->filename, "ab") : NULL;
#endif

  return fc;
}

p4est_file_context_t *
p4est_file_read (p4est_file_context_t * fc, sc_array_t * quadrant_data)
{
#if !defined (P4EST_ENABLE_MPIIO) && defined (P4EST_ENABLE_MPI)
  sc_MPI_Status       status;
#endif
  int                 error_flag;
  size_t              bytes_to_read, num_pad_bytes, array_size, *data_size,
    read_data_size;
  char                array_metadata[NUM_ARRAY_METADATA_BYTES + 1];
  sc_MPI_Offset       size;
  sc_array_t          elem_size;

  P4EST_ASSERT (fc != NULL);

  if (quadrant_data == NULL || quadrant_data->elem_size == 0) {
    /* Nothing to read but we shift our own file pointer */

    sc_array_init (&elem_size, sizeof (size_t));
    p4est_file_info_extra (fc, NULL, NULL, NULL, NULL, &elem_size);     /* TODO: we do not need the whole array */
    data_size = (size_t *) sc_array_index (&elem_size, fc->num_calls);
    /* calculate the padding bytes for this data array */
    array_size = fc->p4est->global_num_quadrants * *data_size;
    get_padding_string (array_size, BYTE_DIV, NULL, &num_pad_bytes);
    fc->accessed_bytes +=
      *data_size * fc->p4est->global_num_quadrants +
      NUM_ARRAY_METADATA_BYTES + num_pad_bytes;
    ++fc->num_calls;
    sc_array_reset (&elem_size);
    return NULL;
  }

  P4EST_ASSERT (quadrant_data->elem_count == fc->p4est->local_num_quadrants);

  /* check how many bytes we read from the disk */
  bytes_to_read = quadrant_data->elem_count * quadrant_data->elem_size;

#ifdef P4EST_ENABLE_MPIIO
  /* check file size; no sync required because the file size does not
   * change in the reading mode.
   */
  sc_mpi_get_file_size (fc->file, &size, "Get file size for read");
#elif defined (P4EST_ENABLE_MPI)
  /* only one process can access the file */
  if (fc->p4est->mpirank == 0) {        /* TODO: really needed? */
    size = get_file_size (fc->file);
  }
  sc_MPI_Bcast (&size, sizeof (long long), sc_MPI_BYTE, 0,
                fc->p4est->mpicomm);
#else
  size = get_file_size (fc->file);
#endif
  if (size - NUM_METADATA_BYTES - fc->header_size < bytes_to_read) {
    /* report wrong file size, collectively close the file and deallocate fc */
    if (fc->p4est->mpirank == 0) {
      P4EST_LERROR (P4EST_STRING
                    "_io: Error reading. File has less bytes than the user wants to read.\n");
    }
    p4est_file_close (fc);
    return NULL;
  }

  /* check the array metadata */
  error_flag = 0;
  if (fc->p4est->mpirank == 0) {
#ifdef P4EST_ENABLE_MPIIO
    sc_mpi_read_at (fc->file,
                    fc->accessed_bytes + NUM_METADATA_BYTES +
                    fc->header_size +
                    fc->p4est->global_first_quadrant[fc->p4est->mpirank] *
                    quadrant_data->elem_size, array_metadata,
                    NUM_ARRAY_METADATA_BYTES, sc_MPI_BYTE,
                    "Reading quadrant-wise metadata");
#else
    fseek (fc->file, fc->accessed_bytes + NUM_METADATA_BYTES + fc->header_size + fc->p4est->global_first_quadrant[fc->p4est->mpirank] * // TODO: Use 0
           quadrant_data->elem_size, SEEK_SET);
    fread (array_metadata, 1, NUM_ARRAY_METADATA_BYTES, fc->file);
#endif
    array_metadata[NUM_ARRAY_METADATA_BYTES] = '\0';
    read_data_size = sc_atol (array_metadata);  // TODO: Check before the reading!
    if (read_data_size != quadrant_data->elem_size) {
      P4EST_LERRORF (P4EST_STRING
                     "_io: Error reading. Wrong array element size (in file = %ld, by parameter = %ld).\n",
                     read_data_size, quadrant_data->elem_size);
      error_flag = 1;
    }
  }
  /* broadcast the error flag to decicde if we continue */
  sc_MPI_Bcast (&error_flag, sizeof (int), sc_MPI_BYTE, 0,
                fc->p4est->mpicomm);
  if (error_flag) {
    p4est_file_close (fc);
    return NULL;
  }

  /* calculate the padding bytes for this data array */
  array_size = fc->p4est->global_num_quadrants * quadrant_data->elem_size;
  get_padding_string (array_size, BYTE_DIV, NULL, &num_pad_bytes);

#ifdef P4EST_ENABLE_MPIIO
  sc_mpi_read_at_all (fc->file,
                      fc->accessed_bytes + NUM_METADATA_BYTES +
                      NUM_ARRAY_METADATA_BYTES + fc->header_size +
                      fc->p4est->global_first_quadrant[fc->p4est->mpirank] *
                      quadrant_data->elem_size, quadrant_data->array,
                      bytes_to_read, sc_MPI_BYTE, "Reading quadrant-wise");
#elif defined (P4EST_ENABLE_MPI)
  if (fc->p4est->mpirank != 0) {
    /* wait until the preceding process finished the I/O operation */
    /* TODO: check the return value */
    /* receive */
    sc_MPI_Recv (&fc->active, 1, sc_MPI_INT,
                 fc->p4est->mpirank - 1, sc_MPI_ANY_TAG, fc->p4est->mpicomm,
                 &status);
  }

  if (fc->active) {
    /* process 0 must not wait */
    if (fc->p4est->mpirank != 0) {
      /* open */
      fc->file = fopen (fc->filename, "rb");
    }

    /* read array data */
    fseek (fc->file, fc->accessed_bytes + NUM_METADATA_BYTES +
           NUM_ARRAY_METADATA_BYTES + fc->header_size +
           fc->p4est->global_first_quadrant[fc->p4est->mpirank] *
           quadrant_data->elem_size, SEEK_SET);
    fread (quadrant_data->array, 1, bytes_to_read, fc->file);

    /* close */
    fclose (fc->file);
    /* only update active process if there are processes left */
    if (fc->p4est->mpirank < fc->p4est->mpisize - 1) {
      /* the current process finished its I/O operations */
      P4EST_ASSERT (fc->active == 1);
      /* send */
      sc_MPI_Send (&fc->active, 1, sc_MPI_INT,
                   fc->p4est->mpirank + 1, 1, fc->p4est->mpicomm);
    }
  }
#else
  fseek (fc->file, fc->accessed_bytes + NUM_METADATA_BYTES +
         NUM_ARRAY_METADATA_BYTES + fc->header_size +
         fc->p4est->global_first_quadrant[fc->p4est->mpirank] *
         quadrant_data->elem_size, SEEK_SET);
  fread (quadrant_data->array, 1, bytes_to_read, fc->file);
#endif

  fc->accessed_bytes +=
    quadrant_data->elem_size * fc->p4est->global_num_quadrants +
    NUM_ARRAY_METADATA_BYTES + num_pad_bytes;
  ++fc->num_calls;
#if !defined (P4EST_ENABLE_MPIIO) && defined (P4EST_ENABLE_MPI)
/* The processes have to wait here because they are not allowed to start
     * other I/O operations.
     */
  sc_MPI_Barrier (fc->p4est->mpicomm);
  /* reset the processor activity values */
  if (fc->p4est->mpirank == 0) {
    fc->active = 1;
  }
  else {
    fc->active = 0;
  }
  fc->file = (fc->p4est->mpirank == 0) ? fopen (fc->filename, "rb") : NULL;
#endif

  return fc;
}

/* always executed on process 0 exclusively */
static void
fill_elem_size (p4est_t * p4est, sc_MPI_File file, size_t header_size,
                sc_array_t * elem_size)
{
  char                array_metadata[NUM_ARRAY_METADATA_BYTES + 1];
  char               *parsing_arg;
  size_t              num_pad_bytes, array_size;
  long               *new_elem;
  sc_MPI_Offset       size, current_position;

  /* read the array metadata */
  P4EST_ASSERT (elem_size->elem_size == sizeof (size_t));
  sc_array_resize (elem_size, 0);

#ifdef P4EST_ENABLE_MPIIO
  /* get the file size */
  sc_mpi_get_file_size (file, &size, "Get file size for info_extra");
#else
  size = get_file_size (file);
#endif
  /* We read the metadata and we skip the user-defined header */
  current_position = NUM_METADATA_BYTES + header_size;
  while (current_position < size) {
#ifdef P4EST_ENABLE_MPIIO
    sc_mpi_read_at (file, current_position, array_metadata,
                    NUM_ARRAY_METADATA_BYTES, sc_MPI_BYTE,
                    "Reading array metadata");
#else
    fseek (file, current_position, SEEK_SET);
    fread (array_metadata, 1, NUM_ARRAY_METADATA_BYTES, file);
#endif

    array_metadata[NUM_ARRAY_METADATA_BYTES] = '\0';
    /* parse and store the element size of the array */
    parsing_arg = strtok (array_metadata, "\n");
    P4EST_ASSERT (parsing_arg != NULL);
    new_elem = (long *) sc_array_push (elem_size);
    *new_elem = sc_atol (parsing_arg);

    /* get padding bytes of the current array */
    array_size = p4est->global_num_quadrants * *new_elem;
    get_padding_string (array_size, BYTE_DIV, NULL, &num_pad_bytes);
    current_position += array_size + NUM_ARRAY_METADATA_BYTES + num_pad_bytes;
  }
}

static void
p4est_file_info_extra (p4est_file_context_t * fc,
                       p4est_gloidx_t * global_num_quads,
                       char p4est_version[16], int *file_io_rev,
                       int *magic_num, sc_array_t * elem_size)
{
  int                 read_file_metadata, count;
  char                metadata[NUM_METADATA_BYTES + 1];
  char               *parsing_arg;
  size_t              current_member;

  P4EST_ASSERT (fc != NULL);

  read_file_metadata = global_num_quads != NULL || p4est_version != NULL
    || file_io_rev != NULL || magic_num != NULL;
  if (fc->p4est->mpirank == 0) {
    if (read_file_metadata) {
#ifdef P4EST_ENABLE_MPIIO
      /* read metadata on rank 0 */
      sc_mpi_read_at (fc->file, 0, metadata, NUM_METADATA_BYTES, sc_MPI_BYTE,
                      "Reading metadata");
#else
      fread (metadata, 1, NUM_METADATA_BYTES, fc->file);
#endif
    }
    metadata[NUM_METADATA_BYTES] = '\0';
    fill_elem_size (fc->p4est, fc->file, fc->header_size, elem_size);
  }
  if (read_file_metadata) {
    /* broadcast to all ranks */
    /* file metadata */
    sc_MPI_Bcast (metadata, NUM_METADATA_BYTES, sc_MPI_BYTE, 0,
                  fc->p4est->mpicomm);
  }

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
                elem_size->elem_count * elem_size->elem_size, sc_MPI_BYTE, 0,
                fc->p4est->mpicomm);

  if (read_file_metadata) {
    /* split the input string */
    count = 0;
    parsing_arg = strtok (metadata, "\n");
    P4EST_ASSERT (parsing_arg != NULL);
    while (parsing_arg != NULL) {
      if (magic_num != NULL && count == 0) {
        *magic_num = sc_atoi (parsing_arg);
      }
      else if (p4est_version != NULL && count == 1) {
        strcpy (p4est_version, parsing_arg);
      }
      else if (file_io_rev != NULL && count == 2) {
        *file_io_rev = sc_atoi (parsing_arg);
      }
      else if (global_num_quads != NULL && count == 3) {
        *global_num_quads = sc_atol (parsing_arg);
      }
      parsing_arg = strtok (NULL, "\n");
      ++count;
    }
  }
}

void
p4est_file_info (p4est_t * p4est, const char *filename,
                 p4est_gloidx_t * global_num_quadrants, size_t * header_size,
                 sc_array_t * elem_size)
{
  sc_MPI_File         file;
  int                 count;
  char                metadata[NUM_METADATA_BYTES];
  char               *parsing_arg;
  size_t              current_member, num_pad_bytes, padded_header;

  P4EST_ASSERT (p4est != NULL && filename != NULL
                && global_num_quadrants != NULL && elem_size != NULL);
  P4EST_ASSERT (elem_size->elem_size == sizeof (size_t));

#ifdef P4EST_ENABLE_MPIIO
  /* Open the file in the reading mode */
  sc_mpi_open (p4est->mpicomm, filename, sc_MPI_MODE_RDONLY, sc_MPI_INFO_NULL,
               &file, "File open file_info");
#else
  file = fopen (filename, "rb");
#endif

  /* read file metadata */
  if (p4est->mpirank == 0) {
#ifdef P4EST_ENABLE_MPIIO
    /* read metadata on rank 0 */
    sc_mpi_read_at (file, 0, metadata, NUM_METADATA_BYTES, sc_MPI_BYTE,
                    "Reading metadata for file_info");
#else
    fread (metadata, 1, NUM_METADATA_BYTES, file);
#endif
  }

  /* broadcast to all ranks */
  /* file metadata */
  sc_MPI_Bcast (metadata, NUM_METADATA_BYTES, sc_MPI_BYTE, 0, p4est->mpicomm);

  /* split the input string */
  count = 0;
  parsing_arg = strtok (metadata, "\n");
  P4EST_ASSERT (parsing_arg != NULL);
  while (parsing_arg != NULL) {
    if (count == 3) {
      *global_num_quadrants = sc_atol (parsing_arg);
    }
    else if (count == 4) {
      *header_size = sc_atol (parsing_arg);
    }
    parsing_arg = strtok (NULL, "\n");
    ++count;
  }

  /* calculate the padding bytes for the user-defined header */
  get_padding_string (*header_size, BYTE_DIV, NULL, &num_pad_bytes);
  padded_header = *header_size + num_pad_bytes;

  if (p4est->mpirank == 0) {
    fill_elem_size (p4est, file, padded_header, elem_size);
  }

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
  sc_mpi_close (&file, "File close");
#else
  fclose (file);
#endif
}

void
p4est_file_close (p4est_file_context_t * fc)
{
  P4EST_ASSERT (fc != NULL);
#ifdef P4EST_ENABLE_MPIIO
  sc_mpi_close (&fc->file, "File close");
#else
  if (fc->file != NULL) {
    fclose (fc->file);
  }
#endif
  P4EST_FREE (fc);
}

//#endif /* !ENABLE_MPIIO */
#endif /* P4_TO_P8 */
