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
#define NUM_METADATA_BYTES 56
#define NUM_ARRAY_METADATA_BYTES 16
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
#ifdef P4EST_ENABLE_MPIIO

struct p4est_file_context
{
  p4est_t            *p4est;
  sc_MPI_File         file;
  size_t              header_size;      /* only the user-defined header */
};

p4est_file_context_t *
p4est_file_open_create (p4est_t * p4est, const char *filename,
                        size_t header_size, void *header_data)
{
  char                metadata[NUM_METADATA_BYTES + 1];
  p4est_file_context_t *file_context = P4EST_ALLOC (p4est_file_context_t, 1);

  /* Open the file and create a new file if necessary */
  sc_mpi_open (p4est->mpicomm, filename,
               sc_MPI_MODE_WRONLY | sc_MPI_MODE_CREATE, sc_MPI_INFO_NULL,
               &file_context->file, "File open create");

  if (p4est->mpirank == 0) {
    /* header_size > 0 must imply header_data != NULL */
    P4EST_ASSERT (header_size <= 0 || header_data != NULL);

    /* write application-defined header */
    snprintf (metadata, NUM_METADATA_BYTES + 1, "%d\n%.15s\n%.15d\n%.15ld\n",
              MAGIC_NUMBER, p4est_version (), FILE_IO_REV,
              p4est->global_num_quadrants);
    sc_mpi_write (file_context->file, metadata,
                  NUM_METADATA_BYTES, sc_MPI_CHAR, "Writing the metadata");

    if (header_data != NULL && header_size != 0) {
      /* Write the user-defined header */
      /* non-collective and blocking */
      sc_mpi_write (file_context->file, header_data,
                    header_size, sc_MPI_CHAR, "Writing the header");
    }
  }
  file_context->p4est = p4est;
  file_context->header_size = header_size;

  return file_context;
}

p4est_file_context_t *
p4est_file_open_read (p4est_t * p4est, const char *filename,
                      size_t header_size, void *header_data)
{
  sc_MPI_Offset       file_size;
  p4est_file_context_t *file_context = P4EST_ALLOC (p4est_file_context_t, 1);

  /* Open the file in the reading mode */
  sc_mpi_open (p4est->mpicomm, filename, sc_MPI_MODE_RDONLY, sc_MPI_INFO_NULL,
               &file_context->file, "File open read");
  file_context->p4est = p4est;
  file_context->header_size = header_size;

  /* check size of the file */
  sc_mpi_get_file_size (file_context->file, &file_size, "Get file size");
  SC_CHECK_ABORT (header_size < file_size,
                  "header_size is bigger than the file size");

  if (p4est->mpirank == 0) {
    /* set file pointer to skip the metadata */
    sc_mpi_file_seek (file_context->file,
                      NUM_METADATA_BYTES, sc_MPI_SEEK_SET,
                      "Set file pointer");

    /* read header on rank 0 */
    sc_mpi_read (file_context->file, header_data, header_size, sc_MPI_CHAR,
                 "Reading header");
  }
  /* broadcast to all ranks */
  sc_MPI_Bcast (header_data, header_size, sc_MPI_CHAR, 0, p4est->mpicomm);

  return file_context;
}

void
p4est_file_write (p4est_file_context_t * fc, sc_array_t * quadrant_data)
{
  size_t              bytes_to_write;
  char                array_metadata[NUM_ARRAY_METADATA_BYTES + 1];

  P4EST_ASSERT (quadrant_data != NULL
                && quadrant_data->elem_count ==
                fc->p4est->local_num_quadrants);

  if (quadrant_data->elem_size == 0) {
    /* nothing to write */
    return;
  }

  /* set file size (collective) */
  sc_mpi_set_file_size (fc->file,
                        NUM_METADATA_BYTES + fc->header_size +
                        fc->p4est->global_num_quadrants *
                        quadrant_data->elem_size, "Set file size");

  /* Check how many bytes we write to the disk */
  bytes_to_write = quadrant_data->elem_count * quadrant_data->elem_size;

  /* set file pointer */
  sc_mpi_file_seek (fc->file,
                    NUM_METADATA_BYTES + fc->header_size +
                    fc->p4est->global_first_quadrant[fc->p4est->mpirank] *
                    quadrant_data->elem_size, sc_MPI_SEEK_SET,
                    "Set file pointer");

  /* array-dependent metadata */
  snprintf (array_metadata, NUM_ARRAY_METADATA_BYTES + 1, "\n%.14ld\n",
            quadrant_data->elem_size);
  sc_mpi_write_all (fc->file, array_metadata,
                    NUM_ARRAY_METADATA_BYTES, sc_MPI_CHAR,
                    "Writing array metadata");

  sc_mpi_write_all (fc->file, quadrant_data->array,
                    bytes_to_write, sc_MPI_CHAR, "Writing quadrant-wise");
}

void
p4est_file_read (p4est_file_context_t * fc, sc_array_t * quadrant_data)
{
  size_t              bytes_to_read;
  sc_MPI_Offset       size;

  P4EST_ASSERT (fc != NULL);

  if (quadrant_data == NULL || quadrant_data->elem_size == 0) {
    /* nothing to read */
    return;
  }

  P4EST_ASSERT (quadrant_data->elem_count == fc->p4est->local_num_quadrants);

  /* Check how many bytes we read from the disk */
  bytes_to_read = quadrant_data->elem_count * quadrant_data->elem_size;

  /* check file size */
  sc_mpi_get_file_size (fc->file, &size, "Get file size");
  SC_CHECK_ABORT (size - NUM_METADATA_BYTES - fc->header_size >=
                  bytes_to_read,
                  "File has less bytes than the user wants to read");

  /* set file pointer */
  sc_mpi_file_seek (fc->file,
                    NUM_METADATA_BYTES + NUM_ARRAY_METADATA_BYTES +
                    fc->header_size +
                    fc->p4est->global_first_quadrant[fc->p4est->mpirank] *
                    quadrant_data->elem_size, sc_MPI_SEEK_SET,
                    "Set file pointer");

  sc_mpi_read_all (fc->file, quadrant_data->array,
                   bytes_to_read, sc_MPI_CHAR, "Reading quadrant-wise");
}

void
p4est_file_info (p4est_file_context_t * fc, p4est_gloidx_t * global_num_quads,
                 char p4est_version[16], int *file_io_rev, int *magic_num,
                 sc_array_t * elem_size)
{
  int                 count;
  char                metadata[NUM_METADATA_BYTES + 1],
    array_metadata[NUM_ARRAY_METADATA_BYTES + 1];
  char               *parsing_arg;
  size_t              current_member;
  long               *new_elem;
  sc_MPI_Offset       size, current_position;

  P4EST_ASSERT (fc != NULL);

  if (fc->p4est->mpirank == 0) {
    /* read metadata on rank 0 */
    sc_mpi_read_at (fc->file, 0, metadata, NUM_METADATA_BYTES, sc_MPI_CHAR,
                    "Reading metadata");

    /* read the array metadata */
    sc_array_init (elem_size, sizeof (long));
    sc_array_resize (elem_size, 0);

    /* get the file size */
    sc_mpi_get_file_size (fc->file, &size, "Get file size");
    /* We read the metadata and we skip the user-defined header */
    current_position = NUM_METADATA_BYTES + fc->header_size;
    while (current_position + NUM_ARRAY_METADATA_BYTES < size) {
      sc_mpi_read_at (fc->file, current_position, array_metadata,
                      NUM_ARRAY_METADATA_BYTES, sc_MPI_CHAR,
                      "Reading array metadata");

      /* parse and store the element size of the array */
      parsing_arg = strtok (array_metadata, "\n");
      new_elem = (long *) sc_array_push (elem_size);
      *new_elem = sc_atol (parsing_arg);

      current_position += fc->p4est->global_num_quadrants * *new_elem;
    }
  }
  /* broadcast to all ranks */
  /* file metadata */
  sc_MPI_Bcast (metadata, NUM_METADATA_BYTES + 1, sc_MPI_CHAR, 0,
                fc->p4est->mpicomm);
  /* array metadata */
  current_member = elem_size->elem_size;
  sc_MPI_Bcast (&current_member, sizeof (size_t), sc_MPI_CHAR, 0,
                fc->p4est->mpicomm);
  if (fc->p4est->mpirank != 0) {
    sc_array_init (elem_size, current_member);
  }
  current_member = elem_size->elem_count;
  sc_MPI_Bcast (&current_member, sizeof (size_t), sc_MPI_CHAR, 0,
                fc->p4est->mpicomm);
  if (fc->p4est->mpirank != 0) {
    sc_array_resize (elem_size, current_member);
  }
  sc_MPI_Bcast (elem_size->array,
                elem_size->elem_count * elem_size->elem_size, sc_MPI_CHAR, 0,
                fc->p4est->mpicomm);
  /* add null termination for atoi */
  metadata[NUM_METADATA_BYTES] = '\0';  /* TODO: Really necessary? */

  /* split the input string */
  count = 0;
  parsing_arg = strtok (metadata, "\n");
  while (parsing_arg != NULL) {
    if (count == 0) {
      *magic_num = sc_atoi (parsing_arg);
    }
    else if (count == 1) {
      strcpy (p4est_version, parsing_arg);
    }
    else if (count == 2) {
      *file_io_rev = sc_atoi (parsing_arg);
    }
    else if (count == 3) {
      *global_num_quads = sc_atol (parsing_arg);
    }
    parsing_arg = strtok (NULL, "\n");
    ++count;
  }

}

void
p4est_file_close (p4est_file_context_t * fc)
{
  P4EST_ASSERT (fc != NULL);
  sc_mpi_close (&fc->file, "File close");
  P4EST_FREE (fc);
}

#endif /* !ENABLE_MPIIO */
#endif /* P4_TO_P8 */
