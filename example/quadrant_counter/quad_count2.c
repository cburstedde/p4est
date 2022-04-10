
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
#include <p4est_build.h>
#include <p4est_extended.h>
#include <p4est_geometry.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <p4est_communication.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_build.h>
#include <p8est_extended.h>
#include <p8est_geometry.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#include <p8est_communication.h>
#endif

static int
count_callback (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                void *point)
{
  if (local_num == -1) {

    /* Decide if rank of first descendant belongs to partition rank */
    p4est_quadrant_t    fd;
    p4est_tree_t       *tree =
      p4est_tree_array_index (p4est->trees, which_tree);
    p4est_quadrant_first_descendant (quadrant, &fd, (int) tree->maxlevel);
    int                 same_rank = 0;
    same_rank = p4est_comm_is_owner (p4est, which_tree, &fd, p4est->mpirank);

    if (same_rank) {
      ++(*(size_t *) (p4est->user_pointer));
    }
    tree = 0;
    return 1;
  }

  /* Leaf always have an unique owner */
  ++(*(size_t *) (p4est->user_pointer));

  return 0;
}

static size_t
count_quadrants (p4est_t * p4est)
{
  const size_t        num_ranks = p4est->mpisize;
  size_t              count_buffer[num_ranks], count, ii;

  count = 0;
  p4est->user_pointer = &count;

  p4est_search_reorder (p4est, 0, NULL, count_callback, NULL, NULL, NULL);
  sc_MPI_Allgather (&count, 1, sc_MPI_UNSIGNED_LONG, count_buffer, 1,
                    sc_MPI_UNSIGNED_LONG, sc_MPI_COMM_WORLD);
  count = 0;
  if (p4est->mpirank == 0) {
    for (ii = 0; ii < (size_t) p4est->mpisize; ++ii)
      count += count_buffer[ii];
  }

  return count;
}

int
main (int argc, char **argv)
{
  int                 mpiret, ntree;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  size_t              total_quadrants;

  /* MPI init */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  /* package init */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* create p4est */
  ntree = 1;

  if (ntree == 1) {
#ifndef P4_TO_P8
    conn = p4est_connectivity_new_unitsquare ();
#else
    conn = p8est_connectivity_new_unitcube ();
#endif
  }

  p4est = p4est_new_ext (mpicomm, conn, 0, 2, 1, 0, NULL, NULL);
  p4est_partition (p4est, 0, NULL);
#ifndef P4_TO_P8
  p4est_vtk_write_file (p4est, NULL, "p4est_vtk");
#else
  p4est_vtk_write_file (p4est, NULL, "p8est_vtk");
#endif

  total_quadrants = count_quadrants (p4est);
  /* quadrant count is only available on master rank */
  if (p4est->mpirank == 0) {
    printf (" total quadrant count = %ld\n ", total_quadrants);
  }

  /* free memory */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  /* close MPI environment */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return EXIT_SUCCESS;
}
