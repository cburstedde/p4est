
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

/** \file p4est_count_quadrants2.c
 *
 * This 2D example program counts all quadrants in a given tree structure
 * using p4est search routines.
 * Both branch quadrants and leaf quadrants are considered.
 */

#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <p4est_communication.h>
#else
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#include <p8est_communication.h>
#endif

/**
 * Data structure to store quadrant statistics.
 */
typedef struct quadrant_stats
{
  p4est_gloidx_t      local_branch_count;
  p4est_gloidx_t      local_leaf_count;
  p4est_gloidx_t      local_aggregate_count;
} quadrant_stats_t;

/**
 * Call-back routine to record quadrant count during the search routine
 * execution.
 * Recursion is allowed to continue in case of a branch quadrant.
 */
static int
count_callback (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                void *point)
{
  if (local_num == -1) {
    /**
     * In the event a branch quadrant is owned by multiple ranks, 
     * p4est_comm_is_owner always returns true for the lowest rank.
     * This property is useful to determine the unique ownership of
     * branch quadrants.
     * A branch quadrant is counted only by its unique owner rank. 
     */
    if (p4est_comm_is_owner (p4est, which_tree, quadrant, p4est->mpirank)) {
      ++((quadrant_stats_t *) p4est->user_pointer)->local_branch_count;
    }

    return 1;
  }

  /* Leaf quadrants always have a unique owner. */
  ++((quadrant_stats_t *) p4est->user_pointer)->local_leaf_count;

  return 0;
}

/**
 * Routine to count all(branch/leaf) quadrants in a forest structure.
 */
static p4est_gloidx_t
count_quadrants (p4est_t * p4est)
{
  p4est_gloidx_t      global_aggregate_count;
  quadrant_stats_t   *qs;
#ifdef P4EST_ENABLE_DEBUG
  p4est_gloidx_t      global_leaf_count;
#endif

  /**
   * Traverse the tree structure of the forest in a way we visit
   * every quadrant only once.
   */
  qs = (quadrant_stats_t *) p4est->user_pointer;
  p4est_search_reorder (p4est, 0, NULL, count_callback, NULL, NULL, NULL);
  qs->local_aggregate_count = qs->local_leaf_count + qs->local_branch_count;

  /**
   * Gather counts from other processors and do summation.
   * The result is broadcasted to every processor.
   */
  sc_MPI_Allreduce (&qs->local_aggregate_count, &global_aggregate_count, 1,
                    sc_MPI_LONG, sc_MPI_SUM, sc_MPI_COMM_WORLD);

  /**
   * Sanity to check to ensure leaf count is consistent with the
   * forest structure.
   */
#ifdef P4EST_ENABLE_DEBUG
  sc_MPI_Allreduce (&qs->local_leaf_count, &global_leaf_count, 1,
                    sc_MPI_LONG, sc_MPI_SUM, sc_MPI_COMM_WORLD);
  P4EST_ASSERT (global_leaf_count == p4est->global_num_quadrants);
#endif

  return global_aggregate_count;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  quadrant_stats_t    qs;
  p4est_gloidx_t      quadrant_count;

  /* MPI init. */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* package init. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* create p4est. */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_unitsquare ();
#else
  conn = p8est_connectivity_new_unitcube ();
#endif
  p4est = p4est_new_ext (sc_MPI_COMM_WORLD, conn, 0, 2, 1, 0, NULL, NULL);
  p4est_partition (p4est, 0, NULL);
#ifndef P4_TO_P8
  p4est_vtk_write_file (p4est, NULL, "p4est_quadrant_count_vtk");
#else
  p4est_vtk_write_file (p4est, NULL, "p8est_quadrant_count_vtk");
#endif

  /* Count quadrants. */
  qs.local_aggregate_count = qs.local_branch_count = qs.local_leaf_count = 0;
  p4est->user_pointer = &qs;
  quadrant_count = count_quadrants (p4est);

  /* Print quadrant count statistics for current rank. */
  P4EST_VERBOSEF ("Local leaf quadrant count = %ld\n", qs.local_leaf_count);
  P4EST_VERBOSEF ("Local branch quadrant count = %ld\n",
                  qs.local_branch_count);

  /* Print total quadrants in the tree structure. */
  P4EST_VERBOSEF ("Global aggregrate quadrant count = %ld\n", quadrant_count);

  /* free memory. */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  /* close MPI environment. */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return EXIT_SUCCESS;
}
