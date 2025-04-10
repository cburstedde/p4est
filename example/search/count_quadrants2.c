
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

/** \file count_quadrants2.c
 *
 * This 2D example program counts all quadrants in a given tree structure
 * using p4est search routine(s).
 * Both branch quadrants and leaf quadrants are added separetely.
 * Leaf quadrant are always local to a given process, so their global number
 * is independent of the partition of any given mesh refinement structure.
 * We count branch quadrants only for the process that owns a first descendant,
 * which means that the number of branch quadrants is also partition independent.
 * This can be verified manually by comparing the output for multiple MPI sizes.
 */

#ifndef P4_TO_P8
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#else
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#endif

/**
 * Data structure to store quadrant statistics.
 */
typedef struct quadrant_stats
{
  p4est_t            *p4est;
  p4est_locidx_t      local_branch_count;
  p4est_locidx_t      local_leaf_count;
  p4est_locidx_t      local_aggregate_count;
  p4est_gloidx_t      global_aggregate_count;
}
quadrant_stats_t;

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
  quadrant_stats_t   *qs;
  qs = (quadrant_stats_t *) p4est->user_pointer;
  P4EST_ASSERT (qs->p4est == p4est);

  if (local_num == -1) {
    /**
     * In the event a branch quadrant is owned by multiple ranks,
     * p4est_comm_is_owner always returns true for the lowest rank.
     * This property is useful to determine the unique ownership of
     * branch quadrants.
     * A branch quadrant is counted only by its unique owner rank.
     */
    if (p4est_comm_is_owner (p4est, which_tree, quadrant, p4est->mpirank)) {
      ++qs->local_branch_count;
    }

    return 1;
  }

  /* Leaf quadrants always have a unique owner. */
  ++qs->local_leaf_count;

  return 0;
}

/**
 * Routine to count all(branch/leaf) quadrants in a forest structure.
 */
static void
count_quadrants (p4est_t * p4est)
{
  quadrant_stats_t   *qs;
  p4est_gloidx_t      gac;
#ifdef P4EST_ENABLE_DEBUG
  p4est_gloidx_t      glc, global_leaf_count;
#endif

  qs = (quadrant_stats_t *) p4est->user_pointer;
  P4EST_ASSERT (qs->p4est == p4est);

  /**
   * Traverse the tree structure of the forest in a way we visit
   * every quadrant only once.
   */
  p4est_search_reorder (p4est, 0, NULL, count_callback, NULL, NULL, NULL);
  qs->local_aggregate_count = qs->local_leaf_count + qs->local_branch_count;

  /**
   * Gather counts from other processors and do summation.
   * The result is broadcasted to every processor.
   */
  gac = qs->local_aggregate_count;
  sc_MPI_Allreduce (&gac, &qs->global_aggregate_count, 1, P4EST_MPI_GLOIDX,
                    sc_MPI_SUM, sc_MPI_COMM_WORLD);

  /**
   * Sanity check to ensure leaf count is consistent with the
   * forest structure.
   */
#ifdef P4EST_ENABLE_DEBUG
  glc = qs->local_leaf_count;
  sc_MPI_Allreduce (&glc, &global_leaf_count, 1, P4EST_MPI_GLOIDX,
                    sc_MPI_SUM, sc_MPI_COMM_WORLD);
  P4EST_ASSERT (global_leaf_count == p4est->global_num_quadrants);
#endif
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 level, wrong_usage;
  const char         *usage;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  p4est_geometry_t   *geom;
  quadrant_stats_t    qs;

  /* MPI initialization. */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* Package init. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* Process command line arguments. */
  /* *INDENT-OFF* */
  usage =
    "Arguments: <connectivity> <level>\n"
    "   Connectivity can be any of\n"
#ifndef P4_TO_P8
    "         unit|three|moebius|\n"
#else
    "         unit|twocubes|rotcubes|sphere\n"
#endif
    "   Level controls the maximum depth of refinement\n";
  /* *INDENT-ON* */
  wrong_usage = 0;

  /* Query basic parameters. */
  if (!wrong_usage && argc != 3) {
    P4EST_GLOBAL_LERROR ("Invalid argument list\n");
    wrong_usage = 1;
  }
  if (!wrong_usage &&
      ((level = atoi (argv[2])) < 0 || level > P4EST_QMAXLEVEL)) {
    P4EST_GLOBAL_LERRORF ("Level out of bounds 0..%d\n", P4EST_QMAXLEVEL);
    wrong_usage = 1;
  }

  /* Create connectivity. */
  geom = NULL;
  if (!wrong_usage) {
    if (!strcmp (argv[1], "unit")) {
#ifndef P4_TO_P8
      conn = p4est_connectivity_new_unitsquare ();
#else
      conn = p8est_connectivity_new_unitcube ();
#endif
    }
#ifndef P4_TO_P8
    else if (!strcmp (argv[1], "three")) {
      conn = p4est_connectivity_new_corner ();
    }
    else if (!strcmp (argv[1], "moebius")) {
      conn = p4est_connectivity_new_moebius ();
    }
#else
    else if (!strcmp (argv[1], "twocubes")) {
      conn = p8est_connectivity_new_twocubes ();
    }
    else if (!strcmp (argv[1], "rotcubes")) {
      conn = p8est_connectivity_new_rotcubes ();
    }
    else if (!strcmp (argv[1], "sphere")) {
      conn = p8est_connectivity_new_sphere ();
      geom = p8est_geometry_new_sphere (conn, 1., 0.191728, 0.039856);
    }
#endif
    else {
      P4EST_GLOBAL_LERROR ("Invalid connectivity\n");
      wrong_usage = 1;
    }
  }
  if (wrong_usage) {
    P4EST_GLOBAL_LERROR (usage);
    sc_abort_collective ("Usage error");
  }

  /* Create p4est. */
  p4est = p4est_new_ext (sc_MPI_COMM_WORLD, conn, 0, level, 1, 0, NULL, NULL);
  p4est_partition (p4est, 1, NULL);
  p4est_vtk_write_file (p4est, geom, P4EST_STRING "_quadrant_count");

  /* Count quadrants. */
  memset (&qs, 0, sizeof (qs));
  qs.p4est = p4est;
  p4est->user_pointer = &qs;
  count_quadrants (p4est);

  /* Print quadrant count statistics for current rank. */
  P4EST_INFOF ("Local leaf quadrant count = %ld\n",
               (long) qs.local_leaf_count);
  P4EST_INFOF ("Local branch quadrant count = %ld\n",
               (long) qs.local_branch_count);

  /* Print total quadrants (leaf + branch) in the tree structure. */
  mpiret = sc_MPI_Barrier (p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_PRODUCTIONF ("Global aggregrate quadrant count = %lld\n",
                            (long long) qs.global_aggregate_count);

  /* Free memory. */
  if (geom != NULL) {
    p4est_geometry_destroy (geom);
  }
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  /* Close MPI environment. */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return EXIT_SUCCESS;
}
