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

#ifdef P4_TO_P8
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_vtk.h>
#else
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_vtk.h>
#endif

/**
 * Runs all tests.
 */
int
main (int argc, char **argv)
{
  const char         *this_fn_name = P4EST_STRING "_test_subcomm";
  /* options */
  const p4est_locidx_t min_quadrants = 15;
  const int           min_level = 4;
  const int           fill_uniform = 0;
  /* parallel environment */
  sc_MPI_Comm         mpicomm = sc_MPI_COMM_WORLD;
  int                 mpisize, submpisize;
  int                 mpiret;
#ifdef P4EST_ENABLE_DEBUG
  int                 rank;
#endif
  /* p4est */
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  p4est_locidx_t     *partition;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* exit if MPI communicator cannot be reduced */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  if (mpisize == 1) {
    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
    return 0;
  }

  /* initialize p4est */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* create connectivity */
#ifdef P4_TO_P8
  connectivity = p8est_connectivity_new_unitcube ();
#else
  connectivity = p4est_connectivity_new_unitsquare ();
#endif

  /* create p4est object */
  p4est = p4est_new_ext (mpicomm, connectivity,
                         min_quadrants, min_level, fill_uniform,
                         0, NULL, NULL);

  /* write vtk: new */
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_subcomm_new");

  /* set variables pertaining to the parallel environment */
#ifdef P4EST_ENABLE_DEBUG
  rank = p4est->mpirank;
#endif
  submpisize = mpisize / 2;
  P4EST_ASSERT (submpisize <= p4est->global_num_quadrants);

  /* construct partitioning with empty ranks */
  {
    p4est_locidx_t      n_quads_per_proc, n_quads_leftover;
    int                 p;

    partition = P4EST_ALLOC (p4est_locidx_t, mpisize);
    n_quads_per_proc = p4est->global_num_quadrants / submpisize;
    n_quads_leftover = p4est->global_num_quadrants -
      (n_quads_per_proc * submpisize);
    for (p = 0; p < mpisize; p++) {
      if (p % 2) {              /* if this rank will get quadrants */
        partition[p] = n_quads_per_proc;
      }
      else {                    /* if this rank will be empty */
        partition[p] = 0;
      }
    }
    partition[1] += n_quads_leftover;

    /* check partitioning */
#ifdef P4EST_ENABLE_DEBUG
    {
      p4est_gloidx_t      sum = 0;

      for (p = 0; p < mpisize; p++) {
        sum += (p4est_gloidx_t) partition[p];
      }
      P4EST_ASSERT (sum == p4est->global_num_quadrants);
    }
#endif
  }

  /*
   * Test 1: Reduce MPI communicator to non-empty ranks
   */

  P4EST_GLOBAL_INFOF ("%s: Into test 1\n", this_fn_name);
  {
    p4est_t            *p4est_subcomm;
    int                 is_nonempty;

    /* create p4est copy and re-partition */
    p4est_subcomm = p4est_copy_ext (p4est, 1, 1);
    (void) p4est_partition_given (p4est_subcomm, partition);

    /* write vtk: partitioned */
    p4est_vtk_write_file (p4est_subcomm, NULL, P4EST_STRING "_subcomm_part");

    /* reduce MPI communicator to non-empty ranks */
    is_nonempty = p4est_comm_parallel_env_reduce (&p4est_subcomm);
    P4EST_ASSERT ((is_nonempty && 0 < partition[rank]) ||
                  (!is_nonempty && 0 == partition[rank]));

    if (is_nonempty) {
      /* write vtk: reduced communicator */
      p4est_vtk_write_file (p4est_subcomm, NULL,
                            P4EST_STRING "_subcomm_sub1");

      /* destroy the p4est that has a reduced MPI communicator */
      p4est_destroy (p4est_subcomm);
    }
  }
  mpiret = sc_MPI_Barrier (mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_INFOF ("%s: Done test 1\n", this_fn_name);

  /*
   * Test 2: Reduce MPI communicator to non-empty ranks, but now the MPI
   * communicator is not owned
   */

  P4EST_GLOBAL_INFOF ("%s: Into test 2\n", this_fn_name);
  {
    p4est_t            *p4est_subcomm;
    int                 is_nonempty;

    /* create p4est copy and re-partition */
    p4est_subcomm = p4est_copy_ext (p4est, 1, 0 /* don't dup. comm. */ );
    (void) p4est_partition_given (p4est_subcomm, partition);

    /* reduce MPI communicator to non-empty ranks */
    is_nonempty = p4est_comm_parallel_env_reduce (&p4est_subcomm);
    P4EST_ASSERT ((is_nonempty && 0 < partition[rank]) ||
                  (!is_nonempty && 0 == partition[rank]));

    if (is_nonempty) {
      /* destroy the p4est that has a reduced MPI communicator */
      p4est_destroy (p4est_subcomm);
    }
  }
  mpiret = sc_MPI_Barrier (mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_INFOF ("%s: Done test 2\n", this_fn_name);

  /*
   * Test 3: Reduce MPI communicator to non-empty ranks, but keep rank 0
   */

  P4EST_GLOBAL_INFOF ("%s: Into test 3\n", this_fn_name);
  {
    p4est_t            *p4est_subcomm;
    int                 sub_exists;
    sc_MPI_Group        group, group_reserve;
    int                 reserve_range[1][3];

    /* create group of full MPI communicator */
    mpiret = sc_MPI_Comm_group (mpicomm, &group);
    SC_CHECK_MPI (mpiret);

    /* create sub-group containing only rank 0 */
    reserve_range[0][0] = 0;
    reserve_range[0][1] = 0;
    reserve_range[0][2] = 1;
    mpiret =
      sc_MPI_Group_range_incl (group, 1, reserve_range, &group_reserve);
    SC_CHECK_MPI (mpiret);

    /* create p4est copy and re-partition */
    p4est_subcomm = p4est_copy_ext (p4est, 1, 1);
    (void) p4est_partition_given (p4est_subcomm, partition);

    /* reduce MPI communicator to non-empty ranks, but keep rank 0 */
    sub_exists = p4est_comm_parallel_env_reduce_ext (&p4est_subcomm,
                                                     group_reserve, 1, NULL);
    P4EST_ASSERT ((sub_exists && (0 < partition[rank] || rank == 0)) ||
                  (!sub_exists && 0 == partition[rank]));

    if (sub_exists) {
      /* write vtk: reduced communicator */
      p4est_vtk_write_file (p4est_subcomm, NULL,
                            P4EST_STRING "_subcomm_sub3");

      /* destroy the p4est that has a reduced MPI communicator */
      p4est_destroy (p4est_subcomm);
    }
  }
  mpiret = sc_MPI_Barrier (mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_INFOF ("%s: Done test 3\n", this_fn_name);

  /*
   * Test 4: Reduce MPI communicator to non-empty ranks, but keep last 2 ranks
   */

  P4EST_GLOBAL_INFOF ("%s: Into test 4\n", this_fn_name);
  {
    p4est_t            *p4est_subcomm;
    int                 sub_exists;
    sc_MPI_Group        group, group_reserve;
    int                 reserve_range[1][3];

    /* create group of full MPI communicator */
    mpiret = sc_MPI_Comm_group (mpicomm, &group);
    SC_CHECK_MPI (mpiret);

    /* create sub-group containing only last 2 ranks */
    reserve_range[0][0] = SC_MAX (0, mpisize - 2);
    reserve_range[0][1] = mpisize - 1;
    reserve_range[0][2] = 1;
    mpiret =
      sc_MPI_Group_range_incl (group, 1, reserve_range, &group_reserve);
    SC_CHECK_MPI (mpiret);

    /* create p4est copy and re-partition */
    p4est_subcomm = p4est_copy_ext (p4est, 1, 1);
    (void) p4est_partition_given (p4est_subcomm, partition);

    /* reduce MPI communicator to non-empty ranks, but keep last 2 ranks */
    sub_exists = p4est_comm_parallel_env_reduce_ext (&p4est_subcomm,
                                                     group_reserve, 0, NULL);
    P4EST_ASSERT ((sub_exists && (0 < partition[rank] || mpisize - 2 <= rank))
                  || (!sub_exists && 0 == partition[rank]));

    if (sub_exists) {
      /* write vtk: reduced communicator */
      p4est_vtk_write_file (p4est_subcomm, NULL,
                            P4EST_STRING "_subcomm_sub4");

      /* destroy the p4est that has a reduced MPI communicator */
      p4est_destroy (p4est_subcomm);
    }
  }
  mpiret = sc_MPI_Barrier (mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_INFOF ("%s: Done test 4\n", this_fn_name);

  /* destroy */
  P4EST_FREE (partition);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  /* finalize */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
