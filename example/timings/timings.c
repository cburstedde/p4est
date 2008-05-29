/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Usage: p4est_timings <configuration> <level>
 *        possible configurations:
 *        o unit      Refinement on the unit square.
 *        o three     Refinement on a forest with three trees.
 *        o moebius   Refinement on a 5-tree Moebius band.
 *        o star      Refinement on a 6-tree star shaped domain.
 *        o periodic  Refinement on the unit square with periodic b.c.
 */

#include <p4est_algorithms.h>
#include <p4est_vtk.h>

enum
{
  P4EST_CONFIG_NULL,
  P4EST_CONFIG_UNIT,
  P4EST_CONFIG_THREE,
  P4EST_CONFIG_MOEBIUS,
  P4EST_CONFIG_STAR,
  P4EST_CONFIG_PERIODIC,
};

typedef struct
{
  MPI_Comm            mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

static int          refine_level = 0;
static int          level_shift = 0;

static int
refine_fractal (p4est_t * p4est, p4est_locidx_t which_tree,
                p4est_quadrant_t * q)
{
  int                 qid;

  if ((int) q->level >= refine_level) {
    return 0;
  }
  if ((int) q->level < refine_level - level_shift) {
    return 1;
  }

  qid = p4est_quadrant_child_id (q);
  return (qid == 0 || qid == 3);
}

static void
abort_fn (void *data)
{
  int                 mpiret;
  mpi_context_t      *mpi = data;

  fprintf (stderr, "[%d] p4est_timings abort handler\n", mpi->mpirank);

  mpiret = MPI_Abort (mpi->mpicomm, 1);
  P4EST_CHECK_MPI (mpiret);
}

int
main (int argc, char **argv)
{
  int                 i;
  int                 mpiret;
  int                 wrongusage, config;
  unsigned            crc;
  p4est_locidx_t     *quadrant_counts;
  p4est_gloidx_t      count_refined, count_balanced;
  p4est_gloidx_t      prev_quadrant, next_quadrant;
  p4est_gloidx_t      global_shipped;
  const char         *config_name, *usage, *errmsg;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  double              start, elapsed_refine;
  double              elapsed_balance, elapsed_rebalance;
  double              elapsed_partition, elapsed_repartition;
  mpi_context_t       mpi_context, *mpi = &mpi_context;

  /* initialize MPI and p4est internals */
  mpiret = MPI_Init (&argc, &argv);
  P4EST_CHECK_MPI (mpiret);
  mpi->mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  P4EST_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  P4EST_CHECK_MPI (mpiret);
  
  sc_init (mpi->mpirank, abort_fn, mpi, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  usage =
    "Arguments: <configuration> <level>\n"
    "   Configuration can be any of\n"
    "      unit|three|moebius|star|periodic\n"
    "   Level controls the maximum depth of refinement\n";
  errmsg = NULL;
  wrongusage = 0;
  config = P4EST_CONFIG_NULL;
  config_name = argv[1];
  if (!wrongusage && argc != 3) {
    wrongusage = 1;
  }
  if (!wrongusage) {
    if (!strcmp (config_name, "unit")) {
      config = P4EST_CONFIG_UNIT;
    }
    else if (!strcmp (config_name, "three")) {
      config = P4EST_CONFIG_THREE;
    }
    else if (!strcmp (config_name, "moebius")) {
      config = P4EST_CONFIG_MOEBIUS;
    }
    else if (!strcmp (config_name, "star")) {
      config = P4EST_CONFIG_STAR;
    }
    else if (!strcmp (config_name, "periodic")) {
      config = P4EST_CONFIG_PERIODIC;
    }
    else {
      wrongusage = 1;
    }
  }
  if (wrongusage) {
    if (mpi->mpirank == 0) {
      fputs ("Usage error\n", stderr);
      fputs (usage, stderr);
      if (errmsg != NULL) {
        fputs (errmsg, stderr);
      }
      sc_abort ();
    }
    mpiret = MPI_Barrier (mpi->mpicomm);
    P4EST_CHECK_MPI (mpiret);
  }

#ifndef P4EST_DEBUG
  sc_set_log_defaults (NULL, SC_LP_STATISTICS, SC_FP_KEEP);
#endif

  /* get command line argument: maximum refinement level */
  refine_level = atoi (argv[2]);
  level_shift = 4;

  /* print general setup information */
  P4EST_GLOBAL_STATISTICSF
    ("Processors %d configuration %s level %d shift %d\n", mpi->mpisize,
     config_name, refine_level, level_shift);

  /* create connectivity and forest structures */
  if (config == P4EST_CONFIG_THREE) {
    connectivity = p4est_connectivity_new_corner ();
  }
  else if (config == P4EST_CONFIG_MOEBIUS) {
    connectivity = p4est_connectivity_new_moebius ();
  }
  else if (config == P4EST_CONFIG_STAR) {
    connectivity = p4est_connectivity_new_star ();
  }
  else if (config == P4EST_CONFIG_PERIODIC) {
    connectivity = p4est_connectivity_new_periodic ();
  }
  else {
    connectivity = p4est_connectivity_new_unitsquare ();
  }
  p4est = p4est_new (mpi->mpicomm, connectivity, 0, NULL);
  quadrant_counts = P4EST_ALLOC (p4est_locidx_t, p4est->mpisize);

  /* time refine */
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  p4est_refine (p4est, refine_fractal, NULL);
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  elapsed_refine = start + MPI_Wtime ();
  if (refine_level <= 8) {
    p4est_vtk_write_file (p4est, "mesh_timings_refined");
  }
  count_refined = p4est->global_num_quadrants;

  /* time balance */
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  p4est_balance (p4est, NULL);
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  elapsed_balance = start + MPI_Wtime ();
  if (refine_level <= 8) {
    p4est_vtk_write_file (p4est, "mesh_timings_balanced");
  }
  count_balanced = p4est->global_num_quadrants;
  crc = p4est_checksum (p4est);

  /* time rebalance - is a noop on the tree */
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  p4est_balance (p4est, NULL);
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  elapsed_rebalance = start + MPI_Wtime ();
  P4EST_ASSERT (count_balanced == p4est->global_num_quadrants);
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* time a uniform partition */
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  p4est_partition (p4est, NULL);
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  elapsed_partition = start + MPI_Wtime ();
  if (refine_level <= 8) {
    p4est_vtk_write_file (p4est, "mesh_timings_partitioned");
  }
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* time a partition with a shift of all elements by one processor */
  for (i = 0, next_quadrant = 0; i < p4est->mpisize; ++i) {
    prev_quadrant = next_quadrant;
    next_quadrant = (p4est->global_num_quadrants * (i + 1)) / p4est->mpisize;
    quadrant_counts[i] = (p4est_locidx_t) (next_quadrant - prev_quadrant);
  }
  if (p4est->mpisize > 1) {
    quadrant_counts[0] += quadrant_counts[p4est->mpisize - 1];
    quadrant_counts[p4est->mpisize - 1] = 0;
  }
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  global_shipped = p4est_partition_given (p4est, quadrant_counts);
  P4EST_GLOBAL_PRODUCTIONF
    ("Done p4est_partition_given shipped %lld quadrants %.3g%%\n",
     (long long) global_shipped,
     global_shipped * 100. / p4est->global_num_quadrants);
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  elapsed_repartition = start + MPI_Wtime ();
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* print checksum and timings */
  P4EST_GLOBAL_STATISTICSF ("Processors %d level %d shift %d"
                            " tree checksum 0x%x\n",
                            mpi->mpisize, refine_level, level_shift, crc);
  P4EST_GLOBAL_STATISTICSF ("Level %d refined to %lld balanced to %lld\n",
                            refine_level, (long long) count_refined,
                            (long long) count_balanced);
  P4EST_GLOBAL_STATISTICSF ("Time for refinement %.3gs\n", elapsed_refine);
  P4EST_GLOBAL_STATISTICSF ("Time for balance %.3gs rebalance %.3gs\n",
                            elapsed_balance, elapsed_rebalance);
  P4EST_GLOBAL_STATISTICSF ("Time for partition %.3gs repartition %.3gs\n",
                            elapsed_partition, elapsed_repartition);

  /* destroy the p4est and its connectivity structure */
  P4EST_FREE (quadrant_counts);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  P4EST_CHECK_MPI (mpiret);

  return 0;
}

/* EOF timings.c */
