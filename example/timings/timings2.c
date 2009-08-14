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
 *        o periodic  Refinement on the unit square with periodic b.c.
 *        o three     Refinement on a forest with three trees.
 *        o moebius   Refinement on a 5-tree Moebius band.
 *        o star      Refinement on a 6-tree star shaped domain.
 *
 * Usage: p8est_timings <configuration> <level>
 *        possible configurations:
 *        o unit      Refinement on the unit cube.
 *        o periodic  Refinement on the unit cube with all-periodic b.c.
 *        o rotwrap   Refinement on the unit cube with weird periodic b.c.
 *        o twocubes  Refinement on a forest with two trees.
 *        o rotcubes  Refinement on a forest with six rotated trees.
 *        o shell     Refinement on a 24-tree spherical shell.
 */

#ifndef P4_TO_P8
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_nodes.h>
#include <p4est_vtk.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_ghost.h>
#include <p8est_nodes.h>
#include <p8est_trilinear.h>
#include <p8est_vtk.h>
#endif
#include <sc_flops.h>
#include <sc_statistics.h>

/* #define P4EST_TIMINGS_VTK */

typedef enum
{
  P4EST_CONFIG_NULL,
  P4EST_CONFIG_UNIT,
  P4EST_CONFIG_PERIODIC,
#ifndef P4_TO_P8
  P4EST_CONFIG_THREE,
  P4EST_CONFIG_MOEBIUS,
  P4EST_CONFIG_STAR,
#else
  P4EST_CONFIG_ROTWRAP,
  P4EST_CONFIG_TWOCUBES,
  P4EST_CONFIG_ROTCUBES,
  P4EST_CONFIG_SHELL,
#endif
}
timings_config_t;

enum
{
  TIMINGS_REFINE,
  TIMINGS_BALANCE,
  TIMINGS_REBALANCE,
  TIMINGS_PARTITION,
  TIMINGS_GHOSTS,
  TIMINGS_NODES,
#ifdef P4_TO_P8
  TIMINGS_TRILINEAR,
#endif
  TIMINGS_REPARTITION,
  TIMINGS_NUM_STATS,
};

typedef struct
{
  timings_config_t    config;
  int                 mpisize;
  int                 level;
  unsigned            checksum;
}
timings_regression_t;

typedef struct
{
  MPI_Comm            mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

static int          refine_level = 0;
static int          level_shift = 0;

/* *INDENT-OFF* */
static const timings_regression_t regression[] =
{
#ifndef P4_TO_P8
  { P4EST_CONFIG_UNIT, 1, 10, 0x6e3e83c4U },
  { P4EST_CONFIG_UNIT, 1, 11, 0x334bc3deU },
  { P4EST_CONFIG_UNIT, 64, 14, 0xad908ce4U },
  { P4EST_CONFIG_UNIT, 256, 15, 0x9e7da646U },
  { P4EST_CONFIG_STAR, 1, 6, 0x14107b57U },
  { P4EST_CONFIG_STAR, 4, 6, 0x14107b57U },
  { P4EST_CONFIG_STAR, 52, 13, 0xc86c74d9U },
  { P4EST_CONFIG_STAR, 64, 13, 0xc86c74d9U },
#else
  { P4EST_CONFIG_UNIT, 1, 5, 0xe1ffa67bU },
  { P4EST_CONFIG_UNIT, 1, 6, 0x2cad814dU },
  { P4EST_CONFIG_UNIT, 3, 8, 0xeb252238U },
  { P4EST_CONFIG_PERIODIC, 1, 5, 0x99874fedU },
  { P4EST_CONFIG_PERIODIC, 2, 5, 0x575af6d5U },
  { P4EST_CONFIG_PERIODIC, 7, 6, 0xbc35524aU },
  { P4EST_CONFIG_ROTWRAP, 2, 6, 0x372f7402U },
  { P4EST_CONFIG_ROTWRAP, 7, 6, 0xa2f1ee48U },
  { P4EST_CONFIG_TWOCUBES, 5, 6, 0xa8b1f54eU },
  { P4EST_CONFIG_TWOCUBES, 8, 5, 0x98d3579dU },
  { P4EST_CONFIG_ROTCUBES, 1, 5, 0x404e4aa8U },
  { P4EST_CONFIG_ROTCUBES, 7, 6, 0x4c381706U },
  { P4EST_CONFIG_SHELL, 1, 4, 0x8c56f159U },
  { P4EST_CONFIG_SHELL, 3, 5, 0xafbc4f8cU },
  { P4EST_CONFIG_SHELL, 5, 6, 0xf6d9efb8U },
#endif
  { P4EST_CONFIG_NULL, 0, 0, 0 }
};
/* *INDENT-ON* */

static int
refine_fractal (p4est_t * p4est, p4est_topidx_t which_tree,
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
  return (qid == 0 || qid == 3
#ifdef P4_TO_P8
          || qid == 5 || qid == 6
#endif
    );
}

int
main (int argc, char **argv)
{
  int                 i;
  int                 mpiret;
  int                 wrongusage;
  int                *ghost_owner;
  unsigned            crc;
  const char         *config_name, *usage;
  p4est_locidx_t     *quadrant_counts;
  p4est_gloidx_t      count_refined, count_balanced;
  p4est_gloidx_t      prev_quadrant, next_quadrant;
  p4est_gloidx_t      global_shipped;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  p4est_nodes_t      *nodes;
  sc_array_t          ghost_layer;
#ifdef P4_TO_P8
  trilinear_mesh_t   *mesh;
#endif
  const timings_regression_t *r;
  timings_config_t    config;
  sc_statinfo_t       stats[TIMINGS_NUM_STATS];
  sc_flopinfo_t       fi, snapshot;
  mpi_context_t       mpi_context, *mpi = &mpi_context;

  /* initialize MPI and p4est internals */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpi->mpicomm, true, true, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  usage =
    "Arguments: <configuration> <level>\n   Configuration can be any of\n"
#ifndef P4_TO_P8
    "      unit|periodic|three|moebius|star\n"
#else
    "      unit|periodic|rotwrap|twocubes|rotcubes|shell\n"
#endif
    "   Level controls the maximum depth of refinement\n";
  wrongusage = 0;
  config = P4EST_CONFIG_NULL;
  config_name = NULL;
  if (!wrongusage && argc != 3) {
    wrongusage = 1;
  }
  if (!wrongusage) {
    config_name = argv[1];
    if (!strcmp (config_name, "unit")) {
      config = P4EST_CONFIG_UNIT;
    }
    else if (!strcmp (config_name, "periodic")) {
      config = P4EST_CONFIG_PERIODIC;
    }
#ifndef P4_TO_P8
    else if (!strcmp (config_name, "three")) {
      config = P4EST_CONFIG_THREE;
    }
    else if (!strcmp (config_name, "moebius")) {
      config = P4EST_CONFIG_MOEBIUS;
    }
    else if (!strcmp (config_name, "star")) {
      config = P4EST_CONFIG_STAR;
    }
#else
    else if (!strcmp (config_name, "rotwrap")) {
      config = P4EST_CONFIG_ROTWRAP;
    }
    else if (!strcmp (config_name, "twocubes")) {
      config = P4EST_CONFIG_TWOCUBES;
    }
    else if (!strcmp (config_name, "rotcubes")) {
      config = P4EST_CONFIG_ROTCUBES;
    }
    else if (!strcmp (config_name, "shell")) {
      config = P4EST_CONFIG_SHELL;
    }
#endif
    else {
      wrongusage = 1;
    }
  }
  if (wrongusage) {
    if (mpi->mpirank == 0) {
      fputs (usage, stderr);
      SC_ABORT ("Usage error");
    }
    mpiret = MPI_Barrier (mpi->mpicomm);
    SC_CHECK_MPI (mpiret);
  }

#ifndef P4EST_DEBUG
  sc_set_log_defaults (NULL, NULL, SC_LP_STATISTICS);
#endif

  /* get command line argument: maximum refinement level */
  refine_level = atoi (argv[2]);
  level_shift = 4;

  /* print general setup information */
  P4EST_GLOBAL_STATISTICSF
    ("Processors %d configuration %s level %d shift %d\n", mpi->mpisize,
     config_name, refine_level, level_shift);

  /* start overall timing */
  sc_flops_start (&fi);

  /* create connectivity and forest structures */
#ifndef P4_TO_P8
  if (config == P4EST_CONFIG_PERIODIC) {
    connectivity = p4est_connectivity_new_periodic ();
  }
  else if (config == P4EST_CONFIG_THREE) {
    connectivity = p4est_connectivity_new_corner ();
  }
  else if (config == P4EST_CONFIG_MOEBIUS) {
    connectivity = p4est_connectivity_new_moebius ();
  }
  else if (config == P4EST_CONFIG_STAR) {
    connectivity = p4est_connectivity_new_star ();
  }
  else {
    connectivity = p4est_connectivity_new_unitsquare ();
  }
#else
  if (config == P4EST_CONFIG_PERIODIC) {
    connectivity = p8est_connectivity_new_periodic ();
  }
  else if (config == P4EST_CONFIG_ROTWRAP) {
    connectivity = p8est_connectivity_new_rotwrap ();
  }
  else if (config == P4EST_CONFIG_TWOCUBES) {
    connectivity = p8est_connectivity_new_twocubes ();
  }
  else if (config == P4EST_CONFIG_ROTCUBES) {
    connectivity = p8est_connectivity_new_rotcubes ();
  }
  else if (config == P4EST_CONFIG_SHELL) {
    connectivity = p8est_connectivity_new_shell ();
  }
  else {
    connectivity = p8est_connectivity_new_unitcube ();
  }
#endif
  p4est = p4est_new (mpi->mpicomm, connectivity, 15, 0, NULL, NULL);
  quadrant_counts = P4EST_ALLOC (p4est_locidx_t, p4est->mpisize);

  /* time refine */
  sc_flops_snap (&fi, &snapshot);
  p4est_refine (p4est, true, refine_fractal, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_REFINE], snapshot.iwtime, "Refine");
#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est, "timings_refined");
#endif
  count_refined = p4est->global_num_quadrants;

  /* time balance */
  sc_flops_snap (&fi, &snapshot);
  p4est_balance (p4est, P4EST_BALANCE_FULL, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_BALANCE], snapshot.iwtime, "Balance");
#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est, "timings_balanced");
#endif
  count_balanced = p4est->global_num_quadrants;
  crc = p4est_checksum (p4est);

  /* time rebalance - is a noop on the tree */
  sc_flops_snap (&fi, &snapshot);
  p4est_balance (p4est, P4EST_BALANCE_FULL, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_REBALANCE], snapshot.iwtime, "Rebalance");
  P4EST_ASSERT (count_balanced == p4est->global_num_quadrants);
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* time a uniform partition */
  sc_flops_snap (&fi, &snapshot);
  p4est_partition (p4est, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_PARTITION], snapshot.iwtime, "Partition");
#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est, "timings_partitioned");
#endif
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* time building the ghost layer */
  sc_array_init (&ghost_layer, sizeof (p4est_quadrant_t));
  sc_flops_snap (&fi, &snapshot);
  p4est_build_ghost_layer (p4est, P4EST_BALANCE_FULL,
                           &ghost_layer, &ghost_owner);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_GHOSTS], snapshot.iwtime, "Ghost layer");
  P4EST_FREE (ghost_owner);

  /* time the node numbering */
  sc_flops_snap (&fi, &snapshot);
  nodes = p4est_nodes_new (p4est, &ghost_layer);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_NODES], snapshot.iwtime, "Nodes");

#ifdef P4_TO_P8
  /* time trilinear mesh extraction */
  sc_flops_snap (&fi, &snapshot);
  mesh = p8est_trilinear_mesh_new (p4est, nodes);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_TRILINEAR], snapshot.iwtime, "Trilinear");

  /* destroy mesh related memory */
  p8est_trilinear_mesh_destroy (mesh);
#endif
  p4est_nodes_destroy (nodes);
  sc_array_reset (&ghost_layer);

  /* time a partition with a shift of all elements by one processor */
  for (i = 0, next_quadrant = 0; i < p4est->mpisize; ++i) {
    prev_quadrant = next_quadrant;
    next_quadrant = (p4est->global_num_quadrants * (i + 1)) / p4est->mpisize;
    quadrant_counts[i] = (p4est_locidx_t) (next_quadrant - prev_quadrant);
  }
  if (p4est->mpisize > 1) {
    quadrant_counts[0] += quadrant_counts[p4est->mpisize - 1];  /* same type */
    quadrant_counts[p4est->mpisize - 1] = 0;
  }

  sc_flops_snap (&fi, &snapshot);
  global_shipped = p4est_partition_given (p4est, quadrant_counts);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_REPARTITION], snapshot.iwtime, "Repartition");

  P4EST_GLOBAL_PRODUCTIONF
    ("Done " P4EST_STRING "_partition_given shipped %lld quadrants %.3g%%\n",
     (long long) global_shipped,
     global_shipped * 100. / p4est->global_num_quadrants);
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* verify forest checksum */
  if (mpi->mpirank == 0) {
    for (r = regression; r->config != P4EST_CONFIG_NULL; ++r) {
      if (r->config != config || r->mpisize != mpi->mpisize
          || r->level != refine_level)
        continue;
      SC_CHECK_ABORT (crc == r->checksum, "Checksum mismatch");
      P4EST_GLOBAL_INFO ("Checksum regression OK\n");
      break;
    }
  }

  /* print status and checksum */
  P4EST_GLOBAL_STATISTICSF ("Processors %d level %d shift %d"
                            " tree checksum 0x%x\n",
                            mpi->mpisize, refine_level, level_shift, crc);
  P4EST_GLOBAL_STATISTICSF ("Level %d refined to %lld balanced to %lld\n",
                            refine_level, (long long) count_refined,
                            (long long) count_balanced);

  /* calculate and print timings */
  sc_stats_compute (mpi->mpicomm, TIMINGS_NUM_STATS, stats);
  sc_stats_print (p4est_package_id, SC_LP_STATISTICS,
                  TIMINGS_NUM_STATS, stats, true, true);

  /* destroy the p4est and its connectivity structure */
  P4EST_FREE (quadrant_counts);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
