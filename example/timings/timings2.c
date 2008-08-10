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
 *        o periodic  Refinement on the unit cube with periodic b.c.
 *        o twocubes  Refinement on a forest with two trees.
 *        o rotcubes  Refinement on a forest with six rotated trees.
 */

#ifndef P4_TO_P8
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_mesh.h>
#include <p4est_vtk.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_ghost.h>
#include <p8est_mesh.h>
#include <p8est_trilinear.h>
#include <p8est_vtk.h>
#endif

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
  P4EST_CONFIG_TWOCUBES,
  P4EST_CONFIG_ROTCUBES,
#endif
}
timings_config_t;

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
  { P4EST_CONFIG_PERIODIC, 2, 6, 0x372f7402U },
  { P4EST_CONFIG_PERIODIC, 7, 6, 0xa2f1ee48U },
  { P4EST_CONFIG_TWOCUBES, 5, 6, 0xa8b1f54eU },
  { P4EST_CONFIG_TWOCUBES, 8, 5, 0x98d3579dU },
  { P4EST_CONFIG_ROTCUBES, 1, 5, 0x404e4aa8U },
  { P4EST_CONFIG_ROTCUBES, 7, 6, 0x4c381706U },
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

static void
abort_fn (void *data)
{
  mpi_context_t      *mpi = data;

  fprintf (stderr, "[%d] " P4EST_STRING "_timings abort handler\n",
           mpi->mpirank);

  /* Don't check the return value */
  MPI_Abort (mpi->mpicomm, 1);
}

int
main (int argc, char **argv)
{
  int                 i;
  int                 mpiret;
  int                 wrongusage;
  int                *ghost_owner;
  unsigned            crc;
  double              start, elapsed_refine;
  double              elapsed_balance, elapsed_rebalance;
  double              elapsed_partition, elapsed_repartition;
  double              elapsed_ghosts, elapsed_nodes;
  const char         *config_name, *usage, *errmsg;
  p4est_locidx_t     *quadrant_counts;
  p4est_gloidx_t      count_refined, count_balanced;
  p4est_gloidx_t      prev_quadrant, next_quadrant;
  p4est_gloidx_t      global_shipped;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  p4est_nodes_t      *nodes;
  sc_array_t          ghost_layer;
#ifdef P4_TO_P8
  double              elapsed_trilinear;
  trilinear_mesh_t   *mesh;
#endif
  const timings_regression_t *r;
  timings_config_t    config;
  mpi_context_t       mpi_context, *mpi = &mpi_context;

  /* initialize MPI and p4est internals */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpi->mpirank, abort_fn, mpi, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  usage =
    "Arguments: <configuration> <level>\n   Configuration can be any of\n"
#ifndef P4_TO_P8
    "      unit|periodic|three|moebius|star\n"
#else
    "      unit|periodic|twocubes|rotcubes\n"
#endif
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
    else if (!strcmp (config_name, "twocubes")) {
      config = P4EST_CONFIG_TWOCUBES;
    }
    else if (!strcmp (config_name, "rotcubes")) {
      config = P4EST_CONFIG_ROTCUBES;
    }
#endif
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
    SC_CHECK_MPI (mpiret);
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
#ifndef P4_TO_P8
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
#else
  if (config == P4EST_CONFIG_TWOCUBES) {
    connectivity = p8est_connectivity_new_twocubes ();
  }
  else if (config == P4EST_CONFIG_ROTCUBES) {
    connectivity = p8est_connectivity_new_rotcubes ();
  }
  else if (config == P4EST_CONFIG_PERIODIC) {
    connectivity = p8est_connectivity_new_periodic ();
  }
  else {
    connectivity = p8est_connectivity_new_unitcube ();
  }
#endif
  p4est = p4est_new (mpi->mpicomm, connectivity, 15, 0, NULL, NULL);
  quadrant_counts = P4EST_ALLOC (p4est_locidx_t, p4est->mpisize);

  /* time refine */
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  p4est_refine (p4est, true, refine_fractal, NULL);
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  elapsed_refine = start + MPI_Wtime ();
#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est, "timings_refined");
#endif
  count_refined = p4est->global_num_quadrants;

  /* time balance */
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  p4est_balance (p4est, NULL);
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  elapsed_balance = start + MPI_Wtime ();
#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est, "timings_balanced");
#endif
  count_balanced = p4est->global_num_quadrants;
  crc = p4est_checksum (p4est);

  /* time rebalance - is a noop on the tree */
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  p4est_balance (p4est, NULL);
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  elapsed_rebalance = start + MPI_Wtime ();
  P4EST_ASSERT (count_balanced == p4est->global_num_quadrants);
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* time a uniform partition */
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  p4est_partition (p4est, NULL);
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  elapsed_partition = start + MPI_Wtime ();
#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est, "timings_partitioned");
#endif
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* time building the ghost layer */
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  sc_array_init (&ghost_layer, sizeof (p4est_quadrant_t));
  p4est_build_ghost_layer (p4est, true, &ghost_layer, &ghost_owner);
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  elapsed_ghosts = start + MPI_Wtime ();
  P4EST_FREE (ghost_owner);

  /* time the node numbering */
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  nodes = p4est_nodes_new (p4est, &ghost_layer);
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  elapsed_nodes = start + MPI_Wtime ();

#ifdef P4_TO_P8
  /* time trilinear mesh extraction */
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  mesh = p8est_trilinear_mesh_new (p4est, nodes);
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  elapsed_trilinear = start + MPI_Wtime ();

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
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  global_shipped = p4est_partition_given (p4est, quadrant_counts);
  P4EST_GLOBAL_PRODUCTIONF
    ("Done " P4EST_STRING "_partition_given shipped %lld quadrants %.3g%%\n",
     (long long) global_shipped,
     global_shipped * 100. / p4est->global_num_quadrants);
  mpiret = MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  elapsed_repartition = start + MPI_Wtime ();
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
#ifndef P4_TO_P8
  P4EST_GLOBAL_STATISTICSF ("Time for ghosts %.3gs nodes %.3gs\n",
                            elapsed_ghosts, elapsed_nodes);
#else
  P4EST_GLOBAL_STATISTICSF
    ("Time for ghosts %.3gs nodes %.3gs trilinear %.3g\n",
     elapsed_ghosts, elapsed_nodes, elapsed_trilinear);
#endif

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

/* EOF timings2.c */
