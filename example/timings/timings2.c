/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
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
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_nodes.h>
#include <p4est_vtk.h>
#include <p4est_lnodes.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_nodes.h>
#include <p8est_vtk.h>
#include <p8est_lnodes.h>
#endif
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

/* #define P4EST_TIMINGS_VTK */

typedef enum
{
  P4EST_CONFIG_NULL,
  P4EST_CONFIG_UNIT,
  P4EST_CONFIG_PERIODIC,
#ifndef P4_TO_P8
  P4EST_CONFIG_THREE,
  P4EST_CONFIG_MOEBIUS,
  P4EST_CONFIG_STAR
#else
  P4EST_CONFIG_ROTWRAP,
  P4EST_CONFIG_TWOCUBES,
  P4EST_CONFIG_ROTCUBES,
  P4EST_CONFIG_SHELL
#endif
}
timings_config_t;

enum
{
  TIMINGS_REFINE,
  TIMINGS_BALANCE,
  TIMINGS_BALANCE_A,
  TIMINGS_BALANCE_COMM,
  TIMINGS_BALANCE_B,
  TIMINGS_BALANCE_A_COUNT_IN,
  TIMINGS_BALANCE_A_COUNT_OUT,
  TIMINGS_BALANCE_COMM_SENT,
  TIMINGS_BALANCE_COMM_NZPEERS,
  TIMINGS_BALANCE_B_COUNT_IN,
  TIMINGS_BALANCE_B_COUNT_OUT,
  TIMINGS_BALANCE_RANGES,
  TIMINGS_BALANCE_NOTIFY,
  TIMINGS_BALANCE_NOTIFY_ALLGATHER,
  TIMINGS_BALANCE_A_ZERO_SENDS,
  TIMINGS_BALANCE_A_ZERO_RECEIVES,
  TIMINGS_BALANCE_B_ZERO_SENDS,
  TIMINGS_BALANCE_B_ZERO_RECEIVES,
  TIMINGS_REBALANCE,
  TIMINGS_REBALANCE_A,
  TIMINGS_REBALANCE_COMM,
  TIMINGS_REBALANCE_B,
  TIMINGS_REBALANCE_A_COUNT_IN,
  TIMINGS_REBALANCE_A_COUNT_OUT,
  TIMINGS_REBALANCE_COMM_SENT,
  TIMINGS_REBALANCE_COMM_NZPEERS,
  TIMINGS_REBALANCE_B_COUNT_IN,
  TIMINGS_REBALANCE_B_COUNT_OUT,
  TIMINGS_PARTITION,
  TIMINGS_GHOSTS,
  TIMINGS_NODES,
  TIMINGS_TRILINEAR_OBSOLETE,
  TIMINGS_REPARTITION,
  TIMINGS_LNODES,
  TIMINGS_LNODES3,
  TIMINGS_LNODES7,
  TIMINGS_NUM_STATS
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
  sc_MPI_Comm         mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

static int          refine_level = 0;
static int          level_shift = 0;

/* *INDENT-OFF* */
static const timings_regression_t regression_oldschool[] =
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

static const timings_regression_t regression_latest[] =
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
  { P4EST_CONFIG_UNIT, 1, 5, 0xc2012e84U },
  { P4EST_CONFIG_UNIT, 1, 6, 0x2cad814dU },
  { P4EST_CONFIG_UNIT, 3, 8, 0xeb252238U },
  { P4EST_CONFIG_PERIODIC, 1, 5, 0x2776c9b7U },
  { P4EST_CONFIG_PERIODIC, 2, 5, 0x2776c9b7U },
  { P4EST_CONFIG_PERIODIC, 7, 6, 0x4f281079U },
  { P4EST_CONFIG_ROTWRAP, 2, 6, 0x372f7402U },
  { P4EST_CONFIG_ROTWRAP, 7, 6, 0x372f7402U },
  { P4EST_CONFIG_TWOCUBES, 5, 6, 0xa8b1f54eU },
  { P4EST_CONFIG_TWOCUBES, 8, 5, 0x0aad11d0U },
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
  unsigned            crc, gcrc;
  const char         *config_name;
  const char         *load_name;
  p4est_locidx_t     *quadrant_counts;
  p4est_gloidx_t      count_refined, count_balanced;
  p4est_gloidx_t      prev_quadrant, next_quadrant;
  p4est_gloidx_t      global_shipped;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  p4est_nodes_t      *nodes = NULL;
  p4est_ghost_t      *ghost;
  p4est_lnodes_t     *lnodes;
  const timings_regression_t *r, *regression;
  timings_config_t    config;
  sc_statinfo_t       stats[TIMINGS_NUM_STATS];
  sc_flopinfo_t       fi, snapshot;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  sc_options_t       *opt;
  int                 overlap;
  int                 subtree;
  int                 borders;
  int                 max_ranges;
  int                 use_ranges, use_ranges_notify, use_balance_verify;
  int                 oldschool, generate;
  int                 first_argc;
  int                 test_multiple_orders;
  int                 skip_nodes, skip_lnodes;

  /* initialize MPI and p4est internals */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpi->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
#ifndef P4EST_ENABLE_DEBUG
  sc_set_log_defaults (NULL, NULL, SC_LP_STATISTICS);
#endif
  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  P4EST_GLOBAL_PRODUCTIONF ("Size of %dtant: %lld bytes\n", P4EST_DIM,
                            (long long) sizeof (p4est_quadrant_t));

  opt = sc_options_new (argv[0]);

  sc_options_add_switch (opt, 'o', "new-balance-overlap", &overlap,
                         "use the new balance overlap algorithm");
  sc_options_add_switch (opt, 's', "new-balance-subtree", &subtree,
                         "use the new balance subtree algorithm");
  sc_options_add_switch (opt, 'b', "new-balance-borders", &borders,
                         "use borders in balance");
  sc_options_add_int (opt, 'm', "max-ranges", &max_ranges, -1,
                      "override p4est_num_ranges");
  sc_options_add_switch (opt, 'r', "ranges", &use_ranges,
                         "use ranges in balance");
  sc_options_add_switch (opt, 't', "ranges-notify", &use_ranges_notify,
                         "use both ranges and notify");
  sc_options_add_switch (opt, 'y', "balance-verify", &use_balance_verify,
                         "use verifications in balance");
  sc_options_add_int (opt, 'l', "level", &refine_level, 0,
                      "initial refine level");
#ifndef P4_TO_P8
  sc_options_add_string (opt, 'c', "configuration", &config_name, "unit",
                         "configuration: unit|periodic|three|moebius|star");
#else
  sc_options_add_string (opt, 'c', "configuration", &config_name, "unit",
                         "configuration: unit|periodic|rotwrap|twocubes|rotcubes|shell");
#endif
  sc_options_add_bool (opt, 0, "oldschool", &oldschool, 0,
                       "Use original p4est_new call");
  sc_options_add_switch (opt, 0, "generate", &generate,
                         "Generate regression commands");
  sc_options_add_string (opt, 'f', "load-forest", &load_name, NULL,
                         "load saved " P4EST_STRING);
  sc_options_add_switch (opt, 0, "multiple-lnodes",
                         &test_multiple_orders,
                         "Also time lnodes for orders 2, 4, and 8");
  sc_options_add_switch (opt, 0, "skip-nodes", &skip_nodes, "Skip nodes");
  sc_options_add_switch (opt, 0, "skip-lnodes", &skip_lnodes, "Skip lnodes");

  first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                 opt, argc, argv);
  if (first_argc < 0 || first_argc != argc) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
    return 1;
  }
  if (max_ranges < -1 || max_ranges == 0) {
    P4EST_GLOBAL_LERROR ("The -m / --max-ranges option must be positive\n");
    return 1;
  }
  sc_options_print_summary (p4est_package_id, SC_LP_PRODUCTION, opt);

  if (skip_lnodes) {
    if (test_multiple_orders) {
      SC_GLOBAL_PRODUCTION
        ("Warning: cannot test multiple lnode orders if --skip-lnodes is given.\n");
      test_multiple_orders = 0;
    }
  }
  wrongusage = 0;
  config = P4EST_CONFIG_NULL;
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
  else if (load_name != NULL) {
    config_name = load_name;
  }
  else {
    wrongusage = 1;
  }
  if (wrongusage) {
    P4EST_GLOBAL_LERRORF ("Wrong configuration name given: %s\n",
                          config_name);
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
    sc_abort_collective ("Usage error");
  }

  /* get command line argument: maximum refinement level */
  level_shift = 4;

  /* print general setup information */
  P4EST_GLOBAL_STATISTICSF
    ("Processors %d configuration %s level %d shift %d\n", mpi->mpisize,
     config_name, refine_level, level_shift);

  /* start overall timing */
  mpiret = sc_MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  sc_flops_start (&fi);

  /* create connectivity and forest structures */
  regression = NULL;
  if (load_name == NULL) {
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

    /* create new p4est from scratch */
    if (oldschool) {
      regression = regression_oldschool;
      p4est = p4est_new_ext (mpi->mpicomm, connectivity,
                             15, 0, 0, 0, NULL, NULL);
    }
    else {
      regression = regression_latest;
      p4est = p4est_new_ext (mpi->mpicomm, connectivity,
                             0, refine_level - level_shift, 1, 0, NULL, NULL);
    }

    /* print all available regression tests */
    if (generate) {
      const char         *config_name = NULL;

      P4EST_GLOBAL_PRODUCTION ("Checksum regression tests available:\n");
      for (r = regression; r->config != P4EST_CONFIG_NULL; ++r) {
        switch (r->config) {
        case P4EST_CONFIG_PERIODIC:
          config_name = "periodic";
          break;
#ifndef P4_TO_P8
        case P4EST_CONFIG_THREE:
          config_name = "three";
          break;
        case P4EST_CONFIG_MOEBIUS:
          config_name = "moebius";
          break;
        case P4EST_CONFIG_STAR:
          config_name = "star";
          break;
#else
        case P4EST_CONFIG_ROTWRAP:
          config_name = "rotwrap";
          break;
        case P4EST_CONFIG_TWOCUBES:
          config_name = "twocubes";
          break;
        case P4EST_CONFIG_ROTCUBES:
          config_name = "rotcubes";
          break;
        case P4EST_CONFIG_SHELL:
          config_name = "shell";
          break;
#endif
        default:
          config_name = "unit";
          break;
        }
        P4EST_GLOBAL_PRODUCTIONF ("mpirun -np %3d %s%s -c %10s -l %2d\n",
                                  r->mpisize, opt->program_path,
                                  oldschool ? " --oldschool" : "",
                                  config_name, r->level);
      }
    }
  }
  else {
    p4est = p4est_load (load_name, mpi->mpicomm, 0, 0, NULL, &connectivity);
  }

  p4est->inspect = P4EST_ALLOC_ZERO (p4est_inspect_t, 1);
  p4est->inspect->use_balance_ranges = use_ranges;
  p4est->inspect->use_balance_ranges_notify = use_ranges_notify;
  p4est->inspect->use_balance_verify = use_balance_verify;
  p4est->inspect->balance_max_ranges = max_ranges;
  P4EST_GLOBAL_STATISTICSF
    ("Balance: new overlap %d new subtree %d borders %d\n", overlap,
     (overlap && subtree), (overlap && borders));
  quadrant_counts = P4EST_ALLOC (p4est_locidx_t, p4est->mpisize);

  /* time refine */
  sc_flops_snap (&fi, &snapshot);
  p4est_refine (p4est, 1, refine_fractal, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_REFINE], snapshot.iwtime, "Refine");
#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est, "timings_refined");
#endif
  count_refined = p4est->global_num_quadrants;

  /* time balance */
  sc_flops_snap (&fi, &snapshot);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_BALANCE], snapshot.iwtime, "Balance");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_A],
                 p4est->inspect->balance_A, "Balance A time");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_COMM],
                 p4est->inspect->balance_comm, "Balance comm time");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_B],
                 p4est->inspect->balance_B, "Balance B time");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_A_COUNT_IN],
                 (double) p4est->inspect->balance_A_count_in,
                 "Balance A count inlist");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_A_COUNT_OUT],
                 (double) p4est->inspect->balance_A_count_out,
                 "Balance A count outlist");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_COMM_SENT],
                 (double) p4est->inspect->balance_comm_sent,
                 "Balance sent second round");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_COMM_NZPEERS],
                 (double) p4est->inspect->balance_comm_nzpeers,
                 "Balance nonzero peers second round");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_B_COUNT_IN],
                 (double) p4est->inspect->balance_B_count_in,
                 "Balance B count inlist");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_B_COUNT_OUT],
                 (double) p4est->inspect->balance_B_count_out,
                 "Balance B count outlist");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_RANGES],
                 p4est->inspect->balance_ranges, "Balance time for ranges");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_NOTIFY],
                 p4est->inspect->balance_notify, "Balance time for notify");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_NOTIFY_ALLGATHER],
                 p4est->inspect->balance_notify_allgather,
                 "Balance time for notify_allgather");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_A_ZERO_RECEIVES],
                 p4est->inspect->balance_zero_receives[0],
                 "Balance A zero receives");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_A_ZERO_SENDS],
                 p4est->inspect->balance_zero_sends[0],
                 "Balance A zero sends");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_B_ZERO_RECEIVES],
                 p4est->inspect->balance_zero_receives[1],
                 "Balance B zero receives");
  sc_stats_set1 (&stats[TIMINGS_BALANCE_B_ZERO_SENDS],
                 p4est->inspect->balance_zero_sends[1],
                 "Balance B zero sends");
#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est, "timings_balanced");
#endif
  count_balanced = p4est->global_num_quadrants;
  crc = p4est_checksum (p4est);

  /* time rebalance - is a noop on the tree */
  sc_flops_snap (&fi, &snapshot);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_REBALANCE], snapshot.iwtime, "Rebalance");
  sc_stats_set1 (&stats[TIMINGS_REBALANCE_A],
                 p4est->inspect->balance_A, "Rebalance A time");
  sc_stats_set1 (&stats[TIMINGS_REBALANCE_COMM],
                 p4est->inspect->balance_comm, "Rebalance comm time");
  sc_stats_set1 (&stats[TIMINGS_REBALANCE_B],
                 p4est->inspect->balance_B, "Rebalance B time");
  sc_stats_set1 (&stats[TIMINGS_REBALANCE_A_COUNT_IN],
                 (double) p4est->inspect->balance_A_count_in,
                 "Rebalance A count inlist");
  sc_stats_set1 (&stats[TIMINGS_REBALANCE_A_COUNT_OUT],
                 (double) p4est->inspect->balance_A_count_out,
                 "Rebalance A count outlist");
  sc_stats_set1 (&stats[TIMINGS_REBALANCE_COMM_SENT],
                 (double) p4est->inspect->balance_comm_sent,
                 "Rebalance sent second round");
  sc_stats_set1 (&stats[TIMINGS_REBALANCE_COMM_NZPEERS],
                 (double) p4est->inspect->balance_comm_nzpeers,
                 "Rebalance nonzero peers second round");
  sc_stats_set1 (&stats[TIMINGS_REBALANCE_B_COUNT_IN],
                 (double) p4est->inspect->balance_B_count_in,
                 "Rebalance B count inlist");
  sc_stats_set1 (&stats[TIMINGS_REBALANCE_B_COUNT_OUT],
                 (double) p4est->inspect->balance_B_count_out,
                 "Rebalance B count outlist");
  P4EST_ASSERT (count_balanced == p4est->global_num_quadrants);
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* time a uniform partition */
  sc_flops_snap (&fi, &snapshot);
  p4est_partition (p4est, 0, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_PARTITION], snapshot.iwtime, "Partition");
#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est, "timings_partitioned");
#endif
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* time building the ghost layer */
  sc_flops_snap (&fi, &snapshot);
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_GHOSTS], snapshot.iwtime, "Ghost layer");
  gcrc = p4est_ghost_checksum (p4est, ghost);

  /* time the node numbering */
  if (!skip_nodes) {
    sc_flops_snap (&fi, &snapshot);
    nodes = p4est_nodes_new (p4est, ghost);
    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[TIMINGS_NODES], snapshot.iwtime, "Nodes");
  }
  else {
    sc_stats_set1 (&stats[TIMINGS_NODES], 0., "Nodes");
  }

  /* set this anyway so the output format is dimension independent */
  sc_stats_set1 (&stats[TIMINGS_TRILINEAR_OBSOLETE], 0., "Unused");

  if (!skip_nodes) {
    p4est_nodes_destroy (nodes);
  }

  /* time the lnode numbering */
  if (!skip_lnodes) {
    sc_flops_snap (&fi, &snapshot);
    lnodes = p4est_lnodes_new (p4est, ghost, 1);
    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[TIMINGS_LNODES], snapshot.iwtime, "L-Nodes");
    p4est_lnodes_destroy (lnodes);
  }
  else {
    sc_stats_set1 (&stats[TIMINGS_LNODES], 0., "L-Nodes");
  }

  if (test_multiple_orders) {
    sc_flops_snap (&fi, &snapshot);
    lnodes = p4est_lnodes_new (p4est, ghost, 3);
    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[TIMINGS_LNODES3], snapshot.iwtime, "L-Nodes 3");
    p4est_lnodes_destroy (lnodes);

    sc_flops_snap (&fi, &snapshot);
    lnodes = p4est_lnodes_new (p4est, ghost, 7);
    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[TIMINGS_LNODES7], snapshot.iwtime, "L-Nodes 7");
    p4est_lnodes_destroy (lnodes);
  }
  else {
    sc_stats_set1 (&stats[TIMINGS_LNODES3], 0., "L-Nodes 3");
    sc_stats_set1 (&stats[TIMINGS_LNODES7], 0., "L-Nodes 7");
  }

  p4est_ghost_destroy (ghost);

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
  if (regression != NULL && mpi->mpirank == 0) {
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
                            " checksums 0x%08x 0x%08x\n",
                            mpi->mpisize, refine_level, level_shift,
                            crc, gcrc);
  P4EST_GLOBAL_STATISTICSF ("Level %d refined to %lld balanced to %lld\n",
                            refine_level, (long long) count_refined,
                            (long long) count_balanced);

  /* calculate and print timings */
  sc_stats_compute (mpi->mpicomm, TIMINGS_NUM_STATS, stats);
  sc_stats_print (p4est_package_id, SC_LP_STATISTICS,
                  TIMINGS_NUM_STATS, stats, 1, 1);

  /* destroy the p4est and its connectivity structure */
  P4EST_FREE (quadrant_counts);
  P4EST_FREE (p4est->inspect);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  sc_options_destroy (opt);

  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
