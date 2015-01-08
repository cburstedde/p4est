/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2014 The University of Texas System
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

#include <p4est_bits.h>
#include <p6est.h>
#include <p6est_extended.h>
#include <p6est_ghost.h>
#include <p6est_vtk.h>
#include <p6est_lnodes.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

char                test_data = 'x';
char               *TEST_USER_POINTER = &test_data;

static int          refine_level = -1;
static int          refine_zlevel = -1;

/* To define a p6est_refine_column_t, all we have to do is take a p4est refine
 * function ... */
static int
p4est_refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * quadrant)
{
  int                 cid;

  if (quadrant->level >= refine_level) {
    return 0;
  }
  if (which_tree == 2 || which_tree == 3) {
    return 0;
  }

  cid = p4est_quadrant_child_id (quadrant);

  if (cid == P4EST_CHILDREN - 1 ||
      (quadrant->x >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2) &&
       quadrant->y >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#ifdef P4_TO_P8
       && quadrant->z >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#endif
      )) {
    return 1;
  }
  if (quadrant->level == 1 && cid == 2) {
    return 1;
  }
  if (quadrant->x == P4EST_QUADRANT_LEN (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->y >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

/* and wrap it.*/
static int
refine_column_fn (p6est_t * p6est, p4est_topidx_t which_tree,
                  p4est_quadrant_t * column)
{
  return p4est_refine_fn (p6est->columns, which_tree, column);
}

static int
refine_layer_fn (p6est_t * p6est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * column, p2est_quadrant_t * layer)
{
  p4est_topidx_t      tohash[4];
  unsigned            hash;

  tohash[0] = (p4est_topidx_t) column->x;
  tohash[1] = (p4est_topidx_t) column->y;
  tohash[2] = (p4est_topidx_t) layer->z;
  tohash[3] = (((p4est_topidx_t) column->level) << 16) |
    ((p4est_topidx_t) layer->level);

  hash = p4est_topidx_hash4 (tohash);

  return (layer->level < refine_zlevel && !((int) hash % 3));
}

void
init_fn (p6est_t * p6est, p4est_topidx_t which_tree,
         p4est_quadrant_t * col, p2est_quadrant_t * layer)
{
  SC_CHECK_ABORT (p6est->user_pointer == TEST_USER_POINTER,
                  "user_pointer corruption\n");
}

static int
coarsen_column_fn (p6est_t * p6est, p4est_topidx_t which_tree,
                   p4est_quadrant_t * column[])
{
  return 1;
}

static int
weight_fn (p6est_t * p6est, p4est_topidx_t which_tree,
           p4est_quadrant_t * col, p2est_quadrant_t * layer)
{
  return 1;
}

static int
coarsen_layer_fn (p6est_t * p6est, p4est_topidx_t which_tree,
                  p4est_quadrant_t * column, p2est_quadrant_t * layers[])
{
  return 1;
}

enum
{
  TIMINGS_CONNECTIVITY,
  TIMINGS_NEW,
  TIMINGS_NEW_EXT,
  TIMINGS_REFINE_COLUMNS_A,
  TIMINGS_REFINE_COLUMNS_B,
  TIMINGS_REFINE_LAYERS,
  TIMINGS_COARSEN_COLUMNS,
  TIMINGS_COARSEN_LAYERS,
  TIMINGS_GHOST_FACE,
  TIMINGS_GHOST_FULL,
  TIMINGS_GHOST_EXPAND_1,
  TIMINGS_GHOST_EXPAND_2,
  TIMINGS_BALANCE_FACE,
  TIMINGS_BALANCE_EDGE,
  TIMINGS_BALANCE_FULL,
  TIMINGS_PARTITION,
  TIMINGS_PARTITION_SAME,
  TIMINGS_LNODES_1,
  TIMINGS_LNODES_2,
  TIMINGS_LNODES_3,
  TIMINGS_SAVE,
  TIMINGS_LOAD,
  TIMINGS_NUM_STATS
};

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm = sc_MPI_COMM_WORLD;
  p4est_connectivity_t *conn4;
  p6est_connectivity_t *conn, *copy_conn;
  p6est_t            *p6est, *copy_p6est;
  p6est_ghost_t      *ghost;
  double              height[3] = { 0., 0., 0.1 };
  int                 i;
  int                 vtk;
  unsigned            crc_computed = 0;
  sc_options_t       *opt;
  int                 first_argc;
  const char         *config_name;
  const char         *save_filename = NULL;
  sc_statinfo_t       stats[TIMINGS_NUM_STATS];
  sc_flopinfo_t       fi, snapshot;
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
#ifndef P4EST_DEBUG
  sc_set_log_defaults (NULL, NULL, SC_LP_STATISTICS);
#endif
  p4est_init (NULL, SC_LP_DEFAULT);

  opt = sc_options_new (argv[0]);

  sc_options_add_int (opt, 'l', "level", &refine_level, 1,
                      "initial refine level");
  sc_options_add_int (opt, 'z', "z-level", &refine_zlevel, 2,
                      "initial refine level");
  sc_options_add_string (opt, 'c', "configuration", &config_name, "unit",
                         "configuration: brick23|corner|cubed|disk|moebius|periodic|pillow|rotwrap|star|unit");
  sc_options_add_string (opt, 'P', "save-file", &save_filename,
                         NULL, "filename for saving");
  sc_options_add_switch (opt, 'w', "write-vtk", &vtk, "write vtk files");

  first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                 opt, argc, argv);

  if (first_argc < 0 || first_argc != argc) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
    return 1;
  }
  sc_options_print_summary (p4est_package_id, SC_LP_PRODUCTION, opt);

  /* start overall timing */
  mpiret = sc_MPI_Barrier (mpicomm);
  SC_CHECK_MPI (mpiret);
  sc_flops_start (&fi);

  conn4 = p4est_connectivity_new_byname (config_name);
  SC_CHECK_ABORTF (conn4 != NULL, "Invalid connectivity name: %s\n",
                   config_name);

  sc_flops_snap (&fi, &snapshot);
  conn = p6est_connectivity_new (conn4, NULL, height);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_CONNECTIVITY], snapshot.iwtime,
                 "Connectivity");

  p4est_connectivity_destroy (conn4);

  sc_flops_snap (&fi, &snapshot);
  p6est = p6est_new (mpicomm, conn, 4, init_fn, TEST_USER_POINTER);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_NEW], snapshot.iwtime, "New");
  p6est_destroy (p6est);

  sc_flops_snap (&fi, &snapshot);
  p6est =
    p6est_new_ext (mpicomm, conn, 0, refine_level, refine_zlevel, 3, 1, 3,
                   init_fn, TEST_USER_POINTER);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_NEW_EXT], snapshot.iwtime, "New extended");

  refine_level += 2;
  sc_flops_snap (&fi, &snapshot);
  p6est_refine_columns (p6est, 1, refine_column_fn, init_fn);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_REFINE_COLUMNS_A], snapshot.iwtime,
                 "Refine columns A");

  refine_zlevel += 2;
  sc_flops_snap (&fi, &snapshot);
  p6est_refine_layers (p6est, 1, refine_layer_fn, init_fn);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_REFINE_LAYERS], snapshot.iwtime,
                 "Refine layers");

  refine_level += 2;
  sc_flops_snap (&fi, &snapshot);
  p6est_refine_columns (p6est, 1, refine_column_fn, init_fn);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_REFINE_COLUMNS_B], snapshot.iwtime,
                 "Refine layers B");

  copy_p6est = p6est_copy (p6est, 1);
  sc_flops_snap (&fi, &snapshot);
  p6est_coarsen_columns (copy_p6est, 1, coarsen_column_fn, init_fn);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_COARSEN_COLUMNS], snapshot.iwtime,
                 "Coarsen columns");
  if (vtk) {
    p6est_vtk_write_file (copy_p6est, "p6est_test_coarsen_columns");
  }
  p6est_destroy (copy_p6est);

  copy_p6est = p6est_copy (p6est, 1);
  sc_flops_snap (&fi, &snapshot);
  p6est_coarsen_layers (copy_p6est, 0, coarsen_layer_fn, init_fn);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_COARSEN_LAYERS], snapshot.iwtime,
                 "Coarsen layers");
  if (vtk) {
    p6est_vtk_write_file (copy_p6est, "p6est_test_coarsen_layers");
  }
  p6est_destroy (copy_p6est);

  if (vtk) {
    p6est_vtk_write_file (p6est, "p6est_test_pre_balance");
  }

  sc_flops_snap (&fi, &snapshot);
  ghost = p6est_ghost_new (p6est, P4EST_CONNECT_FACE);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_GHOST_FACE], snapshot.iwtime, "Ghost face");
  p6est_ghost_destroy (ghost);

  sc_flops_snap (&fi, &snapshot);
  ghost = p6est_ghost_new (p6est, P4EST_CONNECT_FULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_GHOST_FULL], snapshot.iwtime, "Ghost full");

  sc_flops_snap (&fi, &snapshot);
  p6est_ghost_expand (p6est, ghost);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_GHOST_EXPAND_1], snapshot.iwtime,
                 "Ghost expand 1");

  sc_flops_snap (&fi, &snapshot);
  p6est_ghost_expand (p6est, ghost);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_GHOST_EXPAND_2], snapshot.iwtime,
                 "Ghost expand 2");

  p6est_ghost_destroy (ghost);

  sc_flops_snap (&fi, &snapshot);
  p6est_balance (p6est, P8EST_CONNECT_FACE, init_fn);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_BALANCE_FACE], snapshot.iwtime,
                 "Balance face");

  if (vtk) {
    p6est_vtk_write_file (p6est, "p6est_test_balance_face");
  }

  sc_flops_snap (&fi, &snapshot);
  p6est_balance (p6est, P8EST_CONNECT_EDGE, init_fn);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_BALANCE_EDGE], snapshot.iwtime,
                 "Balance edge");

  if (vtk) {
    p6est_vtk_write_file (p6est, "p6est_test_balance_edge");
  }

  sc_flops_snap (&fi, &snapshot);
  p6est_balance (p6est, P8EST_CONNECT_FULL, init_fn);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_BALANCE_FULL], snapshot.iwtime,
                 "Balance full");

  if (vtk) {
    p6est_vtk_write_file (p6est, "p6est_test_balance_full");
  }

  sc_flops_snap (&fi, &snapshot);
  p6est_partition (p6est, weight_fn);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_PARTITION], snapshot.iwtime, "Partition");

  sc_flops_snap (&fi, &snapshot);
  p6est_partition (p6est, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIMINGS_PARTITION_SAME], snapshot.iwtime,
                 "Partition same");

  if (vtk) {
    p6est_vtk_write_file (p6est, "p6est_test_partition");
  }

  for (i = 1; i <= 3; i++) {
    p6est_lnodes_t     *lnodes;

    sc_flops_snap (&fi, &snapshot);
    lnodes = p6est_lnodes_new (p6est, NULL, i);
    sc_flops_shot (&fi, &snapshot);
    switch (i) {
    case 1:
      sc_stats_set1 (&stats[TIMINGS_LNODES_1], snapshot.iwtime, "Lnodes 1");
      break;
    case 2:
      sc_stats_set1 (&stats[TIMINGS_LNODES_2], snapshot.iwtime, "Lnodes 2");
      break;
    case 3:
      sc_stats_set1 (&stats[TIMINGS_LNODES_3], snapshot.iwtime, "Lnodes 3");
      break;
    }

    p6est_lnodes_destroy (lnodes);
  }

#ifdef P4EST_HAVE_ZLIB
  crc_computed = p6est_checksum (p6est);

  P4EST_GLOBAL_PRODUCTIONF ("p6est checksum 0x%08x\n", crc_computed);
#endif

  if (save_filename) {
    sc_flops_snap (&fi, &snapshot);
    p6est_save (save_filename, p6est, 1);
    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[TIMINGS_SAVE], snapshot.iwtime, "Save");

    sc_flops_snap (&fi, &snapshot);
    copy_p6est = p6est_load (save_filename, p6est->mpicomm,
                             p6est->data_size, 1, p6est->user_pointer,
                             &copy_conn);
    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[TIMINGS_LOAD], snapshot.iwtime, "Load");

    p6est_destroy (copy_p6est);

    p6est_connectivity_destroy (copy_conn);
  }
  else {
    sc_stats_set1 (&stats[TIMINGS_SAVE], 0., "Save");
    sc_stats_set1 (&stats[TIMINGS_LOAD], 0., "Load");
  }

  p6est_destroy (p6est);

  p6est_connectivity_destroy (conn);

  /* calculate and print timings */
  sc_stats_compute (mpicomm, TIMINGS_NUM_STATS, stats);
  sc_stats_print (p4est_package_id, SC_LP_STATISTICS,
                  TIMINGS_NUM_STATS, stats, 1, 1);

  sc_options_destroy (opt);

  /* exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
