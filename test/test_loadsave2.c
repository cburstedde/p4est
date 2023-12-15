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
#include <p4est_extended.h>
#include <p4est_io.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_io.h>
#endif
#include <sc_options.h>
#include <sc_statistics.h>

#ifndef P4_TO_P8
#define P4EST_CONN_SUFFIX "p4c"
#define P4EST_FOREST_SUFFIX "p4p"
static const int    default_refine_level = 7;
#else
#define P4EST_CONN_SUFFIX "p8c"
#define P4EST_FOREST_SUFFIX "p8p"
static const int    default_refine_level = 4;
#endif
static int          refine_level = 0;
static int          counter = 0;

static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
  int                *data = (int *) quadrant->p.user_data;

  *data = (counter = counter * 1664525 + 1013904223) + (int) which_tree;
}

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 cid;

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
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && cid == 2) {
    return 1;
  }
  if (quadrant->y == P4EST_QUADRANT_LEN (2) &&
      quadrant->x == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->y >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

enum
{
  STATS_CONN_LOAD,
  STATS_P4EST_SAVE1,
  STATS_P4EST_LOAD1a,
  STATS_P4EST_LOAD1b,
  STATS_P4EST_ELEMS,
  STATS_P4EST_SAVE2,
  STATS_P4EST_LOAD2,
  STATS_P4EST_SAVE3,
  STATS_P4EST_LOAD3,
  STATS_P4EST_LOAD4,
  STATS_COUNT
};

static void
test_deflate (p4est_t * p4est)
{
  p4est_gloidx_t     *pertree;
  p4est_t            *p4est2;
  sc_array_t         *qarr, *darr;

  pertree = P4EST_ALLOC (p4est_gloidx_t, p4est->connectivity->num_trees + 1);
  p4est_comm_count_pertree (p4est, pertree);
  darr = NULL;
  qarr = p4est_deflate_quadrants (p4est, p4est->data_size > 0 ? &darr : NULL);

  /* Data that describes the forest completely
     (a) shared data (identical on all processors):
     p4est->connectivity
     p4est->global_first_quadrant (does not need to be stored away)
     pertree
     (b) per-processor data (partition independent after allgatherv):
     qarr
     darr (if per-quadrant data size is greater 0 and it should be saved)
   */

  /* Create a forest from this information and compare */
  p4est2 = p4est_inflate (p4est->mpicomm, p4est->connectivity,
                          p4est->global_first_quadrant, pertree,
                          qarr, darr, p4est->user_pointer);
  SC_CHECK_ABORT (p4est_is_equal (p4est, p4est2, 1), "de/inflate");
  p4est_destroy (p4est2);

  /* clean up allocated memory */
  P4EST_FREE (pertree);
  sc_array_destroy (qarr);
  if (darr != NULL) {
    sc_array_destroy (darr);
  }
}

static unsigned
test_checksum (p4est_t * p4est, int have_zlib)
{
  return have_zlib ? p4est_checksum (p4est) : 0;
}

static void
test_loadsave (p4est_connectivity_t * connectivity, const char *prefix,
               sc_MPI_Comm mpicomm, int mpirank)
{
  int                 mpiret, retval;
  int                 have_zlib;
  unsigned            csum, csum2;
  double              elapsed, wtime;
  p4est_connectivity_t *conn2;
  p4est_t            *p4est, *p4est2;
  sc_statinfo_t       stats[STATS_COUNT];
  char                conn_name[BUFSIZ];
  char                p4est_name[BUFSIZ];

  snprintf (conn_name, BUFSIZ, "%s.%s", prefix, P4EST_CONN_SUFFIX);
  snprintf (p4est_name, BUFSIZ, "%s.%s", prefix, P4EST_FOREST_SUFFIX);
  P4EST_GLOBAL_INFOF ("Using file names %s and %s\n", conn_name, p4est_name);

  /* check for ZLIB usability */
  if (!(have_zlib = p4est_have_zlib ())) {
    P4EST_GLOBAL_LERROR
      ("Not found a working ZLIB installation: ignoring CRCs\n");
  }

  p4est = p4est_new_ext (mpicomm, connectivity, 0, 0, 0,
                         sizeof (int), init_fn, NULL);
  p4est_refine (p4est, 1, refine_fn, init_fn);
  test_deflate (p4est);

  /* save, synchronize, load connectivity and compare */
  if (mpirank == 0) {
    retval = p4est_connectivity_save (conn_name, connectivity);
    SC_CHECK_ABORT (retval == 0, "connectivity_save failed");
  }
  mpiret = sc_MPI_Barrier (mpicomm);
  SC_CHECK_MPI (mpiret);

  wtime = sc_MPI_Wtime ();
  conn2 = p4est_connectivity_load (conn_name, NULL);
  elapsed = sc_MPI_Wtime () - wtime;
  sc_stats_set1 (stats + STATS_CONN_LOAD, elapsed, "conn load");

  SC_CHECK_ABORT (p4est_connectivity_is_equal (connectivity, conn2),
                  "load/save connectivity mismatch A");
  p4est_connectivity_destroy (conn2);

  /* save, synchronize, load p4est and compare */
  wtime = sc_MPI_Wtime ();
  p4est_save (p4est_name, p4est, 1);
  elapsed = sc_MPI_Wtime () - wtime;
  sc_stats_set1 (stats + STATS_P4EST_SAVE1, elapsed, "p4est save 1");

  wtime = sc_MPI_Wtime ();
  p4est2 = p4est_load (p4est_name, mpicomm, sizeof (int), 1, NULL, &conn2);
  elapsed = sc_MPI_Wtime () - wtime;
  sc_stats_set1 (stats + STATS_P4EST_LOAD1a, elapsed, "p4est load 1a");

  SC_CHECK_ABORT (p4est_connectivity_is_equal (connectivity, conn2),
                  "load/save connectivity mismatch Ba");
  SC_CHECK_ABORT (p4est_is_equal (p4est, p4est2, 1),
                  "load/save p4est mismatch Ba");
  p4est_destroy (p4est2);
  p4est_connectivity_destroy (conn2);

  wtime = sc_MPI_Wtime ();
  p4est2 = p4est_load (p4est_name, mpicomm, 0, 0, NULL, &conn2);
  elapsed = sc_MPI_Wtime () - wtime;
  sc_stats_set1 (stats + STATS_P4EST_LOAD1b, elapsed, "p4est load 1b");

  SC_CHECK_ABORT (p4est_connectivity_is_equal (connectivity, conn2),
                  "load/save connectivity mismatch Bb");
  SC_CHECK_ABORT (p4est_is_equal (p4est, p4est2, 0),
                  "load/save p4est mismatch Bb");
  test_deflate (p4est2);
  p4est_destroy (p4est2);
  p4est_connectivity_destroy (conn2);

  /* partition and balance */
  p4est_partition (p4est, 0, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);
  csum = test_checksum (p4est, have_zlib);
  sc_stats_set1 (stats + STATS_P4EST_ELEMS,
                 (double) p4est->local_num_quadrants, "p4est elements");

  /* save, synchronize, load p4est and compare */
  wtime = sc_MPI_Wtime ();
  p4est_save (p4est_name, p4est, 0);
  elapsed = sc_MPI_Wtime () - wtime;
  sc_stats_set1 (stats + STATS_P4EST_SAVE2, elapsed, "p4est save 2");

  wtime = sc_MPI_Wtime ();
  p4est2 = p4est_load (p4est_name, mpicomm, sizeof (int), 0, NULL, &conn2);
  elapsed = sc_MPI_Wtime () - wtime;
  sc_stats_set1 (stats + STATS_P4EST_LOAD2, elapsed, "p4est load 2");

  SC_CHECK_ABORT (p4est_connectivity_is_equal (connectivity, conn2),
                  "load/save connectivity mismatch C");
  SC_CHECK_ABORT (p4est_is_equal (p4est, p4est2, 0),
                  "load/save p4est mismatch C");
  p4est_destroy (p4est2);
  p4est_connectivity_destroy (conn2);

  /* save, synchronize, load p4est and compare */
  wtime = sc_MPI_Wtime ();
  p4est_save (p4est_name, p4est, 1);
  elapsed = sc_MPI_Wtime () - wtime;
  sc_stats_set1 (stats + STATS_P4EST_SAVE3, elapsed, "p4est save 3");

  wtime = sc_MPI_Wtime ();
  p4est2 = p4est_load (p4est_name, mpicomm, sizeof (int), 0, NULL, &conn2);
  elapsed = sc_MPI_Wtime () - wtime;
  sc_stats_set1 (stats + STATS_P4EST_LOAD3, elapsed, "p4est load 3");

  SC_CHECK_ABORT (p4est_connectivity_is_equal (connectivity, conn2),
                  "load/save connectivity mismatch D");
  SC_CHECK_ABORT (p4est_is_equal (p4est, p4est2, 0),
                  "load/save p4est mismatch D");
  p4est_destroy (p4est2);
  p4est_connectivity_destroy (conn2);

  /* Test autopartition load feature */
  wtime = sc_MPI_Wtime ();
  p4est2 = p4est_load_ext (p4est_name, mpicomm, sizeof (int), 0,
                           1, 0, NULL, &conn2);
  elapsed = sc_MPI_Wtime () - wtime;
  csum2 = test_checksum (p4est2, have_zlib);
  sc_stats_set1 (stats + STATS_P4EST_LOAD4, elapsed, "p4est load 4");

  SC_CHECK_ABORT (p4est_connectivity_is_equal (connectivity, conn2),
                  "load/save connectivity mismatch E");
  SC_CHECK_ABORT (mpirank != 0 || csum == csum2,
                  "load/save p4est mismatch E");
  p4est_destroy (p4est2);
  p4est_connectivity_destroy (conn2);

  /* destroy data structures */
  p4est_destroy (p4est);

  /* compute and print timings */
  sc_stats_compute (mpicomm, STATS_COUNT, stats);
  sc_stats_print (p4est_package_id, SC_LP_STATISTICS,
                  STATS_COUNT, stats, 0, 1);
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpirank;
  int                 first_arg;
  const char         *prefix;
  p4est_connectivity_t *connectivity;
  sc_options_t       *opt;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* initialize libsc and p4est */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* handle command line options */
  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'l', "level", &refine_level,
                      default_refine_level, "Refinement level");
  sc_options_add_string (opt, 'o', "oprefix", &prefix,
                         P4EST_STRING, "Output prefix");
  first_arg = sc_options_parse (p4est_package_id, SC_LP_INFO,
                                opt, argc, argv);
  SC_CHECK_ABORT (first_arg >= 0, "Option error");

  /* create connectivity */
#ifndef P4_TO_P8
  connectivity = p4est_connectivity_new_star ();
#else
  connectivity = p8est_connectivity_new_rotcubes ();
#endif

  /* test with vertex information */
  test_loadsave (connectivity, prefix, mpicomm, mpirank);

  /* test without vertex information */
  connectivity->num_vertices = 0;
  P4EST_FREE (connectivity->vertices);
  connectivity->vertices = NULL;
  P4EST_FREE (connectivity->tree_to_vertex);
  connectivity->tree_to_vertex = NULL;
  p4est_connectivity_set_attr (connectivity, 1);
  memset (connectivity->tree_to_attr, 0,
          connectivity->num_trees * sizeof (int8_t));
  test_loadsave (connectivity, prefix, mpicomm, mpirank);

  /* clean up and exit */
  p4est_connectivity_destroy (connectivity);
  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
