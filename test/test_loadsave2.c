/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2008 Carsten Burstedde, Lucas Wilcox.

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

#ifndef P4_TO_P8
#include <p4est.h>
#include <p4est_bits.h>
#include <p4est_algorithms.h>
#else
#include <p8est.h>
#include <p8est_bits.h>
#include <p8est_algorithms.h>
#endif

#ifndef P4_TO_P8
#define P4EST_CONN_SUFFIX "p4c"
#define P4EST_FOREST_SUFFIX "p4p"
static const int    refine_level = 8;
#else
#define P4EST_CONN_SUFFIX "p8c"
#define P4EST_FOREST_SUFFIX "p8p"
static const int    refine_level = 5;
#endif

static int          counter = 0;

static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
  int                *data = quadrant->p.user_data;

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

void
test_loadsave (p4est_connectivity_t * connectivity,
               MPI_Comm mpicomm, int mpirank)
{
  int                 mpiret;
  p4est_connectivity_t *conn2;
  p4est_t            *p4est, *p4est2;

  p4est = p4est_new (mpicomm, connectivity, 0, sizeof (int), init_fn, NULL);
  p4est_refine (p4est, true, refine_fn, init_fn);

  /* save, synchronize, load connectivity and compare */
  if (mpirank == 0) {
    p4est_connectivity_save (P4EST_STRING "." P4EST_CONN_SUFFIX,
                             connectivity);
  }
  mpiret = MPI_Barrier (mpicomm);
  SC_CHECK_MPI (mpiret);
  conn2 = p4est_connectivity_load (P4EST_STRING "." P4EST_CONN_SUFFIX, NULL);
  SC_CHECK_ABORT (p4est_connectivity_is_equal (connectivity, conn2),
                  "load/save connectivity mismatch A");
  p4est_connectivity_destroy (conn2);

  /* save, synchronize, load p4est and compare */
  p4est_save (P4EST_STRING "." P4EST_FOREST_SUFFIX, p4est, true);
  p4est2 = p4est_load (P4EST_STRING "." P4EST_FOREST_SUFFIX,
                       mpicomm, sizeof (int), true, NULL, &conn2);
  SC_CHECK_ABORT (p4est_connectivity_is_equal (connectivity, conn2),
                  "load/save connectivity mismatch B");
  SC_CHECK_ABORT (p4est_is_equal (p4est, p4est2, true),
                  "load/save p4est mismatch B");
  p4est_destroy (p4est2);
  p4est_connectivity_destroy (conn2);

  /* partition (still not balanced) */
  p4est_partition (p4est, NULL);

  /* save, synchronize, load p4est and compare */
  p4est_save (P4EST_STRING "." P4EST_FOREST_SUFFIX, p4est, false);
  p4est2 = p4est_load (P4EST_STRING "." P4EST_FOREST_SUFFIX,
                       mpicomm, sizeof (int), false, NULL, &conn2);
  SC_CHECK_ABORT (p4est_connectivity_is_equal (connectivity, conn2),
                  "load/save connectivity mismatch C");
  SC_CHECK_ABORT (p4est_is_equal (p4est, p4est2, false),
                  "load/save p4est mismatch C");
  p4est_destroy (p4est2);
  p4est_connectivity_destroy (conn2);

  /* save, synchronize, load p4est and compare */
  p4est_save (P4EST_STRING "." P4EST_FOREST_SUFFIX, p4est, true);
  p4est2 = p4est_load (P4EST_STRING "." P4EST_FOREST_SUFFIX,
                       mpicomm, sizeof (int), false, NULL, &conn2);
  SC_CHECK_ABORT (p4est_connectivity_is_equal (connectivity, conn2),
                  "load/save connectivity mismatch D");
  SC_CHECK_ABORT (p4est_is_equal (p4est, p4est2, false),
                  "load/save p4est mismatch D");
  p4est_destroy (p4est2);
  p4est_connectivity_destroy (conn2);

  /* destroy data structures */
  p4est_destroy (p4est);
}

int
main (int argc, char **argv)
{
  MPI_Comm            mpicomm;
  int                 mpiret;
  int                 mpirank;
  p4est_connectivity_t *connectivity;

  /* initialize MPI */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, true, true, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* create connectivity and p4est (not balanced) */
#ifndef P4_TO_P8
  connectivity = p4est_connectivity_new_star ();
#else
  connectivity = p8est_connectivity_new_rotcubes ();
#endif

  /* test with vertex information */
  test_loadsave (connectivity, mpicomm, mpirank);

#ifdef P4_TO_P8
  /* test without vertex information */
  connectivity->num_vertices = 0;
  P4EST_FREE (connectivity->vertices);
  connectivity->vertices = NULL;
  P4EST_FREE (connectivity->tree_to_vertex);
  connectivity->tree_to_vertex = NULL;
  test_loadsave (connectivity, mpicomm, mpirank);
#endif

  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_loadsave2.c */
