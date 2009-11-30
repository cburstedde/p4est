/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2009 Carsten Burstedde, Lucas Wilcox.

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
#include <p4est_bits.h>
#include <p4est_geometry.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_geometry.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#endif

typedef struct
{
  const char         *name;
  p4est_quadrant_t    quad;
}
test_point_t;

static const int    refine_level = 3;
static int          found_count = -1;

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 cid;

  if ((int) quadrant->level >= refine_level)
    return 0;

  if (which_tree == 2 || which_tree == 5)
    return 0;

  cid = p4est_quadrant_child_id (quadrant);
  if (cid == 0 || cid == 1 || cid == 6)
    return 1;

  if (quadrant->x >= P4EST_LAST_OFFSET (2)
#ifdef P4_TO_P8
      && quadrant->z >= P4EST_LAST_OFFSET (2)
#endif
    ) {
    return 1;
  }

  return 0;
}

static int
search_callback (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * quadrant, int is_leaf, void *point)
{
  test_point_t       *p = (test_point_t *) point;
  int                 is_match;

  P4EST_LDEBUGF ("Tree %lld quadrant %s level %d %d child %d leaf %d\n",
                 (long long) which_tree, p->name,
                 (int) p->quad.level, (int) quadrant->level,
                 p4est_quadrant_child_id (quadrant), is_leaf);

  if (which_tree != p->quad.p.which_tree) {
    P4EST_LDEBUGF ("Rejecting quadrant %s for tree %lld\n",
                   p->name, (long long) which_tree);
    return 0;
  }

  if (quadrant->level < p->quad.level) {
    is_match = p4est_quadrant_is_ancestor (quadrant, &p->quad);
    P4EST_LDEBUGF ("Ancestor for quadrant %s is %d\n", p->name, is_match);
  }
  else {
    is_match = !p4est_quadrant_compare (quadrant, &p->quad);
    P4EST_LDEBUGF ("Tree %lld same size quadrant %s match %d\n",
                   (long long) which_tree, p->name, is_match);
  }

  if (is_match && is_leaf) {
    P4EST_INFOF ("Matched quadrant %s\n", p->name);
    p4est_quadrant_print (SC_LP_INFO, quadrant);
    p4est_quadrant_print (SC_LP_INFO, &p->quad);
    ++found_count;
  }

  return is_match;
}

int
main (int argc, char **argv)
{
  MPI_Comm            mpicomm;
  int                 mpiret;
  int                 found_total;
  p4est_connectivity_t *conn;
  p4est_geometry_t   *geom;
  p4est_t            *p4est;
  sc_array_t         *points;
  test_point_t       *p;
  const char         *vtkname;

  /* Initialize MPI */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;

  /* Initialize packages */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* Create forest */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_star ();
  geom = NULL;
  vtkname = "test_search2";
#else
  conn = p8est_connectivity_new_sphere ();
  geom = p8est_geometry_new_sphere (1., 0.191728, 0.039856);
  vtkname = "test_search3";
#endif
  p4est = p4est_new (mpicomm, conn, 0, 0, NULL, NULL);
  p4est_refine (p4est, 1, refine_fn, NULL);
  p4est_partition (p4est, NULL);
  p4est_vtk_write_file (p4est, geom, vtkname);

  /* Prepare a point search */
  points = sc_array_new (sizeof (test_point_t));

  /* A */
  p = (test_point_t *) sc_array_push (points);
  p->name = "A";
  P4EST_QUADRANT_INIT (&p->quad);
  p4est_quadrant_set_morton (&p->quad, 3, 23);
  p->quad.p.which_tree = 0;

  /* B */
  p = (test_point_t *) sc_array_push (points);
  p->name = "B";
  P4EST_QUADRANT_INIT (&p->quad);
  p4est_quadrant_set_morton (&p->quad, 2, 13);
  p->quad.p.which_tree = conn->num_trees / 2;

  /* Go */
  found_count = 0;
  p4est_search (p4est, search_callback, points);
  mpiret = MPI_Allreduce (&found_count, &found_total,
                          1, MPI_INT, MPI_SUM, mpicomm);
  SC_CHECK_MPI (mpiret);
  SC_CHECK_ABORT (found_total == (int) points->elem_count, "Point search");

  /* Clear memory */
  sc_array_destroy (points);
  p4est_destroy (p4est);
  P4EST_FREE (geom);
  p4est_connectivity_destroy (conn);

  /* Finalize */
  sc_finalize ();
  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
