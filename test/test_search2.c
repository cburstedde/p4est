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

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_geometry.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
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
                 p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                 void *point)
{
  test_point_t       *p = (test_point_t *) point;
  int                 is_leaf;
  int                 is_match;

  is_leaf = local_num >= 0;
  P4EST_ASSERT (!is_leaf || local_num < p4est->local_num_quadrants);
  P4EST_ASSERT (point != NULL);

  P4EST_LDEBUGF ("Tree %lld quadrant %s level %d %d child %d leaf %d\n",
                 (long long) which_tree, p->name,
                 (int) p->quad.level, (int) quadrant->level,
                 p4est_quadrant_child_id (quadrant), is_leaf);

  if (which_tree != p->quad.p.piggy3.which_tree) {
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
    p4est_locidx_t      num = -1;

    if (quadrant->level < p->quad.level) {
      num = p->quad.p.piggy3.local_num = -1;
    }
    else {
      P4EST_ASSERT (local_num >= 0);
      num = p->quad.p.piggy3.local_num = local_num;
    }
    P4EST_INFOF ("Matched quadrant %s at %lld\n", p->name, (long long) num);
    p4est_quadrant_print (SC_LP_INFO, quadrant);
    p4est_quadrant_print (SC_LP_INFO, &p->quad);
    ++found_count;
  }

  return is_match;
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 found_total;
  p4est_locidx_t      jt, Al, Bl;
  p4est_connectivity_t *conn;
  p4est_quadrant_t   *A, *B;
  p4est_geometry_t   *geom;
  p4est_t            *p4est;
  sc_array_t         *points;
  test_point_t       *p;
  const char         *vtkname;

  /* Initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

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
  geom = p8est_geometry_new_sphere (conn, 1., 0.191728, 0.039856);
  vtkname = "test_search3";
#endif
  p4est = p4est_new_ext (mpicomm, conn, 0, 0, 0, 0, NULL, NULL);
  p4est_refine (p4est, 1, refine_fn, NULL);
  p4est_partition (p4est, 0, NULL);
  p4est_vtk_write_file (p4est, geom, vtkname);

  /* Prepare a point search -- fix size so the memory is not relocated */
  points = sc_array_new_size (sizeof (test_point_t), 2);

  /* A */
  p = (test_point_t *) sc_array_index (points, 0);
  p->name = "A";
  A = &p->quad;
  P4EST_QUADRANT_INIT (A);
  p4est_quadrant_set_morton (A, 3, 23);
  A->p.piggy3.which_tree = 0;
  A->p.piggy3.local_num = -1;
  Al = -1;

  /* B */
  p = (test_point_t *) sc_array_index (points, 1);
  p->name = "B";
  B = &p->quad;
  P4EST_QUADRANT_INIT (B);
  p4est_quadrant_set_morton (B, 2, 13);
  B->p.piggy3.which_tree = conn->num_trees / 2;
  B->p.piggy3.local_num = -1;
  Bl = -1;

  /* Find quadrant numbers if existing */
  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    size_t              zz;
    p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, jt);
    p4est_quadrant_t   *quad;
    sc_array_t         *tquadrants = &tree->quadrants;

    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      quad = p4est_quadrant_array_index (tquadrants, zz);
      if (A->p.piggy3.which_tree == jt && !p4est_quadrant_compare (quad, A)) {
        Al = tree->quadrants_offset + (p4est_locidx_t) zz;
        P4EST_VERBOSEF ("Searching for A at %lld\n", (long long) Al);
      }
      if (B->p.piggy3.which_tree == jt && !p4est_quadrant_compare (quad, B)) {
        Bl = tree->quadrants_offset + (p4est_locidx_t) zz;
        P4EST_VERBOSEF ("Searching for B at %lld\n", (long long) Bl);
      }
    }
  }

  /* Go */
  found_count = 0;
  p4est_search (p4est, NULL, search_callback, points);
  mpiret = sc_MPI_Allreduce (&found_count, &found_total,
                             1, sc_MPI_INT, sc_MPI_SUM, mpicomm);
  SC_CHECK_MPI (mpiret);
  SC_CHECK_ABORT (found_total == (int) points->elem_count, "Point search");
  SC_CHECK_ABORT (A->p.piggy3.local_num == Al, "Search A");
  SC_CHECK_ABORT (B->p.piggy3.local_num == Bl, "Search B");

  /* Clear memory */
  sc_array_destroy (points);
  p4est_destroy (p4est);
  if (geom != NULL) {
    p4est_geometry_destroy (geom);
  }
  p4est_connectivity_destroy (conn);

  /* Finalize */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
