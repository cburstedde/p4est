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
#include <p4est_build.h>
#include <p4est_extended.h>
#include <p4est_geometry.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_build.h>
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
count_callback (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                void *point)
{
  P4EST_ASSERT (point == NULL);

  if (local_num == -1) {
    /* keep recursing to reach a leaf eventually */
    return 1;
  }
  else {
    p4est_locidx_t     *local_count = (p4est_locidx_t *) p4est->user_pointer;

    /* return value shall be ignored for leaves */
    P4EST_ASSERT (local_count != NULL);
    SC_CHECK_ABORT (local_num == *local_count, "Count mismatch");
    ++*local_count;
    return 0;
  }
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

typedef struct
{
  int                 maxlevel;
  int                 counter;
  int                 wrapper;
  int                 init_default;
  int                 init_add;
  int                 count_add;
  p4est_topidx_t      last_tree;
  p4est_build_t      *build;
}
test_build_t;

static int
test_build_refine (p4est_t * p4est, p4est_topidx_t which_tree,
                   p4est_quadrant_t * quadrant)
{
  test_build_t       *tb;

  tb = (test_build_t *) p4est->user_pointer;

  if (quadrant->level >= tb->maxlevel) {
    return 0;
  }
  return !(tb->counter = (tb->counter + 1) % tb->wrapper);
}

static int
test_build_coarsen (p4est_t * p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t * quadrants[])
{
  return 1;
}

static int
test_search_local_1 (p4est_t * p4est, p4est_topidx_t which_tree,
                     p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                     void *point)
{
  int                 retval;
  test_build_t       *tb;

  tb = (test_build_t *) p4est->user_pointer;

  /* take all quadrants and try duplicates regularly */
  if (local_num >= 0) {
    P4EST_EXECUTE_ASSERT_TRUE (p4est_build_add
                               (tb->build, which_tree, quadrant));
    if (!(tb->counter = (tb->counter + 1) % tb->wrapper)) {
      /* try to add it twice which should be reported properly */
      retval = p4est_build_add (tb->build, which_tree, quadrant);
      SC_CHECK_ABORT (!retval, "Tried to add a duplicate");
    }
  }

  return 1;
}

static int
test_search_local_2 (p4est_t * p4est, p4est_topidx_t which_tree,
                     p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                     void *point)
{
  return 1;
}

static void
test_search_init_3 (p4est_t * p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t * quadrant)
{
  test_build_t       *tb;

  P4EST_ASSERT (p4est->data_size == 0);

  tb = (test_build_t *) p4est->user_pointer;
  ++tb->init_default;

  quadrant->p.user_int = 1135;
}

static void
test_search_init_add_3 (p4est_t * p4est, p4est_topidx_t which_tree,
                        p4est_quadrant_t * quadrant)
{
  test_build_t       *tb;

  P4EST_ASSERT (p4est->data_size == 0);

  tb = (test_build_t *) p4est->user_pointer;
  ++tb->init_add;

  quadrant->p.user_int = 629;
}

static int
test_search_local_3 (p4est_t * p4est, p4est_topidx_t which_tree,
                     p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                     void *point)
{
  test_build_t       *tb;

  tb = (test_build_t *) p4est->user_pointer;

  /* take every third quadrant or so */
  if (local_num >= 0 && !(tb->counter = (tb->counter + 1) % tb->wrapper)) {
    ++tb->count_add;
    P4EST_EXECUTE_ASSERT_TRUE (p4est_build_add
                               (tb->build, which_tree, quadrant));
  }

  return 1;
}

static void
p4est_build_verify_3 (p4est_t * p4est)
{
  p4est_topidx_t      jt;
  p4est_locidx_t      il, c1, c2;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quadrant;
  test_build_t       *tb;

  tb = (test_build_t *) p4est->user_pointer;

  c1 = c2 = 0;
  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    for (il = 0; il < (p4est_locidx_t) tree->quadrants.elem_count; ++il) {
      quadrant = p4est_quadrant_array_index (&tree->quadrants, il);
      switch (quadrant->p.user_int) {
      case 1135:
        ++c1;
        break;
      case 629:
        ++c2;
        break;
      default:
        SC_ABORT_NOT_REACHED ();
      }
    }
  }
  SC_CHECK_ABORT (c1 + c2 == p4est->local_num_quadrants,
                  "Test 3 count quadrants");
  SC_CHECK_ABORT (c1 + c2 >= (p4est_locidx_t) tb->count_add,
                  "Test 3 count sum");
  SC_CHECK_ABORT (c1 == (p4est_locidx_t) tb->init_default,
                  "Test 3 count default");
  SC_CHECK_ABORT (c2 == (p4est_locidx_t) tb->init_add, "Test 3 count add");
}

static void
test_search_init_4 (p4est_t * p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t * quadrant)
{
  test_build_t       *tb;

  P4EST_ASSERT (p4est->data_size == sizeof (long));

  tb = (test_build_t *) p4est->user_pointer;
  ++tb->init_default;

  *(long *) quadrant->p.user_data = 11321;
}

static void
test_search_init_add_4 (p4est_t * p4est, p4est_topidx_t which_tree,
                        p4est_quadrant_t * quadrant)
{
  test_build_t       *tb;

  P4EST_ASSERT (p4est->data_size == sizeof (long));

  tb = (test_build_t *) p4est->user_pointer;
  ++tb->init_add;

  *(long *) quadrant->p.user_data = -748;
}

static int
test_search_local_4 (p4est_t * p4est, p4est_topidx_t which_tree,
                     p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                     void *point)
{
  test_build_t       *tb;

  tb = (test_build_t *) p4est->user_pointer;

  /* take the first third quadrant or so in every tree */
  if (local_num >= 0 && tb->last_tree != which_tree &&
      !(tb->counter = (tb->counter + 1) % tb->wrapper)) {
    ++tb->count_add;
    P4EST_EXECUTE_ASSERT_TRUE (p4est_build_add
                               (tb->build, which_tree, quadrant));

    tb->last_tree = which_tree;
  }

  return 1;
}

static void
p4est_build_verify_4 (p4est_t * p4est)
{
  p4est_topidx_t      jt;
  p4est_locidx_t      il, c1, c2;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quadrant;
  test_build_t       *tb;

  tb = (test_build_t *) p4est->user_pointer;

  c1 = c2 = 0;
  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    for (il = 0; il < (p4est_locidx_t) tree->quadrants.elem_count; ++il) {
      quadrant = p4est_quadrant_array_index (&tree->quadrants, il);
      switch (*(long *) quadrant->p.user_data) {
      case 11321:
        ++c1;
        break;
      case -748:
        ++c2;
        break;
      default:
        SC_ABORT_NOT_REACHED ();
      }
    }
  }
  SC_CHECK_ABORT (c1 + c2 == p4est->local_num_quadrants,
                  "Test 4 count quadrants");
  SC_CHECK_ABORT (c1 + c2 >= (p4est_locidx_t) tb->count_add,
                  "Test 4 count sum");
  SC_CHECK_ABORT (c1 == (p4est_locidx_t) tb->init_default,
                  "Test 4 count default");
  SC_CHECK_ABORT (c2 == (p4est_locidx_t) tb->init_add, "Test 4 count add");
}

static int
test_search_point_5 (p4est_t * p4est, p4est_topidx_t which_tree,
                     p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                     void *point)
{
  int                 retval;
#ifdef P4EST_ENABLE_DEBUG
  int8_t              ip;
#endif
  test_build_t       *tb;

  P4EST_ASSERT (point != NULL);

#ifdef P4EST_ENABLE_DEBUG
  ip = *(int8_t *) point;
#endif
  P4EST_ASSERT (0 <= ip && ip < 2);

  tb = (test_build_t *) p4est->user_pointer;

  if (!(tb->counter = (tb->counter + 1) % tb->wrapper)) {
    /* rare */
    return 0;
  }
  else {
    /* frequent */

    if (local_num >= 0) {
      /* this is a leaf, add it */
      retval = p4est_build_add (tb->build, which_tree, quadrant);
      if (!retval) {
        /* this is a duplicate leaf */
        ++tb->init_default;
      }
      else {
        /* first addition of this quadrant */
        ++tb->init_add;
      }
    }
    return 1;
  }
}

#if 0

static void
p4est_build_verify_5 (p4est_t * p4est)
{
  test_build_t       *tb;

  tb = (test_build_t *) p4est->user_pointer;

  P4EST_LDEBUGF ("T5 added %d dup %d\n", tb->init_add, tb->init_default);
}

#endif

static void
test_build_local (sc_MPI_Comm mpicomm)
{
  sc_array_t         *points;
  p4est_connectivity_t *conn;
  p4est_t            *p4est, *built, *copy;
  test_build_t        stb, *tb = &stb;

  /* 0. prepare data that we will reuse */
  tb->maxlevel = 7 - P4EST_DIM;
  tb->counter = -1;
  tb->wrapper = 3;
  tb->init_default = -1;
  tb->init_add = -1;
  tb->count_add = -1;
  tb->last_tree = -1;
  tb->build = NULL;
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_moebius ();
#else
  conn = p8est_connectivity_new_rotcubes ();
#endif /* P4_TO_P8 */
  p4est = p4est_new_ext (mpicomm, conn, 0, 0, 2, 0, NULL, tb);
  p4est_refine (p4est, 1, test_build_refine, NULL);
  p4est_partition (p4est, 0, NULL);

  /* TODO: enrich tests with quadrant data */

  /* 1. Create a p4est that shall be identical to the old one. */

  tb->build = p4est_build_new (p4est, 0, NULL, NULL);
  p4est_search_local (p4est, 0, test_search_local_1, NULL, NULL);
  built = p4est_build_complete (tb->build);
  SC_CHECK_ABORT (p4est_is_equal (p4est, built, 0), "Mismatch build_local 1");
  p4est_destroy (built);

  /* 2. Create a p4est that is as coarse as possible.
   *    Coarsen recursively, compare. */

  tb->build = p4est_build_new (p4est, 4, NULL, NULL);
  p4est_search_local (p4est, 0, test_search_local_2, NULL, NULL);
  built = p4est_build_complete (tb->build);
  copy = p4est_copy (p4est, 0);
  p4est_coarsen (copy, 1, test_build_coarsen, NULL);
  SC_CHECK_ABORT (p4est_is_equal (copy, built, 0), "Mismatch build_local 2");
  p4est_destroy (copy);
  p4est_destroy (built);

  /* 3. Create a p4est with some random pattern for demonstration */

  tb->init_default = 0;
  tb->init_add = 0;
  tb->count_add = 0;
  tb->build = p4est_build_new (p4est, 0, test_search_init_3, tb);
  p4est_build_init_add (tb->build, test_search_init_add_3);
  p4est_search_local (p4est, 1, test_search_local_3, NULL, NULL);
  built = p4est_build_complete (tb->build);
  p4est_build_verify_3 (built);
  SC_CHECK_ABORT (p4est_is_valid (built), "Invalid build_local 3");
  p4est_destroy (built);

  /* 4. Create a p4est from a search with one quadrant per tree */

  tb->init_default = 0;
  tb->init_add = 0;
  tb->count_add = 0;
  tb->last_tree = -1;
  tb->build = p4est_build_new (p4est, sizeof (long), test_search_init_4, tb);
  p4est_build_init_add (tb->build, test_search_init_add_4);
  p4est_search_local (p4est, 0, test_search_local_4, NULL, NULL);
  built = p4est_build_complete (tb->build);
  p4est_build_verify_4 (built);
  SC_CHECK_ABORT (p4est_is_valid (built), "Invalid build_local 4");
  p4est_destroy (built);

  /* 5. Create a p4est from a multiple-item search */

  points = sc_array_new_size (sizeof (int8_t), 2);
  *(int8_t *) sc_array_index (points, 0) = 0;
  *(int8_t *) sc_array_index (points, 1) = 1;
  tb->wrapper = 5;
  tb->init_default = 0;
  tb->init_add = 0;
  tb->build = p4est_build_new (p4est, 0, NULL, tb);
  p4est_search_local (p4est, 0, NULL, test_search_point_5, points);
  built = p4est_build_complete (tb->build);
#if 0
  p4est_build_verify_5 (built);
#endif
  SC_CHECK_ABORT (p4est_is_valid (built), "Invalid build_local 5");
  p4est_destroy (built);
  sc_array_destroy (points);

  /* clean up */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 found_total;
  p4est_locidx_t      jt, Al, Bl;
  p4est_locidx_t      local_count;
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
  p4est = p4est_new_ext (mpicomm, conn, 0, 0, 0, 0, NULL, &local_count);
  p4est_refine (p4est, 1, refine_fn, NULL);
  p4est_partition (p4est, 0, NULL);
  p4est_vtk_write_file (p4est, geom, vtkname);

  /* The following code should really be in a separate function. */

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
  p4est_search_local (p4est, 0, NULL, search_callback, points);
  mpiret = sc_MPI_Allreduce (&found_count, &found_total,
                             1, sc_MPI_INT, sc_MPI_SUM, mpicomm);
  SC_CHECK_MPI (mpiret);
  SC_CHECK_ABORT (found_total == (int) points->elem_count, "Point search");
  SC_CHECK_ABORT (A->p.piggy3.local_num == Al, "Search A");
  SC_CHECK_ABORT (B->p.piggy3.local_num == Bl, "Search B");

  /* Use another search to count local quadrants */
  local_count = 0;
  p4est_search_local (p4est, 0, count_callback, NULL, NULL);
  SC_CHECK_ABORT (local_count == p4est->local_num_quadrants, "Count search");

  /* Clear memory */
  sc_array_destroy (points);
  p4est_destroy (p4est);
  if (geom != NULL) {
    p4est_geometry_destroy (geom);
  }
  p4est_connectivity_destroy (conn);

  /* Test the build_local function and friends */
  test_build_local (mpicomm);

  /* Finalize */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
