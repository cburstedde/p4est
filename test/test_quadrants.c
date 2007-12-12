/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007 Carsten Burstedde, Lucas Wilcox.

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

#include <p4est_algorithms.h>
#include <p4est_base.h>

static int
refine_none (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * q)
{
  return 0;
}

static int
refine_some (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * q)
{
  if (q->x < (1 << (P4EST_MAXLEVEL - 2))) {
    return q->level <= 4;
  }
  else if (q->x < (1 << (P4EST_MAXLEVEL - 1))) {
    return q->level <= 3;
  }
  else {
    return q->level <= 2;
  }
}

static int
coarsen_none (p4est_t * p4est, int32_t which_tree,
              p4est_quadrant_t * q0, p4est_quadrant_t * q1,
              p4est_quadrant_t * q2, p4est_quadrant_t * q3)
{
  return 0;
}

static int
coarsen_some (p4est_t * p4est, int32_t which_tree,
              p4est_quadrant_t * q0, p4est_quadrant_t * q1,
              p4est_quadrant_t * q2, p4est_quadrant_t * q3)
{
  if (q0->x < (1 << (P4EST_MAXLEVEL - 2))) {
    return q0->level >= 2;
  }
  else if (q0->x < (1 << (P4EST_MAXLEVEL - 1))) {
    return q0->level >= 3;
  }
  else {
    return q0->level >= 4;
  }
}

static int
coarsen_all (p4est_t * p4est, int32_t which_tree,
             p4est_quadrant_t * q0, p4est_quadrant_t * q1,
             p4est_quadrant_t * q2, p4est_quadrant_t * q3)
{
  return 1;
}

int
main (void)
{
  int                 i, j, incount;
  int8_t              level, mid, cid;
  int64_t             index1, index2;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est1;
  p4est_t            *p4est2;
  p4est_tree_t       *t1, *t2, tree;
  p4est_quadrant_t   *p, *q1, *q2;
  p4est_quadrant_t    r, s;
  p4est_quadrant_t    c0, c1, c2, c3;

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_unitsquare ();
  p4est1 = p4est_new (MPI_COMM_NULL, stdout, connectivity, 0, NULL);
  p4est2 = p4est_new (MPI_COMM_NULL, stdout, connectivity, 8, NULL);

  /* refine the second tree to a uniform level */
  p4est_refine (p4est1, refine_none, NULL);
  p4est_refine (p4est2, refine_some, NULL);
  t1 = p4est_array_index (p4est1->trees, 0);
  t2 = p4est_array_index (p4est2->trees, 0);
  P4EST_CHECK_ABORT (p4est_tree_is_sorted (t1), "is_sorted");
  P4EST_CHECK_ABORT (p4est_tree_is_sorted (t2), "is_sorted");

  /* run a bunch of cross-tests */
  p = NULL;
  for (i = 0; i < t1->quadrants->elem_count; ++i) {
    q1 = p4est_array_index (t1->quadrants, i);

    /* test the index conversion */
    index1 = p4est_quadrant_linear_id (q1, q1->level);
    p4est_quadrant_set_morton (&r, q1->level, index1);
    index2 = p4est_quadrant_linear_id (&r, r.level);
    P4EST_CHECK_ABORT (index1 == index2, "index conversion");
    level = (int8_t) (q1->level - 1);
    if (level >= 0) {
      index1 = p4est_quadrant_linear_id (q1, level);
      p4est_quadrant_set_morton (&r, level, index1);
      index2 = p4est_quadrant_linear_id (&r, level);
      P4EST_CHECK_ABORT (index1 == index2, "index conversion");
    }

    /* test the is_next function */
    if (p != NULL) {
      P4EST_CHECK_ABORT (p4est_quadrant_is_next (p, q1), "is_next");
    }
    p = q1;

    /* test the is_family function */
    p4est_quadrant_children (q1, &c0, &c1, &c2, &c3);
    P4EST_CHECK_ABORT (p4est_quadrant_is_family (&c0, &c1, &c2, &c3),
                       "is_family");
    P4EST_CHECK_ABORT (!p4est_quadrant_is_family (&c1, &c0, &c2, &c3),
                       "is_family");
    P4EST_CHECK_ABORT (!p4est_quadrant_is_family (&c0, &c0, &c1, &c2),
                       "is_family");

    /* test the sibling function */
    mid = p4est_quadrant_child_id (q1);
    for (cid = 0; cid < 4; ++cid) {
      p4est_quadrant_sibling (q1, &r, cid);
      if (cid != mid) {
        P4EST_CHECK_ABORT (p4est_quadrant_is_sibling (q1, &r), "sibling");
      }
      else {
        P4EST_CHECK_ABORT (p4est_quadrant_is_equal (q1, &r), "sibling");
      }
    }

    /* test t1 against itself */
    for (j = 0; j < t1->quadrants->elem_count; ++j) {
      q2 = p4est_array_index (t1->quadrants, j);

      /* test the comparison function */
      P4EST_CHECK_ABORT (p4est_quadrant_compare (q1, q2) ==
                         -p4est_quadrant_compare (q2, q1), "compare");
      P4EST_CHECK_ABORT ((p4est_quadrant_compare (q1, q2) == 0) ==
                         p4est_quadrant_is_equal (q1, q2), "is_equal");

      /* test the descriptive versions of functions */
      P4EST_CHECK_ABORT (p4est_quadrant_is_sibling_D (q1, q2) ==
                         p4est_quadrant_is_sibling (q1, q2), "is_sibling");
      P4EST_CHECK_ABORT (p4est_quadrant_is_parent_D (q1, q2) ==
                         p4est_quadrant_is_parent (q1, q2), "is_parent");
      P4EST_CHECK_ABORT (p4est_quadrant_is_parent_D (q2, q1) ==
                         p4est_quadrant_is_parent (q2, q1), "is_parent");
      P4EST_CHECK_ABORT (p4est_quadrant_is_ancestor_D (q1, q2) ==
                         p4est_quadrant_is_ancestor (q1, q2), "is_ancestor");
      P4EST_CHECK_ABORT (p4est_quadrant_is_ancestor_D (q2, q1) ==
                         p4est_quadrant_is_ancestor (q2, q1), "is_ancestor");
      P4EST_CHECK_ABORT (p4est_quadrant_is_next_D (q1, q2) ==
                         p4est_quadrant_is_next (q1, q2), "is_next");
      P4EST_CHECK_ABORT (p4est_quadrant_is_next_D (q2, q1) ==
                         p4est_quadrant_is_next (q2, q1), "is_next");
      p4est_nearest_common_ancestor_D (q1, q2, &r);
      p4est_nearest_common_ancestor (q1, q2, &s);
      P4EST_CHECK_ABORT (p4est_quadrant_is_equal (&r, &s), "common_ancestor");
      p4est_nearest_common_ancestor_D (q2, q1, &r);
      p4est_nearest_common_ancestor (q2, q1, &s);
      P4EST_CHECK_ABORT (p4est_quadrant_is_equal (&r, &s), "common_ancestor");
    }

    /* test t1 against t2 */
    for (j = 0; j < t2->quadrants->elem_count; ++j) {
      q2 = p4est_array_index (t2->quadrants, j);

      /* test the comparison function */
      P4EST_CHECK_ABORT (p4est_quadrant_compare (q1, q2) ==
                         -p4est_quadrant_compare (q2, q1), "compare");
      P4EST_CHECK_ABORT ((p4est_quadrant_compare (q1, q2) == 0) ==
                         p4est_quadrant_is_equal (q1, q2), "is_equal");

      /* test the descriptive versions of functions */
      P4EST_CHECK_ABORT (p4est_quadrant_is_sibling_D (q1, q2) ==
                         p4est_quadrant_is_sibling (q1, q2), "is_sibling");
      P4EST_CHECK_ABORT (p4est_quadrant_is_parent_D (q1, q2) ==
                         p4est_quadrant_is_parent (q1, q2), "is_parent");
      P4EST_CHECK_ABORT (p4est_quadrant_is_parent_D (q2, q1) ==
                         p4est_quadrant_is_parent (q2, q1), "is_parent");
      P4EST_CHECK_ABORT (p4est_quadrant_is_ancestor_D (q1, q2) ==
                         p4est_quadrant_is_ancestor (q1, q2), "is_ancestor");
      P4EST_CHECK_ABORT (p4est_quadrant_is_ancestor_D (q2, q1) ==
                         p4est_quadrant_is_ancestor (q2, q1), "is_ancestor");
      P4EST_CHECK_ABORT (p4est_quadrant_is_next_D (q1, q2) ==
                         p4est_quadrant_is_next (q1, q2), "is_next");
      P4EST_CHECK_ABORT (p4est_quadrant_is_next_D (q2, q1) ==
                         p4est_quadrant_is_next (q2, q1), "is_next");
      p4est_nearest_common_ancestor_D (q1, q2, &r);
      p4est_nearest_common_ancestor (q1, q2, &s);
      P4EST_CHECK_ABORT (p4est_quadrant_is_equal (&r, &s), "common_ancestor");
      p4est_nearest_common_ancestor_D (q2, q1, &r);
      p4est_nearest_common_ancestor (q2, q1, &s);
      P4EST_CHECK_ABORT (p4est_quadrant_is_equal (&r, &s), "common_ancestor");
    }
  }

  p = NULL;
  for (i = 0; i < t2->quadrants->elem_count; ++i) {
    q1 = p4est_array_index (t2->quadrants, i);

    /* test the is_next function */
    if (p != NULL) {
      P4EST_CHECK_ABORT (p4est_quadrant_is_next (p, q1), "is_next");
    }
    p = q1;
  }

  /* test the coarsen function */
  p4est_coarsen (p4est1, coarsen_none, NULL);
  p4est_coarsen (p4est1, coarsen_all, NULL);
  p4est_coarsen (p4est2, coarsen_some, NULL);

  /* test the linearize algorithm */
  incount = t2->quadrants->elem_count;
  p4est_linearize_subtree (p4est2, t2);
  P4EST_CHECK_ABORT (incount == t2->quadrants->elem_count, "linearize");

  /* this is user_data neutral only when p4est1->data_size == 0 */
  tree.quadrants = p4est_array_new (sizeof (p4est_quadrant_t));
  p4est_array_resize (tree.quadrants, 18);
  q1 = p4est_array_index (tree.quadrants, 0);
  q2 = p4est_array_index (t2->quadrants, 0);
  *q1 = *q2;
  q2 = p4est_array_index (t2->quadrants, 1);
  for (i = 0; i < 3; ++i) {
    q1 = p4est_array_index (tree.quadrants, i + 1);
    *q1 = *q2;
    q1->level = (int8_t) (q1->level + i);
  }
  for (i = 0; i < 10; ++i) {
    q1 = p4est_array_index (tree.quadrants, i + 4);
    q2 = p4est_array_index (t2->quadrants, i + 3);
    *q1 = *q2;
    q1->level = (int8_t) (q1->level + i);
  }
  for (i = 0; i < 4; ++i) {
    q1 = p4est_array_index (tree.quadrants, i + 14);
    q2 = p4est_array_index (t2->quadrants, i + 12);
    *q1 = *q2;
    q1->level = (int8_t) (q1->level + 10 + i);
  }
  tree.maxlevel = 0;
  for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
    tree.quadrants_per_level[i] = 0;
  }
  incount = tree.quadrants->elem_count;
  for (i = 0; i < incount; ++i) {
    q1 = p4est_array_index (tree.quadrants, i);
    ++tree.quadrants_per_level[q1->level];
    tree.maxlevel = (int8_t) P4EST_MAX (tree.maxlevel, q1->level);
  }
  P4EST_CHECK_ABORT (!p4est_tree_is_linear (&tree), "is_linear");
  p4est_linearize_subtree (p4est1, &tree);
  P4EST_CHECK_ABORT (incount - 3 == tree.quadrants->elem_count, "linearize");
  p4est_array_destroy (tree.quadrants);

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est1);
  p4est_destroy (p4est2);
  p4est_connectivity_destroy (connectivity);

  p4est_memory_check ();

  return 0;
}

/* EOF test_quadrants.c */
