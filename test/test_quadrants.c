
#include <p4est_algorithms.h>
#include <p4est_base.h>

static int
refine_fn (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * q)
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

int
main (void)
{
  int                 i, j;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est1;
  p4est_t            *p4est2;
  p4est_tree_t       *t1, *t2;
  p4est_quadrant_t   *q1, *q2;
  p4est_quadrant_t    r, s;

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_unitsquare ();
  p4est1 = p4est_new (MPI_COMM_NULL, stdout, connectivity, 0, NULL);
  p4est2 = p4est_new (MPI_COMM_NULL, stdout, connectivity, 8, NULL);

  /* refine the second tree to a uniform level */
  p4est_refine (p4est2, refine_fn, NULL);
  t1 = p4est_array_index (p4est1->trees, 0);
  t2 = p4est_array_index (p4est2->trees, 0);
  P4EST_CHECK_ABORT (p4est_tree_is_sorted (t1), "is_sorted");
  P4EST_CHECK_ABORT (p4est_tree_is_sorted (t2), "is_sorted");

  /* run a bunch of cross-tests */
  for (i = 0; i < t1->quadrants->elem_count; ++i) {
    q1 = p4est_array_index (t1->quadrants, i);

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
      p4est_nearest_common_ancestor_D (q1, q2, &r);
      p4est_nearest_common_ancestor (q1, q2, &s);
      P4EST_CHECK_ABORT (p4est_quadrant_is_equal (&r, &s), "common_ancestor");
      p4est_nearest_common_ancestor_D (q2, q1, &r);
      p4est_nearest_common_ancestor (q2, q1, &s);
      P4EST_CHECK_ABORT (p4est_quadrant_is_equal (&r, &s), "common_ancestor");
    }
  }

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est1);
  p4est_destroy (p4est2);
  p4est_connectivity_destroy (connectivity);

  p4est_memory_check ();

  return 0;
}

/* EOF test_quadrants.c */
