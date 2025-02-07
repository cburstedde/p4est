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

/* we grab a few static functions from the 2D file */
#include <p4est_to_p8est.h>
#include "test_quadrants2.c"

static void
check_linear_id (const p4est_quadrant_t * q1, const p4est_quadrant_t * q2)
{
  int                 l;
  int                 comp = p4est_quadrant_compare (q1, q2);
  int                 level = (int) SC_MIN (q1->level, q2->level);
  p4est_lid_t         id1;
  p4est_lid_t         id2;
  p4est_quadrant_t    quad, par, anc;

  p4est_quadrant_linear_id_ext128 (q1, level, &id1);
  p4est_quadrant_linear_id_ext128 (q2, level, &id2);

  /* test linear id */
  if (p4est_quadrant_is_ancestor (q1, q2)) {
    SC_CHECK_ABORT (p4est_lid_is_equal (&id1, &id2)
                    && comp < 0, "Ancestor 1");
  }
  else if (p4est_quadrant_is_ancestor (q2, q1)) {
    SC_CHECK_ABORT (p4est_lid_is_equal (&id1, &id2)
                    && comp > 0, "Ancestor 2");
  }
  else {
    SC_CHECK_ABORT ((comp == 0 && p4est_lid_is_equal (&id1, &id2))
                    || (comp < 0 && (p4est_lid_compare (&id1, &id2) < 0))
                    || (comp > 0
                        && (p4est_lid_compare (&id1, &id2) > 0)), "compare");
  }

  /* test ancestor and parent functions */
  par = quad = *q1;
  for (l = quad.level - 1; l >= 0; --l) {
    p4est_quadrant_parent (&par, &par);
    p4est_quadrant_ancestor (&quad, l, &anc);
    SC_CHECK_ABORT (p4est_quadrant_is_equal (&par, &anc), "Ancestor test A");
  }
  par = quad = *q2;
  for (l = quad.level - 1; l >= 0; --l) {
    p4est_quadrant_parent (&par, &par);
    p4est_quadrant_ancestor (&quad, l, &anc);
    SC_CHECK_ABORT (p4est_quadrant_is_equal (&par, &anc), "Ancestor test B");
  }
}

static void
check_successor_predecessor (const p4est_quadrant_t * q)
{
  p4est_quadrant_t    temp1, temp2;
  p4est_lid_t         lid, successor_lid, one;

  p4est_quadrant_linear_id_ext128 (q, q->level, &lid);
  p4est_quadrant_successor (q, &temp1);
  p4est_quadrant_linear_id_ext128 (&temp1, q->level, &successor_lid);
  p4est_lid_set_one (&one);
  p4est_lid_add_inplace (&lid, &one);
  SC_CHECK_ABORT (p4est_lid_is_equal (&successor_lid, &lid), "successor");
  p4est_quadrant_predecessor (&temp1, &temp2);
  /* Check if predecessor inverts successor. */
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&temp2, q), "predecessor");
}

static void
check_predecessor_successor (const p4est_quadrant_t * q)
{
  p4est_quadrant_t    temp1, temp2;
  p4est_lid_t         lid, predecessor_lid, one;

  p4est_quadrant_linear_id_ext128 (q, q->level, &lid);
  p4est_quadrant_predecessor (q, &temp1);
  p4est_quadrant_linear_id_ext128 (&temp1, q->level, &predecessor_lid);
  p4est_lid_set_one (&one);
  p4est_lid_sub_inplace (&lid, &one);
  SC_CHECK_ABORT (p4est_lid_is_equal (&predecessor_lid, &lid), "predecessor");
  p4est_quadrant_successor (&temp1, &temp2);
  /* Check if successor inverts predecessor. */
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&temp2, q), "successor");
}

#define NEG_ONE_MAXL (~((((p4est_qcoord_t) 1) << P8EST_MAXLEVEL) - 1))
#define NEG_ONE_MAXLM1 (~((((p4est_qcoord_t) 1) << (P8EST_MAXLEVEL - 1)) - 1))
#define NEG_ONE_MAXLP1 \
  (NEG_ONE_MAXL & ~(((p4est_qcoord_t) 1) << P8EST_MAXLEVEL))

int
main (int argc, char **argv)
{
  const p4est_qcoord_t qone = 1;
  int                 mpiret;
  int                 k;
  int                 level, mid, cid;
  int                 id0, id1, id2, id3;
  int64_t             index1, index2;
  size_t              iz, jz, incount;
  p4est_qcoord_t      mh = P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est1;
  p4est_t            *p4est2;
  p4est_tree_t       *t1, *t2, tree;
  p4est_quadrant_t   *p, *q1, *q2;
  p4est_quadrant_t    r, s;
  p4est_quadrant_t    c0, c1, c2, c3, c4, c5, c6, c7;
  p4est_quadrant_t    cv[P4EST_CHILDREN], *cp[P4EST_CHILDREN];
  p4est_quadrant_t    A, B, C, D, E, F, G, H, I, P, Q, R, S;
  p4est_quadrant_t    a, f, g, h;
  uint64_t            Aid;
  p4est_lid_t         Fid;
  const int           indices[27] = { 0, 1, 2, 3, 4, 5, 6, 7,
    7, 9, 11, 13, 18, 19, 22, 23, 27, 31,
    36, 37, 38, 39, 45, 47, 54, 55, 63
  };
  p4est_lid_t         id;
  sc_rand_state_t     state;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* create connectivity and forest structures */
  connectivity = p8est_connectivity_new_unitcube ();
  p4est1 = p4est_new_ext (sc_MPI_COMM_SELF, connectivity, 15, 0, 0,
                          0, NULL, NULL);
  p4est2 = p4est_new_ext (sc_MPI_COMM_SELF, connectivity, 15, 0, 0,
                          8, NULL, NULL);

  /* refine the second tree to a uniform level */
  p4est_refine (p4est1, 1, refine_none, NULL);
  p4est_refine (p4est2, 1, refine_some, NULL);
  t1 = p4est_tree_array_index (p4est1->trees, 0);
  t2 = p4est_tree_array_index (p4est2->trees, 0);
  SC_CHECK_ABORT (p4est_tree_is_sorted (t1), "is_sorted");
  SC_CHECK_ABORT (p4est_tree_is_sorted (t2), "is_sorted");

  /* run a bunch of cross-tests */
  p = NULL;
  for (iz = 0; iz < t1->quadrants.elem_count; ++iz) {
    q1 = p4est_quadrant_array_index (&t1->quadrants, iz);

    /* test coordinates of all boundary objects */
    check_coordinates (q1);

    /* test the index conversion */
    index1 = p4est_quadrant_linear_id (q1, (int) q1->level);
    p4est_quadrant_set_morton (&r, (int) q1->level, index1);
    index2 = p4est_quadrant_linear_id (&r, (int) r.level);
    SC_CHECK_ABORT (index1 == index2, "index conversion");
    level = (int) q1->level - 1;
    if (level >= 0) {
      index1 = p4est_quadrant_linear_id (q1, level);
      p4est_quadrant_set_morton (&r, level, index1);
      index2 = p4est_quadrant_linear_id (&r, level);
      SC_CHECK_ABORT (index1 == index2, "index conversion");
    }

    /* test the is_next function */
    if (p != NULL) {
      SC_CHECK_ABORT (p4est_quadrant_is_next (p, q1), "is_next");
    }
    p = q1;

    /* test the is_family function */
    p8est_quadrant_children (q1, &c0, &c1, &c2, &c3, &c4, &c5, &c6, &c7);
    SC_CHECK_ABORT (p8est_quadrant_is_family
                    (&c0, &c1, &c2, &c3, &c4, &c5, &c6, &c7), "is_family");
    SC_CHECK_ABORT (!p8est_quadrant_is_family
                    (&c1, &c0, &c2, &c3, &c4, &c5, &c6, &c7), "is_family");
    SC_CHECK_ABORT (!p8est_quadrant_is_family
                    (&c0, &c1, &c2, &c3, &c4, &c5, &c5, &c7), "is_family");
    p4est_quadrant_childrenv (q1, cv);
    SC_CHECK_ABORT (p4est_quadrant_is_equal (&c0, &cv[0]), "is_family");
    SC_CHECK_ABORT (p4est_quadrant_is_equal (&c1, &cv[1]), "is_family");
    SC_CHECK_ABORT (p4est_quadrant_is_equal (&c2, &cv[2]), "is_family");
    SC_CHECK_ABORT (p4est_quadrant_is_equal (&c3, &cv[3]), "is_family");
    SC_CHECK_ABORT (p4est_quadrant_is_equal (&c4, &cv[4]), "is_family");
    SC_CHECK_ABORT (p4est_quadrant_is_equal (&c5, &cv[5]), "is_family");
    SC_CHECK_ABORT (p4est_quadrant_is_equal (&c6, &cv[6]), "is_family");
    SC_CHECK_ABORT (p4est_quadrant_is_equal (&c7, &cv[7]), "is_family");
    SC_CHECK_ABORT (p8est_quadrant_is_family (&cv[0], &cv[1], &cv[2], &cv[3],
                                              &cv[4], &cv[5], &cv[6], &cv[7]),
                    "is_family");
    cp[0] = &cv[0];
    cp[1] = &cv[1];
    cp[2] = &cv[2];
    cp[3] = &cv[3];
    cp[4] = &cv[4];
    cp[5] = &cv[5];
    cp[6] = &cv[6];
    cp[7] = &cv[7];
    SC_CHECK_ABORT (p4est_quadrant_is_familypv (cp), "is_family");
    cv[1] = cv[0];
    SC_CHECK_ABORT (!p4est_quadrant_is_familyv (cv), "is_family");
    cp[1] = &c1;
    SC_CHECK_ABORT (p4est_quadrant_is_familypv (cp), "is_family");
    cp[6] = &c7;
    SC_CHECK_ABORT (!p4est_quadrant_is_familypv (cp), "is_family");

    /* test the sibling function */
    mid = p4est_quadrant_child_id (q1);
    for (cid = 0; cid < P4EST_CHILDREN; ++cid) {
      p4est_quadrant_sibling (q1, &r, cid);
      if (cid != mid) {
        SC_CHECK_ABORT (p4est_quadrant_is_sibling (q1, &r), "sibling");
      }
      else {
        SC_CHECK_ABORT (p4est_quadrant_is_equal (q1, &r), "sibling");
      }
    }

    /* test t1 against itself */
    for (jz = 0; jz < t1->quadrants.elem_count; ++jz) {
      q2 = p4est_quadrant_array_index (&t1->quadrants, jz);

      /* test the comparison function */
      SC_CHECK_ABORT (p4est_quadrant_compare (q1, q2) ==
                      -p4est_quadrant_compare (q2, q1), "compare");
      SC_CHECK_ABORT ((!p4est_quadrant_compare (q1, q2)) ==
                      p4est_quadrant_is_equal (q1, q2), "is_equal");

      /* test the descriptive versions of functions */
      SC_CHECK_ABORT (p4est_quadrant_is_sibling_D (q1, q2) ==
                      p4est_quadrant_is_sibling (q1, q2), "is_sibling");
      SC_CHECK_ABORT (p4est_quadrant_is_parent_D (q1, q2) ==
                      p4est_quadrant_is_parent (q1, q2), "is_parent");
      SC_CHECK_ABORT (p4est_quadrant_is_parent_D (q2, q1) ==
                      p4est_quadrant_is_parent (q2, q1), "is_parent");
      SC_CHECK_ABORT (p4est_quadrant_is_ancestor_D (q1, q2) ==
                      p4est_quadrant_is_ancestor (q1, q2), "is_ancestor");
      SC_CHECK_ABORT (p4est_quadrant_is_ancestor_D (q2, q1) ==
                      p4est_quadrant_is_ancestor (q2, q1), "is_ancestor");
      SC_CHECK_ABORT (p4est_quadrant_is_next_D (q1, q2) ==
                      p4est_quadrant_is_next (q1, q2), "is_next");
      SC_CHECK_ABORT (p4est_quadrant_is_next_D (q2, q1) ==
                      p4est_quadrant_is_next (q2, q1), "is_next");
      p4est_nearest_common_ancestor_D (q1, q2, &r);
      p4est_nearest_common_ancestor (q1, q2, &s);
      SC_CHECK_ABORT (p4est_quadrant_is_equal (&r, &s), "common_ancestor");
      p4est_nearest_common_ancestor_D (q2, q1, &r);
      p4est_nearest_common_ancestor (q2, q1, &s);
      SC_CHECK_ABORT (p4est_quadrant_is_equal (&r, &s), "common_ancestor");
    }

    /* test t1 against t2 */
    for (jz = 0; jz < t2->quadrants.elem_count; ++jz) {
      q2 = p4est_quadrant_array_index (&t2->quadrants, jz);

      /* test the comparison function */
      SC_CHECK_ABORT (p4est_quadrant_compare (q1, q2) ==
                      -p4est_quadrant_compare (q2, q1), "compare");
      SC_CHECK_ABORT ((!p4est_quadrant_compare (q1, q2)) ==
                      p4est_quadrant_is_equal (q1, q2), "is_equal");

      /* test the descriptive versions of functions */
      SC_CHECK_ABORT (p4est_quadrant_is_sibling_D (q1, q2) ==
                      p4est_quadrant_is_sibling (q1, q2), "is_sibling");
      SC_CHECK_ABORT (p4est_quadrant_is_parent_D (q1, q2) ==
                      p4est_quadrant_is_parent (q1, q2), "is_parent");
      SC_CHECK_ABORT (p4est_quadrant_is_parent_D (q2, q1) ==
                      p4est_quadrant_is_parent (q2, q1), "is_parent");
      SC_CHECK_ABORT (p4est_quadrant_is_ancestor_D (q1, q2) ==
                      p4est_quadrant_is_ancestor (q1, q2), "is_ancestor");
      SC_CHECK_ABORT (p4est_quadrant_is_ancestor_D (q2, q1) ==
                      p4est_quadrant_is_ancestor (q2, q1), "is_ancestor");
      SC_CHECK_ABORT (p4est_quadrant_is_next_D (q1, q2) ==
                      p4est_quadrant_is_next (q1, q2), "is_next");
      SC_CHECK_ABORT (p4est_quadrant_is_next_D (q2, q1) ==
                      p4est_quadrant_is_next (q2, q1), "is_next");
      p4est_nearest_common_ancestor_D (q1, q2, &r);
      p4est_nearest_common_ancestor (q1, q2, &s);
      SC_CHECK_ABORT (p4est_quadrant_is_equal (&r, &s), "common_ancestor");
      p4est_nearest_common_ancestor_D (q2, q1, &r);
      p4est_nearest_common_ancestor (q2, q1, &s);
      SC_CHECK_ABORT (p4est_quadrant_is_equal (&r, &s), "common_ancestor");
    }
  }

  p = NULL;
  for (iz = 0; iz < t2->quadrants.elem_count; ++iz) {
    q1 = p4est_quadrant_array_index (&t2->quadrants, iz);

    /* test the is_next function */
    if (p != NULL) {
      SC_CHECK_ABORT (p4est_quadrant_is_next (p, q1), "is_next");
    }
    p = q1;
  }

  /* test the coarsen function */
  p4est_coarsen (p4est1, 1, coarsen_none, NULL);
  p4est_coarsen (p4est1, 1, coarsen_all, NULL);
  p4est_coarsen (p4est2, 1, coarsen_some, NULL);

  /* test the linearize algorithm */
  incount = t2->quadrants.elem_count;
  (void) p4est_linearize_tree (p4est2, t2);
  SC_CHECK_ABORT (incount == t2->quadrants.elem_count, "linearize");

  /* this is user_data neutral only when p4est1->data_size == 0 */
  sc_array_init (&tree.quadrants, sizeof (p4est_quadrant_t));
  sc_array_resize (&tree.quadrants, 18);
  q1 = p4est_quadrant_array_index (&tree.quadrants, 0);
  q2 = p4est_quadrant_array_index (&t2->quadrants, 0);
  *q1 = *q2;
  q2 = p4est_quadrant_array_index (&t2->quadrants, 1);
  for (k = 0; k < 3; ++k) {
    q1 = p4est_quadrant_array_index (&tree.quadrants, (size_t) (k + 1));
    *q1 = *q2;
    q1->level = (int8_t) (q1->level + k);
  }
  for (k = 0; k < 10; ++k) {
    q1 = p4est_quadrant_array_index (&tree.quadrants, (size_t) (k + 4));
    q2 = p4est_quadrant_array_index (&t2->quadrants, (size_t) (k + 3));
    *q1 = *q2;
    q1->level = (int8_t) (q1->level + k);
  }
  for (k = 0; k < 4; ++k) {
    q1 = p4est_quadrant_array_index (&tree.quadrants, (size_t) (k + 14));
    q2 = p4est_quadrant_array_index (&t2->quadrants, (size_t) (k + 12));
    *q1 = *q2;
    q1->level = (int8_t) (q1->level + 10 + k);
  }
  tree.maxlevel = 0;
  for (k = 0; k <= P4EST_QMAXLEVEL; ++k) {
    tree.quadrants_per_level[k] = 0;
  }
  for (; k <= P4EST_MAXLEVEL; ++k) {
    tree.quadrants_per_level[k] = -1;
  }
  incount = tree.quadrants.elem_count;
  for (iz = 0; iz < incount; ++iz) {
    q1 = p4est_quadrant_array_index (&tree.quadrants, iz);
    ++tree.quadrants_per_level[q1->level];
    tree.maxlevel = (int8_t) SC_MAX (tree.maxlevel, q1->level);
  }
  SC_CHECK_ABORT (!p4est_tree_is_linear (&tree), "is_linear");
  (void) p4est_linearize_tree (p4est1, &tree);
  SC_CHECK_ABORT (incount - 3 == tree.quadrants.elem_count, "linearize");
  sc_array_reset (&tree.quadrants);

  /* create a partial tree and check overlap */
  sc_array_resize (&tree.quadrants, 4);
  q1 = p4est_quadrant_array_index (&tree.quadrants, 0);
  p4est_quadrant_set_morton (q1, 3, 191);
  q1 = p4est_quadrant_array_index (&tree.quadrants, 1);
  p4est_quadrant_set_morton (q1, 1, 3);
  q1 = p4est_quadrant_array_index (&tree.quadrants, 2);
  p4est_quadrant_set_morton (q1, 2, 32);
  q1 = p4est_quadrant_array_index (&tree.quadrants, 3);
  p4est_quadrant_set_morton (q1, 2, 33);
  for (k = 0; k <= P4EST_QMAXLEVEL; ++k) {
    tree.quadrants_per_level[k] = 0;
  }
  for (; k <= P4EST_MAXLEVEL; ++k) {
    tree.quadrants_per_level[k] = -1;
  }
  tree.quadrants_per_level[1] = 1;
  tree.quadrants_per_level[2] = 2;
  tree.quadrants_per_level[3] = 1;
  tree.maxlevel = 3;
  p4est_quadrant_first_descendant (p4est_quadrant_array_index
                                   (&tree.quadrants, 0), &tree.first_desc,
                                   P4EST_QMAXLEVEL);
  p4est_quadrant_last_descendant (p4est_quadrant_array_index
                                  (&tree.quadrants,
                                   tree.quadrants.elem_count - 1),
                                  &tree.last_desc, P4EST_QMAXLEVEL);
  SC_CHECK_ABORT (p4est_tree_is_complete (&tree), "is_complete");

  p4est_quadrant_set_morton (&D, 0, 0);
  SC_CHECK_ABORT (p4est_quadrant_overlaps_tree (&tree, &D), "overlaps 0");

  p4est_quadrant_set_morton (&A, 1, 0);
  SC_CHECK_ABORT (!p4est_quadrant_overlaps_tree (&tree, &A), "overlaps 1");
  p4est_quadrant_set_morton (&A, 1, 2);
  SC_CHECK_ABORT (p4est_quadrant_overlaps_tree (&tree, &A), "overlaps 2");
  p4est_quadrant_set_morton (&A, 1, 3);
  SC_CHECK_ABORT (p4est_quadrant_overlaps_tree (&tree, &A), "overlaps 3");
  p4est_quadrant_set_morton (&A, 1, 4);
  SC_CHECK_ABORT (p4est_quadrant_overlaps_tree (&tree, &A), "overlaps 4");
  p4est_quadrant_set_morton (&A, 1, 5);
  SC_CHECK_ABORT (!p4est_quadrant_overlaps_tree (&tree, &A), "overlaps 5");

  p4est_quadrant_set_morton (&B, 3, 13);
  SC_CHECK_ABORT (!p4est_quadrant_overlaps_tree (&tree, &B), "overlaps 6");
  p4est_quadrant_set_morton (&B, 3, 191);
  SC_CHECK_ABORT (p4est_quadrant_overlaps_tree (&tree, &B), "overlaps 7");
  p4est_quadrant_set_morton (&B, 3, 271);
  SC_CHECK_ABORT (p4est_quadrant_overlaps_tree (&tree, &B), "overlaps 8");
  p4est_quadrant_set_morton (&B, 3, 272);
  SC_CHECK_ABORT (!p4est_quadrant_overlaps_tree (&tree, &B), "overlaps 9");

  p4est_quadrant_set_morton (&C, 4, 2175);
  SC_CHECK_ABORT (p4est_quadrant_overlaps_tree (&tree, &C), "overlaps 10");
  p4est_quadrant_set_morton (&C, 4, 2176);
  SC_CHECK_ABORT (!p4est_quadrant_overlaps_tree (&tree, &C), "overlaps 11");

  sc_array_reset (&tree.quadrants);

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est1);
  p4est_destroy (p4est2);
  p4est_connectivity_destroy (connectivity);

  /* This will test the ability to address negative quadrants */
  P4EST_QUADRANT_INIT (&A);
  P4EST_QUADRANT_INIT (&B);
  P4EST_QUADRANT_INIT (&C);
  P4EST_QUADRANT_INIT (&D);
  P4EST_QUADRANT_INIT (&E);
  P4EST_QUADRANT_INIT (&F);
  P4EST_QUADRANT_INIT (&G);
  P4EST_QUADRANT_INIT (&H);
  P4EST_QUADRANT_INIT (&I);
  P4EST_QUADRANT_INIT (&P);
  P4EST_QUADRANT_INIT (&Q);
  P4EST_QUADRANT_INIT (&R);

  A.x = NEG_ONE_MAXL;
  A.y = NEG_ONE_MAXL;
  A.z = 0;
  A.level = 0;

  B.x = qone << P4EST_MAXLEVEL;
  B.y = NEG_ONE_MAXL;
  B.z = 0;
  B.level = 0;

  C.x = NEG_ONE_MAXL;
  C.y = qone << P4EST_MAXLEVEL;
  C.z = 0;
  C.level = 0;

  D.x = qone << P4EST_MAXLEVEL;
  D.y = qone << P4EST_MAXLEVEL;
  D.z = 0;
  D.level = 0;

  /* this one is outside the 3x3 box */
  E.x = NEG_ONE_MAXLP1;
  E.y = -qone;
  E.z = -qone;
  E.level = 0;

  F.x = P4EST_ROOT_LEN + (P4EST_ROOT_LEN - mh);
  F.y = P4EST_ROOT_LEN + (P4EST_ROOT_LEN - mh);
  F.z = NEG_ONE_MAXL;
  F.level = P4EST_QMAXLEVEL;

  G.x = -mh;
  G.y = -mh;
  G.z = -mh;
  G.level = P4EST_QMAXLEVEL;

  H.x = NEG_ONE_MAXLM1;
  H.y = NEG_ONE_MAXLM1;
  H.z = qone << (P4EST_MAXLEVEL - 1);
  H.level = 1;

  I.x = NEG_ONE_MAXL;
  I.y = NEG_ONE_MAXLM1;
  I.z = P4EST_ROOT_LEN + (P4EST_ROOT_LEN - mh);
  I.level = P4EST_QMAXLEVEL;

  P.x = NEG_ONE_MAXL;
  P.y = NEG_ONE_MAXLM1;
  P.z = qone << (P4EST_MAXLEVEL - 1);
  P.level = 1;

  Q.x = -2 * mh;
  Q.y = -2 * mh;
  Q.z = (qone << P4EST_MAXLEVEL) - 2 * mh;
  Q.level = P4EST_QMAXLEVEL - 1;

  R.x = 2;
  R.y = 2;
  R.z = 2;
  R.level = P4EST_QMAXLEVEL;
  check_coordinates (&R);

  p4est_lid_set_zero (&id);
  p4est_lid_set_bit (&id, 70);
  p4est_quadrant_set_morton_ext128 (&S, P4EST_QMAXLEVEL, &id);

  p4est_quadrant_srand (&R, &state);
  SC_CHECK_ABORT ((uint64_t) state == 7, "quadrant_srand");
  p4est_quadrant_srand (&S, &state);
  SC_CHECK_ABORT ((uint64_t) state == (id.high_bits ^ id.low_bits),
                  "quadrant_srand");

  SC_CHECK_ABORT (p4est_quadrant_compare (&B, &F) < 0, "Comp 1");
  SC_CHECK_ABORT (p4est_quadrant_compare (&A, &G) < 0, "Comp 2");
  SC_CHECK_ABORT (p4est_quadrant_compare (&F, &G) < 0, "Comp 3");
  SC_CHECK_ABORT (p4est_quadrant_compare (&A, &I) < 0, "Comp 4");
  SC_CHECK_ABORT (p4est_quadrant_compare (&D, &C) < 0, "Comp 5");
  SC_CHECK_ABORT (p4est_quadrant_compare (&B, &G) < 0, "Comp 6");
  SC_CHECK_ABORT (p4est_quadrant_compare (&G, &G) == 0, "Comp 7");

  check_linear_id (&A, &A);
  check_linear_id (&A, &B);
  check_linear_id (&A, &C);
  check_linear_id (&A, &D);
  /* check_linear_id (&A, &E); */
  check_linear_id (&A, &F);
  check_linear_id (&A, &G);
  check_linear_id (&A, &H);
  check_linear_id (&A, &I);

  check_linear_id (&B, &A);
  check_linear_id (&B, &B);
  check_linear_id (&B, &C);
  check_linear_id (&B, &D);
  /* check_linear_id (&B, &E); */
  check_linear_id (&B, &F);
  check_linear_id (&B, &G);
  check_linear_id (&B, &H);
  check_linear_id (&B, &I);

  check_linear_id (&D, &A);
  check_linear_id (&D, &B);
  check_linear_id (&D, &C);
  check_linear_id (&D, &D);
  /* check_linear_id (&D, &E); */
  check_linear_id (&D, &F);
  check_linear_id (&D, &G);
  check_linear_id (&D, &H);
  check_linear_id (&D, &I);

  check_linear_id (&G, &A);
  check_linear_id (&G, &B);
  check_linear_id (&G, &C);
  check_linear_id (&G, &D);
  /* check_linear_id (&G, &E); */
  check_linear_id (&G, &F);
  check_linear_id (&G, &G);
  check_linear_id (&G, &H);
  check_linear_id (&G, &I);

  check_linear_id (&I, &A);
  check_linear_id (&I, &B);
  check_linear_id (&I, &C);
  check_linear_id (&I, &D);
  /* check_linear_id (&I, &E); */
  check_linear_id (&I, &F);
  check_linear_id (&I, &G);
  check_linear_id (&I, &H);
  check_linear_id (&I, &I);

  check_linear_id (&P, &F);
  check_linear_id (&P, &G);
  check_linear_id (&P, &H);
  check_linear_id (&P, &Q);

  check_linear_id (&Q, &F);
  check_linear_id (&Q, &B);
  check_linear_id (&Q, &H);
  check_linear_id (&Q, &I);

  SC_CHECK_ABORT (p4est_quadrant_is_extended (&A), "is_extended A");
  SC_CHECK_ABORT (p4est_quadrant_is_extended (&B), "is_extended B");
  SC_CHECK_ABORT (p4est_quadrant_is_extended (&C), "is_extended C");
  SC_CHECK_ABORT (p4est_quadrant_is_extended (&D), "is_extended D");
  SC_CHECK_ABORT (!p4est_quadrant_is_extended (&E), "!is_extended E");
  SC_CHECK_ABORT (p4est_quadrant_is_extended (&F), "is_extended F");
  SC_CHECK_ABORT (p4est_quadrant_is_extended (&G), "is_extended G");

  SC_CHECK_ABORT (!p4est_quadrant_compare (&A, &A), "compare");
  SC_CHECK_ABORT (p4est_quadrant_compare (&A, &B) > 0, "compare");
  SC_CHECK_ABORT (p4est_quadrant_compare (&B, &A) < 0, "compare");

  SC_CHECK_ABORT (!p4est_quadrant_compare (&F, &F), "compare");
  SC_CHECK_ABORT (p4est_quadrant_compare (&G, &F) > 0, "compare");
  SC_CHECK_ABORT (p4est_quadrant_compare (&F, &G) < 0, "compare");

  A.p.piggy1.which_tree = 0;
  B.p.piggy2.which_tree = 0;
  SC_CHECK_ABORT (!p4est_quadrant_compare_piggy (&A, &A), "compare_piggy");
  SC_CHECK_ABORT (p4est_quadrant_compare_piggy (&A, &B) > 0, "compare_piggy");
  SC_CHECK_ABORT (p4est_quadrant_compare_piggy (&B, &A) < 0, "compare_piggy");

  F.p.which_tree = 0;
  G.p.piggy1.which_tree = 0;
  SC_CHECK_ABORT (!p4est_quadrant_compare_piggy (&F, &F), "compare_piggy");
  SC_CHECK_ABORT (p4est_quadrant_compare_piggy (&G, &F) > 0, "compare_piggy");
  SC_CHECK_ABORT (p4est_quadrant_compare_piggy (&F, &G) < 0, "compare_piggy");

  F.p.piggy2.which_tree = (p4est_topidx_t) P4EST_TOPIDX_MAX - 3;
  G.p.which_tree = (p4est_topidx_t) P4EST_TOPIDX_MAX / 2;
  SC_CHECK_ABORT (!p4est_quadrant_compare_piggy (&F, &F), "compare_piggy");
  SC_CHECK_ABORT (p4est_quadrant_compare_piggy (&G, &F) < 0, "compare_piggy");
  SC_CHECK_ABORT (p4est_quadrant_compare_piggy (&F, &G) > 0, "compare_piggy");

  SC_CHECK_ABORT (p4est_quadrant_is_equal (&A, &A), "is_equal");
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&F, &F), "is_equal");
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&G, &G), "is_equal");

  SC_CHECK_ABORT (p4est_quadrant_is_sibling (&P, &H), "is_sibling");
  SC_CHECK_ABORT (!p4est_quadrant_is_sibling (&A, &H), "is_sibling");
  SC_CHECK_ABORT (p4est_quadrant_is_sibling_D (&P, &H), "is_sibling_D");
  SC_CHECK_ABORT (!p4est_quadrant_is_sibling_D (&A, &H), "is_sibling_D");

  SC_CHECK_ABORT (p4est_quadrant_is_parent (&A, &H), "is_parent");
  SC_CHECK_ABORT (!p4est_quadrant_is_parent (&H, &A), "is_parent");
  SC_CHECK_ABORT (!p4est_quadrant_is_parent (&A, &Q), "is_parent");
  SC_CHECK_ABORT (p4est_quadrant_is_parent_D (&A, &H), "is_parent_D");

  SC_CHECK_ABORT (p4est_quadrant_is_ancestor (&A, &Q), "is_ancestor");
  SC_CHECK_ABORT (!p4est_quadrant_is_ancestor (&A, &A), "is_ancestor");

  SC_CHECK_ABORT (p4est_quadrant_is_ancestor_D (&A, &P), "is_ancestor_D");
  SC_CHECK_ABORT (!p4est_quadrant_is_ancestor_D (&G, &G), "is_ancestor_D");

  /* SC_CHECK_ABORT (p4est_quadrant_is_next (&F, &E), "is_next"); */
  SC_CHECK_ABORT (!p4est_quadrant_is_next (&A, &H), "is_next");
  /* SC_CHECK_ABORT (p4est_quadrant_is_next_D (&F, &E), "is_next_D"); */
  SC_CHECK_ABORT (!p4est_quadrant_is_next_D (&A, &H), "is_next_D");

  p4est_quadrant_parent (&H, &a);
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&A, &a), "parent");

  p4est_quadrant_sibling (&P, &h, 7);
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&H, &h), "sibling");

  p8est_quadrant_children (&A, &c0, &c1, &c2, &c3, &c4, &c5, &c6, &c7);
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&c6, &P), "children");
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&c7, &H), "children");
  SC_CHECK_ABORT (!p4est_quadrant_is_equal (&c7, &Q), "children");

  SC_CHECK_ABORT (p8est_quadrant_is_family (&c0, &c1, &c2, &c3,
                                            &c4, &c5, &c6, &c7), "is_family");
  id0 = p4est_quadrant_child_id (&c0);
  id1 = p4est_quadrant_child_id (&c1);
  id2 = p4est_quadrant_child_id (&c2);
  id3 = p4est_quadrant_child_id (&c6);
  SC_CHECK_ABORT (id0 == 0 && id1 == 1 && id2 == 2 && id3 == 6, "child_id");
  SC_CHECK_ABORT (p4est_quadrant_child_id (&G) == 7, "child_id");

  p4est_quadrant_first_descendant (&A, &c1, 1);
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&c0, &c1), "first_descendant");

  p4est_quadrant_last_descendant (&A, &g, P4EST_QMAXLEVEL - 1);
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&Q, &g), "last_descendant");

  p4est_quadrant_linear_id_ext128 (&F, P4EST_QMAXLEVEL, &Fid);
  p4est_quadrant_set_morton_ext128 (&f, P4EST_QMAXLEVEL, &Fid);
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&F, &f), "set_morton/linear_id");

  Aid = p4est_quadrant_linear_id (&A, 0);
  p4est_quadrant_set_morton (&a, 0, Aid);
  SC_CHECK_ABORT (Aid == 27, "linear_id");
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&A, &a), "set_morton/linear_id");

  check_successor_predecessor (&F);
  check_predecessor_successor (&F);
  check_predecessor_successor (&G);
  check_predecessor_successor (&H);
  check_successor_predecessor (&I);
  check_predecessor_successor (&I);
  check_successor_predecessor (&P);
  check_predecessor_successor (&P);
  check_predecessor_successor (&Q);

  p4est_nearest_common_ancestor (&P, &H, &a);
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&A, &a), "ancestor");

  p4est_nearest_common_ancestor_D (&P, &Q, &a);
  SC_CHECK_ABORT (p4est_quadrant_is_equal (&A, &a), "ancestor_D");

  for (k = 0; k < 27; ++k) {
    p4est_quadrant_set_morton (&E, 0, (uint64_t) indices[k]);
  }
  p4est_quadrant_set_morton (&P, 0, 54);
  p4est_quadrant_set_morton (&Q, 0, 55);
  SC_CHECK_ABORT (p4est_quadrant_is_next (&P, &Q), "is_next");
  SC_CHECK_ABORT (!p4est_quadrant_is_next (&A, &Q), "is_next");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
