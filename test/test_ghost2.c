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
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#else
#include <p8est_bits.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#endif

#ifndef P4_TO_P8
static int          refine_level = 5;
#else
static int          refine_level = 4;
#endif

/* #undef P4EST_TEST_CHATTY */

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
  if (quadrant->x == P4EST_QUADRANT_LEN (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->y >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

#define TEST_EXCHANGE_MAGIC 427482.18e-13

typedef struct test_exchange
{
  p4est_gloidx_t      gi;
  long long           ll;
  double              magic;
}
test_exchange_t;

static p4est_ghost_exchange_t *
test_exchange_begin (p4est_t * p4est, p4est_ghost_t * ghost)
{
  size_t              zz;
  p4est_topidx_t      nt;
  p4est_gloidx_t      gnum;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  void              **ghost_void_data;

  /* Test begin/end: p4est data size is 0, transfer what's in the user_data void* */

  p4est_reset_data (p4est, 0, NULL, NULL);
  gnum = p4est->global_first_quadrant[p4est->mpirank];
  for (nt = p4est->first_local_tree; nt <= p4est->last_local_tree; ++nt) {
    tree = p4est_tree_array_index (p4est->trees, nt);
    for (zz = 0; zz < tree->quadrants.elem_count; ++gnum, ++zz) {
      q = p4est_quadrant_array_index (&tree->quadrants, zz);
      q->p.user_long = (long) (3 * gnum + 17);
    }
  }
  P4EST_ASSERT (gnum == p4est->global_first_quadrant[p4est->mpirank + 1]);

  /* allocate data to receive the messages */
  ghost_void_data = P4EST_ALLOC (void *, ghost->ghosts.elem_count);

  /* the ghost data are safely packed into the return type */
  return p4est_ghost_exchange_data_begin (p4est, ghost, ghost_void_data);
}

static void
test_exchange_end (p4est_ghost_exchange_t * exc)
{
  void              **ghost_void_data = (void **) exc->ghost_data;
  p4est_t            *p4est = exc->p4est;
  p4est_ghost_t      *ghost = exc->ghost;
  int                 p;
  p4est_locidx_t      gexcl, gincl, gl;
  p4est_gloidx_t      gnum;
  p4est_quadrant_t   *q;

  /* free this structure on the inside */
  p4est_ghost_exchange_data_end (exc);

  /* verify results of ghost exchange */
  gexcl = 0;
  for (p = 0; p < p4est->mpisize; ++p) {
    gincl = ghost->proc_offsets[p + 1];
    gnum = p4est->global_first_quadrant[p];
#ifdef P4EST_TEST_CHATTY
    P4EST_LDEBUGF ("In test begin/end for %d with %d %d\n", p, gexcl, gincl);
#endif
    for (gl = gexcl; gl < gincl; ++gl) {
      q = p4est_quadrant_array_index (&ghost->ghosts, gl);
      SC_CHECK_ABORT (3 * (gnum + (p4est_gloidx_t) q->p.piggy3.local_num) + 17
                      == (p4est_gloidx_t) ghost_void_data[gl],
                      "Ghost exchange mismatch begin/end");
    }
    gexcl = gincl;
  }
  P4EST_ASSERT (gexcl == (p4est_locidx_t) ghost->ghosts.elem_count);
  P4EST_FREE (ghost_void_data);
}

static void
test_exchange_A (p4est_t * p4est, p4est_ghost_t * ghost)
{
  int                 p;
  size_t              zz;
  p4est_topidx_t      nt;
  p4est_locidx_t      gexcl, gincl, gl;
  p4est_gloidx_t      gnum;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  void              **ghost_void_data;

  /* Test A: p4est data size is 0, transfer what's in the user_data void* */

  p4est_reset_data (p4est, 0, NULL, NULL);
  gnum = p4est->global_first_quadrant[p4est->mpirank];
  for (nt = p4est->first_local_tree; nt <= p4est->last_local_tree; ++nt) {
    tree = p4est_tree_array_index (p4est->trees, nt);
    for (zz = 0; zz < tree->quadrants.elem_count; ++gnum, ++zz) {
      q = p4est_quadrant_array_index (&tree->quadrants, zz);
      q->p.user_long = (long) gnum;
    }
  }
  P4EST_ASSERT (gnum == p4est->global_first_quadrant[p4est->mpirank + 1]);

  ghost_void_data = P4EST_ALLOC (void *, ghost->ghosts.elem_count);
  p4est_ghost_exchange_data (p4est, ghost, ghost_void_data);

  gexcl = 0;
  for (p = 0; p < p4est->mpisize; ++p) {
    gincl = ghost->proc_offsets[p + 1];
    gnum = p4est->global_first_quadrant[p];
#ifdef P4EST_TEST_CHATTY
    P4EST_LDEBUGF ("In test A for %d with %d %d\n", p, gexcl, gincl);
#endif
    for (gl = gexcl; gl < gincl; ++gl) {
      q = p4est_quadrant_array_index (&ghost->ghosts, gl);
      SC_CHECK_ABORT (gnum + (p4est_gloidx_t) q->p.piggy3.local_num ==
                      (p4est_gloidx_t) ghost_void_data[gl],
                      "Ghost exchange mismatch A");
    }
    gexcl = gincl;
  }
  P4EST_ASSERT (gexcl == (p4est_locidx_t) ghost->ghosts.elem_count);
  P4EST_FREE (ghost_void_data);
}

static void
test_exchange_B (p4est_t * p4est, p4est_ghost_t * ghost)
{
  int                 p;
  size_t              zz;
  p4est_topidx_t      nt;
  p4est_locidx_t      gexcl, gincl, gl;
  p4est_gloidx_t      gnum;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  test_exchange_t    *ghost_struct_data, *e;

  /* Test B: p4est data size is > 0, transfer what's in user_data */

  p4est_reset_data (p4est, sizeof (test_exchange_t), NULL, NULL);
  gnum = p4est->global_first_quadrant[p4est->mpirank];
  for (nt = p4est->first_local_tree; nt <= p4est->last_local_tree; ++nt) {
    tree = p4est_tree_array_index (p4est->trees, nt);
    for (zz = 0; zz < tree->quadrants.elem_count; ++gnum, ++zz) {
      q = p4est_quadrant_array_index (&tree->quadrants, zz);
      e = (test_exchange_t *) q->p.user_data;
      e->gi = gnum;
      e->ll = (long) gnum;
      e->magic = TEST_EXCHANGE_MAGIC;
    }
  }
  P4EST_ASSERT (gnum == p4est->global_first_quadrant[p4est->mpirank + 1]);

  ghost_struct_data = P4EST_ALLOC (test_exchange_t, ghost->ghosts.elem_count);
  p4est_ghost_exchange_data (p4est, ghost, ghost_struct_data);

  gexcl = 0;
  for (p = 0; p < p4est->mpisize; ++p) {
    gincl = ghost->proc_offsets[p + 1];
    gnum = p4est->global_first_quadrant[p];
#ifdef P4EST_TEST_CHATTY
    P4EST_LDEBUGF ("In test B for %d with %d %d\n", p, gexcl, gincl);
#endif
    for (gl = gexcl; gl < gincl; ++gl) {
      q = p4est_quadrant_array_index (&ghost->ghosts, gl);
      e = ghost_struct_data + gl;
      SC_CHECK_ABORT (gnum + (p4est_gloidx_t) q->p.piggy3.local_num ==
                      e->gi, "Ghost exchange mismatch B1");
      SC_CHECK_ABORT (gnum + (p4est_gloidx_t) q->p.piggy3.local_num ==
                      (p4est_gloidx_t) e->ll, "Ghost exchange mismatch B2");
      SC_CHECK_ABORT (e->magic == TEST_EXCHANGE_MAGIC,
                      "Ghost exchange mismatch B3");
    }
    gexcl = gincl;
  }
  P4EST_ASSERT (gexcl == (p4est_locidx_t) ghost->ghosts.elem_count);
  P4EST_FREE (ghost_struct_data);
}

static void
test_exchange_C (p4est_t * p4est, p4est_ghost_t * ghost)
{
  int                 p;
  size_t              zz;
  p4est_locidx_t      gexcl, gincl, gl;
  p4est_gloidx_t      gnum;
  p4est_quadrant_t   *q;
  void              **mirror_data;
  test_exchange_t    *mirror_struct_data;
  test_exchange_t    *ghost_struct_data, *e;

  /* Test C: don't use p4est user_data at all */

  mirror_struct_data =
    P4EST_ALLOC (test_exchange_t, ghost->mirrors.elem_count);
  mirror_data = P4EST_ALLOC (void *, ghost->mirrors.elem_count);
  for (zz = 0; zz < ghost->mirrors.elem_count; ++zz) {
    q = p4est_quadrant_array_index (&ghost->mirrors, zz);
    gnum = p4est->global_first_quadrant[p4est->mpirank] +
      (p4est_gloidx_t) q->p.piggy3.local_num;
    mirror_data[zz] = e = mirror_struct_data + zz;
    e->gi = gnum;
    e->ll = (long) gnum;
    e->magic = TEST_EXCHANGE_MAGIC;
  }

  ghost_struct_data = P4EST_ALLOC (test_exchange_t, ghost->ghosts.elem_count);
  p4est_ghost_exchange_custom (p4est, ghost, sizeof (test_exchange_t),
                               mirror_data, ghost_struct_data);

  P4EST_FREE (mirror_data);
  P4EST_FREE (mirror_struct_data);

  gexcl = 0;
  for (p = 0; p < p4est->mpisize; ++p) {
    gincl = ghost->proc_offsets[p + 1];
    gnum = p4est->global_first_quadrant[p];
#ifdef P4EST_TEST_CHATTY
    P4EST_LDEBUGF ("In test C for %d with %d %d\n", p, gexcl, gincl);
#endif
    for (gl = gexcl; gl < gincl; ++gl) {
      q = p4est_quadrant_array_index (&ghost->ghosts, gl);
      e = ghost_struct_data + gl;
      SC_CHECK_ABORT (gnum + (p4est_gloidx_t) q->p.piggy3.local_num ==
                      e->gi, "Ghost exchange mismatch C1");
      SC_CHECK_ABORT (gnum + (p4est_gloidx_t) q->p.piggy3.local_num ==
                      (p4est_gloidx_t) e->ll, "Ghost exchange mismatch C2");
      SC_CHECK_ABORT (e->magic == TEST_EXCHANGE_MAGIC,
                      "Ghost exchange mismatch C3");
    }
    gexcl = gincl;
  }
  P4EST_ASSERT (gexcl == (p4est_locidx_t) ghost->ghosts.elem_count);
  P4EST_FREE (ghost_struct_data);
}

static void
test_exchange_D (p4est_t * p4est, p4est_ghost_t * ghost)
{
  const int           exchange_minlevel = 1;
  const int           exchange_maxlevel = refine_level - 1;
  int                 p;
  size_t              zz;
  p4est_locidx_t      gexcl, gincl, gl;
  p4est_gloidx_t      gnum;
  p4est_quadrant_t   *q;
  void              **mirror_data;
  test_exchange_t    *mirror_struct_data;
  test_exchange_t    *ghost_struct_data, *e;

  /* Test C: don't use p4est user_data at all */

  mirror_struct_data =
    P4EST_ALLOC (test_exchange_t, ghost->mirrors.elem_count);
  mirror_data = P4EST_ALLOC (void *, ghost->mirrors.elem_count);
  for (zz = 0; zz < ghost->mirrors.elem_count; ++zz) {
    q = p4est_quadrant_array_index (&ghost->mirrors, zz);
    gnum = p4est->global_first_quadrant[p4est->mpirank] +
      (p4est_gloidx_t) q->p.piggy3.local_num;
    mirror_data[zz] = e = mirror_struct_data + zz;
    e->gi = gnum;
    e->ll = (long) gnum;
    e->magic = TEST_EXCHANGE_MAGIC;
  }

  ghost_struct_data = P4EST_ALLOC (test_exchange_t, ghost->ghosts.elem_count);
  p4est_ghost_exchange_custom_levels (p4est, ghost,
                                      exchange_minlevel, exchange_maxlevel,
                                      sizeof (test_exchange_t),
                                      mirror_data, ghost_struct_data);

  P4EST_FREE (mirror_data);
  P4EST_FREE (mirror_struct_data);

  gexcl = 0;
  for (p = 0; p < p4est->mpisize; ++p) {
    gincl = ghost->proc_offsets[p + 1];
    gnum = p4est->global_first_quadrant[p];
#ifdef P4EST_TEST_CHATTY
    P4EST_LDEBUGF ("In test D for %d with %d %d\n", p, gexcl, gincl);
#endif
    for (gl = gexcl; gl < gincl; ++gl) {
      q = p4est_quadrant_array_index (&ghost->ghosts, gl);
      if (exchange_minlevel <= (int) q->level &&
          (int) q->level <= exchange_maxlevel) {
        e = ghost_struct_data + gl;
        SC_CHECK_ABORT (gnum + (p4est_gloidx_t) q->p.piggy3.local_num ==
                        e->gi, "Ghost exchange mismatch D1");
        SC_CHECK_ABORT (gnum + (p4est_gloidx_t) q->p.piggy3.local_num ==
                        (p4est_gloidx_t) e->ll, "Ghost exchange mismatch D2");
        SC_CHECK_ABORT (e->magic == TEST_EXCHANGE_MAGIC,
                        "Ghost exchange mismatch D3");
      }
    }
    gexcl = gincl;
  }
  P4EST_ASSERT (gexcl == (p4est_locidx_t) ghost->ghosts.elem_count);
  P4EST_FREE (ghost_struct_data);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t      *ghost;
  p4est_ghost_exchange_t *exc;
  int                 num_cycles = 2;
  int                 i;
  p4est_lnodes_t     *lnodes;
  int                 type;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

#ifndef P4_TO_P8
  conn = p4est_connectivity_new_moebius ();
#else
  conn = p8est_connectivity_new_rotcubes ();
#endif

  p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

  /* refine to make the number of elements interesting */
  p4est_refine (p4est, 1, refine_fn, NULL);

  /* balance the forest */
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  /* do a uniform partition */
  p4est_partition (p4est, 0, NULL);

  /* create the ghost layer */
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

  /* test ghost data exchange */
  test_exchange_A (p4est, ghost);
  test_exchange_B (p4est, ghost);
  test_exchange_C (p4est, ghost);
  test_exchange_D (p4est, ghost);

  for (i = 0; i < num_cycles; i++) {
    /* expand and test that the ghost layer can still exchange data properly
     * */
    p4est_ghost_expand (p4est, ghost);
    test_exchange_A (p4est, ghost);
    test_exchange_B (p4est, ghost);
    test_exchange_C (p4est, ghost);
    test_exchange_D (p4est, ghost);
  }

  p4est_ghost_destroy (ghost);
  /* repeate the cyle, but with lnodes */
  /* create the ghost layer */
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  type = p4est_connect_type_int (ghost->btype);
  lnodes = p4est_lnodes_new (p4est, ghost, -type);
  p4est_ghost_support_lnodes (p4est, lnodes, ghost);
  /* test ghost data exchange */
  test_exchange_A (p4est, ghost);
  test_exchange_B (p4est, ghost);
  test_exchange_C (p4est, ghost);
  test_exchange_D (p4est, ghost);

  for (i = 0; i < num_cycles; i++) {
    /* expand and test that the ghost layer can still exchange data properly
     * */
    p4est_ghost_expand_by_lnodes (p4est, lnodes, ghost);
    exc = test_exchange_begin (p4est, ghost);
    test_exchange_A (p4est, ghost);
    test_exchange_B (p4est, ghost);
    test_exchange_C (p4est, ghost);
    test_exchange_D (p4est, ghost);
    test_exchange_end (exc);
  }

  /* clean up */
  p4est_lnodes_destroy (lnodes);
  p4est_ghost_destroy (ghost);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  /* exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
