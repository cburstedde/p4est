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
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#else
#include <p8est_algorithms.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#endif

typedef struct
{
  p4est_topidx_t      a;
  int64_t             sum;
}
user_data_t;

static int          weight_counter;
static int          weight_index;

static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
  user_data_t        *data = (user_data_t *) quadrant->p.user_data;

  /* prevent uninitialized bytes due to compiler padding */
  memset (data, -1, sizeof (user_data_t));

  data->a = which_tree;
  data->sum = quadrant->x + quadrant->y + quadrant->level;
}

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  if (quadrant->level >= 6) {
    return 0;
  }
#ifdef P4_TO_P8
  if (quadrant->level >= 5 && quadrant->z <= P4EST_QUADRANT_LEN (3)) {
    return 0;
  }
#endif

  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->x >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

static int
weight_one (p4est_t * p4est, p4est_topidx_t which_tree,
            p4est_quadrant_t * quadrant)
{
  return 1;
}

static int
weight_once (p4est_t * p4est, p4est_topidx_t which_tree,
             p4est_quadrant_t * quadrant)
{
  if (weight_counter++ == weight_index) {
    return 1;
  }

  return 0;
}

static int
traverse_fn (p4est_t * p4est, p4est_topidx_t which_tree,
             p4est_quadrant_t * quadrant, int pfirst, int plast, void *point)
{
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (0 <= which_tree &&
                which_tree < p4est->connectivity->num_trees);
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (0 <= pfirst && pfirst <= plast && plast < p4est->mpisize);

  return 1;
}

static int          circle_count;

static void
circle_init (p4est_t * p4est, p4est_topidx_t which_tree,
             p4est_quadrant_t * quadrant)
{
  int                *idata = (int *) quadrant->p.user_data;

  *idata = ++circle_count;
}

typedef struct test_transfer
{
  p4est_t            *p4est;
  p4est_t            *back;
}
test_transfer_t;

static test_transfer_t *
test_transfer_pre (p4est_t * p4est)
{
  test_transfer_t    *tt;

  tt = P4EST_ALLOC (test_transfer_t, 1);
  tt->p4est = p4est;
  tt->back = p4est_copy (tt->p4est, 1);

  return tt;
}

static void
test_transfer_post (test_transfer_t * tt, p4est_t * p4est)
{
  size_t              pds, gds, data_size;
  int                 i;
  int                *dest_sizes;
  int                *dest_vdata;
  int                *src_sizes;
  int                *src_vdata;
  int                *ti;
  char               *dest_data;
  char               *src_data;
  char               *td, *cmp;
  size_t              zz;
  size_t              vz, vcountd, vcounts;
  p4est_topidx_t      tid;
  p4est_locidx_t      li;
  p4est_gloidx_t      tog;
  p4est_t            *back;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  p4est_transfer_context_t *tf;

  P4EST_ASSERT (tt != NULL);
  P4EST_ASSERT (tt->p4est == p4est);
  back = tt->back;
  P4EST_ASSERT (p4est->data_size == back->data_size);

  /* now back is a copy of the p4est before partition */
  /* p4est has been partitioned once */

  /* put together some buffers */
  gds = sizeof (p4est_gloidx_t);
  pds = SC_MAX (p4est->data_size, gds);
  data_size = gds + pds;
  cmp = P4EST_ALLOC (char, data_size);
  dest_data = P4EST_ALLOC (char, data_size * p4est->local_num_quadrants);
  dest_sizes = P4EST_ALLOC (int, p4est->local_num_quadrants);
  src_data = P4EST_ALLOC_ZERO (char, data_size * back->local_num_quadrants);
  src_sizes = P4EST_ALLOC (int, back->local_num_quadrants);

  /* assemble data to send that includes the user_data */
  td = src_data;
  li = 0;
  vcounts = 0;
  tog = back->global_first_quadrant[back->mpirank];
  for (tid = back->first_local_tree; tid <= back->last_local_tree; ++tid) {
    tree = p4est_tree_array_index (back->trees, tid);
    for (zz = 0; zz < tree->quadrants.elem_count; ++zz) {
      vcounts += vz = tog % 4;
      quad = p4est_quadrant_array_index (&tree->quadrants, zz);
      *(p4est_gloidx_t *) td = tog++;
      td += gds;
      memcpy (td, quad->p.user_data, p4est->data_size);
      td += pds;
      src_sizes[li++] = vz * sizeof (int);
    }
  }
  P4EST_ASSERT (li == back->local_num_quadrants);
  P4EST_ASSERT (td - src_data ==
                (ptrdiff_t) (data_size * back->local_num_quadrants));
  P4EST_ASSERT (tog == back->global_first_quadrant[back->mpirank + 1]);

  /* do data transfer part I */
  tf = p4est_transfer_fixed_begin (p4est->global_first_quadrant,
                                   back->global_first_quadrant,
                                   p4est->mpicomm, 0,
                                   dest_data, src_data, data_size);
  p4est_transfer_fixed (p4est->global_first_quadrant,
                        back->global_first_quadrant, p4est->mpicomm, 1,
                        dest_sizes, src_sizes, sizeof (int));
  p4est_transfer_fixed_end (tf);

  /* we verify the fixed data we have sent */
  td = dest_data;
  li = 0;
  vcountd = 0;
  tog = p4est->global_first_quadrant[p4est->mpirank];
  for (tid = p4est->first_local_tree; tid <= p4est->last_local_tree; ++tid) {
    tree = p4est_tree_array_index (p4est->trees, tid);
    for (zz = 0; zz < tree->quadrants.elem_count; ++zz) {
      vcountd += vz = tog % 4;
      quad = p4est_quadrant_array_index (&tree->quadrants, zz);
      SC_CHECK_ABORT (*(p4est_gloidx_t *) td == tog,
                      "Transfer index mismatch");
      td += gds;
      SC_CHECK_ABORT (!memcmp (td, quad->p.user_data, p4est->data_size),
                      "Transfer data mismatch");
      td += pds;
      ++tog;
      SC_CHECK_ABORT (dest_sizes[li++] == (int) (vz * sizeof (int)),
                      "Transfer size mismatch");
    }
  }
  P4EST_ASSERT (li == p4est->local_num_quadrants);
  P4EST_ASSERT (td - dest_data ==
                (ptrdiff_t) (data_size * p4est->local_num_quadrants));
  P4EST_ASSERT (tog == p4est->global_first_quadrant[p4est->mpirank + 1]);

  /* allocate space for variable data sizes */
  ti = src_vdata = P4EST_ALLOC (int, vcounts * sizeof (int));
  for (li = 0; li < back->local_num_quadrants; ++li) {
    for (i = 0; i < src_sizes[li] / (int) sizeof (int); ++i) {
      *ti++ = i;
    }
  }
  P4EST_ASSERT (ti - src_vdata == (ptrdiff_t) vcounts);
  dest_vdata = P4EST_ALLOC (int, vcountd * sizeof (int));

  /* do data transfer part II */
  p4est_transfer_custom (p4est->global_first_quadrant,
                         back->global_first_quadrant, p4est->mpicomm, 1,
                         dest_vdata, dest_sizes, src_vdata, src_sizes);

  /* we verify the variable data we have sent */
  ti = dest_vdata;
  for (li = 0; li < p4est->local_num_quadrants; ++li) {
    for (i = 0; i < dest_sizes[li] / (int) sizeof (int); ++i) {
      SC_CHECK_ABORT (*ti == i, "Transfer variable mismatch");
      ++ti;
    }
  }
  P4EST_ASSERT (ti - dest_vdata == (ptrdiff_t) vcountd);

  /* cleanup memory */
  P4EST_FREE (dest_data);
  P4EST_FREE (dest_vdata);
  P4EST_FREE (dest_sizes);
  P4EST_FREE (src_data);
  P4EST_FREE (src_vdata);
  P4EST_FREE (src_sizes);
  P4EST_FREE (cmp);

  /* cleanup context */
  p4est_destroy (tt->back);
  P4EST_FREE (tt);
}

static void
test_pertree (p4est_t * p4est, const p4est_gloidx_t * prev_pertree,
              p4est_gloidx_t * new_pertree)
{
  const p4est_topidx_t num_trees = p4est->connectivity->num_trees;
  p4est_gloidx_t     *pertree;

  /* test counting of quadrants in individual trees */
  P4EST_ASSERT ((size_t) num_trees == p4est->trees->elem_count);
  if (new_pertree == NULL) {
    pertree = P4EST_ALLOC (p4est_gloidx_t, num_trees + 1);
  }
  else {
    pertree = new_pertree;
  }
  p4est_comm_count_pertree (p4est, pertree);
  SC_CHECK_ABORT (pertree[num_trees] == p4est->global_num_quadrants,
                  "pertree check failed");
  if (prev_pertree != NULL) {
    SC_CHECK_ABORT (!memcmp (pertree, prev_pertree,
                             sizeof (p4est_gloidx_t) * (num_trees + 1)),
                    "pertree now different");
  }
  if (new_pertree == NULL) {
    P4EST_FREE (pertree);
  }

  /* test traversal routine */
  p4est_search_partition (p4est, 1, traverse_fn, NULL, NULL);
}

static unsigned
test_checksum (p4est_t * p4est, int have_zlib)
{
  return have_zlib ? p4est_checksum (p4est) : 0;
}

static void
test_partition_circle (sc_MPI_Comm mpicomm,
                       p4est_connectivity_t * connectivity,
                       p4est_gloidx_t * pertree1, p4est_gloidx_t * pertree2,
                       int have_zlib)
{
  int                 i, j;
  int                 num_procs;
  int                 empty_proc1, empty_proc2;
  unsigned            crc1, crc2;
  p4est_gloidx_t      global_num;
  p4est_locidx_t     *new_counts;
  p4est_t            *p4est, *copy;
  test_transfer_t    *tt;

  /* Create a forest and make a copy */

  circle_count = 0;
  p4est = p4est_new_ext (mpicomm, connectivity, 0, 3, 1,
                         sizeof (int), circle_init, NULL);
  num_procs = p4est->mpisize;
  test_pertree (p4est, NULL, pertree1);

  global_num = p4est->global_num_quadrants;
  crc1 = test_checksum (p4est, have_zlib);
  copy = p4est_copy (p4est, 1);
  P4EST_ASSERT (test_checksum (copy, have_zlib) == crc1);

  new_counts = P4EST_ALLOC (p4est_locidx_t, num_procs);

  /* Partition with one empty processor */
  if (num_procs > 1) {
    P4EST_GLOBAL_INFO ("First circle partition\n");
    empty_proc1 = num_procs / 3;
    j = 0;
    for (i = 0; i < num_procs; ++i) {
      if (i == empty_proc1) {
        new_counts[i] = 0;
      }
      else {
        new_counts[i] =
          p4est_partition_cut_gloidx (global_num, j + 1, num_procs - 1) -
          p4est_partition_cut_gloidx (global_num, j, num_procs - 1);
        P4EST_ASSERT (new_counts[i] >= 0);
        ++j;
      }
    }
    P4EST_ASSERT (j == num_procs - 1);
    tt = test_transfer_pre (p4est);
    p4est_partition_given (p4est, new_counts);
    test_transfer_post (tt, p4est);
    test_pertree (p4est, pertree1, pertree2);
    crc2 = test_checksum (p4est, have_zlib);
    SC_CHECK_ABORT (crc1 == crc2, "First checksum mismatch");
  }

  /* Partition with two empty processors */
  if (num_procs > 2) {
    P4EST_GLOBAL_INFO ("Second circle partition\n");
    empty_proc1 = (2 * num_procs) / 3 - 2;
    empty_proc2 = (2 * num_procs) / 3;
    j = 0;
    for (i = 0; i < num_procs; ++i) {
      if (i == empty_proc1 || i == empty_proc2) {
        new_counts[i] = 0;
      }
      else {
        new_counts[i] =
          p4est_partition_cut_gloidx (global_num, j + 1, num_procs - 2) -
          p4est_partition_cut_gloidx (global_num, j, num_procs - 2);
        P4EST_ASSERT (new_counts[i] >= 0);
        ++j;
      }
    }
    P4EST_ASSERT (j == num_procs - 2);
    tt = test_transfer_pre (p4est);
    p4est_partition_given (p4est, new_counts);
    test_transfer_post (tt, p4est);
    test_pertree (p4est, pertree1, pertree2);
    crc2 = test_checksum (p4est, have_zlib);
    SC_CHECK_ABORT (crc1 == crc2, "Second checksum mismatch");
  }

  /* Uniform partition */
  P4EST_GLOBAL_INFO ("Third circle partition\n");
  tt = test_transfer_pre (p4est);
  p4est_partition (p4est, 0, NULL);
  test_transfer_post (tt, p4est);
  test_pertree (p4est, pertree1, pertree2);
  crc2 = test_checksum (p4est, have_zlib);
  SC_CHECK_ABORT (crc1 == crc2, "Third checksum mismatch");
  SC_CHECK_ABORT (p4est_is_equal (p4est, copy, 1), "Forest mismatch");

  P4EST_FREE (new_counts);
  p4est_destroy (copy);
  p4est_destroy (p4est);
}

int
main (int argc, char **argv)
{
  int                 rank;
  int                 num_procs;
  int                 mpiret;
  int                 have_zlib;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est, *copy;
  p4est_connectivity_t *connectivity;
  int                 i;
  p4est_topidx_t      t;
  size_t              qz;
  p4est_locidx_t      num_quadrants_on_last;
  p4est_locidx_t     *num_quadrants_in_proc;
  p4est_gloidx_t     *pertree1, *pertree2;
  p4est_quadrant_t   *quad;
  p4est_tree_t       *tree;
  user_data_t        *user_data;
  int64_t             sum;
  unsigned            crc;
  test_transfer_t    *tt;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* establish parallel logging */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* check for ZLIB usability */
  if (!(have_zlib = p4est_have_zlib ())) {
    P4EST_GLOBAL_LERROR
      ("Not found a working ZLIB installation: ignoring CRCs\n");
  }

  /* create connectivity and forest structures */
#ifdef P4_TO_P8
  connectivity = p8est_connectivity_new_twocubes ();
#else
  connectivity = p4est_connectivity_new_corner ();
#endif
  p4est = p4est_new_ext (mpicomm, connectivity, 15, 0, 0,
                         sizeof (user_data_t), init_fn, NULL);

  pertree1 = P4EST_ALLOC (p4est_gloidx_t, p4est->connectivity->num_trees + 1);
  pertree2 = P4EST_ALLOC (p4est_gloidx_t, p4est->connectivity->num_trees + 1);
  num_procs = p4est->mpisize;
  num_quadrants_in_proc = P4EST_ALLOC (p4est_locidx_t, num_procs);

  /* refine and balance to make the number of elements interesting */
  test_pertree (p4est, NULL, pertree1);
  p4est_refine (p4est, 1, refine_fn, init_fn);
  test_pertree (p4est, NULL, pertree1);

  /* Set an arbitrary partition.
   *
   * Since this is just a test we assume the global number of
   * quadrants will fit in an int32_t
   */
  num_quadrants_on_last = (p4est_locidx_t) p4est->global_num_quadrants;
  for (i = 0; i < num_procs - 1; ++i) {
    num_quadrants_in_proc[i] = (p4est_locidx_t) i + 1;  /* type ok */
    num_quadrants_on_last -= (p4est_locidx_t) i + 1;    /* type ok */
  }
  num_quadrants_in_proc[num_procs - 1] = num_quadrants_on_last;
  SC_CHECK_ABORT (num_quadrants_on_last > 0,
                  "Negative number of quadrants on the last processor");

  /* Save a checksum of the original forest */
  crc = test_checksum (p4est, have_zlib);

  /* partition the forest */
  tt = test_transfer_pre (p4est);
  (void) p4est_partition_given (p4est, num_quadrants_in_proc);
  test_transfer_post (tt, p4est);
  test_pertree (p4est, pertree1, pertree2);

  /* Double check that we didn't loose any quads */
  SC_CHECK_ABORT (crc == test_checksum (p4est, have_zlib),
                  "bad checksum, missing a quad");

  /* count the actual number of quadrants per proc */
  SC_CHECK_ABORT (num_quadrants_in_proc[rank]
                  == p4est->local_num_quadrants,
                  "partition failed, wrong number of quadrants");

  /* check user data content */
  for (t = p4est->first_local_tree; t <= p4est->last_local_tree; ++t) {
    tree = p4est_tree_array_index (p4est->trees, t);
    for (qz = 0; qz < tree->quadrants.elem_count; ++qz) {
      quad = p4est_quadrant_array_index (&tree->quadrants, qz);
      user_data = (user_data_t *) quad->p.user_data;
      sum = quad->x + quad->y + quad->level;

      SC_CHECK_ABORT (user_data->a == t, "bad user_data, a");
      SC_CHECK_ABORT (user_data->sum == sum, "bad user_data, sum");
    }
  }

  /* do a weighted partition with uniform weights */
  tt = test_transfer_pre (p4est);
  p4est_partition (p4est, 0, weight_one);
  test_transfer_post (tt, p4est);
  test_pertree (p4est, pertree1, pertree2);
  SC_CHECK_ABORT (crc == test_checksum (p4est, have_zlib),
                  "bad checksum after uniformly weighted partition");

  /* copy the p4est */
  copy = p4est_copy (p4est, 1);
  SC_CHECK_ABORT (crc == test_checksum (copy, have_zlib),
                  "bad checksum after copy");

  /* do a weighted partition with many zero weights */
  weight_counter = 0;
  weight_index = (rank == 1) ? 1342 : 0;
  tt = test_transfer_pre (copy);
  p4est_partition (copy, 0, weight_once);
  test_transfer_post (tt, copy);
  test_pertree (copy, pertree1, pertree2);
  SC_CHECK_ABORT (crc == test_checksum (copy, have_zlib),
                  "bad checksum after unevenly weighted partition 1");

  /* do a weighted partition with many zero weights */
  weight_counter = 0;
  weight_index = 0;
  tt = test_transfer_pre (copy);
  p4est_partition (copy, 0, weight_once);
  test_transfer_post (tt, copy);
  test_pertree (copy, pertree1, pertree2);
  SC_CHECK_ABORT (crc == test_checksum (copy, have_zlib),
                  "bad checksum after unevenly weighted partition 2");

  /* do a weighted partition with many zero weights
   *
   * Since this is just a test we assume the local number of
   * quadrants will fit in an int
   */
  weight_counter = 0;
  weight_index =
    (rank == num_procs - 1) ? ((int) copy->local_num_quadrants - 1) : 0;
  tt = test_transfer_pre (copy);
  p4est_partition (copy, 0, weight_once);
  test_transfer_post (tt, copy);
  test_pertree (copy, pertree1, pertree2);
  SC_CHECK_ABORT (crc == test_checksum (copy, have_zlib),
                  "bad checksum after unevenly weighted partition 3");

  /* check user data content */
  for (t = copy->first_local_tree; t <= copy->last_local_tree; ++t) {
    tree = p4est_tree_array_index (copy->trees, t);
    for (qz = 0; qz < tree->quadrants.elem_count; ++qz) {
      quad = p4est_quadrant_array_index (&tree->quadrants, qz);
      user_data = (user_data_t *) quad->p.user_data;
      sum = quad->x + quad->y + quad->level;

      SC_CHECK_ABORT (user_data->a == t, "bad user_data, a");
      SC_CHECK_ABORT (user_data->sum == sum, "bad user_data, sum");
    }
  }

  /* Add another test.  Overwrites pertree1, pertree2 */
  test_partition_circle (mpicomm, connectivity, pertree1, pertree2,
                         have_zlib);

  /* clean up and exit */
  P4EST_FREE (pertree1);
  P4EST_FREE (pertree2);
  P4EST_FREE (num_quadrants_in_proc);
  p4est_destroy (p4est);
  p4est_destroy (copy);
  p4est_connectivity_destroy (connectivity);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
