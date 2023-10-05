/*
  This file is part of p4est, version 3.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2019 individual authors
  Originally written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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
#include <p4est3_quadrant_yx.h>
#include <p4est3_p4est.h>
#else
#include <p4est3_quadrant_zyx.h>
#include <p4est3_p8est.h>
#endif

#define N_QUADS_2_TEST 6000000

typedef struct timeavx2
{
  int                 mpirank;
  p4est3_locidx       n_quads;

  sc3_allocator_t    *alloc;
  sc3_array_t        *qarr;
  sc3_array_t        *qarr_avx;
  const p4est3_quadrant_vtable_t *qvt;
  const p4est3_quadrant_vtable_t *qvt_avx;
}
timeavx2_t;

static sc3_error_t *
testavx2_prepare (timeavx2_t * t, int *avx_is_success)
{
  void               *p;
  t->n_quads = N_QUADS_2_TEST;

  SC3E_RETVAL (avx_is_success, 0);
  SC3A_CHECK (t != NULL);
  SC3A_CHECK (t->n_quads > 0);

  /* the standard p4est2 virtual table always exists */
  SC3E (p4est3_quadrant_vtable_p4est (&t->qvt));

  /* the AVX virtual table can only be set with hardware support */
  SC3E (p4est3_quadrant_yx_vtable (&t->qvt_avx));
  if (t->qvt_avx == NULL) {
    return NULL;
  }
  SC3E_DEMAND (t->qvt_avx != NULL, "AVX is not supported by hardware "
              "or p4est is not build neither in 2D nor 3D");

  /* create a toplevel allocator */
  SC3E (sc3_allocator_new (sc3_allocator_nocount (), &t->alloc));
  SC3E (sc3_allocator_setup (t->alloc));

  /* allocate quadrant arrays */
  SC3E (p4est3_quadrant_array_new (t->alloc, t->qvt, t->n_quads, &t->qarr));
  SC3E (p4est3_quadrant_array_new (t->alloc,
                                   t->qvt_avx, t->n_quads, &t->qarr_avx));

  /* initialize first element */
  SC3E (sc3_array_index (t->qarr, 0, &p));
  SC3E (p4est3_quadrant_root (t->qvt, p));
  SC3E (sc3_array_index (t->qarr_avx, 0, &p));
  SC3E (p4est3_quadrant_root (t->qvt_avx, p));

  /* clean and successful return */
  *avx_is_success = 1;
  return NULL;
}

sc3_error_t        *
create_quadrants (p4est3_quadrant_vtable_t * qvt_avx,
                  p4est3_quadrant_vtable_t * qvt, sc3_array_t * v,
                  sc3_array_t * q, int n_quad)
{
  int                 quad, put_ind, child;
  void               *in, *out;

  for (quad = 0, put_ind = 1; P4EST_CHILDREN * (quad + 1) < n_quad;
       ++quad, put_ind += P4EST_CHILDREN) {
    for (child = 0; child < P4EST_CHILDREN; ++child) {
      SC3E (sc3_array_index (v, quad, &in));
      SC3E (sc3_array_index (v, put_ind + child, &out));
      SC3E (p4est3_quadrant_child (qvt_avx, in, child, out));

      SC3E (sc3_array_index (q, quad, &in));
      SC3E (sc3_array_index (q, put_ind + child, &out));
      SC3E (p4est3_quadrant_child (qvt, in, child, out));
    }
  }

  for (child = 0; child < n_quad - put_ind; ++child) {
    SC3E (sc3_array_index (v, quad, &in));
    SC3E (sc3_array_index (v, put_ind + child, &out));
    SC3E (p4est3_quadrant_child (qvt_avx, in, child, out));

    SC3E (sc3_array_index (q, quad, &in));
    SC3E (sc3_array_index (q, put_ind + child, &out));
    SC3E (p4est3_quadrant_child (qvt, in, child, out));
  }
  return NULL;
}

sc3_error_t        *
is_equal_child (p4est3_quadrant_vtable_t * qvt_avx,
                p4est3_quadrant_vtable_t * qvt, sc3_array_t * v,
                sc3_array_t * q, int32_t n_quad)
{
  int32_t             is_equal, i;
  int                 xy_avx[P4EST_DIM], xy[P4EST_DIM], l;
  void               *p;
  void               *v_p;
  p4est_quadrant_t   *q_row;

  SC3E (create_quadrants (qvt_avx, qvt, v, q, n_quad));

  SC3E (sc3_array_index (q, 0, &p));
  q_row = (p4est_quadrant_t *) p;

  for (i = 0; i < n_quad; ++i) {
    SC3E (sc3_array_index (v, i, &v_p));
    SC3E (p4est3_quadrant_level (qvt_avx, v_p, &l));
    SC3E (p4est3_quadrant_coordinates (qvt_avx, v_p, xy_avx));
    SC3E (p4est3_quadrant_coordinates (qvt, &q_row[i], xy));
    is_equal =
      (int32_t) q_row[i].level == l && xy[0] == xy_avx[0] &&
      xy[1] == xy_avx[1] &&
#ifdef P4_TO_P8
      xy[2] == xy_avx[2] &&
#endif
      1;
    SC3E_DEMAND (is_equal, "Comparing of quadrant_child results");
  }
  return NULL;
}

sc3_error_t        *
is_equal_parent (p4est3_quadrant_vtable_t * qvt_avx,
                 p4est3_quadrant_vtable_t * qvt, sc3_array_t * v,
                 sc3_array_t * q, int n_quad)
{
  int32_t             is_equal, i;
  int                 xy_avx[P4EST_DIM], xy[P4EST_DIM], l;
  void               *v_p;
  p4est_quadrant_t   *q_p;
  void               *in, *p;

  SC3E (sc3_array_index (v, 0, &v_p));

  SC3E (sc3_array_index (q, 0, &p));
  q_p = (p4est_quadrant_t *) p;

  for (i = 1; i < n_quad; ++i) {
    SC3E (sc3_array_index (v, i, &in));
    SC3E (p4est3_quadrant_parent (qvt_avx, in, v_p));

    SC3E (sc3_array_index (q, i, &in));
    SC3E (p4est3_quadrant_parent (qvt, in, q_p));

    SC3E (p4est3_quadrant_level (qvt_avx, v_p, &l));
    SC3E (p4est3_quadrant_coordinates (qvt_avx, v_p, xy_avx));
    SC3E (p4est3_quadrant_coordinates (qvt, q_p, xy));
    is_equal =
      (int32_t) q_p->level == l && xy[0] == xy_avx[0] && xy[1] == xy_avx[1] &&
#ifdef P4_TO_P8
      xy[2] == xy_avx[2] &&
#endif
      1;
    SC3E_DEMAND (is_equal, "Comparing of quadrant_parent results");
  }
  return NULL;
}

sc3_error_t        *
is_equal_compare (p4est3_quadrant_vtable_t * qvt_avx,
                  p4est3_quadrant_vtable_t * qvt, sc3_array_t * v,
                  sc3_array_t * q, int n_quad)
{
  int                 res, res_avx, i;
  void               *lhs, *rhs;
  for (i = 0; i < n_quad; ++i) {
    SC3E (sc3_array_index (v, i, &lhs));
    SC3E (sc3_array_index (v, n_quad - i - 1, &rhs));
    SC3E (p4est3_quadrant_compare (qvt_avx, lhs, rhs, &res_avx));

    SC3E (sc3_array_index (q, i, &lhs));
    SC3E (sc3_array_index (q, n_quad - i - 1, &rhs));
    SC3E (p4est3_quadrant_compare (qvt, lhs, rhs, &res));
    SC3E_DEMAND (res_avx == res, "Comparing of quadrant_compare results");
  }
  return NULL;
}

sc3_error_t        *
is_equal_successor (p4est3_quadrant_vtable_t * qvt_avx,
                    p4est3_quadrant_vtable_t * qvt, sc3_array_t * v,
                    sc3_array_t * q, int n_quad)
{
  int32_t             is_equal, i, j;
  int                 xy_avx[P4EST_DIM], xy[P4EST_DIM], l;
  void               *v_p;
  p4est_quadrant_t   *q_row;
  void               *in, *p;

  SC3E (sc3_array_index (v, 0, &v_p));

  SC3E (sc3_array_index (q, 0, &p));
  q_row = (p4est_quadrant_t *) p;

  for (i = 1; i < n_quad; i += P4EST_CHILDREN) {
    for (j = 0; j < P4EST_CHILDREN - 1 && i + j < n_quad; ++j) {
      SC3E (sc3_array_index (v, i + j, &in));
      SC3E (p4est3_quadrant_successor (qvt_avx, in, v_p));

      SC3E (sc3_array_index (q, i + j, &in));
      SC3E (p4est3_quadrant_successor (qvt, in, q_row));

      SC3E (p4est3_quadrant_level (qvt_avx, v_p, &l));
      SC3E (p4est3_quadrant_coordinates (qvt_avx, v_p, xy_avx));
      SC3E (p4est3_quadrant_coordinates (qvt, &q_row[i], xy));
      is_equal =
        (int32_t) q_row[i].level == l && xy[0] == xy_avx[0] &&
        xy[1] == xy_avx[1] &&
#ifdef P4_TO_P8
        xy[2] == xy_avx[2] &&
#endif
        1;
      if (!is_equal) {
        return NULL;
      }
      SC3E_DEMAND (is_equal, "Comparing of quadrant_successor results");
    }
  }
  return NULL;
}

static sc3_error_t *
testavx2_compare (timeavx2_t * t)
{
  SC3A_CHECK (t != NULL);
  SC3A_CHECK (t->n_quads > 0);

  SC3X (is_equal_child
        (t->qvt_avx, t->qvt, t->qarr_avx, t->qarr, t->n_quads));
  SC3X (is_equal_parent
        (t->qvt_avx, t->qvt, t->qarr_avx, t->qarr, t->n_quads));
  SC3X (is_equal_compare
        (t->qvt_avx, t->qvt, t->qarr_avx, t->qarr, t->n_quads));
  SC3X (is_equal_successor
        (t->qvt_avx, t->qvt, t->qarr_avx, t->qarr, t->n_quads));
  return NULL;
}

static sc3_error_t *
timeavx2_cleanup (timeavx2_t * t)
{
  SC3A_CHECK (t != NULL);
  SC3A_CHECK (t->n_quads > 0);

  SC3E (sc3_array_destroy (&t->qarr_avx));
  SC3E (sc3_array_destroy (&t->qarr));
  SC3E (sc3_allocator_destroy (&t->alloc));
  return NULL;
}

int
main (int argc, char **argv)
{
  int                 avx_is_success = 0;
  timeavx2_t          st, *t = &st;

  SC3X (sc3_MPI_Init (&argc, &argv));
  SC3X (sc3_MPI_Comm_rank (SC3_MPI_COMM_WORLD, &t->mpirank));

  testavx2_prepare (t, &avx_is_success);

  if (avx_is_success) {
    if (t->mpirank == 0) {
      SC3X (testavx2_compare (t));
    }
    SC3X (timeavx2_cleanup (t));
  }

  SC3X (sc3_MPI_Finalize ());
  return 0;
}
