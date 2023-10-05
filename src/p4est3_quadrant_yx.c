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
#include <p4est.h>
#include <p4est3_quadrant_yx.h>
#else
#include <p8est.h>
#include <p4est3_quadrant_zyx.h>
#endif /* !P4_TO_P8 */

/* These have the same values for the 2D and 3D implementations */
#define P4EST3_YX_MAXLEVEL 30
#define P4EST3_YX_QMAXLEVEL 30
/** The length of a side of the root quadrant */
#define P4EST3_YX_ROOT_LEN ((int32_t) 1 << P4EST3_YX_MAXLEVEL)
/** The length of a quadrant of level l */
#define P4EST3_YX_QUADRANT_LEN(l) ((int32_t) 1 << (P4EST3_YX_MAXLEVEL - (l)))
#ifdef P4EST_ENABLE_AVX2

#include <immintrin.h>
#include <smmintrin.h>
#include <emmintrin.h>

static              p4est3_gloidx
p4est3_quadrant_zyx_num_uniform (int level)
{
  if (level < 0) {
    return -1;
  }
  return p4est3_glopow (P4EST_CHILDREN, level);
}

static int
p4est3_quadrant_zyx_is_valid (const __m128i * q, char *reason)
{
  int32_t             level = _mm_extract_epi32 (*q, 0);
  SC3E_TEST (level >= 0, reason);
  SC3E_TEST (level <= P4EST3_YX_MAXLEVEL, reason);
/* *INDENT-OFF* */
  SC3E_TEST (
    _mm_testz_si128 (*q, _mm_set_epi32 (P4EST3_YX_QUADRANT_LEN (level) - 1
                                      , P4EST3_YX_QUADRANT_LEN (level) - 1
                                      , P4EST3_YX_QUADRANT_LEN (level) - 1
                                      , 0)) == 1
  , reason
  );
/* *INDENT-ON* */
  SC3E_YES (reason);
}

static sc3_error_t *
p4est3_quadrant_zyx_is_parent (const __m128i * q, const __m128i * r, int *j)
{
  int32_t             r_level = _mm_extract_epi32 (*r, 0);

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  SC3E_RETVAL (j, 0);

  if (_mm_extract_epi32 (*q, 0) + 1 != r_level) {
    return NULL;
  }

/* *INDENT-OFF* */
  __m128i             lhs = _mm_add_epi32 (*q, _mm_set_epi32 (0, 0, 0, 1));
  __m128i             rhs = _mm_and_si128 (*r
                            , _mm_set_epi32 ( /* Optimized by compiler? */
                                             ~P4EST3_YX_QUADRANT_LEN (r_level)
                                           , ~P4EST3_YX_QUADRANT_LEN (r_level)
#ifdef P4_TO_P8
                                           , ~P4EST3_YX_QUADRANT_LEN (r_level)
#else
                                           , 0
#endif
                                           , 0xFFFFFFFF)
                            );
/* *INDENT-ON* */
  if (_mm_testc_si128 (rhs, lhs) == 1) {
    *j = 1;
  }
  return NULL;
}

#ifdef P4EST_ENABLE_DEBUG

static int
p4est3_quadrant_zyx_is_parent_internal (const __m128i * q, const __m128i * r,
                                        char *reason)
{
  int             j;
  SC3E_DO (p4est3_quadrant_zyx_is_parent (q, r, &j), reason);
  SC3E_TEST (j, reason);
  SC3E_YES (reason);
}

#endif

static sc3_error_t *
p4est3_quadrant_zyx_is_ancestor (const __m128i * q, const __m128i * r, int *j)
{
  int                 lq, lr;
  __m128i             exclor;

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);

  lq = _mm_extract_epi32 (*q, 0);
  lr = _mm_extract_epi32 (*r, 0);
  if (lq >= lr) {
    *j = 0;
    return NULL;
  }
/* *INDENT-OFF* */
  exclor = _mm_srlv_epi32 (
             _mm_xor_si128 (*q, *r)
           , _mm_set_epi32 (P4EST3_YX_MAXLEVEL - lq
                          , P4EST3_YX_MAXLEVEL - lq
                          , P4EST3_YX_MAXLEVEL - lq
                          , 32)
           );
/* *INDENT-ON* */
  *j = _mm_test_all_zeros (exclor, _mm_set1_epi32 (0xFFFFFFFF));
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_get_tree_boundary (const __m128i * q, int face, int *j)
{
  const int32_t       l = _mm_extract_epi32 (*q, 0);
  const int           direction = face / 2;
  int32_t             coord = 0, bound;
  switch (direction) {
  case 0:
    coord = _mm_extract_epi32 (*q, 3);
    break;
  case 1:
    coord = _mm_extract_epi32 (*q, 2);
    break;
#ifdef P4_TO_P8
  case 2:
    coord = _mm_extract_epi32 (*q, 1);
    break;
#endif
  default:
    break;
  }

  bound = face % 2 == 0 ? 0 : P4EST3_YX_ROOT_LEN - P4EST3_YX_QUADRANT_LEN (l);
  *j = (coord == bound);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_tree_boundaries (const __m128i * q, sc3_array_t * nf)
{
  SC3A_CHECK (q != NULL);
  SC3A_CHECK (nf != NULL);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);

  const int           level = _mm_extract_epi32 (*q, 0);
  const int           upper_bound =
    P4EST3_YX_ROOT_LEN - P4EST3_YX_QUADRANT_LEN (level);
  __m128i             r;
  int                *_mem_addr;
  SC3E (sc3_array_index (nf, 0, &_mem_addr));

/* *INDENT-OFF* */
  if (level == 0) {
    _mem_addr[0] = -2;
    _mem_addr[1] = -2;
#ifdef P4_TO_P8
    _mem_addr[2] = -2;
#endif
    return NULL;
  }

  r =
  _mm_or_si128 (
    _mm_and_si128 (
      _mm_cmpeq_epi32 (
        _mm_setzero_si128 ()
      , *q)
    , _mm_set_epi32 (1, 3, 5, 0)
    )
  , _mm_and_si128 (
      _mm_cmpeq_epi32 (
        _mm_set1_epi32 (upper_bound)
      , *q)
    , _mm_set_epi32 (2, 4, 6, 0)
    )
  );

  r =
  _mm_sub_epi32 (
    r
  , _mm_set1_epi32 (1)
  );
  _mem_addr[0] = _mm_extract_epi32 (r, 3);
  _mem_addr[1] = _mm_extract_epi32 (r, 2);
#ifdef P4_TO_P8
  _mem_addr[2] = _mm_extract_epi32 (r, 1);
#endif
/* *INDENT-ON* */

  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_child (const __m128i * q, int child_id, __m128i * r)
{
  const p4est_qcoord_t level = _mm_extract_epi32 (*q, 0);

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (level < P4EST3_YX_QMAXLEVEL);
  SC3A_CHECK (0 <= child_id && child_id < P4EST_CHILDREN);

/* *INDENT-OFF* */
  *r =
   _mm_or_si128 (
     _mm_sllv_epi32 (
       _mm_srlv_epi32 (
         _mm_and_si128 (_mm_set_epi32 (child_id, child_id, child_id, 0x00)
                      , _mm_set_epi32 (0x01, 0x02, 0x04, 0x00))
       , _mm_set_epi32 (0, 1, 2, 0)
       )
     , _mm_set1_epi32 (P4EST3_YX_MAXLEVEL - (level + 1))
     )
   , *q
   );
/* *INDENT-ON* */
  *r = _mm_insert_epi32 (*r, level + 1, 0);
#ifndef P4_TO_P8
  *r = _mm_insert_epi32 (*r, 0, 1);
#endif
  SC3A_IS2 (p4est3_quadrant_zyx_is_parent_internal, q, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_parent (const __m128i * q, __m128i * r)
{
  int32_t             level = _mm_extract_epi32 (*q, 0);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (level > 0);

/* *INDENT-OFF* */
  *r =
   _mm_sub_epi32 (
     _mm_and_si128 (*q, _mm_set_epi32 (~P4EST3_YX_QUADRANT_LEN (level)
                                     , ~P4EST3_YX_QUADRANT_LEN (level)
                                     , ~P4EST3_YX_QUADRANT_LEN (level)
                                     , 0xFFFFFFFF))
   , _mm_set_epi32 (0, 0, 0, 1)
   );
/* *INDENT-ON* */

  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  SC3A_IS2 (p4est3_quadrant_zyx_is_parent_internal, r, q);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_face_neighbor (const __m128i * q, int i, __m128i * r)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  int32_t             shift =
    P4EST3_YX_QUADRANT_LEN (_mm_extract_epi32 (*q, 0));
  shift = shift * (i & 0x01 ? 1 : -1);
/* *INDENT-OFF* */
  switch (i / 2)
  {
  case 0:
  *r =
    _mm_add_epi32 (
      _mm_insert_epi32 (
        _mm_setzero_si128 (), shift, 3
      ),
      *q);
    break;
  case 1:
  *r =
    _mm_add_epi32 (
      _mm_insert_epi32 (
        _mm_setzero_si128 (), shift, 2
      ),
      *q);
    break;
  case 2:
  *r =
    _mm_add_epi32 (
      _mm_insert_epi32 (
        _mm_setzero_si128 (), shift, 1
      ),
      *q);
    break;
  default:
    break;
  }
/* *INDENT-ON* */

  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_tree_face_neighbor (const __m128i * q,
                                        sc3_array_t * transform,
                                        int face, __m128i * r)
{
  int                *my_axis;
  int                *target_axis;
  int                *edge_reverse;
  int32_t             target_xyzl[3] = { 0, 0, 0 };
  int32_t             q_coords[3];
  int32_t             level, mh, Rmh;

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (q != r);

  SC3E (sc3_array_index (transform, 0, &my_axis));
  SC3E (sc3_array_index (transform, 3, &target_axis));
  SC3E (sc3_array_index (transform, 6, &edge_reverse));

  q_coords[0] = _mm_extract_epi32 (*q, 3);
  q_coords[1] = _mm_extract_epi32 (*q, 2);
  q_coords[2] = _mm_extract_epi32 (*q, 1);
  level = _mm_extract_epi32 (*q, 0);

#ifdef P4EST_ENABLE_DEBUG
  int                 i;
  for (i = 0; i < 3; ++i) {
    SC3A_CHECK (0 <= my_axis[i] && my_axis[i] < P4EST_DIM);
    SC3A_CHECK (0 <= target_axis[i] && target_axis[i] < P4EST_DIM);
  }
  SC3A_CHECK (my_axis[0] != my_axis[2]);
  SC3A_CHECK (target_axis[0] != target_axis[2]);
  SC3A_CHECK (0 <= edge_reverse[0] && edge_reverse[0] < 2);
  SC3A_CHECK (0 <= edge_reverse[2] && edge_reverse[2] < 4);

#ifdef P4_TO_P8
  *r = _mm_set_epi32 (INT32_MIN, INT32_MIN, INT32_MIN, 0);
  SC3A_CHECK (my_axis[0] != my_axis[1] && my_axis[1] != my_axis[2]);
  SC3A_CHECK (target_axis[0] != target_axis[1] &&
              target_axis[1] != target_axis[2]);
  SC3A_CHECK (0 <= edge_reverse[1] && edge_reverse[1] < 2);
#else
  *r = _mm_set_epi32 (INT32_MIN, INT32_MIN, 0, 0);
  SC3A_CHECK (my_axis[1] == 0 && target_axis[1] == 0);
  SC3A_CHECK (edge_reverse[1] == 0);
#endif /* P4_TO_P8 */
#endif /* P4EST_ENABLE_DEBUG */

  if (level == P4EST3_YX_MAXLEVEL) {
    /* If (P4EST3_YX_MAXLEVEL == 31) Rmh will overflow.
       P4EST3_YX_MAXLEVEL == 30 so far,
       but will be increased in long term perspective. */
    mh = 0;
  }
  else {
    mh = P4EST3_YX_QUADRANT_LEN (level);
  }
  Rmh = P4EST3_YX_ROOT_LEN - mh;

  if (!edge_reverse[0]) {
    target_xyzl[target_axis[0]] = q_coords[my_axis[0]];
  }
  else {
    target_xyzl[target_axis[0]] = Rmh - q_coords[my_axis[0]];
  }
#ifdef P4_TO_P8
  if (!edge_reverse[1]) {
    target_xyzl[target_axis[1]] = q_coords[my_axis[1]];
  }
  else {
    target_xyzl[target_axis[1]] = Rmh - q_coords[my_axis[1]];
  }
#endif

  if (edge_reverse[2] == 1 || edge_reverse[2] == 3) {
    target_xyzl[target_axis[2]] = Rmh;
  }
  *r = _mm_set_epi32 (target_xyzl[0], target_xyzl[1], target_xyzl[2], level);

  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_copy (const __m128i * q, __m128i * copy)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  *copy = *q;
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_compare (const __m128i * q1, const __m128i * q2, int *j)
{
  int64_t             diff;
  __m128i             p, q;

/* *INDENT-OFF* */
  { /* Waiting for the corresponding SC3A_.. macro */
    SC3A_CHECK (p4est3_quadrant_zyx_is_valid (q1, NULL));
    SC3A_CHECK (p4est3_quadrant_zyx_is_valid (q2, NULL));
  }
  __m128i            exclor_coords = _mm_xor_si128 (*q1, *q2);
  __m128i            exclor =
                      _mm_set1_epi32 (_mm_extract_epi32 (exclor_coords, 3)
                                    | _mm_extract_epi32 (exclor_coords, 2));
#ifdef P4_TO_P8
  exclor = _mm_or_si128 (
             _mm_set_epi32 (0, 0, _mm_extract_epi32 (exclor_coords, 1), 0)
          , exclor);
#endif

  if (!_mm_extract_epi32 (exclor, 1)) {
    *j = _mm_cvtsi128_si32 (_mm_sub_epi32 (*q1, *q2));
    return NULL;
  }

  __m128i             cond = _mm_cmpgt_epi32 (
                               exclor_coords,
                              _mm_xor_si128 (exclor, exclor_coords)
                             );
#ifdef P4_TO_P8
  if (_mm_extract_epi32 (cond, 1)) {
    q = _mm_set_epi64x ((int64_t) _mm_extract_epi32 (*q2, 1)
                      , (int64_t) _mm_extract_epi32 (*q1, 1));
  } else
#if 0
    ;
#endif
#endif
  if (_mm_extract_epi32 (cond, 2)) {
    q = _mm_set_epi64x ((int64_t) _mm_extract_epi32 (*q2, 2)
                      , (int64_t) _mm_extract_epi32 (*q1, 2));
  } else {
    q = _mm_set_epi64x ((int64_t) _mm_extract_epi32 (*q2, 3)
                      , (int64_t) _mm_extract_epi32 (*q1, 3));
  }

  p = _mm_or_si128 (_mm_cmpeq_epi64 (q, _mm_setzero_si128 ())
                  , _mm_cmpgt_epi64 (q, _mm_setzero_si128 ()));

  q = _mm_add_epi64 (q
      , _mm_set_epi64x (
            _mm_extract_epi64 (p, 0) ? 0 : ((int64_t) 1 << (P4EST3_YX_MAXLEVEL + 2))
          , _mm_extract_epi64 (p, 1) ? 0 : ((int64_t) 1 << (P4EST3_YX_MAXLEVEL + 2))
        )
      );
/* *INDENT-ON* */
  diff = _mm_extract_epi64 (q, 0) - _mm_extract_epi64 (q, 1);
  *j = ((diff == 0) ? 0 : ((diff < 0) ? -1 : 1));
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_ancestor_id (const __m128i * q, int level, int *j)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (0 <= level && level <= P4EST3_YX_MAXLEVEL);
  SC3A_CHECK (_mm_extract_epi32 (*q, 0) >= level);
  *j = 0;

  if (level == 0) {
    return NULL;
  }

  __m128i             t =
    _mm_and_si128 (*q, _mm_set1_epi32 (P4EST3_YX_QUADRANT_LEN (level)));
  *j |= _mm_extract_epi32 (t, 3) ? 0x01 : 0;
  *j |= _mm_extract_epi32 (t, 2) ? 0x02 : 0;
#ifdef P4_TO_P8
  *j |= _mm_extract_epi32 (t, 1) ? 0x04 : 0;
#endif

  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_child_id (const __m128i * q, int *j)
{
  int                 level = _mm_extract_epi32 (*q, 0);
  SC3E (p4est3_quadrant_zyx_ancestor_id (q, level, j));
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_coordinates (const __m128i * q, int *j)
{
  __m128i             r;
  int                 d = P4EST3_REF_MAXLEVEL - P4EST3_YX_MAXLEVEL;
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (d >= 0);

  r = _mm_slli_epi32 (*q, d);
  j[0] = _mm_extract_epi32 (r, 3);
  j[1] = _mm_extract_epi32 (r, 2);
#ifdef P4_TO_P8
  j[2] = _mm_extract_epi32 (r, 1);
#endif
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_quadrant (const p4est_qcoord_t * c, int l, __m128i * q)
{
  SC3A_CHECK (c != NULL);
  SC3A_CHECK (q != NULL);
  SC3A_CHECK (0 <= l && l <= P4EST3_YX_MAXLEVEL);
  const int           d = P4EST3_REF_MAXLEVEL - P4EST3_YX_MAXLEVEL;

  /* Since the coordinates are normalized by the P4EST3_REF_MAXLEVEL,
     we shift them according to qvt maxlevel */
/* *INDENT-OFF* */
  *q =
  _mm_srli_epi32 (
#ifdef P4_TO_P8
    _mm_set_epi32 (c[0], c[1], c[2], 0)
#else
    _mm_set_epi32 (c[0], c[1], 0, 0)
#endif
  , d
  );
/* *INDENT-ON* */
  *q = _mm_insert_epi32 (*q, l, 0);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_level (const __m128i * q, int *l)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  *l = _mm_extract_epi32 (*q, 0);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_ancestor (const __m128i * q, int level, __m128i * r)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (_mm_extract_epi32 (*q, 0) > level && level >= 0);
/* *INDENT-OFF* */
  *r =
    _mm_insert_epi32 (
     _mm_and_si128 (*q, _mm_set1_epi32 (~(P4EST3_YX_QUADRANT_LEN (level) - 1))),
     level, 0);
/* *INDENT-ON* */
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_sibling (const __m128i * q, int sibling_id, __m128i * r)
{
  const p4est_qcoord_t q_level = _mm_extract_epi32 (*q, 0);
  const p4est_qcoord_t shift = P4EST3_YX_QUADRANT_LEN (q_level);
/* *INDENT-OFF* */
  const __m128i       add =
                      _mm_slli_epi32 (
                       _mm_srlv_epi32 (
                        _mm_and_si128 (_mm_set1_epi32 (sibling_id)
                                     , _mm_set_epi32 (0x01, 0x02, 0x04, 0))
                      , _mm_set_epi32 (0, 1, 2, 0))
                     , P4EST3_YX_MAXLEVEL - q_level);
/* *INDENT-ON* */
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (q_level > 0);
  SC3A_CHECK (sibling_id >= 0 && sibling_id < P4EST_CHILDREN);

  *r = _mm_or_si128 (add, _mm_andnot_si128 (_mm_set1_epi32 (shift), *q));
  *r = _mm_insert_epi32 (*r, q_level, 0);
#ifndef P4_TO_P8
  *r = _mm_insert_epi32 (*r, 0, 1);
#endif
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_first_descendant (const __m128i * q, int level,
                                      __m128i * fd)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (_mm_extract_epi32 (*q, 0) <= level &&
              level <= P4EST3_YX_MAXLEVEL);

  *fd = _mm_insert_epi32 (*q, level, 0);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_last_descendant (const __m128i * q, int level,
                                     __m128i * ld)
{
  p4est_qcoord_t      shift;

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (_mm_extract_epi32 (*q, 0) <= level &&
              level <= P4EST3_YX_QMAXLEVEL);

  shift = P4EST3_YX_QUADRANT_LEN (_mm_extract_epi32 (*q, 0))
    - P4EST3_YX_QUADRANT_LEN (level);

  *ld = _mm_insert_epi32 (_mm_add_epi32 (*q, _mm_set1_epi32 (shift)),
                          level, 0);
  return NULL;
}

static sc3_error_t *
p4est3_zyx_nearest_common_ancestor (const __m128i * q1, const __m128i * q2,
                                    __m128i * r)
{
  int                 maxlevel, min;
  __m128i             exclor;
  int32_t             maxclor;

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q1);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q2);

  exclor = _mm_xor_si128 (*q1, *q2);
#ifdef P4_TO_P8
/* *INDENT-OFF* */
  maxclor = _mm_extract_epi32 (exclor, 3)
          | _mm_extract_epi32 (exclor, 2)
          | _mm_extract_epi32 (exclor, 1);
/* *INDENT-ON* */
#else
/* *INDENT-OFF* */
  maxclor = _mm_extract_epi32 (exclor, 3)
          | _mm_extract_epi32 (exclor, 2);
/* *INDENT-ON* */
#endif
  maxlevel = SC_LOG2_32 (maxclor) + 1;

  SC3A_CHECK (maxlevel <= P4EST3_YX_MAXLEVEL);
  *r = _mm_and_si128 (*q1, _mm_set1_epi32 (~((1 << maxlevel) - 1)));
/* *INDENT-OFF* */
  min = (int) SC_MIN (P4EST3_YX_MAXLEVEL - maxlevel,
                       SC_MIN (_mm_extract_epi32 (*q1, 0)
                             , _mm_extract_epi32 (*q2, 0)));
/* *INDENT-ON* */
  *r = _mm_insert_epi32 (*r, min, 0);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_successor (const __m128i * q, __m128i * r)
{
  int                 level, q_level;
  int                 successor_id = 0;
  q_level = level = _mm_extract_epi32 (*q, 0);

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (q_level > 0);

  SC3E (p4est3_quadrant_zyx_ancestor_id (q, level, &successor_id));
  successor_id++;

  /* iterate until it is possible to increment the child/ancestor_id */
  while (successor_id == P4EST_CHILDREN) {
    SC3E (p4est3_quadrant_zyx_ancestor_id (q, --level, &successor_id));
    successor_id++;
    SC3A_CHECK (level > 0);
  }
  SC3A_CHECK (0 < successor_id && successor_id < P4EST_CHILDREN);

  /* compute result */
  if (level < q_level) {
    /* coarsen to level - 1 and add shifts according to the successor_id */
/* *INDENT-OFF* */
    *r =
    _mm_add_epi32 (
      _mm_sllv_epi32 (
        _mm_srlv_epi32 (
          _mm_and_si128 (_mm_set1_epi32 (successor_id)
                       , _mm_set_epi32  (0x01, 0x02, 0x04, 0x00))
        , _mm_set_epi32 (0, 1, 2, 0)
        )
      , _mm_set1_epi32 (P4EST3_YX_MAXLEVEL - level)
      )
    , _mm_and_si128 (*q, _mm_set_epi32 (~(P4EST3_YX_QUADRANT_LEN (level - 1) - 1)
                                      , ~(P4EST3_YX_QUADRANT_LEN (level - 1) - 1)
                                      , ~(P4EST3_YX_QUADRANT_LEN (level - 1) - 1)
                                      , 0xFFFFFFFF))
    );
/* *INDENT-ON* */
  }
  else {
    SC3E (p4est3_quadrant_zyx_sibling (q, successor_id, r));
  }
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_predecessor (const __m128i * q, __m128i * r)
{
  int                 level, q_level;
  int                 predecessor_id;
  q_level = level = _mm_extract_epi32 (*q, 0);

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (q_level > 0);

  SC3E (p4est3_quadrant_zyx_ancestor_id (q, level, &predecessor_id));
  predecessor_id--;

  /* iterate until it is possible to decrement the child/ancestor_id */
  while (predecessor_id == -1) {
    SC3E (p4est3_quadrant_zyx_ancestor_id (q, --level, &predecessor_id));
    predecessor_id--;
    SC3A_CHECK (level > 0);
  }
  SC3A_CHECK (0 <= predecessor_id && predecessor_id < P4EST_CHILDREN - 1);

  /* compute result */
  if (level < q_level) {
    /* coarsen to level - 1 and add shifts according to the predecessor_id */
/* *INDENT-OFF* */
    *r =
    _mm_add_epi32 (
      _mm_sllv_epi32 (
        _mm_srlv_epi32 (
          _mm_and_si128 (_mm_set1_epi32 (predecessor_id)
                       , _mm_set_epi32  (0x01, 0x02, 0x04, 0x00))
        , _mm_set_epi32 (0, 1, 2, 0)
        )
      , _mm_set1_epi32 (P4EST3_YX_MAXLEVEL - level)
      )
    , _mm_and_si128 (*q, _mm_set_epi32 (~(P4EST3_YX_QUADRANT_LEN (level - 1) - 1)
                                      , ~(P4EST3_YX_QUADRANT_LEN (level - 1) - 1)
                                      , ~(P4EST3_YX_QUADRANT_LEN (level - 1) - 1)
                                      , 0xFFFFFFFF))
    );
/* *INDENT-ON* */
  }
  else {
    SC3E (p4est3_quadrant_zyx_sibling (q, predecessor_id, r));
  }
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_linear_id (const __m128i * quadrant, int level,
                               p4est3_gloidx * id)
{
  int                 i;
  __m128i             shifted, res;

  SC3A_IS (p4est3_quadrant_zyx_is_valid, quadrant);
  SC3A_CHECK (0 <= level && level <= P4EST3_YX_MAXLEVEL);

  /* this preserves the high bits from negative numbers */
  shifted = _mm_srli_epi32 (*quadrant, P4EST3_YX_MAXLEVEL - level);

  *id = (p4est3_gloidx) 0;
  for (i = 0; i < level; ++i) {
    res = _mm_and_si128 (shifted, _mm_set1_epi32 ((uint32_t) 1 << i));
    *id |= ((uint64_t) _mm_extract_epi32 (res, 3)) << ((P4EST_DIM - 1) * i);
    *id |=
      ((uint64_t) _mm_extract_epi32 (res, 2)) << ((P4EST_DIM - 1) * i + 1);
    *id |=
      ((uint64_t) _mm_extract_epi32 (res, 1)) << ((P4EST_DIM - 1) * i + 2);
  }
  SC3A_CHECK (*id >= 0L && *id < ((int64_t) 1 << P4EST_DIM * level));
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_morton (int level, p4est3_gloidx id, __m128i * quadrant)
{
  int                 i;

  SC3A_CHECK (0 <= level && level <= P4EST3_YX_QMAXLEVEL);
  if (level < P4EST3_YX_QMAXLEVEL) {
    SC3A_CHECK (id < (((p4est3_gloidx) 1) << P4EST_DIM * level));
  }

  *quadrant = _mm_setzero_si128 ();

  /* this may set the sign bit to create negative numbers */
  for (i = 0; i < level; ++i) {
/* *INDENT-OFF* */
    __m128i xy_coord_id =
    _mm_srlv_epi64 (
      _mm_and_si128 (
        _mm_set1_epi64x (id)
      , _mm_sllv_epi64 (
          _mm_set1_epi64x (1ULL)
        , _mm_set_epi64x (P4EST_DIM * i
                        , P4EST_DIM * i + 1)
        )
      )
    , _mm_set_epi64x ((P4EST_DIM - 1) * i
                    , (P4EST_DIM - 1) * i + 1)
    );
    *quadrant
    = _mm_or_si128 (*quadrant
                  , _mm_set_epi32 (
                      (p4est_qcoord_t) _mm_extract_epi64 (xy_coord_id, 1)
                    , (p4est_qcoord_t) _mm_extract_epi64 (xy_coord_id, 0)
#ifdef P4_TO_P8
                    , (p4est_qcoord_t) ((id & (1ULL << (P4EST_DIM * i + 2)))
                        >> ((P4EST_DIM - 1) * i + 2))
#else
                    , 0
#endif /* P4_TO_P8 */
                    , 0
                  ));
/* *INDENT-ON* */
  }

  *quadrant = _mm_slli_epi32 (*quadrant, P4EST3_YX_MAXLEVEL - level);
  *quadrant = _mm_insert_epi32 (*quadrant, level, 0);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, quadrant);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_root (__m128i * r)
{
  SC3E (p4est3_quadrant_zyx_morton (0, 0, r));
  return NULL;
}

static const p4est3_quadrant_vtable_t quadrant_vtable_yx =
{
  P4EST_STRING "_quadrant_vtable_zyx",
  P4EST_DIM,
  P4EST3_YX_MAXLEVEL,
  sizeof (__m128i),

  (p4est3_quadrant_num_uniform_t) p4est3_quadrant_zyx_num_uniform,
  (p4est3_quadrant_root_t) p4est3_quadrant_zyx_root,
  (p4est3_quadrant_morton_t) p4est3_quadrant_zyx_morton,
  (p4est3_quadrant_quadrant_t) p4est3_quadrant_zyx_quadrant,
  (p4est3_quadrant_is_t) p4est3_quadrant_zyx_is_valid,
  (p4est3_quadrant_level_t) p4est3_quadrant_zyx_level,
  (p4est3_quadrant_child_id_t) p4est3_quadrant_zyx_child_id,
  (p4est3_quadrant_ancestor_id_t) p4est3_quadrant_zyx_ancestor_id,
  (p4est3_quadrant_in_out_t) p4est3_quadrant_zyx_coordinates,
  (p4est3_quadrant_linear_id_t) p4est3_quadrant_zyx_linear_id,
  (p4est3_quadrant_get_tree_boundary_t) p4est3_quadrant_zyx_get_tree_boundary,
  (p4est3_quadrant_tree_boundaries_t) p4est3_quadrant_zyx_tree_boundaries,
  (p4est3_quadrant_copy_t) p4est3_quadrant_zyx_copy,
  (p4est3_quadrant_child_t) p4est3_quadrant_zyx_child,
  (p4est3_quadrant_sibling_t) p4est3_quadrant_zyx_sibling,
  (p4est3_quadrant_parent_t) p4est3_quadrant_zyx_parent,
  (p4est3_quadrant_ancestor_t) p4est3_quadrant_zyx_ancestor,
  (p4est3_quadrant_predecessor_t) p4est3_quadrant_zyx_predecessor,
  (p4est3_quadrant_successor_t) p4est3_quadrant_zyx_successor,
  (p4est3_quadrant_first_descendant_t) p4est3_quadrant_zyx_first_descendant,
  (p4est3_quadrant_last_descendant_t) p4est3_quadrant_zyx_last_descendant,
  (p4est3_quadrant_face_neighbor_t) p4est3_quadrant_zyx_face_neighbor,
  (p4est3_quadrant_tree_face_neighbor_t) p4est3_quadrant_zyx_tree_face_neighbor,
  (p4est3_quadrant_compare_t) p4est3_quadrant_zyx_compare,
  NULL, /*< use generic implementation */
  (p4est3_quadrant_is_ancestor_t) p4est3_quadrant_zyx_is_ancestor,
  (p4est3_quadrant_is_parent_t) p4est3_quadrant_zyx_is_parent,
  (p4est3_nearest_common_ancestor_t) p4est3_zyx_nearest_common_ancestor
};
#endif /* P4EST_ENABLE_AVX2 */

static const p4est3_quadrant_vtable_t * qvt_yx =
#if ((P4EST_DIM == 2 && defined(P4EST_ENABLE_BUILD_2D)) \
 || (P4EST_DIM == 3 && defined(P4EST_ENABLE_BUILD_3D))) \
 && defined(P4EST_ENABLE_AVX2)
  &quadrant_vtable_yx
#else
  NULL
#endif
;

sc3_error_t        *
p4est3_quadrant_yx_vtable (const p4est3_quadrant_vtable_t ** qvt)
{
  SC3A_CHECK (qvt != NULL);

  if (qvt_yx != NULL) {
    /* verify correctness */
    SC3A_IS (p4est3_quadrant_vtable_is_valid, qvt_yx);

    /* pass internally constructed virtual table to the outside */
    *qvt = qvt_yx;
  }
  else {
    *qvt = NULL;
  }
  return NULL;
}
