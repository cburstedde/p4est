/*
  This file is part of p4est, version 3
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
#include <p4est_bits.h>
#include <p4est_algorithms.h>
#include <p4est3_p4est.h>
#else
#include <p8est_bits.h>
#include <p8est_algorithms.h>
#include <p4est3_p8est.h>
#endif

static              p4est3_gloidx
p4est_quadrant_vtable_num_uniform (int level)
{
  if (level < 0) {
    return -1;
  }
  return p4est3_glopow (P4EST_CHILDREN, level);
}

static int
p4est_quadrant_vtable_is_valid (const void *q, char *reason)
{
  SC3E_TEST (p4est_quadrant_is_valid ((const p4est_quadrant_t *) q), reason);
  SC3E_YES (reason);
}

static int
p4est_quadrant_vtable_is_equal (const void *q1, const void *q2, char *reason)
{
  SC3E_TEST (p4est_quadrant_is_equal ((const p4est_quadrant_t *) q1,
                                      (const p4est_quadrant_t *) q2), reason);
  SC3E_YES (reason);
}

static sc3_error_t *
p4est_quadrant_vtable_get_tree_boundary (const void *q, int face, int *j)
{
  const p4est_quadrant_t *quad = (const p4est_quadrant_t *) q;
  int                 direction = face / 2;
  p4est_qcoord_t      coord;
  int                 bound;

#ifdef P4_TO_P8
  coord = (direction == 0) ? quad->x : (direction == 1) ? quad->y : quad->z;
#else
  coord = (direction == 0) ? quad->x : quad->y;
#endif
  bound =
    face % 2 == 0 ? 0 : P4EST_ROOT_LEN - P4EST_QUADRANT_LEN (quad->level);

  *j = (coord == bound);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_tree_boundaries (const void *q, sc3_array_t * nf)
{
  SC3A_CHECK (q != NULL);
  SC3A_CHECK (nf != NULL);
  SC3A_IS (p4est_quadrant_vtable_is_valid, q);
  const p4est_quadrant_t *quad = (const p4est_quadrant_t *) q;
  const int           upper_bound =
    P4EST_ROOT_LEN - P4EST_QUADRANT_LEN (quad->level);
  int                *x, *y;
#ifdef P4_TO_P8
  int                *z;
#endif

  SC3E (sc3_array_index (nf, 0, &x));
  SC3E (sc3_array_index (nf, 1, &y));
#ifdef P4_TO_P8
  SC3E (sc3_array_index (nf, 2, &z));
#endif

  if (quad->level == 0) {
    *x = *y = -2;
#ifdef P4_TO_P8
    *z = -2;
#endif
    return NULL;
  }

  *x = quad->x == 0 ? 0 : (quad->x == upper_bound) ? 1 : -1;
  *y = quad->y == 0 ? 2 : (quad->y == upper_bound) ? 3 : -1;
#ifdef P4_TO_P8
  *z = quad->z == 0 ? 4 : (quad->z == upper_bound) ? 5 : -1;
#endif

  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_level (const void *q, int *l)
{
  SC3A_CHECK (q != NULL);
  SC3A_CHECK (l != NULL);

  *l = ((const p4est_quadrant_t *) q)->level;
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_child_id (const void *q, int *j)
{
  SC3A_CHECK (j != NULL);

  *j = p4est_quadrant_child_id ((const p4est_quadrant_t *) q);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_ancestor_id (const void *q, int i, int *j)
{
  SC3A_CHECK (j != NULL);

  *j = p4est_quadrant_ancestor_id ((const p4est_quadrant_t *) q, i);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_coordinates (const void *q, void *j)
{
  const p4est_quadrant_t *quad = (const p4est_quadrant_t *) q;
  p4est_qcoord_t     *coords = (p4est_qcoord_t *) j;
  int                 d = P4EST3_REF_MAXLEVEL - P4EST_MAXLEVEL;
  SC3A_CHECK (coords != NULL);
  SC3A_CHECK (d >= 0);
  coords[0] = quad->x << d;
  coords[1] = quad->y << d;
#ifdef P4_TO_P8
  coords[2] = quad->z << d;
#endif
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_quadrant (const void *c, int l, void *q)
{
  const p4est_qcoord_t *coords = (const p4est_qcoord_t *) c;
  p4est_quadrant_t   *quad = (p4est_quadrant_t *) q;
  const int           d = P4EST3_REF_MAXLEVEL - P4EST_MAXLEVEL;
  SC3A_CHECK (coords != NULL);
  SC3A_CHECK (q != NULL);
  SC3A_CHECK (0 <= l && l <= P4EST_MAXLEVEL);

  /* Since the coordinates are normalized by the P4EST3_REF_MAXLEVEL,
     we shift them according to qvt maxlevel */
  quad->x = coords[0] >> d;
  quad->y = coords[1] >> d;
#ifdef P4_TO_P8
  quad->z = coords[2] >> d;
#endif
  quad->level = l;
  SC3A_IS (p4est_quadrant_vtable_is_valid, q);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_compare (const void *q1, const void *q2, int *j)
{
  SC3A_CHECK (j != NULL);

  *j = p4est_quadrant_compare ((const p4est_quadrant_t *) q1,
                               (const p4est_quadrant_t *) q2);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_root (void *r)
{
  p4est_quadrant_set_morton ((p4est_quadrant_t *) r, 0, 0);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_copy (const void *q, void *r)
{
  p4est_quadrant_copy ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_parent (const void *q, void *r)
{
  p4est_quadrant_parent
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_sibling (const void *q, int i, void *r)
{
  p4est_quadrant_sibling
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r, i);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_face_neighbor (const void *q, int face, void *r)
{
  p4est_quadrant_face_neighbor
    ((const p4est_quadrant_t *) q, face, (p4est_quadrant_t *) r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_vtable_tree_face_neighbor (const void *q,
                                           sc3_array_t * transform,
                                           int face, void *r)
{
  p4est_quadrant_t    temp;
  int                *idx;

  p4est_quadrant_face_neighbor
    ((const p4est_quadrant_t *) q, face, (p4est_quadrant_t *) r);

  /* Input and output pointing on the same memory are forbidden.
     See the documentation for p4est_quadrant_transform_face */
  temp = *((p4est_quadrant_t *) r);
  SC3E (sc3_array_index (transform, 0, &idx));
  p4est_quadrant_transform_face (&temp, (p4est_quadrant_t *) r, idx);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_predecessor (const void *q, void *r)
{
  p4est_quadrant_predecessor
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_successor (const void *q, void *r)
{
  p4est_quadrant_successor
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_child (const void *q, int i, void *r)
{
  p4est_quadrant_child
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r, i);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_ancestor (const void *q, int l, void *r)
{
  p4est_quadrant_ancestor
    ((const p4est_quadrant_t *) q, l, (p4est_quadrant_t *) r);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_first_descendant (const void *q, int l, void *r)
{
  p4est_quadrant_first_descendant
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r, l);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_last_descendant (const void *q, int l, void *r)
{
  p4est_quadrant_last_descendant
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r, l);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_morton (int level, p4est_gloidx_t id, void *r)
{
  p4est_quadrant_set_morton ((p4est_quadrant_t *) r, level, (uint64_t) id);
  return NULL;
}

static sc3_error_t *
p4est_vtable_nearest_common_ancestor (const void *q1, const void *q2, void *r)
{
  p4est_nearest_common_ancestor ((const p4est_quadrant_t *) q1,
                                 (const p4est_quadrant_t *) q2,
                                 (p4est_quadrant_t *) r);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_linear_id (const void *q, int level,
                                 p4est_gloidx_t * id)
{
  *id = p4est_quadrant_linear_id ((const p4est_quadrant_t *) q, level);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_is_ancestor (const void *q1, const void *q2, int *j)
{
  *j = p4est_quadrant_is_ancestor ((const p4est_quadrant_t *) q1,
                                   (const p4est_quadrant_t *) q2);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_is_parent (const void *q, const void *r, int *j)
{
  *j = p4est_quadrant_is_parent ((const p4est_quadrant_t *) q,
                                 (const p4est_quadrant_t *) r);
  return NULL;
}

static const p4est3_quadrant_vtable_t quadrant_vtable_p4est =
{
  P4EST_STRING "quadrant_vtable_p4est",
  P4EST_DIM,
  P4EST_QMAXLEVEL,
  sizeof (p4est_quadrant_t),

  (p4est3_quadrant_num_uniform_t) p4est_quadrant_vtable_num_uniform,
  (p4est3_quadrant_root_t) p4est_quadrant_vtable_root,
  (p4est3_quadrant_morton_t) p4est_quadrant_vtable_morton,
  (p4est3_quadrant_quadrant_t) p4est_quadrant_vtable_quadrant,
  (p4est3_quadrant_is_t) p4est_quadrant_vtable_is_valid,
  (p4est3_quadrant_level_t) p4est_quadrant_vtable_level,
  (p4est3_quadrant_child_id_t) p4est_quadrant_vtable_child_id,
  (p4est3_quadrant_ancestor_id_t) p4est_quadrant_vtable_ancestor_id,
  (p4est3_quadrant_in_out_t) p4est_quadrant_vtable_coordinates,
  (p4est3_quadrant_linear_id_t) p4est_quadrant_vtable_linear_id,
  (p4est3_quadrant_get_tree_boundary_t) p4est_quadrant_vtable_get_tree_boundary,
  (p4est3_quadrant_tree_boundaries_t) p4est_quadrant_vtable_tree_boundaries,
  (p4est3_quadrant_copy_t) p4est_quadrant_vtable_copy,
  (p4est3_quadrant_child_t) p4est_quadrant_vtable_child,
  (p4est3_quadrant_sibling_t) p4est_quadrant_vtable_sibling,
  (p4est3_quadrant_parent_t) p4est_quadrant_vtable_parent,
  (p4est3_quadrant_ancestor_t) p4est_quadrant_vtable_ancestor,
  (p4est3_quadrant_predecessor_t) p4est_quadrant_vtable_predecessor,
  (p4est3_quadrant_successor_t) p4est_quadrant_vtable_successor,
  (p4est3_quadrant_first_descendant_t) p4est_quadrant_vtable_first_descendant,
  (p4est3_quadrant_last_descendant_t) p4est_quadrant_vtable_last_descendant,
  (p4est3_quadrant_face_neighbor_t) p4est_quadrant_vtable_face_neighbor,
  (p4est3_quadrant_tree_face_neighbor_t) p4est3_quadrant_vtable_tree_face_neighbor,
  (p4est3_quadrant_compare_t) p4est_quadrant_vtable_compare,
  (p4est3_quadrant_is2_t) p4est_quadrant_vtable_is_equal,
  (p4est3_quadrant_is_ancestor_t) p4est_quadrant_vtable_is_ancestor,
  (p4est3_quadrant_is_parent_t) p4est_quadrant_vtable_is_parent,
  (p4est3_nearest_common_ancestor_t) p4est_vtable_nearest_common_ancestor
};

static const p4est3_quadrant_vtable_t * qvt_p4est =
#if (P4EST_DIM == 2 && defined(P4EST_ENABLE_BUILD_2D)) \
 || (P4EST_DIM == 3 && defined(P4EST_ENABLE_BUILD_3D))
  &quadrant_vtable_p4est
#else
  NULL
#endif
;

sc3_error_t        *
p4est3_quadrant_vtable_p4est (const p4est3_quadrant_vtable_t ** qvt)
{
  SC3A_CHECK (qvt != NULL);

  if (qvt_p4est != NULL) {
    /* verify correctness */
    SC3A_IS (p4est3_quadrant_vtable_is_valid, qvt_p4est);

    /* pass internally constructed virtual table to the outside */
    *qvt = qvt_p4est;
  }
  else {
    *qvt = NULL;
  }
  return NULL;
}
