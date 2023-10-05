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

#include <p4est3_quadrant_vtable.h>

int
p4est3_quadrant_vtable_is_valid (const p4est3_quadrant_vtable_t * qvt, char *reason)
{
  /* arguments */
  SC3E_TEST (qvt != NULL, reason);

  /*** test member variables ***/
  SC3E_TEST (0 < qvt->dim && qvt->dim <= 3, reason);
  SC3E_TEST (0 < qvt->max_level, reason);

  /*** test member functions ***/
  SC3E_TEST (qvt->quadrant_tree_boundaries != NULL, reason);
  SC3E_TEST (qvt->quadrant_get_tree_boundary != NULL, reason);
  SC3E_TEST (qvt->quadrant_num_uniform != NULL, reason);
  SC3E_TEST (qvt->quadrant_level != NULL, reason);
  SC3E_TEST (qvt->quadrant_child_id != NULL, reason);
  SC3E_TEST (qvt->quadrant_ancestor_id != NULL, reason);
  SC3E_TEST (qvt->quadrant_coordinates != NULL, reason);
  SC3E_TEST (qvt->quadrant_quadrant != NULL, reason);
  SC3E_TEST (qvt->quadrant_compare != NULL, reason);
  SC3E_TEST (qvt->quadrant_root != NULL, reason);
  SC3E_TEST (qvt->quadrant_copy != NULL, reason);
  SC3E_TEST (qvt->quadrant_parent != NULL || qvt->quadrant_ancestor != NULL,
             reason);
  SC3E_TEST (qvt->quadrant_face_neighbor != NULL, reason);
  SC3E_TEST (qvt->quadrant_tree_face_neighbor != NULL, reason);
  SC3E_TEST (qvt->quadrant_predecessor != NULL, reason);
  SC3E_TEST (qvt->quadrant_successor != NULL, reason);
  SC3E_TEST (qvt->quadrant_child != NULL, reason);
  SC3E_TEST (qvt->quadrant_first_descendant != NULL, reason);
  SC3E_TEST (qvt->quadrant_last_descendant != NULL, reason);
  SC3E_TEST (qvt->quadrant_morton != NULL, reason);

  /* this is it */
  SC3E_YES (reason);
}

const char *
p4est3_quadrant_name (const p4est3_quadrant_vtable_t * qvt)
{
  if (qvt == NULL) {
    return NULL;
  }
  return qvt->name;
}

int
p4est3_quadrant_dim (const p4est3_quadrant_vtable_t * qvt)
{
  if (qvt == NULL || qvt->dim <= 0 || qvt->dim > 3) {
    return -1;
  }
  return qvt->dim;
}

int
p4est3_quadrant_max_level (const p4est3_quadrant_vtable_t * qvt)
{
  if (qvt == NULL || qvt->max_level < 0) {
    return -1;
  }
  return qvt->max_level;
}

size_t
p4est3_quadrant_size (const p4est3_quadrant_vtable_t * qvt)
{
  if (qvt == NULL) {
    return 0;
  }
  return qvt->quadrant_size;
}

int
p4est3_quadrant_num_children (const p4est3_quadrant_vtable_t * qvt)
{
  if (qvt == NULL || (qvt->dim <= 1 || qvt->dim > 3)) {
    return -1;
  }
  return 1 << qvt->dim;
}

p4est3_gloidx
p4est3_quadrant_num_uniform (const p4est3_quadrant_vtable_t * qvt, int level)
{
  if (qvt == NULL || qvt->quadrant_num_uniform == NULL) {
    return -1;
  }
  return qvt->quadrant_num_uniform (level);
}

int
p4est3_quadrant_vtable_is2_valid (const p4est3_quadrant_vtable_t * qvt,
                                  const void *q, char *reason)
{
  SC3E_TEST (qvt != NULL, reason);
  if (qvt->quadrant_is_valid != NULL) {
    SC3E_IS (qvt->quadrant_is_valid, q, reason);
  }
  SC3E_YES (reason);
}

int
p4est3_quadrant_is3_equal (const p4est3_quadrant_vtable_t * qvt,
                           const void *q1, const void *q2, char *reason)
{
  SC3E_TEST (qvt != NULL, reason);
  if (qvt->quadrant_is_equal != NULL) {
    SC3E_IS2 (qvt->quadrant_is_equal, q1, q2, reason);
  }
  else {
    SC3E_TEST (!memcmp (q1, q2, qvt->quadrant_size), reason);
  }
  SC3E_YES (reason);
}

sc3_error_t        *
p4est3_quadrant_get_tree_boundary (const p4est3_quadrant_vtable_t * qvt,
                                   const void *q, int i, int *j)
{
  /* TODO: either require both get_tree_boundary and tree_boundaries or none */
  SC3A_CHECK (qvt != NULL);
  if (qvt->quadrant_get_tree_boundary != NULL) {
    SC3E (qvt->quadrant_get_tree_boundary (q, i, j));
  }
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_tree_boundaries (const p4est3_quadrant_vtable_t * qvt,
                                 const void *q, sc3_array_t * nf)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_tree_boundaries != NULL);
  SC3E (qvt->quadrant_tree_boundaries (q, nf));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_level (const p4est3_quadrant_vtable_t * qvt, const void *q, int *l)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_level != NULL);
  SC3E (qvt->quadrant_level (q, l));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_child_id (const p4est3_quadrant_vtable_t * qvt,
                          const void *q, int *j)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_child_id != NULL);
  SC3E (qvt->quadrant_child_id (q, j));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_ancestor_id (const p4est3_quadrant_vtable_t * qvt,
                             const void *q, int l, int *j)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_ancestor_id != NULL);
  SC3E (qvt->quadrant_ancestor_id (q, l, j));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_coordinates (const p4est3_quadrant_vtable_t * qvt,
                             const void *q, void *j)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_coordinates != NULL);
  SC3E (qvt->quadrant_coordinates (q, j));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_quadrant (const p4est3_quadrant_vtable_t * qvt,
                          const void *c, int l, void *q)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_quadrant != NULL);
  SC3E (qvt->quadrant_quadrant (c, l, q));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_translate (const p4est3_quadrant_vtable_t * vtold,
                           const void *qin,
                           const p4est3_quadrant_vtable_t * vtnew,
                           void *qout)
{
  int                 level;
  int32_t             c[3] = {-1, -1, -1};

  SC3A_CHECK (vtold != NULL && vtnew != NULL);
  SC3A_CHECK (vtold->quadrant_is_valid != NULL);
  SC3A_IS (vtold->quadrant_is_valid, qin);
  SC3A_CHECK (vtold->dim == vtnew->dim);

  if (vtold == vtnew) {
    /* just hardcopy the quadrant */
    SC3A_CHECK (vtold->quadrant_copy != NULL);
    SC3E (vtold->quadrant_copy(qin, qout));
  }
  else {
    SC3A_CHECK (vtold->quadrant_level != NULL);
    SC3E (vtold->quadrant_level(qin, &level));
    SC3A_CHECK (level <= vtnew->max_level);
    SC3A_CHECK (vtold->quadrant_coordinates != NULL
                && vtnew->quadrant_quadrant != NULL);
    SC3E (vtold->quadrant_coordinates(qin, c));
    SC3E (vtnew->quadrant_quadrant(c, level, qout));
  }
  SC3A_IS (vtnew->quadrant_is_valid, qout);
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_compare (const p4est3_quadrant_vtable_t * qvt,
                         const void *q1, const void *q2, int *j)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_compare != NULL);
  SC3E (qvt->quadrant_compare (q1, q2, j));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_root (const p4est3_quadrant_vtable_t * qvt, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_root != NULL);
  SC3E (qvt->quadrant_root (r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_copy (const p4est3_quadrant_vtable_t * qvt, const void *q, void *r)
{
  SC3A_CHECK (qvt != NULL);

  if (qvt->quadrant_copy != NULL) {
    SC3E (qvt->quadrant_copy (q, r));
  }
  else {
    memcpy (r, q, qvt->quadrant_size);
  }
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_parent (const p4est3_quadrant_vtable_t * qvt,
                        const void *q, void *r)
{
  SC3A_CHECK (qvt != NULL);

  if (qvt->quadrant_parent != NULL) {
    SC3E (qvt->quadrant_parent (q, r));
  }
  else {
    int                 level;

    SC3A_CHECK (qvt->quadrant_level != NULL);
    SC3E (qvt->quadrant_level (q, &level));
    SC3A_CHECK (level > 0);

    SC3A_CHECK (qvt->quadrant_ancestor != NULL);
    SC3E (qvt->quadrant_ancestor (q, level - 1, r));
  }
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_sibling (const p4est3_quadrant_vtable_t * qvt,
                         const void *q, int i, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_sibling != NULL);
  SC3E (qvt->quadrant_sibling (q, i, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_face_neighbor (const p4est3_quadrant_vtable_t * qvt,
                               const void *q, int i, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_face_neighbor != NULL);
  SC3A_CHECK (qvt->quadrant_get_tree_boundary != NULL);
#ifdef P4EST_ENABLE_DEBUG
  int                 j;
  SC3E (qvt->quadrant_get_tree_boundary (q, i, &j));
  SC3A_CHECK (!j);
#endif

  SC3E (qvt->quadrant_face_neighbor (q, i, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_tree_face_neighbor (const p4est3_quadrant_vtable_t * qvt,
                                    const void *q, sc3_array_t * transform,
                                    int i, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_tree_face_neighbor != NULL);
  SC3A_CHECK (qvt->quadrant_get_tree_boundary != NULL);
#ifdef P4EST_ENABLE_DEBUG
  int                 j;
  SC3E (qvt->quadrant_get_tree_boundary (q, i, &j));
  SC3A_CHECK (j);
#endif

  SC3E (qvt->quadrant_tree_face_neighbor (q, transform, i, r));

  return NULL;
}

sc3_error_t        *
p4est3_quadrant_predecessor (const p4est3_quadrant_vtable_t * qvt,
                             const void *q, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_predecessor != NULL);
  SC3E (qvt->quadrant_predecessor (q, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_successor (const p4est3_quadrant_vtable_t * qvt,
                           const void *q, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_successor != NULL);
  SC3E (qvt->quadrant_successor (q, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_child (const p4est3_quadrant_vtable_t * qvt,
                       const void *q, int i, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_child != NULL);
  SC3E (qvt->quadrant_child (q, i, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_ancestor (const p4est3_quadrant_vtable_t * qvt,
                          const void *q, int l, void *r)
{
  SC3A_CHECK (qvt != NULL);

  if (qvt->quadrant_ancestor != NULL) {
    SC3E (qvt->quadrant_ancestor (q, l, r));
  }
  else {
    int                 level;

    SC3A_CHECK (l >= 0);
    SC3E (qvt->quadrant_level (q, &level));
    SC3A_CHECK (level >= l);

    SC3E (p4est3_quadrant_copy (qvt, q, r));

    SC3A_CHECK (qvt->quadrant_parent != NULL);
    while (level > l) {
      SC3E (qvt->quadrant_parent (r, r));
      SC3E (qvt->quadrant_level (r, &level));
    }
  }
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_first_descendant (const p4est3_quadrant_vtable_t * qvt,
                                  const void *q, int l, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_first_descendant != NULL);
  SC3E (qvt->quadrant_first_descendant (q, l, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_last_descendant (const p4est3_quadrant_vtable_t * qvt,
                                 const void *q, int l, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_last_descendant != NULL);
  SC3E (qvt->quadrant_last_descendant (q, l, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_morton (const p4est3_quadrant_vtable_t * qvt,
                        int level, p4est3_gloidx id, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_morton != NULL);
  SC3E (qvt->quadrant_morton (level, id, r));
  return NULL;
}

sc3_error_t        *
p4est3_nearest_common_ancestor (const p4est3_quadrant_vtable_t * qvt,
                                const void *q1, const void *q2, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->nearest_common_ancestor != NULL);
  SC3E (qvt->nearest_common_ancestor (q1, q2, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_linear_id (const p4est3_quadrant_vtable_t * qvt,
                           const void *q, int l, p4est3_gloidx * id)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_linear_id != NULL);
  SC3E (qvt->quadrant_linear_id (q, l, id));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_is_ancestor (const p4est3_quadrant_vtable_t * qvt,
                             const void *q1, const void *q2, int *j)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_is_ancestor != NULL);
  SC3E (qvt->quadrant_is_ancestor (q1, q2, j));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_is_parent (const p4est3_quadrant_vtable_t * qvt,
                           const void *q1, const void *q2, int *j)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_is_parent != NULL);
  SC3E (qvt->quadrant_is_parent (q1, q2, j));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_array_new (sc3_allocator_t * alloc,
                           const p4est3_quadrant_vtable_t * qvt,
                           p4est3_locidx n, sc3_array_t ** arr)
{
  SC3E_RETVAL (arr, NULL);
  SC3A_IS (sc3_allocator_is_setup, alloc);
  SC3A_CHECK (qvt != NULL);
  SC3A_CHECK (n >= 0);

  SC3E (sc3_array_new (alloc, arr));
  SC3E (sc3_array_set_elem_size (*arr, qvt->quadrant_size));
  SC3E (sc3_array_set_elem_count (*arr, n));
  SC3E (sc3_array_setup (*arr));

  return NULL;
}
