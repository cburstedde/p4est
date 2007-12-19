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

static const int    hash_table_minsize = 1361;
static const int    hash_table_maxsize = 99133;

/* *INDENT-OFF* */

static const int8_t log_lookup_table[256] =
{ -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
   4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
   5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
   5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
   7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
   7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
   7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
   7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
   7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
   7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
   7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
   7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
};

/** The offsets of the 3 indirect neighbors in units of h.
 * Indexing [cid][neighbor][xy] where cid is the child id.
 * Neighbors are indexed in z-order.
 */
static const int32_t indirect_neighbors[4][3][2] =
{{{-1, -1}, { 1, -1}, {-1, 1}},
 {{ 0, -1}, { 2, -1}, { 1, 0}},
 {{-1,  0}, {-2,  1}, { 0, 1}},
 {{ 1, -1}, {-1,  1}, { 1, 1}}
};

/** Indicate which neighbor to omit if edges are balanced, not corners
 * Indexing [cid] where cid is the child id.
 */
static const int corners_omitted[4] =
{ 0, 1, 1, 2 };

/* *INDENT-ON* */

#define P4EST_LOG2_8(x) (log_lookup_table[(x)])
#define P4EST_LOG2_16(x) (((x) > 0xff) ? \
                          (P4EST_LOG2_8 ((x) >> 8) + 8) : P4EST_LOG2_8 (x))
#define P4EST_LOG2_32(x) (((x) > 0xffff) ? \
                          (P4EST_LOG2_16 ((x) >> 16)) + 16 : P4EST_LOG2_16 (x))

/* here come small auxiliary functions */

int
p4est_quadrant_compare (const void *v1, const void *v2)
{
  const p4est_quadrant_t *q1 = v1;
  const p4est_quadrant_t *q2 = v2;

  int32_t             exclorx, exclory;

  P4EST_ASSERT (p4est_quadrant_is_valid (q1));
  P4EST_ASSERT (p4est_quadrant_is_valid (q2));

  exclorx = q1->x ^ q2->x;
  exclory = q1->y ^ q2->y;

  if (exclory == 0 && exclorx == 0) {
    return q1->level - q2->level;
  }
  else if (P4EST_LOG2_32 (exclory) >= P4EST_LOG2_32 (exclorx)) {
    return q1->y - q2->y;
  }
  else {
    return q1->x - q2->x;
  }
}

int
p4est_quadrant_is_equal (const void *v1, const void *v2)
{
  const p4est_quadrant_t *q1 = v1;
  const p4est_quadrant_t *q2 = v2;

  P4EST_ASSERT (p4est_quadrant_is_valid (q1));
  P4EST_ASSERT (p4est_quadrant_is_valid (q2));

  return (q1->level == q2->level && q1->x == q2->x && q1->y == q2->y);
}

int
p4est_quadrant_hash (const void *v)
{
  const p4est_quadrant_t *q = v;

  return p4est_quadrant_linear_id (q, q->level) % (1LL << 30);
}

int8_t
p4est_quadrant_child_id (const p4est_quadrant_t * q)
{
  int                 id = 0;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  id |= ((q->x & (1 << (P4EST_MAXLEVEL - q->level))) ? 0x01 : 0);
  id |= ((q->y & (1 << (P4EST_MAXLEVEL - q->level))) ? 0x02 : 0);

  return (int8_t) id;
}

int
p4est_quadrant_is_valid (const p4est_quadrant_t * q)
{
  return
    (q->level >= 0 && q->level <= P4EST_MAXLEVEL) &&
    (q->x >= 0 && q->x < (1 << P4EST_MAXLEVEL)) &&
    (q->y >= 0 && q->y < (1 << P4EST_MAXLEVEL)) &&
    ((q->x & ((1 << (P4EST_MAXLEVEL - q->level)) - 1)) == 0) &&
    ((q->y & ((1 << (P4EST_MAXLEVEL - q->level)) - 1)) == 0);
}

int
p4est_quadrant_is_sibling (const p4est_quadrant_t * q1,
                           const p4est_quadrant_t * q2)
{
  int32_t             exclorx, exclory;

  P4EST_ASSERT (p4est_quadrant_is_valid (q1));
  P4EST_ASSERT (p4est_quadrant_is_valid (q2));

  exclorx = q1->x ^ q2->x;
  exclory = q1->y ^ q2->y;
  if (exclorx == 0 && exclory == 0) {
    return 0;
  }

  return
    (q1->level == q2->level) &&
    ((exclorx & ~(1 << (P4EST_MAXLEVEL - q1->level))) == 0) &&
    ((exclory & ~(1 << (P4EST_MAXLEVEL - q1->level))) == 0);
}

int
p4est_quadrant_is_sibling_D (const p4est_quadrant_t * q1,
                             const p4est_quadrant_t * q2)
{
  p4est_quadrant_t    p1, p2;

  /* validity of q1 and q2 is asserted in p4est_quadrant_is_equal */
  if (p4est_quadrant_is_equal (q1, q2)) {
    return 0;
  }

  p4est_quadrant_parent (q1, &p1);
  p4est_quadrant_parent (q2, &p2);

  return p4est_quadrant_is_equal (&p1, &p2);
}

int
p4est_quadrant_is_family (const p4est_quadrant_t * q0,
                          const p4est_quadrant_t * q1,
                          const p4est_quadrant_t * q2,
                          const p4est_quadrant_t * q3)
{
  int32_t             inc;

  P4EST_ASSERT (p4est_quadrant_is_valid (q0));
  P4EST_ASSERT (p4est_quadrant_is_valid (q1));
  P4EST_ASSERT (p4est_quadrant_is_valid (q2));
  P4EST_ASSERT (p4est_quadrant_is_valid (q3));

  if (q0->level != q1->level ||
      q0->level != q2->level || q0->level != q3->level) {
    return 0;
  }

  inc = (1 << (P4EST_MAXLEVEL - q0->level));
  return ((q0->x + inc == q1->x && q0->y == q1->y) &&
          (q0->x == q2->x && q0->y + inc == q2->y) &&
          (q1->x == q3->x && q2->y == q3->y));
}

int
p4est_quadrant_is_parent (const p4est_quadrant_t * q,
                          const p4est_quadrant_t * r)
{
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_valid (r));

  return
    (q->level + 1 == r->level) &&
    (q->x == (r->x & ~(1 << (P4EST_MAXLEVEL - r->level)))) &&
    (q->y == (r->y & ~(1 << (P4EST_MAXLEVEL - r->level))));
}

int
p4est_quadrant_is_parent_D (const p4est_quadrant_t * q,
                            const p4est_quadrant_t * r)
{
  p4est_quadrant_t    p;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  /* validity of r is asserted in p4est_quadrant_parent */
  p4est_quadrant_parent (r, &p);

  return p4est_quadrant_is_equal (q, &p);
}

int
p4est_quadrant_is_ancestor (const p4est_quadrant_t * q,
                            const p4est_quadrant_t * r)
{
  int32_t             exclorx;
  int32_t             exclory;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_valid (r));

  if (q->level >= r->level) {
    return 0;
  }

  exclorx = (q->x ^ r->x) >> (P4EST_MAXLEVEL - q->level);
  exclory = (q->y ^ r->y) >> (P4EST_MAXLEVEL - q->level);

  return (exclorx == 0 && exclory == 0);
}

int
p4est_quadrant_is_ancestor_D (const p4est_quadrant_t * q,
                              const p4est_quadrant_t * r)
{
  p4est_quadrant_t    s;

  /* validity of q and r is asserted in p4est_quadrant_is_equal */
  if (p4est_quadrant_is_equal (q, r)) {
    return 0;
  }

  p4est_nearest_common_ancestor_D (q, r, &s);

  return p4est_quadrant_is_equal (q, &s);
}

int
p4est_quadrant_is_next (const p4est_quadrant_t * q,
                        const p4est_quadrant_t * r)
{
  int8_t              minlevel;
  int32_t             mask;
  int64_t             i1, i2;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_valid (r));

  /* the condition q < r is checked implicitly */

  if (q->level > r->level) {
    /* check if q is the last child up to the common level */
    mask =
      (1 << (P4EST_MAXLEVEL - r->level)) - (1 << (P4EST_MAXLEVEL - q->level));
    if ((q->x & mask) != mask || (q->y & mask) != mask) {
      return 0;
    }
    minlevel = r->level;
  }
  else {
    minlevel = q->level;
  }
  i1 = p4est_quadrant_linear_id (q, minlevel);
  i2 = p4est_quadrant_linear_id (r, minlevel);

  return (i1 + 1 == i2);
}

int
p4est_quadrant_is_next_D (const p4est_quadrant_t * q,
                          const p4est_quadrant_t * r)
{
  int64_t             i1, i2;
  p4est_quadrant_t    a, b;

  /* validity of q and r is asserted in p4est_quadrant_compare */
  if (p4est_quadrant_compare (q, r) >= 0) {
    return 0;
  }

  a = *q;
  b = *r;
  while (a.level > b.level) {
    if (p4est_quadrant_child_id (&a) != 3) {
      return 0;
    }
    p4est_quadrant_parent (&a, &a);
  }
  i1 = p4est_quadrant_linear_id (&a, a.level);
  i2 = p4est_quadrant_linear_id (&b, a.level);

  return (i1 + 1 == i2);
}

void
p4est_quadrant_parent (const p4est_quadrant_t * q, p4est_quadrant_t * r)
{
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (q->level > 0);

  r->x = q->x & ~(1 << (P4EST_MAXLEVEL - q->level));
  r->y = q->y & ~(1 << (P4EST_MAXLEVEL - q->level));
  r->level = (int8_t) (q->level - 1);

  P4EST_ASSERT (p4est_quadrant_is_valid (r));
}

void
p4est_quadrant_sibling (const p4est_quadrant_t * q, p4est_quadrant_t * r,
                        int8_t sibling_id)
{
  int                 addx = (sibling_id & 0x01);
  int                 addy = (sibling_id & 0x02) >> 1;
  int                 shift = (1 << (P4EST_MAXLEVEL - q->level));

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (q->level > 0);
  P4EST_ASSERT (sibling_id >= 0 && sibling_id < 4);

  r->x = (addx ? (q->x | shift) : (q->x & ~shift));
  r->y = (addy ? (q->y | shift) : (q->y & ~shift));
  r->level = q->level;
}

void
p4est_quadrant_children (const p4est_quadrant_t * q,
                         p4est_quadrant_t * c0, p4est_quadrant_t * c1,
                         p4est_quadrant_t * c2, p4est_quadrant_t * c3)
{
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (q->level < P4EST_MAXLEVEL);

  c0->x = q->x;
  c0->y = q->y;
  c0->level = (int8_t) (q->level + 1);

  c1->x = c0->x | (1 << (P4EST_MAXLEVEL - c0->level));
  c1->y = c0->y;
  c1->level = c0->level;

  c2->x = c0->x;
  c2->y = c0->y | (1 << (P4EST_MAXLEVEL - c0->level));
  c2->level = c0->level;

  c3->x = c1->x;
  c3->y = c2->y;
  c3->level = c0->level;

  P4EST_ASSERT (p4est_quadrant_is_valid (c0));
  P4EST_ASSERT (p4est_quadrant_is_valid (c1));
  P4EST_ASSERT (p4est_quadrant_is_valid (c2));
  P4EST_ASSERT (p4est_quadrant_is_valid (c3));
}

void
p4est_quadrant_first_descendent (const p4est_quadrant_t * q,
                                 p4est_quadrant_t * fd, int8_t level)
{
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (q->level <= level && level <= P4EST_MAXLEVEL);

  fd->x = q->x;
  fd->y = q->y;
  fd->level = level;
}

void
p4est_quadrant_last_descendent (const p4est_quadrant_t * q,
                                p4est_quadrant_t * ld, int8_t level)
{
  int32_t             shift;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (q->level <= level && level <= P4EST_MAXLEVEL);

  shift =
    (1 << (P4EST_MAXLEVEL - q->level)) - (1 << (P4EST_MAXLEVEL - level));

  ld->x = q->x + shift;
  ld->y = q->y + shift;
  ld->level = level;
}

void
p4est_nearest_common_ancestor (const p4est_quadrant_t * q1,
                               const p4est_quadrant_t * q2,
                               p4est_quadrant_t * r)
{
  int32_t             exclorx, exclory;
  int32_t             maxclor, maxlevel;

  P4EST_ASSERT (p4est_quadrant_is_valid (q1));
  P4EST_ASSERT (p4est_quadrant_is_valid (q2));

  exclorx = q1->x ^ q2->x;
  exclory = q1->y ^ q2->y;
  maxclor = exclorx | exclory;
  maxlevel = P4EST_LOG2_32 (maxclor) + 1;

  r->x = q1->x & ~((1 << maxlevel) - 1);
  r->y = q1->y & ~((1 << maxlevel) - 1);
  r->level = (int8_t) P4EST_MIN (P4EST_MAXLEVEL - maxlevel,
                                 P4EST_MIN (q1->level, q2->level));

  P4EST_ASSERT (p4est_quadrant_is_valid (r));
}

void
p4est_nearest_common_ancestor_D (const p4est_quadrant_t * q1,
                                 const p4est_quadrant_t * q2,
                                 p4est_quadrant_t * r)
{
  p4est_quadrant_t    s1 = *q1;
  p4est_quadrant_t    s2 = *q2;

  P4EST_ASSERT (p4est_quadrant_is_valid (q1));
  P4EST_ASSERT (p4est_quadrant_is_valid (q2));

  /* first stage: promote the deepest one to the same level */
  while (s1.level > s2.level) {
    p4est_quadrant_parent (&s1, &s1);
  }
  while (s1.level < s2.level) {
    p4est_quadrant_parent (&s2, &s2);
  }

  /* second stage: simultaneously go through their parents */
  while (!p4est_quadrant_is_equal (&s1, &s2)) {
    p4est_quadrant_parent (&s1, &s1);
    p4est_quadrant_parent (&s2, &s2);
  }

  /* don't overwrite r's user_data */
  r->x = s1.x;
  r->y = s1.y;
  r->level = s1.level;

  P4EST_ASSERT (p4est_quadrant_is_valid (r));
}

int64_t
p4est_quadrant_linear_id (const p4est_quadrant_t * quadrant, int8_t level)
{
  int8_t              i;
  int64_t             x, y;
  int64_t             id;

  P4EST_ASSERT (p4est_quadrant_is_valid (quadrant));
  P4EST_ASSERT (quadrant->level >= level && level >= 0);

  x = quadrant->x >> (P4EST_MAXLEVEL - level);
  y = quadrant->y >> (P4EST_MAXLEVEL - level);

  id = 0;
  for (i = 0; i < level; ++i) {
    id |= ((x & (1 << i)) << i);
    id |= ((y & (1 << i)) << (i + 1));
  }

  return id;
}

void
p4est_quadrant_set_morton (p4est_quadrant_t * quadrant,
                           int8_t level, int64_t id)
{
  int8_t              i;

  P4EST_ASSERT (0 <= level && level <= P4EST_MAXLEVEL);
  P4EST_ASSERT (id < (1LL << (2 * level)));

  quadrant->level = level;
  quadrant->x = 0;
  quadrant->y = 0;

  for (i = 0; i < level; ++i) {
    quadrant->x |= (int32_t) ((id & (1LL << (2 * i))) >> i);
    quadrant->y |= (int32_t) ((id & (1LL << (2 * i + 1))) >> (i + 1));
  }

  quadrant->x <<= (P4EST_MAXLEVEL - level);
  quadrant->y <<= (P4EST_MAXLEVEL - level);

  P4EST_ASSERT (p4est_quadrant_is_valid (quadrant));
}

void
p4est_quadrant_init_data (p4est_t * p4est, int32_t which_tree,
                          p4est_quadrant_t * quad, p4est_init_t init_fn)
{
  P4EST_ASSERT (p4est_quadrant_is_valid (quad));

  if (p4est->data_size > 0) {
    quad->user_data = p4est_mempool_alloc (p4est->user_data_pool);
  }
  else {
    quad->user_data = NULL;
  }
  if (init_fn != NULL) {
    init_fn (p4est, which_tree, quad);
  }
}

void
p4est_quadrant_free_data (p4est_t * p4est, p4est_quadrant_t * quad)
{
  P4EST_ASSERT (p4est_quadrant_is_valid (quad));

  if (p4est->data_size > 0) {
    p4est_mempool_free (p4est->user_data_pool, quad->user_data);
  }
  quad->user_data = NULL;
}

void
p4est_quadrant_print (const p4est_quadrant_t * q, int identifier, FILE * nout)
{
  char                prefix[BUFSIZ];

  if (nout == NULL) {
    return;
  }

  if (identifier >= 0) {
    snprintf (prefix, BUFSIZ, "[%d] ", identifier);
  }
  else {
    prefix[0] = '\0';
  }

  fprintf (nout, "%sx 0x%x y 0x%x level %d\n", prefix, q->x, q->y, q->level);
}

int
p4est_tree_is_sorted (p4est_tree_t * tree)
{
  int                 i;
  p4est_quadrant_t   *q1, *q2;

  if (tree->quadrants->elem_count <= 1) {
    return 1;
  }

  q1 = p4est_array_index (tree->quadrants, 0);
  for (i = 1; i < tree->quadrants->elem_count; ++i) {
    q2 = p4est_array_index (tree->quadrants, i);
    if (p4est_quadrant_compare (q1, q2) >= 0) {
      return 0;
    }
    q1 = q2;
  }

  return 1;
}

int
p4est_tree_is_linear (p4est_tree_t * tree)
{
  int                 i;
  p4est_quadrant_t   *q1, *q2;

  if (tree->quadrants->elem_count <= 1) {
    return 1;
  }

  q1 = p4est_array_index (tree->quadrants, 0);
  for (i = 1; i < tree->quadrants->elem_count; ++i) {
    q2 = p4est_array_index (tree->quadrants, i);
    if (p4est_quadrant_compare (q1, q2) >= 0) {
      return 0;
    }
    if (p4est_quadrant_is_ancestor (q1, q2)) {
      return 0;
    }
    q1 = q2;
  }

  return 1;
}

int
p4est_tree_is_complete (p4est_tree_t * tree)
{
  int                 i;
  p4est_quadrant_t   *q1, *q2;

  if (tree->quadrants->elem_count <= 1) {
    return 1;
  }

  q1 = p4est_array_index (tree->quadrants, 0);
  for (i = 1; i < tree->quadrants->elem_count; ++i) {
    q2 = p4est_array_index (tree->quadrants, i);
    if (!p4est_quadrant_is_next (q1, q2)) {
      return 0;
    }
    q1 = q2;
  }

  return 1;
}

void
p4est_tree_print (p4est_tree_t * tree, int identifier, FILE * nout)
{
  int                 j, childid, comp;
  char                prefix[BUFSIZ];
  p4est_quadrant_t   *q1, *q2;

  if (nout == NULL) {
    return;
  }

  if (identifier >= 0) {
    snprintf (prefix, BUFSIZ, "[%d] ", identifier);
  }
  else {
    prefix[0] = '\0';
  }

  q1 = NULL;
  for (j = 0; j < tree->quadrants->elem_count; ++j) {
    q2 = p4est_array_index (tree->quadrants, j);
    childid = p4est_quadrant_child_id (q2);
    fprintf (nout, "%s0x%x 0x%x %d", prefix, q2->x, q2->y, q2->level);
    if (j > 0) {
      comp = p4est_quadrant_compare (q1, q2);
      if (comp > 0) {
        fputs (" R", nout);
      }
      else if (comp == 0) {
        fputs (" I", nout);
      }
      else {
        if (p4est_quadrant_is_sibling (q1, q2)) {
          fprintf (nout, " S%d", childid);
        }
        else if (p4est_quadrant_is_parent (q1, q2)) {
          fprintf (nout, " C%d", childid);
        }
        else if (p4est_quadrant_is_ancestor (q1, q2)) {
          fputs (" D", nout);
        }
        else if (p4est_quadrant_is_next (q1, q2)) {
          fprintf (nout, " N%d", childid);
        }
        else {
          fprintf (nout, " q%d", childid);
        }
      }
    }
    else {
      fprintf (nout, " F%d", childid);
    }
    fputs ("\n", nout);
    q1 = q2;
  }
}

int
p4est_is_valid (p4est_t * p4est)
{
  int                 i;
  int8_t              maxlevel;
  int32_t             j, nquadrants, lquadrants, perlevel;
  p4est_tree_t       *tree;

  lquadrants = 0;
  for (j = 0; j < p4est->trees->elem_count; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    if (!p4est_tree_is_complete (tree)) {
      printf ("not complete\n");
      return 0;
    }
    if ((j < p4est->first_local_tree || j > p4est->last_local_tree) &&
        tree->quadrants->elem_count > 0) {
      printf ("outside count\n");
      return 0;
    }

    maxlevel = 0;
    nquadrants = 0;
    for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
      perlevel = tree->quadrants_per_level[i];

      P4EST_ASSERT (perlevel >= 0);
      nquadrants += perlevel;
      if (perlevel > 0) {
        maxlevel = (int8_t) i;
      }
    }
    lquadrants += nquadrants;

    if (maxlevel != tree->maxlevel) {
      printf ("wrong maxlevel\n");
      return 0;
    }
    if (nquadrants != tree->quadrants->elem_count) {
      printf ("wrong tree quadrant count\n");
      return 0;
    }
  }

  if (lquadrants != p4est->local_num_quadrants) {
    printf ("wrong local quadrant count\n");
    return 0;
  }

  return 1;
}

/* here come the heavyweight algorithms */

/** Find the lowest position tq in a quadrant array such that tq >= q.
 * \return  Returns the guess id if found, -1 otherwise.
 */
static int
p4est_find_lower_bound (p4est_array_t * array,
                        const p4est_quadrant_t * q, int guess)
{
  int                 count, comp;
  int                 quad_low, quad_high;
  p4est_quadrant_t   *cur;

  count = array->elem_count;
  quad_low = 0;
  quad_high = count - 1;

  for (;;) {
    P4EST_ASSERT (quad_low <= quad_high);
    P4EST_ASSERT (0 <= quad_low && quad_low < count);
    P4EST_ASSERT (0 <= quad_high && quad_high < count);
    P4EST_ASSERT (quad_low <= guess && guess <= quad_high);

    /* compare two quadrants */
    cur = p4est_array_index (array, guess);
    comp = p4est_quadrant_compare (q, cur);

    /* check if guess is higher or equal q and there's room below it */
    if (comp <= 0 && (guess > 0 && p4est_quadrant_compare (q, cur - 1) <= 0)) {
      quad_high = guess - 1;
      guess = (quad_low + quad_high + 1) / 2;
      continue;
    }

    /* check if guess is lower than q */
    if (comp > 0) {
      quad_low = guess + 1;
      if (quad_low > quad_high) {
        return -1;
      }
      guess = (quad_low + quad_high) / 2;
      continue;
    }

    /* otherwise guess is the correct quadrant */
    break;
  }

  return guess;
}

/** Find the highest position tq in a quadrant array such that tq <= q.
 * \return  Returns the id of the matching quadrant, or -1 if not found.
 */
static int
p4est_find_higher_bound (p4est_array_t * array,
                         const p4est_quadrant_t * q, int guess)
{
  int                 count, comp;
  int                 quad_low, quad_high;
  p4est_quadrant_t   *cur;

  count = array->elem_count;
  quad_low = 0;
  quad_high = count - 1;

  for (;;) {
    P4EST_ASSERT (quad_low <= quad_high);
    P4EST_ASSERT (0 <= quad_low && quad_low < count);
    P4EST_ASSERT (0 <= quad_high && quad_high < count);
    P4EST_ASSERT (quad_low <= guess && guess <= quad_high);

    /* compare two quadrants */
    cur = p4est_array_index (array, guess);
    comp = p4est_quadrant_compare (cur, q);

    /* check if guess is lower or equal q and there's room above it */
    if (comp <= 0 &&
        (guess < count - 1 && p4est_quadrant_compare (cur + 1, q) <= 0)) {
      quad_low = guess + 1;
      guess = (quad_low + quad_high) / 2;
      continue;
    }

    /* check if guess is higher than q */
    if (comp > 0) {
      quad_high = guess - 1;
      if (quad_high < quad_low) {
        return -1;
      }
      guess = (quad_low + quad_high + 1) / 2;
      continue;
    }

    /* otherwise guess is the correct quadrant */
    break;
  }

  return guess;
}

void
p4est_tree_compute_overlap (p4est_tree_t * tree, p4est_array_t * in,
                            p4est_array_t * out)
{
  int                 i, j, guess;
  int                 treecount, incount, outcount;
  int                 first_index, last_index;
  int32_t             qh, rh;
  p4est_quadrant_t    treefd, treeld;
  p4est_quadrant_t    low_ins, high_ins, fd, ld;
  p4est_quadrant_t   *tq;
  p4est_quadrant_t   *inq, *outq;

  P4EST_ASSERT (p4est_tree_is_complete (tree));
  P4EST_ASSERT (out->elem_count == 0);

  /* assign some numbers */
  treecount = tree->quadrants->elem_count;
  incount = in->elem_count;
  outcount = 0;
  rh = (1 << P4EST_MAXLEVEL);

  /* return if there is nothing to do */
  if (treecount == 0 || incount == 0) {
    return;
  }

  printf ("Into compute overlap with %d %d\n", treecount, incount);

  /* compute first and last descendants in the tree */
  tq = p4est_array_index (tree->quadrants, 0);
  p4est_quadrant_first_descendent (tq, &treefd, P4EST_MAXLEVEL);
  tq = p4est_array_index (tree->quadrants, treecount - 1);
  p4est_quadrant_last_descendent (tq, &treeld, P4EST_MAXLEVEL);

  /* loop over input list of quadrants */
  for (i = 0; i < incount; ++i) {
    inq = p4est_array_index (in, i);
    qh = (1 << (P4EST_MAXLEVEL - inq->level));

    /* compute first descendent of first insulation quadrant */
    low_ins = *inq;
    if (low_ins.x > 0) {
      low_ins.x -= qh;
      P4EST_ASSERT (low_ins.x >= 0);
    }
    if (low_ins.y > 0) {
      low_ins.y -= qh;
      P4EST_ASSERT (low_ins.y >= 0);
    }
    p4est_quadrant_first_descendent (&low_ins, &fd, P4EST_MAXLEVEL);

    /* compute last descendent of last insulation quadrant */
    high_ins = *inq;
    if (high_ins.x + qh < rh) {
      high_ins.x += qh;
      P4EST_ASSERT (high_ins.x + qh <= rh);
    }
    if (high_ins.y + qh < rh) {
      high_ins.y += qh;
      P4EST_ASSERT (high_ins.y + qh <= rh);
    }
    p4est_quadrant_last_descendent (&high_ins, &ld, P4EST_MAXLEVEL);

    /* skip this insulation layer if there is no overlap */
    if (p4est_quadrant_compare (&ld, &treefd) < 0 ||
        p4est_quadrant_compare (&treeld, &fd) < 0) {
      continue;
    }

    /* find first quadrant in tree that fits between fd and ld */
    guess = treecount / 2;
    if (p4est_quadrant_compare (&fd, &treefd) <= 0) {
      /* the first tree quadrant is contained in insulation layer */
      first_index = 0;
    }
    else {
      /* do a binary search for the lowest tree quadrant >= low_ins */
      first_index = p4est_find_lower_bound (tree->quadrants, &low_ins, guess);
      if (first_index < 0) {
        continue;
      }
      guess = first_index;
    }

    /* find last quadrant in tree that fits between fd and ld */
    if (p4est_quadrant_compare (&treeld, &ld) <= 0) {
      /* the last tree quadrant is contained in insulation layer */
      last_index = treecount - 1;
    }
    else {
      /* do a binary search for the highest tree quadrant <= ld */
      last_index = p4est_find_higher_bound (tree->quadrants, &ld, guess);
      if (last_index < 0) {
        printf ("Last index < 0\n");
        continue;
      }
    }
    P4EST_ASSERT (first_index <= last_index);

    /* copy all overlapping quadrants that are small enough into out */
    for (j = first_index; j <= last_index; ++j) {
      tq = p4est_array_index (tree->quadrants, j);
      if (tq->level > inq->level + 1) {
        p4est_array_resize (out, outcount + 1);
        outq = p4est_array_index (out, outcount);
        *outq = *tq;
        ++outcount;
      }
    }
  }
  if (outcount == 0) {
    return;
  }

  /* sort array and remove duplicates */
  p4est_array_sort (out, p4est_quadrant_compare);
  i = 0;                        /* read counter */
  j = 0;                        /* write counter */
  inq = p4est_array_index (out, i);
  while (i < outcount) {
    tq = ((i < outcount - 1) ? p4est_array_index (out, i + 1) : NULL);
    if (i < outcount - 1 && p4est_quadrant_is_equal (inq, tq)) {
      ++i;
    }
    else {
      if (i > j) {
        outq = p4est_array_index (out, j);
        *outq = *inq;
      }
      ++i;
      ++j;
    }
    inq = tq;
  }
  P4EST_ASSERT (i == outcount);
  P4EST_ASSERT (j <= outcount);
  p4est_array_resize (out, j);
}

void
p4est_complete_region (p4est_t * p4est,
                       const p4est_quadrant_t * q1,
                       int include_q1,
                       const p4est_quadrant_t * q2,
                       int include_q2,
                       p4est_tree_t * tree,
                       int32_t which_tree, p4est_init_t init_fn)
{
  p4est_tree_t       *R;
  p4est_list_t       *W;

  p4est_quadrant_t    a = *q1;
  p4est_quadrant_t    b = *q2;

  p4est_quadrant_t    Afinest;
  p4est_quadrant_t   *c0, *c1, *c2, *c3;

  p4est_array_t      *quadrants;
  p4est_mempool_t    *quadrant_pool = p4est->quadrant_pool;

  p4est_quadrant_t   *w;
  p4est_quadrant_t   *r;

  int                 comp;
  int                 quadrant_pool_size, data_pool_size;
  int8_t              level;
  int8_t              maxlevel = 0;
  int32_t            *quadrants_per_level;
  int32_t             num_quadrants = 0;

  W = p4est_list_new (NULL);
  R = tree;

  /* needed for sanity check */
  quadrant_pool_size = p4est->quadrant_pool->elem_count;
  if (p4est->user_data_pool != NULL) {
    data_pool_size = p4est->user_data_pool->elem_count;
  }

  quadrants = R->quadrants;
  quadrants_per_level = R->quadrants_per_level;

  /* Assert that we have an empty tree */
  P4EST_ASSERT (quadrants->elem_count == 0);

  comp = p4est_quadrant_compare (&a, &b);
  /* Assert that a<b */
  P4EST_ASSERT (comp < 0);

  /* R <- R + a */
  if (include_q1) {
    p4est_array_resize (quadrants, 1);
    r = p4est_array_index (quadrants, 0);
    *r = a;
    maxlevel = (int8_t) P4EST_MAX (a.level, maxlevel);
    ++quadrants_per_level[a.level];
    ++num_quadrants;
  }

  if (comp < 0) {
    /* W <- C(A_{finest}(a,b)) */
    p4est_nearest_common_ancestor (&a, &b, &Afinest);

    c0 = p4est_mempool_alloc (quadrant_pool);
    c1 = p4est_mempool_alloc (quadrant_pool);
    c2 = p4est_mempool_alloc (quadrant_pool);
    c3 = p4est_mempool_alloc (quadrant_pool);

    p4est_quadrant_children (&Afinest, c0, c1, c2, c3);

    p4est_list_append (W, c0);
    p4est_list_append (W, c1);
    p4est_list_append (W, c2);
    p4est_list_append (W, c3);

    /* for each w in W */
    while (W->elem_count > 0) {
      w = p4est_list_pop (W);
      level = w->level;

      /* if (a < w < b) and (w not in {A(b)}) */
      if (((p4est_quadrant_compare (&a, w) < 0) &&
           (p4est_quadrant_compare (w, &b) < 0)
          ) && !p4est_quadrant_is_ancestor (w, &b)
        ) {
        /* R <- R + w */
        p4est_array_resize (quadrants, num_quadrants + 1);
        r = p4est_array_index (quadrants, num_quadrants);
        *r = *w;
        p4est_quadrant_init_data (p4est, which_tree, r, init_fn);
        maxlevel = (int8_t) P4EST_MAX (level, maxlevel);
        ++quadrants_per_level[level];
        ++num_quadrants;
      }
      /* else if (w in {{A(a)}, {A(b)}}) */
      else if (p4est_quadrant_is_ancestor (w, &a)
               || p4est_quadrant_is_ancestor (w, &b)) {
        /* W <- W + C(w) */
        c0 = p4est_mempool_alloc (quadrant_pool);
        c1 = p4est_mempool_alloc (quadrant_pool);
        c2 = p4est_mempool_alloc (quadrant_pool);
        c3 = p4est_mempool_alloc (quadrant_pool);

        p4est_quadrant_children (w, c0, c1, c2, c3);

        p4est_list_prepend (W, c3);
        p4est_list_prepend (W, c2);
        p4est_list_prepend (W, c1);
        p4est_list_prepend (W, c0);
      }

      /* W <- W - w */
      p4est_mempool_free (quadrant_pool, w);
    }                           /* end for */

    /* R <- R + b */
    if (include_q2) {
      p4est_array_resize (quadrants, num_quadrants + 1);
      r = p4est_array_index (quadrants, num_quadrants);
      *r = b;
      maxlevel = (int8_t) P4EST_MAX (b.level, maxlevel);
      ++quadrants_per_level[b.level];
      ++num_quadrants;
    }
  }

  R->maxlevel = maxlevel;

  P4EST_ASSERT (W->first == NULL && W->last == NULL);
  p4est_list_destroy (W);

  P4EST_ASSERT (p4est_tree_is_complete (R));
  P4EST_ASSERT (quadrant_pool_size == p4est->quadrant_pool->elem_count);
  P4EST_ASSERT (num_quadrants == quadrants->elem_count);
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size + quadrants->elem_count ==
                  p4est->user_data_pool->elem_count + (include_q1 ? 1 : 0)
                  + (include_q2 ? 1 : 0));
  }
}

/** Internal function to realize local completion / balancing.
 * \param [in] balance  can be 0: no balancing
 *                             1: balance across edges
 *                             2: balance across edges and corners
 */
static void
p4est_complete_or_balance (p4est_t * p4est, p4est_tree_t * tree, int balance,
                           int32_t which_tree, p4est_init_t init_fn)
{
  int                 i, j;
  int                 incount, curcount, ocount;
  int                 comp, lookup, inserted, isfamily;
  int                 quadrant_pool_size, data_pool_size;
  int                 hash_size;
  int                 count_outside_root, count_outside_tree;
  int                 count_already_inlist, count_already_outlist;
  int                *key, *parent_key;
  int8_t              l, inmaxl, bbound;
  int8_t              qid, sid, pid;
  int32_t             ph, rh;
  void               *vlookup;
  p4est_quadrant_t   *family[4];
  p4est_quadrant_t   *q, *r;
  p4est_quadrant_t   *qalloc, *qlookup, **qpointer;
  p4est_quadrant_t    ld, tree_first, tree_last, parent;
  p4est_array_t      *inlist, *olist;
  p4est_mempool_t    *list_alloc, *qpool;
  p4est_hash_t       *hash[P4EST_MAXLEVEL + 1];
  p4est_array_t      *outlist[P4EST_MAXLEVEL + 1];

  P4EST_ASSERT (p4est_tree_is_sorted (tree));

  /*
   * Algorithm works with these data structures
   * inlist  --  sorted list of input quadrants
   * hash    --  hash table to hold additional quadrants not in inlist
   *             this is filled bottom-up to ensure balance condition
   * outlist --  filled simultaneously with hash, holding pointers
   *             don't rely on addresses of elements, it is resized
   * In the end, the elements of hash are appended to inlist
   * and inlist is sorted and linearized. This can be optimized later.
   */

  /* assign some shortcut variables */
  bbound = (int8_t) ((balance == 0) ? 5 : 8);
  inlist = tree->quadrants;
  incount = inlist->elem_count;
  inmaxl = tree->maxlevel;
  qpool = p4est->quadrant_pool;
  key = &comp;                  /* unique user_data pointer */
  parent_key = &lookup;

  /* needed for sanity check */
  quadrant_pool_size = qpool->elem_count;
  if (p4est->user_data_pool != NULL) {
    data_pool_size = p4est->user_data_pool->elem_count;
  }

  /* if tree is empty or a single block, there is nothing to do */
  if (incount <= 1) {
    return;
  }

  /* determine the first and last small quadrants contained in the tree */
  q = p4est_array_index (inlist, 0);
  p4est_quadrant_first_descendent (q, &tree_first, inmaxl);
  p4est_quadrant_last_descendent (q, &tree_last, inmaxl);
  for (i = 1; i < incount; ++i) {
    q = p4est_array_index (inlist, i);
    p4est_quadrant_last_descendent (q, &ld, inmaxl);
    comp = p4est_quadrant_compare (&tree_last, &ld);
    if (comp < 0) {
      tree_last = ld;
    }
  }
  /*
     if (p4est->nout != NULL) {
     fprintf (p4est->nout, "[%d] First descendent 0x%x 0x%x %d\n",
     p4est->mpirank, tree_first.x, tree_first.y, tree_first.level);
     fprintf (p4est->nout, "[%d] Last descendent 0x%x 0x%x %d\n",
     p4est->mpirank, tree_last.x, tree_last.y, tree_last.level);
     }
   */

  /* initialize some counters */
  count_outside_root = count_outside_tree = 0;
  count_already_inlist = count_already_outlist = 0;

  /* initialize temporary storage */
  list_alloc = p4est_mempool_new (sizeof (p4est_link_t));
  for (l = 0; l <= inmaxl; ++l) {
    hash_size = P4EST_MAX (incount / (P4EST_MAXLEVEL * 2),
                           tree->quadrants_per_level[l] * 2 - 1);
    hash_size = P4EST_MIN (hash_table_maxsize, hash_size);
    hash_size = P4EST_MAX (hash_table_minsize, hash_size);
    hash[l] = p4est_hash_new (hash_size, p4est_quadrant_hash,
                              p4est_quadrant_is_equal, list_alloc);
    outlist[l] = p4est_array_new (sizeof (p4est_quadrant_t *));
  }
  for (l = (int8_t) (inmaxl + 1); l <= P4EST_MAXLEVEL; ++l) {
    hash[l] = NULL;
    outlist[l] = NULL;
  }

  /* walk through the input tree bottom-up */
  ph = 0;
  pid = -1;
  parent.x = -1;
  parent.y = -1;
  parent.level = -1;
  qalloc = p4est_mempool_alloc (qpool);
  qalloc->user_data = key;
  rh = (1 << P4EST_MAXLEVEL);
  for (l = inmaxl; l > 0; --l) {
    ocount = outlist[l]->elem_count;    /* fix ocount here, it is growing */
    for (i = 0; i < incount + ocount; ++i) {
      isfamily = 0;
      if (i < incount) {
        q = p4est_array_index (inlist, i);
        if (q->level != l) {
          continue;
        }
        /* this is an optimization to catch adjacent siblings */
        if (i + 4 <= incount) {
          family[0] = q;
          for (j = 1; j < 4; ++j) {
            family[j] = p4est_array_index (inlist, i + j);
          }
          if (p4est_quadrant_is_family (family[0], family[1],
                                        family[2], family[3])) {
            isfamily = 1;
            i += 3;             /* skip siblings */
          }
        }
      }
      else {
        qpointer = p4est_array_index (outlist[l], i - incount);
        q = *qpointer;
        P4EST_ASSERT (q->level == l);
      }

      /*
       * check for q and its siblings,
       * then for q's parent and parent's indirect relevant neighbors
       * sid == 0..3  siblings including q
       *        4     parent of q
       *        5..7  relevant indirect neighbors of parent
       */
      qid = p4est_quadrant_child_id (q);        /* 0 <= qid < 4 */
      for (sid = 0; sid < bbound; ++sid) {
        /* stage 1: determine candidate qalloc */
        if (sid < 4) {
          if (qid == sid || isfamily) {
            /* q (or its family) is included in inlist */
            continue;
          }
          p4est_quadrant_sibling (q, qalloc, sid);
        }
        else if (sid == 4) {
          /* compute the parent */
          p4est_quadrant_parent (q, qalloc);
          if (bbound > 5) {
            parent = *qalloc;   /* copy parent for cases 5..7 */
            ph = (1 << (P4EST_MAXLEVEL - parent.level));        /* its size */
            pid = p4est_quadrant_child_id (&parent);    /* and position */
          }
        }
        else {
          /* determine the 3 parent's relevant indirect neighbors */
          P4EST_ASSERT (sid >= 5 && sid < 8);
          if (balance < 2 && sid - 5 == corners_omitted[pid]) {
            /* this quadrant would only be needed for corner balance */
            continue;
          }
          qalloc->x = parent.x + indirect_neighbors[pid][sid - 5][0] * ph;
          qalloc->y = parent.y + indirect_neighbors[pid][sid - 5][1] * ph;
          qalloc->level = parent.level;
          if ((qalloc->x < 0 || qalloc->x >= rh) ||
              (qalloc->y < 0 || qalloc->y >= rh)) {
            /* quadrant is outside the root */
            ++count_outside_root;
            continue;
          }
        }
        /*
           printf ("Candidate level %d qxy 0x%x 0x%x at sid %d\n",
           qalloc->level, qalloc->x, qalloc->y, sid);
         */

        /* stage 2: include qalloc if necessary */
        p4est_quadrant_last_descendent (qalloc, &ld, inmaxl);
        if ((p4est_quadrant_compare (&tree_first, qalloc) > 0 &&
             (qalloc->x != tree_first.x || qalloc->y != tree_first.y)) ||
            p4est_quadrant_compare (&ld, &tree_last) > 0) {
          /* qalloc is outside the tree */
          ++count_outside_tree;
          continue;
        }
        lookup = p4est_hash_lookup (hash[qalloc->level], qalloc, &vlookup);
        if (lookup) {
          /* qalloc is already included in output list, this catches most */
          ++count_already_outlist;
          qlookup = vlookup;
          if (sid == 4 && qlookup->user_data == parent_key) {
            break;              /* this parent has been triggered before */
          }
          continue;
        }
        r = p4est_array_bsearch (inlist, qalloc, p4est_quadrant_compare);
        if (r != NULL) {
          /* qalloc is included in inlist, this is more expensive to test */
          ++count_already_inlist;
          continue;
        }
        /* insert qalloc into the output list as well */
        if (sid == 4) {
          qalloc->user_data = parent_key;
        }
        inserted = p4est_hash_insert_unique (hash[qalloc->level], qalloc,
                                             NULL);
        P4EST_ASSERT (inserted);
        olist = outlist[qalloc->level];
        p4est_array_resize (olist, olist->elem_count + 1);
        qpointer = p4est_array_index (olist, olist->elem_count - 1);
        *qpointer = qalloc;
        /* we need a new quadrant now, the old one is stored away */
        qalloc = p4est_mempool_alloc (qpool);
        qalloc->user_data = key;
      }
    }
  }
  p4est_mempool_free (qpool, qalloc);

  /* merge outlist into input list and free temporary storage */
  curcount = 0;
  for (l = 0; l <= inmaxl; ++l) {
    /* print statistics and free hash tables */
#ifdef P4EST_HAVE_DEBUG
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Tree %d Level %d ",
               p4est->mpirank, which_tree, l);
    }
    p4est_hash_print_statistics (hash[l], p4est->nout);
#endif
    p4est_hash_unlink_destroy (hash[l]);        /* performance optimization */

    /* merge outlist into inlist */
    curcount = inlist->elem_count;
    ocount = outlist[l]->elem_count;
    p4est_array_resize (inlist, curcount + ocount);
    for (i = 0; i < ocount; ++i) {
      /* copy and free temporary quadrants */
      q = p4est_array_index (inlist, curcount + i);
      qpointer = p4est_array_index (outlist[l], i);
      qalloc = *qpointer;
      P4EST_ASSERT (qalloc->level == l);
      P4EST_ASSERT (qalloc->user_data == key ||
                    qalloc->user_data == parent_key);
      *q = *qalloc;
      p4est_mempool_free (qpool, qalloc);

      /* complete quadrant initialization */
      p4est_quadrant_init_data (p4est, which_tree, q, init_fn);
    }
    tree->quadrants_per_level[l] += ocount;
    if (ocount > 0 && l > tree->maxlevel) {
      tree->maxlevel = l;
    }
    p4est_array_destroy (outlist[l]);
  }
  p4est_mempool_destroy (list_alloc);

  /* print more statistics */
  if (p4est->nout != NULL) {
    fprintf (p4est->nout, "[%d] Tree %d Outside root %d tree %d\n",
             p4est->mpirank, which_tree,
             count_outside_root, count_outside_tree);
    fprintf (p4est->nout, "[%d] Tree %d Already in inlist %d outlist %d\n",
             p4est->mpirank, which_tree,
             count_already_inlist, count_already_outlist);
    fprintf (p4est->nout, "[%d] Tree %d Insertions %d\n",
             p4est->mpirank, which_tree, curcount - incount);
  }

  /* sort and linearize tree */
  p4est_array_sort (inlist, p4est_quadrant_compare);
  P4EST_ASSERT (p4est_tree_is_sorted (tree));
  p4est_linearize_subtree (p4est, tree);

  /* run sanity checks */
  P4EST_ASSERT (quadrant_pool_size == qpool->elem_count);
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size + inlist->elem_count ==
                  p4est->user_data_pool->elem_count + incount);
  }
  P4EST_ASSERT (p4est_tree_is_linear (tree));
  P4EST_ASSERT (p4est_tree_is_complete (tree));
}

void
p4est_complete_subtree (p4est_t * p4est, p4est_tree_t * tree,
                        int32_t which_tree, p4est_init_t init_fn)
{
  p4est_complete_or_balance (p4est, tree, 0, which_tree, init_fn);
}

void
p4est_balance_subtree (p4est_t * p4est, p4est_tree_t * tree,
                       int32_t which_tree, p4est_init_t init_fn)
{
  p4est_complete_or_balance (p4est, tree, 2, which_tree, init_fn);
}

void
p4est_linearize_subtree (p4est_t * p4est, p4est_tree_t * tree)
{
  int                 data_pool_size;
  int                 incount, removed;
  int                 current, rest, num_quadrants;
  int8_t              i, maxlevel;
  p4est_quadrant_t   *q1, *q2;

  P4EST_ASSERT (p4est_tree_is_sorted (tree));

  incount = tree->quadrants->elem_count;
  if (incount <= 1) {
    return;
  }
  if (p4est->user_data_pool != NULL) {
    data_pool_size = p4est->user_data_pool->elem_count;
  }
  removed = 0;

  /* run through the array and remove ancestors */
  current = 0;
  rest = current + 1;
  q1 = p4est_array_index (tree->quadrants, current);
  while (rest < incount) {
    q2 = p4est_array_index (tree->quadrants, rest);
    if (p4est_quadrant_is_ancestor (q1, q2)) {
      --tree->quadrants_per_level[q1->level];
      p4est_quadrant_free_data (p4est, q1);
      *q1 = *q2;
      ++removed;
      ++rest;
    }
    else {
      ++current;
      if (current < rest) {
        q1 = p4est_array_index (tree->quadrants, current);
        *q1 = *q2;
      }
      else {
        q1 = q2;
      }
      ++rest;
    }
  }

  /* resize array */
  p4est_array_resize (tree->quadrants, current + 1);

  /* update level counters */
  maxlevel = 0;
  num_quadrants = 0;
  for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
    P4EST_ASSERT (tree->quadrants_per_level[i] >= 0);
    num_quadrants += tree->quadrants_per_level[i];
    if (tree->quadrants_per_level[i] > 0) {
      maxlevel = i;
    }
  }
  tree->maxlevel = maxlevel;

  /* sanity checks */
  P4EST_ASSERT (num_quadrants == tree->quadrants->elem_count);
  P4EST_ASSERT (tree->quadrants->elem_count == incount - removed);
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size - removed ==
                  p4est->user_data_pool->elem_count);
  }
  P4EST_ASSERT (p4est_tree_is_sorted (tree));
  P4EST_ASSERT (p4est_tree_is_linear (tree));
}

/* EOF p4est_algorithms.c */
