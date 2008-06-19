/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

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

#ifdef P4_TO_P8
#include <p8est_bits.h>
#else
#include <p4est_bits.h>
#endif /* !P4_TO_P8 */

int
p4est_quadrant_compare (const void *v1, const void *v2)
{
  const p4est_quadrant_t *q1 = v1;
  const p4est_quadrant_t *q2 = v2;

  uint32_t            exclorx, exclory;
  int                 log2x, log2y;
#ifdef P4_TO_P8
  uint32_t            exclorz;
  int                 log2z;
#endif
  int64_t             p1, p2, diff;

  P4EST_ASSERT (p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_extended (q2));

  /* these are unsigned variables that inherit the sign bits */
  exclorx = q1->x ^ q2->x;
  exclory = q1->y ^ q2->y;
#ifdef P4_TO_P8
  exclorz = q1->z ^ q2->z;
#endif

  if (exclory == 0 && exclorx == 0 &&
#ifdef P4_TO_P8
      exclorz == 0 &&
#endif
      true) {
    return (int) q1->level - (int) q2->level;
  }
  else {
    log2x = SC_LOG2_32 (exclorx);
    log2y = SC_LOG2_32 (exclory);
#ifdef P4_TO_P8
    log2z = SC_LOG2_32 (exclorz);

    if (log2z >= log2x && log2z >= log2y) {
      p1 = q1->z + ((q1->z >= 0) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2)));
      p2 = q2->z + ((q2->z >= 0) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2)));
    }
    else
#if 0
      ;                         /* let indent survive */
#endif
#endif
    if (log2y >= log2x) {
      p1 = q1->y + ((q1->y >= 0) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2)));
      p2 = q2->y + ((q2->y >= 0) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2)));
    }
    else {
      p1 = q1->x + ((q1->x >= 0) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2)));
      p2 = q2->x + ((q2->x >= 0) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2)));
    }
    diff = p1 - p2;
    return (diff == 0) ? 0 : ((diff < 0) ? -1 : 1);
  }
}

int
p4est_quadrant_compare_piggy (const void *v1, const void *v2)
{
  const p4est_quadrant_t *q1 = v1;
  const p4est_quadrant_t *q2 = v2;

  /* expect non-negative tree information */
  /* *INDENT-OFF* horrible indent bug */
  const p4est_topidx_t diff =
    q1->p.piggy.which_tree - q2->p.piggy.which_tree;  /* same type */
  /* *INDENT-ON* */

  return (diff == 0) ?
    p4est_quadrant_compare (v1, v2) : ((diff < 0) ? -1 : 1);
}

bool
p4est_quadrant_is_equal (const void *v1, const void *v2)
{
  const p4est_quadrant_t *q1 = v1;
  const p4est_quadrant_t *q2 = v2;

  P4EST_ASSERT (p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_extended (q2));

  return (q1->level == q2->level && q1->x == q2->x && q1->y == q2->y)
#ifdef P4_TO_P8
    && (q1->z == q2->z)
#endif
    ;
}

unsigned
p4est_quadrant_hash (const void *v)
{
  const p4est_quadrant_t *q = v;

  return (unsigned) (p4est_quadrant_linear_id (q, (int) q->level) %
                     ((uint64_t) 1 << 30));
}

int
p4est_quadrant_child_id (const p4est_quadrant_t * q)
{
  int                 id = 0;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));

  if (q->level == 0) {
    return 0;
  }

  id |= ((q->x & P4EST_QUADRANT_LEN (q->level)) ? 0x01 : 0);
  id |= ((q->y & P4EST_QUADRANT_LEN (q->level)) ? 0x02 : 0);
#ifdef P4_TO_P8
  id |= ((q->z & P4EST_QUADRANT_LEN (q->level)) ? 0x04 : 0);
#endif

  return id;
}

bool
p4est_quadrant_is_inside_root (const p4est_quadrant_t * q)
{
  return
    (q->x >= 0 && q->x < P4EST_ROOT_LEN) &&
    (q->y >= 0 && q->y < P4EST_ROOT_LEN) &&
#ifdef P4_TO_P8
    (q->z >= 0 && q->z < P4EST_ROOT_LEN) &&
#endif
    true;
}

bool
p4est_quadrant_is_inside_3x3 (const p4est_quadrant_t * q)
{
  return
    (q->x >= -P4EST_ROOT_LEN &&
     q->x <= P4EST_ROOT_LEN + (P4EST_ROOT_LEN - 1)) &&
    (q->y >= -P4EST_ROOT_LEN &&
     q->y <= P4EST_ROOT_LEN + (P4EST_ROOT_LEN - 1)) &&
#ifdef P4_TO_P8
    (q->z >= -P4EST_ROOT_LEN &&
     q->z <= P4EST_ROOT_LEN + (P4EST_ROOT_LEN - 1)) &&
#endif
    true;
}

bool
p4est_quadrant_is_valid (const p4est_quadrant_t * q)
{
  return
    (q->level >= 0 && q->level <= P4EST_MAXLEVEL) &&
    ((q->x & (P4EST_QUADRANT_LEN (q->level) - 1)) == 0) &&
    ((q->y & (P4EST_QUADRANT_LEN (q->level) - 1)) == 0) &&
#ifdef P4_TO_P8
    ((q->z & (P4EST_QUADRANT_LEN (q->level) - 1)) == 0) &&
#endif
    p4est_quadrant_is_inside_root (q);
}

bool
p4est_quadrant_is_extended (const p4est_quadrant_t * q)
{
  return
    (q->level >= 0 && q->level <= P4EST_MAXLEVEL) &&
    ((q->x & (P4EST_QUADRANT_LEN (q->level) - 1)) == 0) &&
    ((q->y & (P4EST_QUADRANT_LEN (q->level) - 1)) == 0) &&
#ifdef P4_TO_P8
    ((q->z & (P4EST_QUADRANT_LEN (q->level) - 1)) == 0) &&
#endif
    p4est_quadrant_is_inside_3x3 (q);
}

bool
p4est_quadrant_is_sibling (const p4est_quadrant_t * q1,
                           const p4est_quadrant_t * q2)
{
  p4est_qcoord_t      exclorx, exclory;
#ifdef P4_TO_P8
  p4est_qcoord_t      exclorz;
#endif

  P4EST_ASSERT (p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_extended (q2));

  if (q1->level == 0) {
    return false;
  }

  exclorx = q1->x ^ q2->x;
  exclory = q1->y ^ q2->y;
#ifdef P4_TO_P8
  exclorz = q1->z ^ q2->z;

  if (exclorx == 0 && exclory == 0 && exclorz == 0) {
    return false;
  }
#else
  if (exclorx == 0 && exclory == 0) {
    return false;
  }
#endif

  return
    (q1->level == q2->level) &&
    ((exclorx & ~P4EST_QUADRANT_LEN (q1->level)) == 0) &&
    ((exclory & ~P4EST_QUADRANT_LEN (q1->level)) == 0) &&
#ifdef P4_To_P8
    ((exclorz & ~P4EST_QUADRANT_LEN (q1->level)) == 0) &&
#endif
    true;
}

bool
p4est_quadrant_is_sibling_D (const p4est_quadrant_t * q1,
                             const p4est_quadrant_t * q2)
{
  p4est_quadrant_t    p1, p2;

  /* make sure the quadrant_parent functions below don't abort */
  if (q1->level == 0) {
    return false;
  }

  /* validity of q1 and q2 is asserted in p4est_quadrant_is_equal */
  if (p4est_quadrant_is_equal (q1, q2)) {
    return false;
  }

  p4est_quadrant_parent (q1, &p1);
  p4est_quadrant_parent (q2, &p2);

  return p4est_quadrant_is_equal (&p1, &p2);
}

#ifndef P4_TO_P8

bool
p4est_quadrant_is_family (const p4est_quadrant_t * q0,
                          const p4est_quadrant_t * q1,
                          const p4est_quadrant_t * q2,
                          const p4est_quadrant_t * q3)
{
  const int8_t        level = q0->level;
  p4est_qcoord_t      inc;

  P4EST_ASSERT (p4est_quadrant_is_extended (q0));
  P4EST_ASSERT (p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_extended (q2));
  P4EST_ASSERT (p4est_quadrant_is_extended (q3));

  if (level == 0 || level != q1->level ||
      level != q2->level || level != q3->level) {
    return false;
  }

  inc = P4EST_QUADRANT_LEN (level);
  return ((q0->x + inc == q1->x && q0->y == q1->y) &&
          (q0->x == q2->x && q0->y + inc == q2->y) &&
          (q1->x == q3->x && q2->y == q3->y));
}

#endif /* !P4_TO_P8 */

bool
p4est_quadrant_is_parent (const p4est_quadrant_t * q,
                          const p4est_quadrant_t * r)
{
  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (r));

  return
    (q->level + 1 == r->level) &&
    (q->x == (r->x & ~P4EST_QUADRANT_LEN (r->level))) &&
    (q->y == (r->y & ~P4EST_QUADRANT_LEN (r->level))) &&
#ifdef P4_TO_P8
    (q->z == (r->z & ~P4EST_QUADRANT_LEN (r->level))) &&
#endif
    true;
}

bool
p4est_quadrant_is_parent_D (const p4est_quadrant_t * q,
                            const p4est_quadrant_t * r)
{
  p4est_quadrant_t    p;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));

  /* make sure the quadrant_parent function below doesn't abort */
  if (r->level == 0) {
    return false;
  }

  /* validity of r is asserted in p4est_quadrant_parent */
  p4est_quadrant_parent (r, &p);

  return p4est_quadrant_is_equal (q, &p);
}

bool
p4est_quadrant_is_ancestor (const p4est_quadrant_t * q,
                            const p4est_quadrant_t * r)
{
  p4est_qcoord_t      exclorx;
  p4est_qcoord_t      exclory;
#ifdef P4_TO_P8
  p4est_qcoord_t      exclorz;
#endif

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (r));

  if (q->level >= r->level) {
    return false;
  }

  exclorx = (q->x ^ r->x) >> (P4EST_MAXLEVEL - q->level);
  exclory = (q->y ^ r->y) >> (P4EST_MAXLEVEL - q->level);
#ifdef P4_TO_P8
  exclorz = (q->z ^ r->z) >> (P4EST_MAXLEVEL - q->level);

  return (exclorx == 0 && exclory == 0 && exclorz == 0);
#else
  return (exclorx == 0 && exclory == 0);
#endif
}

bool
p4est_quadrant_is_ancestor_D (const p4est_quadrant_t * q,
                              const p4est_quadrant_t * r)
{
  p4est_quadrant_t    s;

  /* validity of q and r is asserted in p4est_quadrant_is_equal */
  if (p4est_quadrant_is_equal (q, r)) {
    return false;
  }

  /* this will abort if q and r are in different trees */
  p4est_nearest_common_ancestor_D (q, r, &s);

  return p4est_quadrant_is_equal (q, &s);
}

bool
p4est_quadrant_is_next (const p4est_quadrant_t * q,
                        const p4est_quadrant_t * r)
{
  int                 minlevel;
  uint64_t            i1, i2;
  p4est_qcoord_t      mask;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (r));

  /* the condition q < r is checked implicitly */

  if (q->level > r->level) {
    /* check if q is the last child up to the common level */
    mask = P4EST_QUADRANT_LEN (r->level) - P4EST_QUADRANT_LEN (q->level);
    if ((q->x & mask) != mask || (q->y & mask) != mask ||
#ifdef P4_TO_P8
        (q->z & mask) != mask ||
#endif
        false) {
      return false;
    }
    minlevel = (int) r->level;
  }
  else {
    minlevel = (int) q->level;
  }
  i1 = p4est_quadrant_linear_id (q, minlevel);
  i2 = p4est_quadrant_linear_id (r, minlevel);

  return (i1 + 1 == i2);
}

bool
p4est_quadrant_is_next_D (const p4est_quadrant_t * q,
                          const p4est_quadrant_t * r)
{
  uint64_t            i1, i2;
  p4est_quadrant_t    a, b;

  /* validity of q and r is asserted in p4est_quadrant_compare */
  if (p4est_quadrant_compare (q, r) >= 0) {
    return false;
  }

  a = *q;
  b = *r;
  while (a.level > b.level) {
    if (p4est_quadrant_child_id (&a) != P4EST_CHILDREN - 1) {
      return false;
    }
    p4est_quadrant_parent (&a, &a);
  }
  i1 = p4est_quadrant_linear_id (&a, (int) a.level);
  i2 = p4est_quadrant_linear_id (&b, (int) a.level);

  return (i1 + 1 == i2);
}

bool
p4est_quadrant_overlaps_tree (p4est_tree_t * tree, const p4est_quadrant_t * q)
{
  int                 maxl;
  p4est_quadrant_t    desc;
  p4est_quadrant_t   *treeq;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (tree->quadrants.elem_count == 0)
    return false;

  maxl = SC_MAX (tree->maxlevel, q->level);

  /* check if q is before the first tree quadrant */
  treeq = sc_array_index (&tree->quadrants, 0);
  P4EST_ASSERT (p4est_quadrant_is_inside_root (treeq));
  p4est_quadrant_last_descendent (q, &desc, maxl);
  if (p4est_quadrant_compare (&desc, treeq) < 0)
    return false;

  /* check if q is after the last tree quadrant */
  treeq = sc_array_index (&tree->quadrants, tree->quadrants.elem_count - 1);
  P4EST_ASSERT (p4est_quadrant_is_inside_root (treeq));
  p4est_quadrant_last_descendent (treeq, &desc, maxl);
  if (p4est_quadrant_compare (&desc, q) < 0)
    return false;

  return true;
}

void
p4est_quadrant_parent (const p4est_quadrant_t * q, p4est_quadrant_t * r)
{
  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (q->level > 0);

  r->x = q->x & ~P4EST_QUADRANT_LEN (q->level);
  r->y = q->y & ~P4EST_QUADRANT_LEN (q->level);
#ifdef P4_TO_P8
  r->z = q->z & ~P4EST_QUADRANT_LEN (q->level);
#endif
  r->level = (int8_t) (q->level - 1);

  P4EST_ASSERT (p4est_quadrant_is_extended (r));
}

void
p4est_quadrant_sibling (const p4est_quadrant_t * q, p4est_quadrant_t * r,
                        int sibling_id)
{
  const int           addx = (sibling_id & 0x01);
  const int           addy = (sibling_id & 0x02) >> 1;
#ifdef P4_TO_P8
  const int           addz = (sibling_id & 0x04) >> 2;
#endif
  const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (q->level);

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (q->level > 0);
  P4EST_ASSERT (sibling_id >= 0 && sibling_id < P4EST_CHILDREN);

  r->x = (addx ? (q->x | shift) : (q->x & ~shift));
  r->y = (addy ? (q->y | shift) : (q->y & ~shift));
#ifdef P4_TO_P8
  r->z = (addz ? (q->z | shift) : (q->z & ~shift));
#endif
  r->level = q->level;
}

#ifndef P4_TO_P8

void
p4est_quadrant_children (const p4est_quadrant_t * q,
                         p4est_quadrant_t * c0, p4est_quadrant_t * c1,
                         p4est_quadrant_t * c2, p4est_quadrant_t * c3)
{
  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (q->level < P4EST_MAXLEVEL);

  c0->x = q->x;
  c0->y = q->y;
  c0->level = (int8_t) (q->level + 1);

  c1->x = c0->x | P4EST_QUADRANT_LEN (c0->level);
  c1->y = c0->y;
  c1->level = c0->level;

  c2->x = c0->x;
  c2->y = c0->y | P4EST_QUADRANT_LEN (c0->level);
  c2->level = c0->level;

  c3->x = c1->x;
  c3->y = c2->y;
  c3->level = c0->level;

  P4EST_ASSERT (p4est_quadrant_is_extended (c0));
  P4EST_ASSERT (p4est_quadrant_is_extended (c1));
  P4EST_ASSERT (p4est_quadrant_is_extended (c2));
  P4EST_ASSERT (p4est_quadrant_is_extended (c3));
  P4EST_ASSERT (p4est_quadrant_is_family (c0, c1, c2, c3));
}

#endif /* !P4_TO_P8 */

void
p4est_quadrant_first_descendent (const p4est_quadrant_t * q,
                                 p4est_quadrant_t * fd, int level)
{
  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT ((int) q->level <= level && level <= P4EST_MAXLEVEL);

  fd->x = q->x;
  fd->y = q->y;
#ifdef P4_TO_P8
  fd->z = q->z;
#endif
  fd->level = (int8_t) level;
}

void
p4est_quadrant_last_descendent (const p4est_quadrant_t * q,
                                p4est_quadrant_t * ld, int level)
{
  p4est_qcoord_t      shift;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT ((int) q->level <= level && level <= P4EST_MAXLEVEL);

  shift = P4EST_QUADRANT_LEN (q->level) - P4EST_QUADRANT_LEN (level);

  ld->x = q->x + shift;
  ld->y = q->y + shift;
#ifdef P4_TO_P8
  ld->z = q->z + shift;
#endif
  ld->level = (int8_t) level;
}

void
p4est_nearest_common_ancestor (const p4est_quadrant_t * q1,
                               const p4est_quadrant_t * q2,
                               p4est_quadrant_t * r)
{
  int                 maxlevel;
  uint32_t            exclorx, exclory;
#ifdef P4_TO_P8
  uint32_t            exclorz;
#endif
  uint32_t            maxclor;

  P4EST_ASSERT (p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_extended (q2));

  exclorx = q1->x ^ q2->x;
  exclory = q1->y ^ q2->y;
#ifdef P4_TO_P8
  exclorz = q1->z ^ q2->z;

  maxclor = exclorx | exclory | exclorz;
#else
  maxclor = exclorx | exclory;
#endif
  maxlevel = SC_LOG2_32 (maxclor) + 1;

  P4EST_ASSERT (maxlevel <= P4EST_MAXLEVEL);

  r->x = q1->x & ~((1 << maxlevel) - 1);
  r->y = q1->y & ~((1 << maxlevel) - 1);
#ifdef P4_TO_P8
  r->z = q1->z & ~((1 << maxlevel) - 1);
#endif
  r->level = (int8_t) SC_MIN (P4EST_MAXLEVEL - maxlevel,
                              (int) SC_MIN (q1->level, q2->level));

  P4EST_ASSERT (p4est_quadrant_is_extended (r));
}

void
p4est_nearest_common_ancestor_D (const p4est_quadrant_t * q1,
                                 const p4est_quadrant_t * q2,
                                 p4est_quadrant_t * r)
{
  p4est_quadrant_t    s1 = *q1;
  p4est_quadrant_t    s2 = *q2;

  P4EST_ASSERT (p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_extended (q2));

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
#ifdef P4_TO_P8
  r->z = s1.z;
#endif
  r->level = s1.level;

  P4EST_ASSERT (p4est_quadrant_is_extended (r));
}

#ifndef P4_TO_P8

int
p4est_quadrant_corner_level (const p4est_quadrant_t * q,
                             int zcorner, int level)
{
  p4est_qcoord_t      th, stepx, stepy;
  p4est_quadrant_t    quad, sibling;
  const p4est_qcoord_t zcorner_steps[4][2] =
    { {-1, -1,}, {1, -1}, {-1, 1}, {1, 1} };

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (zcorner >= 0 && zcorner < 4);
  P4EST_ASSERT (level >= 0 && level <= P4EST_MAXLEVEL);

  P4EST_QUADRANT_INIT (&quad);
  P4EST_QUADRANT_INIT (&sibling);
  stepx = zcorner_steps[zcorner][0];
  stepy = zcorner_steps[zcorner][1];

  quad = *q;
  while ((int) quad.level > level) {
    th = P4EST_LAST_OFFSET (quad.level);
    p4est_quadrant_sibling (&quad, &sibling, zcorner);
    if ((zcorner == 0 && sibling.x <= 0 && sibling.y <= 0) ||
        (zcorner == 1 && sibling.x >= th && sibling.y <= 0) ||
        (zcorner == 2 && sibling.x <= 0 && sibling.y >= th) ||
        (zcorner == 3 && sibling.x >= th && sibling.y >= th)) {
      return (int) quad.level;
    }
    p4est_quadrant_parent (&quad, &quad);
    quad.x += stepx * P4EST_QUADRANT_LEN (quad.level);
    quad.y += stepy * P4EST_QUADRANT_LEN (quad.level);
    P4EST_ASSERT (p4est_quadrant_is_extended (&quad));
  }
  P4EST_ASSERT ((int) quad.level == level);

  return level;
}

void
p4est_quadrant_corner (p4est_quadrant_t * q, int zcorner, int inside)
{
  p4est_qcoord_t      lshift, rshift;

  P4EST_ASSERT (q->level >= 0 && q->level <= P4EST_MAXLEVEL);

  lshift = (inside ? 0 : -P4EST_QUADRANT_LEN (q->level));
  rshift = (inside ? P4EST_LAST_OFFSET (q->level) : P4EST_ROOT_LEN);

  switch (zcorner) {
  case 0:
    q->x = q->y = lshift;
    break;
  case 1:
    q->x = rshift;
    q->y = lshift;
    break;
  case 2:
    q->x = lshift;
    q->y = rshift;
    break;
  case 3:
    q->x = q->y = rshift;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }

  P4EST_ASSERT ((inside && p4est_quadrant_is_valid (q)) ||
                (!inside && p4est_quadrant_is_extended (q)));
}

void
p4est_quadrant_translate (p4est_quadrant_t * q, int face)
{
  P4EST_ASSERT (p4est_quadrant_is_extended (q));

  switch (face) {
  case 0:
    q->y += P4EST_ROOT_LEN;
    break;
  case 1:
    q->x -= P4EST_ROOT_LEN;
    break;
  case 2:
    q->y -= P4EST_ROOT_LEN;
    break;
  case 3:
    q->x += P4EST_ROOT_LEN;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
}

int
p4est_node_transform (int node, int transform_type)
{
  int                 trans_node;

  P4EST_ASSERT (0 <= node && node < 4);

  switch (transform_type) {
  case 0:                      /* identity */
    trans_node = node;
    break;
  case 1:                      /* rotate -90 degrees */
    trans_node =
      p4est_corner_to_zorder[(p4est_corner_to_zorder[node] + 1) % 4];
    break;
  case 2:                      /* rotate 180 degrees */
    trans_node = 3 - node;
    break;
  case 3:                      /* rotate 90 degrees */
    trans_node =
      p4est_corner_to_zorder[(p4est_corner_to_zorder[node] + 3) % 4];
    break;
  case 4:                      /* mirror across 0 degree axis */
    switch (node) {
    case 0:
      trans_node = 2;
      break;
    case 1:
      trans_node = 3;
      break;
    case 2:
      trans_node = 0;
      break;
    case 3:
      trans_node = 1;
      break;
    default:
      SC_CHECK_NOT_REACHED ();
      break;
    }
    break;
  case 5:                      /* mirror across 45 degree axis */
    switch (node) {
    case 0:
      trans_node = 0;
      break;
    case 1:
      trans_node = 2;
      break;
    case 2:
      trans_node = 1;
      break;
    case 3:
      trans_node = 3;
      break;
    default:
      SC_CHECK_NOT_REACHED ();
      break;
    }
    break;
  case 6:                      /* mirror across 90 degree axis */
    switch (node) {
    case 0:
      trans_node = 1;
      break;
    case 1:
      trans_node = 0;
      break;
    case 2:
      trans_node = 3;
      break;
    case 3:
      trans_node = 2;
      break;
    default:
      SC_CHECK_NOT_REACHED ();
      break;
    }
    break;
  case 7:                      /* mirror across 135 degree axis */
    switch (node) {
    case 0:
      trans_node = 3;
      break;
    case 1:
      trans_node = 1;
      break;
    case 2:
      trans_node = 2;
      break;
    case 3:
      trans_node = 0;
      break;
    default:
      SC_CHECK_NOT_REACHED ();
      break;
    }
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }
  return trans_node;
}

void
p4est_quadrant_transform (const p4est_quadrant_t * q,
                          p4est_quadrant_t * r, int transform_type)
{
  p4est_qcoord_t      th;

  P4EST_ASSERT (q != r);
  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (0 <= transform_type && transform_type < 8);

  th = P4EST_LAST_OFFSET (q->level);

  switch (transform_type) {
  case 0:                      /* identity */
    r->x = q->x;
    r->y = q->y;
    break;
  case 1:                      /* rotate -90 degrees */
    r->x = th - q->y;
    r->y = q->x;
    break;
  case 2:                      /* rotate 180 degrees */
    r->x = th - q->x;
    r->y = th - q->y;
    break;
  case 3:                      /* rotate 90 degrees */
    r->x = q->y;
    r->y = th - q->x;
    break;
  case 4:                      /* mirror across 0 degree axis */
    r->x = q->x;
    r->y = th - q->y;
    break;
  case 5:                      /* mirror across 45 degree axis */
    r->x = q->y;
    r->y = q->x;
    break;
  case 6:                      /* mirror across 90 degree axis */
    r->x = th - q->x;
    r->y = q->y;
    break;
  case 7:                      /* mirror across 135 degree axis */
    r->x = th - q->y;
    r->y = th - q->x;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
    break;
  }
  r->level = q->level;

  P4EST_ASSERT (p4est_quadrant_is_extended (r));
}

#endif /* !P4_TO_P8 */

uint64_t
p4est_quadrant_linear_id (const p4est_quadrant_t * quadrant, int level)
{
  int                 i;
  uint64_t            id;
  uint64_t            x, y;
#ifdef P4_TO_P8
  uint64_t            z;
#endif

  P4EST_ASSERT (p4est_quadrant_is_extended (quadrant));
  P4EST_ASSERT ((int) quadrant->level >= level && level >= 0);

  /* this preserves the high bits from negative numbers */
  x = quadrant->x >> (P4EST_MAXLEVEL - level);
  y = quadrant->y >> (P4EST_MAXLEVEL - level);
#ifdef P4_TO_P8
  z = quadrant->z >> (P4EST_MAXLEVEL - level);
#endif

  id = 0;
  for (i = 0; i < level + 2; ++i) {
    id |= ((x & ((uint64_t) 1 << i)) << ((P4EST_DIM - 1) * i));
    id |= ((y & ((uint64_t) 1 << i)) << ((P4EST_DIM - 1) * i + 1));
#ifdef P4_TO_P8
    id |= ((z & ((uint64_t) 1 << i)) << ((P4EST_DIM - 1) * i + 2));
#endif
  }

  return id;
}

void
p4est_quadrant_set_morton (p4est_quadrant_t * quadrant,
                           int level, uint64_t id)
{
  int                 i;

  P4EST_ASSERT (0 <= level && level <= P4EST_MAXLEVEL);
  if (level < P4EST_MAXLEVEL) {
    P4EST_ASSERT (id < ((uint64_t) 1 << P4EST_DIM * (level + 2)));
  }

  quadrant->level = (int8_t) level;
  quadrant->x = 0;
  quadrant->y = 0;
#ifdef P4_TO_P8
  quadrant->z = 0;
#endif

  /* this may set the sign bit to create negative numbers */
  for (i = 0; i < level + 2; ++i) {
    quadrant->x |= (p4est_qcoord_t) ((id & (1ULL << (P4EST_DIM * i)))
                                     >> ((P4EST_DIM - 1) * i));
    quadrant->y |= (p4est_qcoord_t) ((id & (1ULL << (P4EST_DIM * i + 1)))
                                     >> ((P4EST_DIM - 1) * i + 1));
#ifdef P4_TO_P8
    quadrant->z |= (p4est_qcoord_t) ((id & (1ULL << (P4EST_DIM * i + 2)))
                                     >> ((P4EST_DIM - 1) * i + 2));
#endif
  }

  quadrant->x <<= (P4EST_MAXLEVEL - level);
  quadrant->y <<= (P4EST_MAXLEVEL - level);
#ifdef P4_TO_P8
  quadrant->z <<= (P4EST_MAXLEVEL - level);
#endif

  P4EST_ASSERT (p4est_quadrant_is_extended (quadrant));
}

/* EOF p4est_bits.c */
