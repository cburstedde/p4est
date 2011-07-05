/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2011 The University of Texas System
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

#ifdef P4_TO_P8
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_search.h>
#else
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_search.h>
#endif /* !P4_TO_P8 */

/* We have a location, and a \a level quadrant that must be shifted by
 * \a distance (>= 0) to be at the location.  This returns the largest quadrant
 * that can exist at \a location and be in balance with the \a level quadrant.
 */
static inline       int8_t
p4est_balance_kernel_1d (p4est_qcoord_t distance, int8_t level)
{
  int                 shift = P4EST_MAXLEVEL - (int) level;
  P4EST_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);
  /* the distance only makes sense if it is an integer number of \a level
   * distances */
  P4EST_ASSERT (distance >= 0);
  P4EST_ASSERT (!(distance & (~(((p4est_qcoord_t) - 1) << shift))));
  distance >>= shift;
  /* The theory says we should use ((distance + 1)&(~1) + 1), but
   * using distance + 1 is equivalent for all distance >= 0 */
  distance++;

  return SC_MAX (0, level - (int8_t) SC_LOG2_32 (distance));
}

/* This is the kernel for 2D balance with face-only balancing */
static inline       int8_t
p4est_balance_kernel_2d (p4est_qcoord_t dx, p4est_qcoord_t dy, int8_t level)
{
  int                 shift = P4EST_MAXLEVEL - (int) level;
  p4est_qcoord_t      distance;

  P4EST_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);
  /* the distance only makes sense if it is an integer number of \a level
   * distances */
  P4EST_ASSERT (dx >= 0);
  P4EST_ASSERT (!(dx & (~(((p4est_qcoord_t) - 1) << shift))));
  P4EST_ASSERT (dy >= 0);
  P4EST_ASSERT (!(dy & (~(((p4est_qcoord_t) - 1) << shift))));

  dx >>= shift;
  /* get the smallest even number greater than or equal to dx */
  dx = (dx + 1) & (~((p4est_qcoord_t) 0x1));

  dy >>= shift;
  /* get the smallest even number greater than or equal to dy */
  dy = (dy + 1) & (~((p4est_qcoord_t) 0x1));

  /* the + 1 guarantees the correct answer even for (0, 0) */
  distance = dx + dy + 1;

  return SC_MAX (0, level - (int8_t) SC_LOG2_32 (distance));
}

#ifdef P4_TO_P8
/* This is the kernel for 3d balance with face and edge balancing */
static inline       int8_t
p8est_balance_kernel_3d_edge (p4est_qcoord_t dx, p4est_qcoord_t dy,
                              p4est_qcoord_t dz, int8_t level)
{
  int                 shift = P4EST_MAXLEVEL - (int) level;
  int                 xbit, ybit, zbit;
  int                 maxbit, thisbit;
  int                 count;

  P4EST_ASSERT (dx >= 0);
  P4EST_ASSERT (!(dx & (~(((p4est_qcoord_t) - 1) << shift))));
  P4EST_ASSERT (dy >= 0);
  P4EST_ASSERT (!(dy & (~(((p4est_qcoord_t) - 1) << shift))));
  P4EST_ASSERT (dz >= 0);
  P4EST_ASSERT (!(dz & (~(((p4est_qcoord_t) - 1) << shift))));

  if (!dx && !dy && !dz) {
    return level;
  }

  dx >>= shift;
  /* get the smallest even number greater than or equal to dx */
  dx = (dx + 1) & (~((p4est_qcoord_t) 0x1));

  dy >>= shift;
  /* get the smallest even number greater than or equal to dy */
  dy = (dy + 1) & (~((p4est_qcoord_t) 0x1));

  dz >>= shift;
  /* get the smallest even number greater than or equal to dz */
  dz = (dz + 1) & (~((p4est_qcoord_t) 0x1));

  xbit = SC_LOG2_32 (dx);
  maxbit = xbit;
  ybit = SC_LOG2_32 (dy);
  maxbit = SC_MAX (maxbit, ybit);
  zbit = SC_LOG2_32 (dz);
  maxbit = SC_MAX (maxbit, zbit);

  P4EST_ASSERT (maxbit >= 1);

  count = (xbit == maxbit);
  count += (ybit == maxbit);
  count += (zbit == maxbit);

  switch (count) {
  case 1:
    /* There is always a path where at most one of the other two dimensions
     * adds a bit in this position, so there is always a path where we don't
     * create a more significant bit */
    return SC_MAX (0, level - (int8_t) maxbit);
  case 2:
    /* This is the start of a chain of 2s.  If this chain ends in 0 or 1,
     * there is always a path where we don't add a more significant bit; if
     * this chain ends in a 3, there is always a path where only one dimension
     * adds a more significant bit */
    thisbit = maxbit - 1;
    do {
      P4EST_ASSERT (thisbit >= 0);
      count = ((dx & (1 << thisbit)) != 0);
      count += ((dy & (1 << thisbit)) != 0);
      count += ((dz & (1 << thisbit)) != 0);
      switch (count) {
      case 0:
      case 1:
        return SC_MAX (0, level - (int8_t) maxbit);
      case 2:
        break;
      case 3:
        return SC_MAX (0, level - (int8_t) (maxbit + 1));
      default:
        SC_ABORT_NOT_REACHED ();
      }
      P4EST_ASSERT (count == 2);
      thisbit--;
    } while (count == 2);
    SC_ABORT_NOT_REACHED ();
    return -1;
  case 3:
    /* There is always a path where only one dimension adds a more signifcant
     * bit */
    return SC_MAX (0, level - (int8_t) (maxbit + 1));
  default:
    SC_ABORT_NOT_REACHED ();
    return -1;
  }
}

static inline int   p8est_balance_2chain (p4est_qcoord_t dx,
                                          p4est_qcoord_t dy,
                                          p4est_qcoord_t dz, int thisbit,
                                          int xbit, int ybit, int zbit);

/* This is the kernel for 3d balance with face balancing only */
static inline       int8_t
p8est_balance_kernel_3d_face (p4est_qcoord_t dx, p4est_qcoord_t dy,
                              p4est_qcoord_t dz, int8_t level)
{
  int                 shift = P4EST_MAXLEVEL - (int) level;
  int                 xbit, ybit, zbit;
  int                 maxbit, thisbit;
  int                 count;

  P4EST_ASSERT (dx >= 0);
  P4EST_ASSERT (!(dx & (~(((p4est_qcoord_t) - 1) << shift))));
  P4EST_ASSERT (dy >= 0);
  P4EST_ASSERT (!(dy & (~(((p4est_qcoord_t) - 1) << shift))));
  P4EST_ASSERT (dz >= 0);
  P4EST_ASSERT (!(dz & (~(((p4est_qcoord_t) - 1) << shift))));

  if (!dx && !dy && !dz) {
    return level;
  }

  dx >>= shift;
  /* get the smallest even number greater than or equal to dx */
  dx = (dx + 1) & (~((p4est_qcoord_t) 0x1));

  dy >>= shift;
  /* get the smallest even number greater than or equal to dy */
  dy = (dy + 1) & (~((p4est_qcoord_t) 0x1));

  dz >>= shift;
  /* get the smallest even number greater than or equal to dz */
  dz = (dz + 1) & (~((p4est_qcoord_t) 0x1));

  xbit = SC_LOG2_32 (dx);
  maxbit = xbit;
  ybit = SC_LOG2_32 (dy);
  maxbit = SC_MAX (maxbit, ybit);
  zbit = SC_LOG2_32 (dz);
  maxbit = SC_MAX (maxbit, zbit);

  P4EST_ASSERT (maxbit >= 1);

  count = (xbit == maxbit);
  count += (ybit == maxbit);
  count += (zbit == maxbit);

  switch (count) {
  case 1:
    /* This is the start of a chain of 1s.  If this chain ends in a 2, 3, or
     * 03, * then every path results in adding one more significant bit.  If
     * this chain ends in 00 or 01, then there is a path that results in no
     * more significant bit.  If this chain ends in 02, we need more
     * information */
    thisbit = maxbit - 1;
    do {
      P4EST_ASSERT (thisbit >= 0);
      count = ((dx & (1 << thisbit)) != 0);
      count += ((dy & (1 << thisbit)) != 0);
      count += ((dz & (1 << thisbit)) != 0);
      switch (count) {
      case 0:
        if (thisbit-- == 0) {
          return SC_MAX (0, level - (int8_t) (maxbit + 1));
        }
        xbit = (dx & (1 << thisbit));
        ybit = (dy & (1 << thisbit));
        zbit = (dz & (1 << thisbit));
        count = (xbit != 0);
        count += (ybit != 0);
        count += (zbit != 0);
        switch (count) {
        case 0:
        case 1:
          return SC_MAX (0, level - (int8_t) maxbit);
        case 2:
          /* If there is a path where this 2 advances only 1 bit, then the
           * initial one does not advance a bit; otherwise this 2 advances 2
           * bits, and the initial 1 advances a bit */
          if (p8est_balance_2chain (dx, dy, dz, thisbit, xbit, ybit, zbit)) {
            return SC_MAX (0, level - (int8_t) (maxbit + 1));
          }
          else {
            return SC_MAX (0, level - (int8_t) maxbit);
          }
        case 3:
          return SC_MAX (0, level - (int8_t) (maxbit + 1));
        default:
          SC_ABORT_NOT_REACHED ();
          return -1;
        }
      case 1:
        break;
      case 2:
      case 3:
        return SC_MAX (0, level - (int8_t) (maxbit + 1));
      default:
        SC_ABORT_NOT_REACHED ();
      }
      P4EST_ASSERT (count == 1);
      thisbit--;
    } while (count == 1);
    SC_ABORT_NOT_REACHED ();
    return -1;
  case 2:
    /* Check to see if there is a path where this 2 advances only 1 bit:
     * otherwise it advances 2 bits */
    if (p8est_balance_2chain (dx, dy, dz, maxbit, xbit, ybit, zbit)) {
      return SC_MAX (0, level - (int8_t) (maxbit + 2));
    }
    else {
      return SC_MAX (0, level - (int8_t) (maxbit + 1));
    }
  case 3:
    /* Every path results in a bit twice more significant */
    return SC_MAX (0, level - (int8_t) (maxbit + 2));
  default:
    SC_ABORT_NOT_REACHED ();
    return -1;
  }
}

static inline int
p8est_balance_2chain (p4est_qcoord_t dx, p4est_qcoord_t dy, p4est_qcoord_t dz,
                      int thisbit, int xbit, int ybit, int zbit)
{
  int                 zfx = xbit;
  int                 zfy = ybit;
  int                 zfz = zbit;
  int                 count;

  thisbit--;

  P4EST_ASSERT ((zfx != 0) + (zfy != 0) + (zfz != 0) == 2);

  do {
    P4EST_ASSERT (thisbit >= 0);

    xbit = (dx & (1 << thisbit));
    ybit = (dy & (1 << thisbit));
    zbit = (dz & (1 << thisbit));
    count = (xbit != 0);
    count += (ybit != 0);
    count += (zbit != 0);

    switch (count) {
    case 0:
      /* if the chain of 2's with at least one zero-free dimension ends in a
       * 0, then there is always a path where the most significant 2 advances
       * only one bit */
      return 0;
    case 1:
      /* if the chain of 2's with at least one zero-free dimension ends in a
       * 1, and that one is in a zero-free dimension, then there is always a
       * path where the most significant 2 advances only one bit; if the one
       * is not in a zero-free dimension, we can treat it as though it is a
       * two that maintains the zero-free dimensions and continue */
      if ((xbit != 0) && zfx) {
        return 0;
      }
      if ((ybit != 0) && zfy) {
        return 0;
      }
      if ((zbit != 0) && zfy) {
        return 0;
      }
      break;
    case 2:
      /* if we no longer have a zero free dimension, then the most significant
       * 2 advances 2 bits; otherwise we continue */
      zfx = (zfx && xbit);
      zfy = (zfy && ybit);
      zfz = (zfz && zbit);
      if (!(zfx || zfy || zfz)) {
        return 1;
      }
      break;
    case 3:
      /* the most significant 2 advances 2 bits */
      return 1;
    default:
      SC_ABORT_NOT_REACHED ();
    }

    thisbit--;
  } while (thisbit >= 0);

  SC_ABORT_NOT_REACHED ();
  return -1;
}
#endif

/* a quadrant q is outside of quadrant p, as a descendent of p's neighbor n
 * outside of \a corner.  We know the \a child of n that contains q.
 * Depending on the \a balance type, q is either relevant (return true) or
 * irrelevant (return false) to the balance of p. */
#ifndef P4_TO_P8
static inline int
p4est_balance_corner_relevant (int balance, int corner, int child)
{
  /* if balance == 1 (i.e. face and corner balance), then every descendent of
   * n is relevant; if balance == 0 (face balance), the child of n farthest
   * from p's corner (i.e. that has the same number as p's corner) is the only
   * child that is not relevant) */

  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= child && child < P4EST_CHILDREN);
  return (balance || corner != child);
}
#else

#define p4est_balance_corner_relevant p8est_balance_corner_relevant
/* assuming face only balance \a corner == 0, only the closest child of n
 * (child 7) and its face neighbors (3, 5, 6) are relevant to p's balance */
static const int    p8est_balance_corner_face_relevant[8] =
  { 0, 0, 0, 1, 0, 1, 1, 1 };

static inline int
p8est_balance_corner_relevant (int balance, int corner, int child)
{
  /* if balance == 2 (i.e. face, edge and corner balance), then every
   * descendent of n is relevant; if balance == 1 (face and edge balance),
   * the child of n farthest from p's corner (i.e. that has the same number as
   * p's corner) is the only child that is not relevant); if balance == 0
   * (face only balance), we transform \a child as though \a corner == 0 and
   * use the table above */

  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= child && child < P4EST_CHILDREN);
  return (balance == 2 ||
          (balance == 1 && child != corner) ||
          p8est_balance_corner_face_relevant[child ^ corner]);
}

/* a quadrant q is outside of quadrant p, as a descendent of p's neighbor n
 * outside of \a edge.  We know the \a child of n that contains q.
 * Depending on the \a balance type, q is either relevant (return true) or
 * irrelevant (return false) to the balance of p. */

static inline int
p8est_balance_edge_relevant (int balance, int edge, int child)
{
  /* if there is edge balancing (balance > 0), then every descendent of n is
   * relevant; if there is only face balancing, \a child is not relevant if it
   * touches the farthest edge of n away from p */

  P4EST_ASSERT (0 <= edge && edge < P8EST_EDGES);
  P4EST_ASSERT (0 <= child && child < P4EST_CHILDREN);
  return (balance || p8est_corner_edges[child][edge >> 2] != edge);
}
#endif

static void         p4est_bal_corner_con_internal (p4est_quadrant_t const
                                                   *restrict q,
                                                   p4est_quadrant_t *
                                                   restrict p, int corner,
                                                   int balance,
                                                   int *consistent,
                                                   int recurse);

static void         p4est_bal_face_con_internal (p4est_quadrant_t const
                                                 *restrict q,
                                                 p4est_quadrant_t *
                                                 restrict p, int face,
                                                 int balance, int *consistent,
                                                 int *conextra, int recurse);

#ifdef P4_TO_P8
static void         p8est_bal_edge_con_internal (p4est_quadrant_t const
                                                 *restrict q,
                                                 p4est_quadrant_t *
                                                 restrict p, int edge,
                                                 int balance, int *consistent,
                                                 int *conextra, int recurse);
#endif

/* the recursive call of p4est_bal_corner_con_internal: we have determined
 * that q causes p to split: this means that they are not consistent, so if \a
 * consistent is given, it is set to false.  If we are to \a recurse, we set p
 * to be its own child in \a corner and recur */
#ifndef P4_TO_P8
#define p4est_bal_corner_recurse(q, p, c, b, cn, r)                    \
do {                                                                   \
  if ((cn) != NULL) {                                                  \
    *(cn) = 0;                                                         \
  }                                                                    \
  if (r) {                                                             \
    (p)->level++;                                                      \
    (p)->x += (((c) & 1) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    (p)->y += (((c) & 2) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    P4EST_ASSERT (p4est_quadrant_is_extended (p));                     \
    p4est_bal_corner_con_internal (q, p, c, b, NULL, r);               \
  }                                                                    \
} while (0)
#else
#define p4est_bal_corner_recurse(q, p, c, b, cn, r)                    \
do {                                                                   \
  if ((cn) != NULL) {                                                  \
    *(cn) = 0;                                                         \
  }                                                                    \
  if (r) {                                                             \
    (p)->level++;                                                      \
    (p)->x += (((c) & 1) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    (p)->y += (((c) & 2) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    (p)->z += (((c) & 4) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    P4EST_ASSERT (p4est_quadrant_is_extended (p));                     \
    p4est_bal_corner_con_internal (q, p, c, b, NULL, r);               \
  }                                                                    \
} while (0)
#endif

#ifndef P4_TO_P8
#define p4est_bal_corner_no_recurse(q, p, c, b, cn, r)                 \
do {                                                                   \
  if ((cn) != NULL) {                                                  \
    *(cn) = 0;                                                         \
  }                                                                    \
  if (r) {                                                             \
    (p)->level++;                                                      \
    (p)->x += (((c) & 1) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    (p)->y += (((c) & 2) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    P4EST_ASSERT (p4est_quadrant_is_extended (p));                     \
  }                                                                    \
} while (0)
#else
#define p4est_bal_corner_no_recurse(q, p, c, b, cn, r)                 \
do {                                                                   \
  if ((cn) != NULL) {                                                  \
    *(cn) = 0;                                                         \
  }                                                                    \
  if (r) {                                                             \
    (p)->level++;                                                      \
    (p)->x += (((c) & 1) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    (p)->y += (((c) & 2) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    (p)->z += (((c) & 4) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    P4EST_ASSERT (p4est_quadrant_is_extended (p));                     \
  }                                                                    \
} while (0)
#endif

/* \a corner is the corner of \a p closest to \a q: the corner of \a q closest
 * to \a p is the dual of \a corner */
static void
p4est_bal_corner_con_internal (p4est_quadrant_t const *restrict q,
                               p4est_quadrant_t * restrict p,
                               int corner, int balance,
                               int *consistent, int recurse)
{
  int                 level = p->level + 1;     /* level is now the size of a
                                                   child */
  int                 child;
  int                 dual = corner ^ (P4EST_CHILDREN - 1);
  int                 relative;
  int                 face;
#ifdef P4_TO_P8
  int                 edge;
  int                 grandchild;
#endif
  int                 recon;
  p4est_quadrant_t    temp;
  p4est_quadrant_t    a;

#ifdef P4EST_DEBUG
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (p));
  P4EST_ASSERT (q->level >= p->level);
  /* get the ancestor of q at the same level as p */
  a = *q;
  a.x &= (-1 << (P4EST_MAXLEVEL - p->level));
  a.y &= (-1 << (P4EST_MAXLEVEL - p->level));
#ifdef P4_TO_P8
  a.z &= (-1 << (P4EST_MAXLEVEL - p->level));
#endif
  a.level = p->level;
  /* check that q really is a descendent of the correct neighbor of p */
  p4est_quadrant_corner_neighbor (&a, dual, &temp);
  P4EST_ASSERT (p4est_quadrant_is_equal (&temp, p));
#endif

  /* start with the assumption that the two are consistent */
  if (consistent != NULL) {
    *consistent = 1;
  }

  /* if q is larger than or equal to half of p,
   * or p is a smallest quadrant, then of
   * course q and p are balance consistent */
  if (q->level <= level || level > P4EST_QMAXLEVEL) {
    return;
  }

  /* early exit if q is in an irrelevant descendent of n */
  child = p4est_quadrant_ancestor_id (q, level);
  if (!p4est_balance_corner_relevant (balance, corner, child)) {
    return;
  }

  /* if q is in the child of n closest to p, i.e. directly outside \a corner
   * */
  if (child == dual) {
    if (balance == (P4EST_DIM - 1)) {
      /* q is smaller than a child of n, which means that a quadrant smaller
       * than a child of p is directly outside of \a corner, which means
       * that q and p are not balance compatible */
      p4est_bal_corner_recurse (q, p, corner, balance, consistent, recurse);
      return;
    }
    if (++level == q->level) {  /* level is now the size of a grandchild */
      /* in all of the remaining cases, a grandchild is consistent */
      return;
    }
#ifdef P4_TO_P8
    if (balance == 1) {         /* edge balance */
#endif
      do {
        child = p4est_quadrant_ancestor_id (q, level);
        if (child != corner) {
          p4est_bal_corner_recurse (q, p, corner, balance, consistent,
                                    recurse);
          return;
        }
      } while (++level < q->level);
      return;
#ifdef P4_TO_P8
    }

    P4EST_ASSERT (balance == 0);
    P4EST_ASSERT (q->level > level);
    /* this is the trickiest case: outside a corner, 3D, face only balance */
    child = p4est_quadrant_ancestor_id (q, level);
    relative = child ^ dual;
    switch (relative) {
    case 7:
      /* if we are in the (grand)child farthest from p, we are consistent */
      return;
    case 0:
      /* if we are in the (grand)child closest to p, it looks like the edge
       * balance case */
      while (++level < q->level) {
        child = p4est_quadrant_ancestor_id (q, level);
        if (child != corner) {
          p4est_bal_corner_recurse (q, p, corner, balance, consistent,
                                    recurse);
          return;
        }
      }
      return;
    case 1:
    case 2:
    case 4:
      /* we are in a (grand)child that is face adjacent to the closest
       * (grand)child to p, which looks like trying to balance
       * outside an edge */
      edge = p8est_corner_edges[corner][relative >> 1];
      while (++level < q->level) {
        child = p4est_quadrant_ancestor_id (q, level);
        if (p8est_balance_edge_relevant (balance, edge, child)) {
          p4est_bal_corner_recurse (q, p, corner, balance, consistent,
                                    recurse);
          return;
        }
      }
      return;
    case 3:
    case 5:
    case 6:
      /* we are in a (grand)child that is edge adjacent to the closest
       * (grand)child to p.  this is the bizarre fractal case */
      while (++level < q->level) {
        grandchild = p4est_quadrant_ancestor_id (q, level);
        if (grandchild == child || grandchild == (child ^ 7)) {
          /* this is the self-similar case */
          continue;
        }
        relative = grandchild ^ dual;
        switch (relative) {
        case 7:
          /* consistent */
          return;
        case 0:
          while (++level < q->level) {
            /* this looks like the outside corner, balance edge case */
            grandchild = p4est_quadrant_ancestor_id (q, level);
            if (grandchild != corner) {
              p4est_bal_corner_recurse (q, p, corner, balance, consistent,
                                        recurse);
              return;
            }
          }
          return;
        case 1:
        case 2:
        case 4:
          /* this case looks like the outside edge, balance face case */
          edge = p8est_corner_edges[corner][relative >> 1];
          while (++level < q->level) {
            grandchild = p4est_quadrant_ancestor_id (q, level);
            if (p8est_balance_edge_relevant (balance, edge, grandchild)) {
              p4est_bal_corner_recurse (q, p, corner, balance, consistent,
                                        recurse);
              return;
            }
          }
          return;
        case 3:
        case 5:
        case 6:
          /* this case looks sort of like the outside edge,
           * balance face case */
          edge = p8est_corner_edges[corner][((child ^ grandchild) ^ 7) >> 1];
          while (++level < q->level) {
            grandchild = p4est_quadrant_ancestor_id (q, level);
            if (!p8est_balance_edge_relevant (balance, edge, grandchild)) {
              /* consistent */
              return;
            }
            if (!p8est_balance_edge_relevant (balance, edge ^ 3, grandchild)) {
              while (++level < q->level) {
                grandchild = p4est_quadrant_ancestor_id (q, level);
                if (p8est_balance_edge_relevant (balance, edge, grandchild)) {
                  p4est_bal_corner_recurse (q, p, corner, balance, consistent,
                                            recurse);
                  return;
                }
              }
              return;
            }
          }
          return;
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }
      return;
    default:
      SC_ABORT_NOT_REACHED ();
    }
#endif
  }

  /* set a to be the child of p's neighbor that contains q */
  a.x = (q->x & (-1 << (P4EST_MAXLEVEL - level)));
  a.y = (q->y & (-1 << (P4EST_MAXLEVEL - level)));
#ifdef P4_TO_P8
  a.z = (q->z & (-1 << (P4EST_MAXLEVEL - level)));
#endif
  a.level = level;
  P4EST_ASSERT (p4est_quadrant_is_ancestor (&a, q));

  if (balance == (P4EST_DIM - 1)) {

    relative = child ^ dual;
    switch (relative) {
    case (P4EST_CHILDREN - 1):
      p4est_quadrant_corner_neighbor (&a, dual, &temp);
      p4est_bal_corner_con_internal (q, &temp, corner, balance, &recon, 0);
      if (!recon) {
        p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                     recurse);
        return;
      }
      return;
    case 1:
    case 2:
#ifdef P4_TO_P8
    case 4:
#endif
      face = p4est_corner_faces[corner][relative >> 1];
      p4est_quadrant_face_neighbor (&a, face ^ 1, &temp);
      p4est_bal_face_con_internal (q, &temp, face, balance, &recon, NULL, 0);
      if (!recon) {
        p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                     recurse);
        return;
      }
      return;
#ifdef P4_TO_P8
    case 3:
    case 5:
    case 6:
      edge = p8est_corner_edges[corner][(relative ^ 7) >> 1];
      p8est_quadrant_edge_neighbor (&a, edge ^ 3, &temp);
      p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon, NULL, 0);
      if (!recon) {
        p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                     recurse);
        return;
      }
      return;
#endif
    default:
      SC_ABORT_NOT_REACHED ();
    }
  }

#ifndef P4_TO_P8
  P4EST_ASSERT (!balance);
#else
  if (!balance) {
#endif
    p4est_quadrant_corner_neighbor (&a, dual, &temp);
    p4est_bal_corner_con_internal (q, &temp, corner, balance, &recon, 0);
    if (!recon) {
      p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                   recurse);
      return;
    }
    return;
#ifdef P4_TO_P8
  }

  P4EST_ASSERT (balance == 1);

  relative = child ^ dual;
  switch (relative) {
  case 1:
  case 2:
  case 4:
    relative >>= 1;
    edge = p8est_corner_edges[corner][(relative + 1) % 3];
    p8est_quadrant_edge_neighbor (&a, edge ^ 3, &temp);
    p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon, NULL, 0);
    if (!recon) {
      p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                   recurse);
      return;
    }
    edge = p8est_corner_edges[corner][(relative + 2) % 3];
    p8est_quadrant_edge_neighbor (&a, edge ^ 3, &temp);
    p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon, NULL, 0);
    if (!recon) {
      p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                   recurse);
      return;
    }
    return;
  case 3:
  case 5:
  case 6:
    p4est_quadrant_corner_neighbor (&a, dual, &temp);
    p4est_bal_corner_con_internal (q, &temp, corner, balance, &recon, 0);
    if (!recon) {
      p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                   recurse);
      return;
    }
    return;
  default:
    SC_ABORT_NOT_REACHED ();
  }
#endif
}

static void
p4est_bal_face_con_internal (p4est_quadrant_t const *restrict q,
                             p4est_quadrant_t * restrict p, int face,
                             int balance, int *consistent, int *conextra,
                             int recurse)
{
  int                 level = p->level + 1;     /* level is now the size of a
                                                   child */
  int                 child;
  int                 dual = face ^ 1;
#ifdef P4_TO_P8
  int                 edge;
#endif
  int                 recon;
  p4est_quadrant_t    temp;
  p4est_quadrant_t    a;
  int                 cfc;
  int                 cfcextra;
  int                 i;
#ifndef P4_TO_P8
  int                 nconextra = 3;
#else
  int                 nconextra = 9;
  int                 reconextra[3];
  int                 j;
#endif
  int                 extraid;
  double              distance;
  int                 blevel;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (p));

  if (q->level <= p->level) {
    if (consistent != NULL) {
      *consistent = 1;
    }
    if (conextra != NULL) {
      conextra[nconextra / 2] = p->level;
    }
    return;
  }

  switch (face) {
  case 0:
    distance = p->x - q->x;
    break;
  case 1:
    distance = (q->x + P4EST_QUADRANT_LEN (q->level)) -
      (p->x + P4EST_QUADRANT_LEN (p->level));
    break;
  case 2:
    distance = p->y - q->y;
    break;
  case 3:
    distance = (q->y + P4EST_QUADRANT_LEN (q->level)) -
      (p->y + P4EST_QUADRANT_LEN (p->level));
    break;
#ifdef P4_TO_P8
  case 4:
    distance = p->z - q->z;
    break;
  case 5:
    distance = (q->z + P4EST_QUADRANT_LEN (q->level)) -
      (p->z + P4EST_QUADRANT_LEN (p->level));
    break;
#endif
  default:
    SC_ABORT_NOT_REACHED ();
  }

  blevel = p4est_balance_kernel_1d (distance, q->level);

  if (blevel <= p->level) {
    if (consistent != NULL) {
      *consistent = 1;
    }
    if (conextra != NULL) {
      conextra[nconextra / 2] = p->level;
    }
    return;
  }

  if (consistent != NULL) {
    *consistent = 0;
  }
  if (conextra != NULL) {
    conextra[nconextra / 2] = blevel;

    if (blevel <= p->level + 1) {
      return;
    }

    a = *q;

    /* shift a until it is inside p */
    switch (face) {
    case 0:
      a.x = q->x + distance;
      break;
    case 1:
      a.x = q->x - distance;
      break;
    case 2:
      a.y = q->y + distance;
      break;
    case 3:
      a.y = q->y - distance;
      break;
#ifdef P4_TO_P8
    case 4:
      a.z = q->z + distance;
      break;
    case 5:
      a.z = q->z - distance;
      break;
#endif
    default:
      SC_ABORT_NOT_REACHED ();
    }

    /* get a's ancestor on blevel - 1 */
    a.x &= (-1 << (P4EST_MAXLEVEL - (blevel - 1)));
    a.y &= (-1 << (P4EST_MAXLEVEL - (blevel - 1)));
#ifdef P4_TO_P8
    a.z &= (-1 << (P4EST_MAXLEVEL - (blevel - 1)));
#endif
    a.level = blevel - 1;

#ifndef P4_TO_P8
    for (i = -1; i <= 1; i += 2) {
      temp = a;
      /* a is in a neighboring family groupone family group over from temp */
      if (face / 2 == 0) {
        temp.y += i * P4EST_QUADRANT_LEN (blevel - 1);
      }
      else {
        temp.x += i * P4EST_QUADRANT_LEN (blevel - 1);
      }

      if ((temp.x & (-1 << (P4EST_MAXLEVEL - (p->level)))) != p->x ||
          (temp.y & (-1 << (P4EST_MAXLEVEL - (p->level)))) != p->y) {
        /* only test other descendents of p */
        continue;
      }

      child = p4est_face_corners[face][(1 - i) / 2];

      p4est_bal_corner_con_internal (q, &temp, child, balance, &recon,
                                     recurse);

      conextra[1 + i] = recon ? blevel - 1 : blevel;
    }
#else
    for (j = -1; j <= 1; j++) {
      for (i = -1; i <= 1; i++) {
        if (!i & !j) {
          continue;
        }
        temp = a;
        switch (face / 2) {
        case 0:
          temp.y += i * P4EST_QUADRANT_LEN (blevel - 1);
          temp.z += j * P4EST_QUADRANT_LEN (blevel - 1);
          break;
        case 1:
          temp.x += i * P4EST_QUADRANT_LEN (blevel - 1);
          temp.z += j * P4EST_QUADRANT_LEN (blevel - 1);
          break;
        case 2:
          temp.x += i * P4EST_QUADRANT_LEN (blevel - 1);
          temp.y += j * P4EST_QUADRANT_LEN (blevel - 1);
          break;
        default:
          SC_ABORT_NOT_REACHED ();
        }

        if ((temp.x & (-1 << (P4EST_MAXLEVEL - (p->level)))) != p->x ||
            (temp.y & (-1 << (P4EST_MAXLEVEL - (p->level)))) != p->y ||
            (temp.z & (-1 << (P4EST_MAXLEVEL - (p->level)))) != p->z) {
          /* only test other descendents of p */
          continue;
        }

        if (i & j) {
          child = p4est_face_corners[face][(1 - j) + (1 - i) / 2];

          p4est_bal_corner_con_internal (q, &temp, child, balance, &recon,
                                         recurse);

          conextra[4 + 3 * j + i] = recon ? blevel - 1 : blevel;
        }
        else {
          if (!i) {
            edge = p8est_face_edges[face][(1 - j) / 2];
          }
          else {
            edge = p8est_face_edges[face][2 + (1 - i) / 2];
          }

          p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon,
                                       reconextra, recurse);

          conextra[4 + 3 * j + i] = recon ? blevel - 1 : blevel;
        }

      }
    }
#endif
  }
  return;

#ifdef P4EST_DEBUG
  P4EST_ASSERT (q->level >= p->level);
  /* get the ancestor of q at the same level as p */
  a = *q;
  a.x &= (-1 << (P4EST_MAXLEVEL - p->level));
  a.y &= (-1 << (P4EST_MAXLEVEL - p->level));
#ifdef P4_TO_P8
  a.z &= (-1 << (P4EST_MAXLEVEL - p->level));
#endif
  a.level = p->level;

  /* check that q really is a descendent of the correct neighbor of p */
  p4est_quadrant_face_neighbor (&a, dual, &temp);
  P4EST_ASSERT (p4est_quadrant_is_equal (&temp, p));
#endif

  /* start with the assumption that the two are consistent */
  if (consistent != NULL) {
    *consistent = 1;
    if (conextra != NULL) {
      conextra[nconextra / 2] = p->level;
    }
  }

  /* if q is larger than half of p, or p is a smallest quadrant, then of
   * course q and p are balance consistent */
  if (q->level <= level || level > P4EST_QMAXLEVEL) {
    return;
  }

  child = p4est_quadrant_ancestor_id (q, level);
  cfc = p4est_corner_face_corners[child][dual];
  if (cfc >= 0) {
    if (consistent != NULL) {
      *consistent = 0;
    }
    if (recurse) {
      child = p4est_face_corners[face][cfc];
      p->level++;
      p->x += ((child & 1) ? P4EST_QUADRANT_LEN (p->level) : 0);
      p->y += ((child & 2) ? P4EST_QUADRANT_LEN (p->level) : 0);
#ifdef P4_TO_P8
      p->z += ((child & 4) ? P4EST_QUADRANT_LEN (p->level) : 0);
#endif
      P4EST_ASSERT (p4est_quadrant_is_extended (p));

      if (conextra != NULL) {
        conextra[nconextra / 2] = p->level;
#ifndef P4_TO_P8
        cfcextra = cfc ^ 1;
        extraid = (cfcextra < cfc) ? 0 : 2;
#else
        cfcextra = cfc ^ 3;
        extraid = ((cfcextra & 1) < (cfc & 1)) ? 0 : 2;
        extraid += ((cfcextra & 2) < (cfc & 2)) ? 0 : 6;
#endif

        p4est_quadrant_sibling (p, &temp, p4est_face_corners[face][cfcextra]);
        p4est_bal_corner_con_internal (q, &temp, child, balance, &recon,
                                       recurse);
        if (temp.level > p->level) {
          conextra[extraid] = temp.level;
        }

#ifdef P4_TO_P8
        cfcextra = cfc ^ 1;
        extraid = 3 + ((cfcextra < cfc) ? 0 : 2);

        p4est_quadrant_sibling (p, &temp, p4est_face_corners[face][cfcextra]);
        edge = p8est_face_edges[face][2 + (cfc & 1)];
        reconextra[0] = reconextra[1] = reconextra[2] = p->level;
        p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon,
                                     reconextra, recurse);

        for (i = 0; i < 3; i++) {
          if (reconextra[i] > p->level) {
            conextra[extraid + (i - 1) * 3] = reconextra[i];
          }
        }

        cfcextra = cfc ^ 2;
        extraid = 1 + ((cfcextra < cfc) ? 0 : 6);

        p4est_quadrant_sibling (p, &temp, p4est_face_corners[face][cfcextra]);
        edge = p8est_face_edges[face][((cfc & 2) >> 1)];
        reconextra[0] = reconextra[1] = reconextra[2] = p->level;
        p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon,
                                     reconextra, recurse);

        for (i = 0; i < 3; i++) {
          if (reconextra[i] > p->level) {
            conextra[extraid + (i - 1)] = reconextra[i];
          }
        }
#endif
      }

      p4est_bal_face_con_internal (q, p, face, balance, NULL, conextra,
                                   recurse);
    }
    return;
  }

  /* set a to be the child of p's neighbor that contains q */
  a.x = (q->x & (-1 << (P4EST_MAXLEVEL - level)));
  a.y = (q->y & (-1 << (P4EST_MAXLEVEL - level)));
#ifdef P4_TO_P8
  a.z = (q->z & (-1 << (P4EST_MAXLEVEL - level)));
#endif
  a.level = level;
  P4EST_ASSERT (p4est_quadrant_is_ancestor (&a, q));

  p4est_quadrant_face_neighbor (&a, dual, &temp);
  p4est_bal_face_con_internal (q, &temp, face, balance, &recon, NULL, 0);
  if (!recon) {
    if (consistent != NULL) {
      *consistent = 0;
    }
    if (recurse) {
      p->level++;
      p->x += ((child & 1) ? P4EST_QUADRANT_LEN (p->level) : 0);
      p->y += ((child & 2) ? P4EST_QUADRANT_LEN (p->level) : 0);
#ifdef P4_TO_P8
      p->z += ((child & 4) ? P4EST_QUADRANT_LEN (p->level) : 0);
#endif
      P4EST_ASSERT (p4est_quadrant_is_extended (p));

      if (conextra != NULL) {
        conextra[nconextra / 2] = p->level;
      }
    }
    return;
  }
}

#ifdef P4_TO_P8
static void
p8est_bal_edge_con_internal (p4est_quadrant_t const *restrict q,
                             p4est_quadrant_t * restrict p, int edge,
                             int balance, int *consistent, int *conextra,
                             int recurse)
{
  int                 level = p->level + 1;     /* level is now the size of a
                                                   child */
  int                 child;
  int                 child2;
  int                 cec;
  int                 cecextra;
  int                 dual = edge ^ 3;
  int                 recon;
  p4est_quadrant_t    temp;
  p4est_quadrant_t    a;
  int                 extraid;
  int                 face;
  int                 relative;

#ifdef P4EST_DEBUG
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (p));
  P4EST_ASSERT (q->level >= p->level);
  /* get the ancestor of q at the same level as p */
  a = *q;
  a.x &= (-1 << (P4EST_MAXLEVEL - p->level));
  a.y &= (-1 << (P4EST_MAXLEVEL - p->level));
#ifdef P4_TO_P8
  a.z &= (-1 << (P4EST_MAXLEVEL - p->level));
#endif
  a.level = p->level;

  /* check that q really is a descendent of the correct neighbor of p */
  p8est_quadrant_edge_neighbor (&a, dual, &temp);
  P4EST_ASSERT (p4est_quadrant_is_equal (&temp, p));
#endif

  /* start with the assumption that the two are consistent */
  if (consistent != NULL) {
    *consistent = 1;
    if (conextra != NULL) {
      conextra[1] = p->level;
    }
  }

  /* if q is larger than half of p, or p is a smallest quadrant, then of
   * course q and p are balance consistent */
  if (q->level <= level || level > P4EST_QMAXLEVEL) {
    return;
  }

  child = p4est_quadrant_ancestor_id (q, level);
  if (!p8est_balance_edge_relevant (balance, edge, child)) {
    return;
  }

  cec = -1;
  if (p8est_edge_corners[dual][0] == child) {
    cec = 0;
  }
  else if (p8est_edge_corners[dual][1] == child) {
    cec = 1;
  }

  if (cec >= 0) {
    if (balance) {
      if (consistent != NULL) {
        *consistent = 0;
      }
      if (recurse) {
        child = p8est_edge_corners[edge][cec];
        p->level++;
        p->x += ((child & 1) ? P4EST_QUADRANT_LEN ((p)->level) : 0);
        p->y += ((child & 2) ? P4EST_QUADRANT_LEN ((p)->level) : 0);
        p->z += ((child & 4) ? P4EST_QUADRANT_LEN ((p)->level) : 0);
        P4EST_ASSERT (p4est_quadrant_is_extended (p));

        if (conextra != NULL) {
          conextra[1] = p->level;

          cecextra = cec ^ 1;
          extraid = (cecextra < cec) ? 0 : 2;

          p4est_quadrant_sibling (p, &temp,
                                  p8est_edge_corners[edge][cecextra]);
          p4est_bal_corner_con_internal (q, &temp, child, balance, &recon,
                                         recurse);

          if (temp.level > p->level) {
            conextra[extraid] = temp.level;
          }
        }

        p8est_bal_edge_con_internal (q, p, edge, balance, NULL, conextra,
                                     recurse);
      }
      return;
    }

    while (++level < q->level) {
      child = p4est_quadrant_ancestor_id (q, level);
      if (p8est_balance_edge_relevant (balance, edge, child)) {
        if (consistent != NULL) {
          *consistent = 0;
        }
        if (recurse) {
          child = p8est_edge_corners[edge][cec];
          p->level++;
          p->x += ((child & 1) ? P4EST_QUADRANT_LEN ((p)->level) : 0);
          p->y += ((child & 2) ? P4EST_QUADRANT_LEN ((p)->level) : 0);
          p->z += ((child & 4) ? P4EST_QUADRANT_LEN ((p)->level) : 0);
          P4EST_ASSERT (p4est_quadrant_is_extended (p));

          if (conextra != NULL) {
            conextra[1] = p->level;

            cecextra = cec ^ 1;
            extraid = (cecextra < cec) ? 0 : 2;

            p4est_quadrant_sibling (p, &temp,
                                    p8est_edge_corners[edge][cecextra]);
            p4est_bal_corner_con_internal (q, &temp, child, balance, &recon,
                                           recurse);

            if (temp.level > p->level) {
              conextra[extraid] = temp.level;
            }
          }

          p8est_bal_edge_con_internal (q, p, edge, balance, NULL, conextra,
                                       recurse);
        }
        return;
      }
    }
    return;
  }

  /* set a to be the child of p's neighbor that contains q */
  a.x = (q->x & (-1 << (P4EST_MAXLEVEL - level)));
  a.y = (q->y & (-1 << (P4EST_MAXLEVEL - level)));
  a.z = (q->z & (-1 << (P4EST_MAXLEVEL - level)));
  a.level = level;
  P4EST_ASSERT (p4est_quadrant_is_ancestor (&a, q));

  if (!balance) {
    p8est_quadrant_edge_neighbor (&a, dual, &temp);
    p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon, NULL, 0);
    if (!recon) {
      if (consistent != NULL) {
        *consistent = 0;
      }
      if (recurse) {
        p->level++;
        p->x += ((child & 1) ? P4EST_QUADRANT_LEN (p->level) : 0);
        p->y += ((child & 2) ? P4EST_QUADRANT_LEN (p->level) : 0);
        p->z += ((child & 4) ? P4EST_QUADRANT_LEN (p->level) : 0);
        P4EST_ASSERT (p4est_quadrant_is_extended (p));

        if (conextra != NULL) {
          conextra[1] = p->level;
        }
      }
      return;
    }
    return;
  }

  if (p8est_edge_corners[p8est_corner_edges[child][edge / 4]][0] == child) {
    cec = 0;
  }
  else {
    cec = 1;
  }
  child2 = p8est_edge_corners[dual][cec];
  face = p8est_corner_faces[child][edge / 4];
  P4EST_ASSERT (face == p8est_corner_faces[child2][edge / 4]);

  relative = (p8est_corner_face_corners[child][face] ^
              p8est_corner_face_corners[child2][face]);

  switch (relative) {
  case 3:
    p8est_quadrant_edge_neighbor (&a, dual, &temp);
    p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon, NULL, 0);
    if (!recon) {
      if (consistent != NULL) {
        *consistent = 0;
      }
      if (recurse) {
        p->level++;
        p->x += ((child & 1) ? P4EST_QUADRANT_LEN (p->level) : 0);
        p->y += ((child & 2) ? P4EST_QUADRANT_LEN (p->level) : 0);
        p->z += ((child & 4) ? P4EST_QUADRANT_LEN (p->level) : 0);
        P4EST_ASSERT (p4est_quadrant_is_extended (p));

        if (conextra != NULL) {
          conextra[1] = p->level;
        }
      }
      return;
    }
    return;
  case 1:
  case 2:
    face = p8est_edge_faces[edge][relative >> 1];
    p8est_quadrant_face_neighbor (&a, face ^ 1, &temp);
    p4est_bal_face_con_internal (q, &temp, face, balance, &recon, NULL, 0);
    if (!recon) {
      if (consistent != NULL) {
        *consistent = 0;
      }
      if (recurse) {
        p->level++;
        p->x += ((child & 1) ? P4EST_QUADRANT_LEN (p->level) : 0);
        p->y += ((child & 2) ? P4EST_QUADRANT_LEN (p->level) : 0);
        p->z += ((child & 4) ? P4EST_QUADRANT_LEN (p->level) : 0);
        P4EST_ASSERT (p4est_quadrant_is_extended (p));

        if (conextra != NULL) {
          conextra[1] = p->level;
        }
      }
      return;
    }
    return;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}
#endif

int
p4est_balance_face_test (p4est_quadrant_t * restrict q,
                         p4est_quadrant_t * restrict p,
                         int face, p4est_connect_type_t balance,
                         sc_array_t * seeds)
{
  p4est_quadrant_t    temp = *p;
  p4est_quadrant_t   *s;
  int                 ibalance;
  int                 consistent;
#ifndef P4_TO_P8
  int                 nextra = 3;
  int                 conextra[3];
#else
  int                 nextra = 9;
  int                 conextra[9];
#endif
  int                 i;

  P4EST_ASSERT (seeds == NULL ||
                seeds->elem_size == sizeof (p4est_quadrant_t));

  if (balance == P4EST_CONNECT_FULL) {
    ibalance = P4EST_DIM - 1;
  }
#ifdef P4_TO_P8
  else if (balance == P8EST_CONNECT_EDGE) {
    ibalance = 1;
  }
#endif
  else {
    ibalance = 0;
  }

  if (seeds == NULL) {
    p4est_bal_face_con_internal (q, &temp, face, ibalance, &consistent, NULL,
                                 0);
    return consistent;
  }
  else {
    for (i = 0; i < nextra; i++) {
      conextra[i] = p->level;
    }
    p4est_bal_face_con_internal (q, &temp, face, ibalance, &consistent,
                                 conextra, 1);

    sc_array_resize (seeds, 0);
    if (!consistent) {
      for (i = 0; i < nextra; i++) {
        if (conextra[i] == temp.level) {
          sc_array_resize (seeds, seeds->elem_count + 1);
          s = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);

          *s = temp;

          switch (i % 3) {
          case 0:
            switch (face / 2) {
            case 0:
              s->y += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 1:
#ifdef P4_TO_P8
            case 2:
#endif
              s->x += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
            break;
          case 1:
            break;
          case 2:
            switch (face / 2) {
            case 0:
              s->y += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 1:
#ifdef P4_TO_P8
            case 2:
#endif
              s->x += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
            break;
          default:
            SC_ABORT_NOT_REACHED ();
          }

#ifdef P4_TO_P8
          switch ((i / 3) % 3) {
          case 0:
            switch (face / 2) {
            case 0:
            case 1:
              s->z += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 2:
              s->y += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
            break;
          case 1:
            break;
          case 2:
            switch (face / 2) {
            case 0:
            case 1:
              s->z += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 2:
              s->y += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
            break;
          default:
            SC_ABORT_NOT_REACHED ();
          }
#endif
        }
      }
    }

    return consistent;
  }
}

int
p4est_balance_corner_test (p4est_quadrant_t * restrict q,
                           p4est_quadrant_t * restrict p,
                           int corner, p4est_connect_type_t balance,
                           sc_array_t * seeds)
{
  p4est_quadrant_t    temp = *p;
  p4est_quadrant_t   *s;
  int                 ibalance;
  int                 consistent;

  P4EST_ASSERT (seeds == NULL ||
                seeds->elem_size == sizeof (p4est_quadrant_t));

  if (balance == P4EST_CONNECT_FULL) {
    ibalance = P4EST_DIM - 1;
  }
#ifdef P4_TO_P8
  else if (balance == P8EST_CONNECT_EDGE) {
    ibalance = 1;
  }
#endif
  else {
    ibalance = 0;
  }

  if (seeds == NULL) {
    p4est_bal_corner_con_internal (q, &temp, corner, ibalance, &consistent,
                                   0);
    return consistent;
  }
  else {
    p4est_bal_corner_con_internal (q, &temp, corner, ibalance, &consistent,
                                   1);

    sc_array_resize (seeds, 0);
    if (!consistent) {
      sc_array_resize (seeds, seeds->elem_count + 1);
      s = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);
      *s = temp;
    }

    return consistent;
  }
}

#ifdef P4_TO_P8
int
p8est_balance_edge_test (p4est_quadrant_t * restrict q,
                         p4est_quadrant_t * restrict p,
                         int edge, p4est_connect_type_t balance,
                         sc_array_t * seeds)
{
  p4est_quadrant_t    temp = *p;
  p4est_quadrant_t   *s;
  int                 ibalance;
  int                 consistent;
  int                 nextra = 3;
  int                 conextra[3];
  int                 i;

  P4EST_ASSERT (seeds == NULL ||
                seeds->elem_size == sizeof (p4est_quadrant_t));

  if (balance == P4EST_CONNECT_FULL) {
    ibalance = P4EST_DIM - 1;
  }
#ifdef P4_TO_P8
  else if (balance == P8EST_CONNECT_EDGE) {
    ibalance = 1;
  }
#endif
  else {
    ibalance = 0;
  }

  if (seeds == NULL) {
    p8est_bal_edge_con_internal (q, &temp, edge, ibalance, &consistent, NULL,
                                 0);
    return consistent;
  }
  else {
    for (i = 0; i < nextra; i++) {
      conextra[i] = p->level;
    }
    p8est_bal_edge_con_internal (q, &temp, edge, ibalance, &consistent,
                                 conextra, 1);

    sc_array_resize (seeds, 0);
    if (!consistent) {
      for (i = 0; i < nextra; i++) {
        if (conextra[i] == temp.level) {
          sc_array_resize (seeds, seeds->elem_count + 1);
          s = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);

          *s = temp;

          switch (i % 3) {
          case 0:
            switch (edge / 4) {
            case 0:
              s->x += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 1:
              s->y += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 2:
              s->z += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
            break;
          case 1:
            break;
          case 2:
            switch (edge / 4) {
            case 0:
              s->x += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 1:
              s->y += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 2:
              s->z += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
            break;
          default:
            SC_ABORT_NOT_REACHED ();
          }
        }
      }
    }

    return consistent;
  }
}
#endif

int
p4est_balance_test (p4est_quadrant_t * restrict q,
                    p4est_quadrant_t * restrict p,
                    p4est_connect_type_t balance, sc_array_t * seeds)
{
  int                 outside[P4EST_DIM];
  int                 i;
  int                 type = 0;
  p4est_qcoord_t      diff;
  p4est_qcoord_t      qc, pc;
  p4est_qcoord_t      pdist = P4EST_QUADRANT_LEN (p->level);
  p4est_qcoord_t      qdist = P4EST_QUADRANT_LEN (q->level);
  p4est_quadrant_t   *s;
  int                 f, c;
#ifdef P4_TO_P8
  int                 e;
#endif

  if (seeds != NULL) {
    sc_array_resize (seeds, 0);
  }

  if (q->level <= p->level + 1) {
    return 1;
  }

  for (i = 0; i < P4EST_DIM; i++) {
    switch (i) {
    case 0:
      qc = q->x;
      pc = p->x;
      break;
    case 1:
      qc = q->y;
      pc = p->y;
      break;
#ifdef P4_TO_P8
    case 2:
      qc = q->z;
      pc = p->z;
      break;
#endif
    default:
      SC_ABORT_NOT_REACHED ();
      break;
    }
    outside[i] = 0;
    if (qc < pc) {
      diff = pc - qc;
      if (diff > pdist) {
        return 1;
      }
      outside[i] = -1;
    }
    else {
      diff = (qc + qdist) - (pc + pdist);
      if (diff > pdist) {
        return 1;
      }
      if (diff > 0) {
        outside[i] = 1;
      }
    }
    type += (outside[i] ? 1 : 0);
  }

  switch (type) {
  case 0:
    sc_array_resize (seeds, seeds->elem_count + 1);
    s = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);
    *s = *q;
    return 0;
  case 1:
    for (i = 0; i < P4EST_DIM; i++) {
      if (outside[i]) {
        f = 2 * i + (outside[i] > 0 ? 1 : 0);
        return p4est_balance_face_test (q, p, f, balance, seeds);
      }
    }
    SC_ABORT_NOT_REACHED ();
    return -1;
  case P4EST_DIM:
    c = 0;
    for (i = 0; i < P4EST_DIM; i++) {
      c += (outside[i] > 0 ? (1 << i) : 0);
    }
    return p4est_balance_corner_test (q, p, c, balance, seeds);
#ifdef P4_TO_P8
  case 2:
    e = 0;
    c = 0;
    for (i = 2; i >= 0; i--) {
      if (outside[i]) {
        c <<= 1;
        c |= (outside[i] > 0 ? 1 : 0);
      }
      else {
        e |= (i << 2);
      }
    }
    e |= c;
    return p8est_balance_edge_test (q, p, e, balance, seeds);
#endif
  default:
    SC_ABORT_NOT_REACHED ();
    return -1;
  }
}
