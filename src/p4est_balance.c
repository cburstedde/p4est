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

  count = (xbit == maxbit) + (ybit == maxbit) + (zbit == maxbit);

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
        xbit = ((dx & (1 << thisbit)) != 0);
        ybit = ((dy & (1 << thisbit)) != 0);
        zbit = ((dz & (1 << thisbit)) != 0);
        count = xbit + ybit + zbit;
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
    if (p8est_balance_2chain (dx, dy, dz, maxbit, (xbit == maxbit),
                              (ybit == maxbit), (zbit == maxbit))) {
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

  P4EST_ASSERT (zfx + zfy + zfz == 2);

  do {
    P4EST_ASSERT (thisbit >= 0);

    xbit = ((dx & (1 << thisbit)) != 0);
    ybit = ((dy & (1 << thisbit)) != 0);
    zbit = ((dz & (1 << thisbit)) != 0);
    count = xbit + ybit + zbit;

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
      if (xbit && zfx) {
        return 0;
      }
      if (ybit && zfy) {
        return 0;
      }
      if (zbit && zfy) {
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

static void         p4est_bal_corner_con_internal (p4est_quadrant_t const
                                                   *restrict q,
                                                   p4est_quadrant_t *
                                                   restrict p, int corner,
                                                   int balance,
                                                   int *consisent);

static void         p4est_bal_face_con_internal (p4est_quadrant_t const
                                                 *restrict q,
                                                 p4est_quadrant_t *
                                                 restrict p, int face,
                                                 int balance,
                                                 int *consisent,
                                                 p4est_quadrant_t * add);

#ifdef P4_TO_P8
static void         p8est_bal_edge_con_internal (p4est_quadrant_t const
                                                 *restrict q,
                                                 p4est_quadrant_t *
                                                 restrict p, int edge,
                                                 int balance,
                                                 int *consistent,
                                                 p4est_quadrant_t * add);
#endif

/* \a corner is the corner of \a p closest to \a q: the corner of \a q closest
 * to \a p is the dual of \a corner */
static void
p4est_bal_corner_con_internal (p4est_quadrant_t const *restrict q,
                               p4est_quadrant_t * restrict p,
                               int corner, int balance, int *consistent)
{
  int                 qlevel = q->level;
  int                 plevel = p->level;
  int                 blevel;
  p4est_qcoord_t      dx, dy, dist;
#ifdef P4_TO_P8
  p4est_qcoord_t      dz;
#endif

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (p));

  if (qlevel <= plevel) {
    if (consistent != NULL) {
      *consistent = 1;
    }
    return;
  }

  dx = (corner & 1) ? ((q->x + P4EST_QUADRANT_LEN (qlevel)) -
                       (p->x + P4EST_QUADRANT_LEN (plevel))) : p->x - q->x;
  P4EST_ASSERT (dx >= 0);
  dy = (corner & 2) ? ((q->y + P4EST_QUADRANT_LEN (qlevel)) -
                       (p->y + P4EST_QUADRANT_LEN (plevel))) : p->y - q->y;
  P4EST_ASSERT (dy >= 0);
#ifdef P4_TO_P8
  dz = (corner & 4) ? ((q->z + P4EST_QUADRANT_LEN (qlevel)) -
                       (p->z + P4EST_QUADRANT_LEN (plevel))) : p->z - q->z;
  P4EST_ASSERT (dz >= 0);
#endif

#ifndef P4_TO_P8
  if (balance) {
    dist = SC_MAX (dx, dy);
    blevel = p4est_balance_kernel_1d (dist, qlevel);
  }
  else {
    blevel = p4est_balance_kernel_2d (dx, dy, qlevel);
  }
#else
  switch (balance) {
  case 0:
    blevel = p8est_balance_kernel_3d_face (dx, dy, dz, qlevel);
    break;
  case 1:
    blevel = p8est_balance_kernel_3d_edge (dx, dy, dz, qlevel);
    break;
  case 2:
    dist = SC_MAX (dx, dy);
    dist = SC_MAX (dist, dz);
    blevel = p4est_balance_kernel_1d (dist, qlevel);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
#endif

  if (blevel <= plevel) {
    if (consistent != NULL) {
      *consistent = 1;
    }
    return;
  }

  if (consistent != NULL) {
    *consistent = 0;
  }

  p->x = q->x + ((corner & 1) ? -dx : dx);
  p->x &= (-1 << (P4EST_MAXLEVEL - blevel));
  p->y = q->y + ((corner & 2) ? -dy : dy);
  p->y &= (-1 << (P4EST_MAXLEVEL - blevel));
#ifdef P4_TO_P8
  p->z = q->z + ((corner & 4) ? -dz : dz);
  p->z &= (-1 << (P4EST_MAXLEVEL - blevel));
#endif
  p->level = blevel;
  P4EST_ASSERT (p4est_quadrant_is_extended (p));
}

static void
p4est_bal_face_con_internal (p4est_quadrant_t const *restrict q,
                             p4est_quadrant_t * restrict p, int face,
                             int balance, int *consistent,
                             p4est_quadrant_t * add)
{
  int                 qlevel = q->level;
  int                 plevel = p->level;
  int                 blevel;

  int                 child;
#ifdef P4_TO_P8
  int                 edge;
#endif
  int                 recon;
  p4est_quadrant_t    porig;
  p4est_quadrant_t    temp;
  p4est_quadrant_t    a;
  int                 i;
#ifndef P4_TO_P8
  int                 nconextra = 3;
#else
  int                 nconextra = 9;
  int                 j;
#endif
  double              distance;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (p));

  if (qlevel <= plevel) {
    if (consistent != NULL) {
      *consistent = 1;
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

  P4EST_ASSERT (distance >= 0);

  blevel = p4est_balance_kernel_1d (distance, q->level);

  if (blevel <= plevel) {
    if (consistent != NULL) {
      *consistent = 1;
    }
    return;
  }

  if (consistent != NULL) {
    *consistent = 0;
  }

  porig = *p;

  *p = *q;

  /* shift a until it is inside p */
  switch (face) {
  case 0:
    p->x += distance;
    break;
  case 1:
    p->x -= distance;
    break;
  case 2:
    p->y += distance;
    break;
  case 3:
    p->y -= distance;
    break;
#ifdef P4_TO_P8
  case 4:
    p->z += distance;
    break;
  case 5:
    p->z -= distance;
    break;
#endif
  default:
    SC_ABORT_NOT_REACHED ();
  }

  p->x &= (-1 << (P4EST_MAXLEVEL - blevel));
  p->y &= (-1 << (P4EST_MAXLEVEL - blevel));
#ifdef P4_TO_P8
  p->z &= (-1 << (P4EST_MAXLEVEL - blevel));
#endif
  p->level = blevel;
  P4EST_ASSERT (p4est_quadrant_is_extended (p));

  if (add != NULL) {

    add[nconextra / 2] = *p;

    a = *p;
    a.x &= (-1 << (P4EST_MAXLEVEL - (blevel - 1)));
    a.y &= (-1 << (P4EST_MAXLEVEL - (blevel - 1)));
#ifdef P4_TO_P8
    a.z &= (-1 << (P4EST_MAXLEVEL - (blevel - 1)));
#endif
    a.level = blevel - 1;

#ifndef P4_TO_P8
    for (i = -1; i <= 1; i += 2) {
      temp = a;
      /* temp is in a family group one family group over from temp */
      if (face / 2 == 0) {
        temp.y += i * P4EST_QUADRANT_LEN (blevel - 1);
      }
      else {
        temp.x += i * P4EST_QUADRANT_LEN (blevel - 1);
      }

      if ((temp.x & (-1 << (P4EST_MAXLEVEL - (plevel)))) != porig.x ||
          (temp.y & (-1 << (P4EST_MAXLEVEL - (plevel)))) != porig.y) {
        /* only test other descendents of p */
        continue;
      }

      child = p4est_face_corners[face][(1 - i) / 2];

      p4est_bal_corner_con_internal (q, &temp, child, balance, &recon);

      if (!recon) {
        add[1 + i] = temp;
      }
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

        if ((temp.x & (-1 << (P4EST_MAXLEVEL - (plevel)))) != porig.x ||
            (temp.y & (-1 << (P4EST_MAXLEVEL - (plevel)))) != porig.y ||
            (temp.z & (-1 << (P4EST_MAXLEVEL - (plevel)))) != porig.z) {
          /* only test other descendents of p */
          continue;
        }

        if (i && j) {

          /* for face only balance, we need to check a larger neighbor in one
           * instance */
          if (!balance) {
            switch (face / 2) {
            case 0:
              if ((temp.y & (-1 << (P4EST_MAXLEVEL - (blevel - 2)))) !=
                  (a.y & (-1 << (P4EST_MAXLEVEL - (blevel - 2)))) &&
                  (temp.z & (-1 << (P4EST_MAXLEVEL - (blevel - 2)))) !=
                  (a.z & (-1 << (P4EST_MAXLEVEL - (blevel - 2))))) {
                temp.x &= (-1 << (P4EST_MAXLEVEL - (blevel - 2)));
                temp.y &= (-1 << (P4EST_MAXLEVEL - (blevel - 2)));
                temp.z &= (-1 << (P4EST_MAXLEVEL - (blevel - 2)));
                temp.level = blevel - 2;
              }
              break;
            case 1:
              if ((temp.x & (-1 << (P4EST_MAXLEVEL - (blevel - 2)))) !=
                  (a.x & (-1 << (P4EST_MAXLEVEL - (blevel - 2)))) &&
                  (temp.z & (-1 << (P4EST_MAXLEVEL - (blevel - 2)))) !=
                  (a.z & (-1 << (P4EST_MAXLEVEL - (blevel - 2))))) {
                temp.x &= (-1 << (P4EST_MAXLEVEL - (blevel - 2)));
                temp.y &= (-1 << (P4EST_MAXLEVEL - (blevel - 2)));
                temp.z &= (-1 << (P4EST_MAXLEVEL - (blevel - 2)));
                temp.level = blevel - 2;
              }
              break;
            case 2:
              if ((temp.x & (-1 << (P4EST_MAXLEVEL - (blevel - 2)))) !=
                  (a.x & (-1 << (P4EST_MAXLEVEL - (blevel - 2)))) &&
                  (temp.y & (-1 << (P4EST_MAXLEVEL - (blevel - 2)))) !=
                  (a.y & (-1 << (P4EST_MAXLEVEL - (blevel - 2))))) {
                temp.x &= (-1 << (P4EST_MAXLEVEL - (blevel - 2)));
                temp.y &= (-1 << (P4EST_MAXLEVEL - (blevel - 2)));
                temp.z &= (-1 << (P4EST_MAXLEVEL - (blevel - 2)));
                temp.level = blevel - 2;
              }
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
          }

          child = p4est_face_corners[face][(1 - j) + (1 - i) / 2];

          p4est_bal_corner_con_internal (q, &temp, child, balance, &recon);

          if (!recon) {
            add[4 + 3 * j + i] = temp;
          }
        }
        else {
          if (!i) {
            edge = p8est_face_edges[face][(1 - j) / 2];
          }
          else {
            edge = p8est_face_edges[face][2 + (1 - i) / 2];
          }

          p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon, NULL);

          if (!recon) {
            add[4 + 3 * j + i] = temp;
          }
        }

      }
    }
#endif
  }
}

#ifdef P4_TO_P8
static void
p8est_bal_edge_con_internal (p4est_quadrant_t const *restrict q,
                             p4est_quadrant_t * restrict p, int edge,
                             int balance, int *consistent,
                             p4est_quadrant_t * add)
{
  int                 plevel = p->level;
  int                 qlevel = q->level;
  int                 blevel;
  int                 child;
  int                 recon;
  p4est_quadrant_t    porig;
  p4est_quadrant_t    temp;
  p4est_quadrant_t    a;
  p4est_qcoord_t      dx, dy;
  p4est_qcoord_t      dist;
  int                 i;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (p));

  if (qlevel <= plevel) {
    if (consistent != NULL) {
      *consistent = 1;
    }
    return;
  }

  switch (edge / 4) {
  case 0:
    dx = (edge & 1) ? (q->y + P4EST_QUADRANT_LEN (qlevel)) -
      (p->y + P4EST_QUADRANT_LEN (plevel)) : p->y - q->y;
    dy = (edge & 2) ? (q->z + P4EST_QUADRANT_LEN (qlevel)) -
      (p->z + P4EST_QUADRANT_LEN (plevel)) : p->z - q->z;
    break;
  case 1:
    dx = (edge & 1) ? (q->x + P4EST_QUADRANT_LEN (qlevel)) -
      (p->x + P4EST_QUADRANT_LEN (plevel)) : p->x - q->x;
    dy = (edge & 2) ? (q->z + P4EST_QUADRANT_LEN (qlevel)) -
      (p->z + P4EST_QUADRANT_LEN (plevel)) : p->z - q->z;
    break;
  case 2:
    dx = (edge & 1) ? (q->x + P4EST_QUADRANT_LEN (qlevel)) -
      (p->x + P4EST_QUADRANT_LEN (plevel)) : p->x - q->x;
    dy = (edge & 2) ? (q->y + P4EST_QUADRANT_LEN (qlevel)) -
      (p->y + P4EST_QUADRANT_LEN (plevel)) : p->y - q->y;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  P4EST_ASSERT (dx >= 0);
  P4EST_ASSERT (dy >= 0);

  if (balance) {
    dist = SC_MAX (dx, dy);
    blevel = p4est_balance_kernel_1d (dist, qlevel);
  }
  else {
    blevel = p4est_balance_kernel_2d (dx, dy, qlevel);
  }

  if (blevel <= plevel) {
    if (consistent != NULL) {
      *consistent = 1;
    }
    return;
  }

  if (consistent != NULL) {
    *consistent = 0;
  }

  porig = *p;
  *p = *q;

  switch (edge / 4) {
  case 0:
    p->y += (edge & 1) ? -dx : dx;
    p->z += (edge & 2) ? -dy : dy;
    break;
  case 1:
    p->x += (edge & 1) ? -dx : dx;
    p->z += (edge & 2) ? -dy : dy;
    break;
  case 2:
    p->x += (edge & 1) ? -dx : dx;
    p->y += (edge & 2) ? -dy : dy;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  p->x &= (-1 << (P4EST_MAXLEVEL - blevel));
  p->y &= (-1 << (P4EST_MAXLEVEL - blevel));
  p->z &= (-1 << (P4EST_MAXLEVEL - blevel));
  p->level = blevel;
  P4EST_ASSERT (p4est_quadrant_is_extended (p));

  if (add != NULL) {
    add[1] = *p;

    a = *p;
    a.x &= (-1 << (P4EST_MAXLEVEL - (blevel - 1)));
    a.y &= (-1 << (P4EST_MAXLEVEL - (blevel - 1)));
    a.z &= (-1 << (P4EST_MAXLEVEL - (blevel - 1)));
    a.level = blevel - 1;

    for (i = -1; i <= 1; i += 2) {
      temp = a;
      /* temp is in a family group one family group over from temp */
      switch (edge / 4) {
      case 0:
        temp.x += i * P4EST_QUADRANT_LEN (blevel - 1);
        break;
      case 1:
        temp.y += i * P4EST_QUADRANT_LEN (blevel - 1);
        break;
      case 2:
        temp.z += i * P4EST_QUADRANT_LEN (blevel - 1);
        break;
      default:
        SC_ABORT_NOT_REACHED ();
      }

      if ((temp.x & (-1 << (P4EST_MAXLEVEL - (plevel)))) != porig.x ||
          (temp.y & (-1 << (P4EST_MAXLEVEL - (plevel)))) != porig.y ||
          (temp.z & (-1 << (P4EST_MAXLEVEL - (plevel)))) != porig.z) {
        /* only test other descendents of p */
        continue;
      }

      child = p8est_edge_corners[edge][(1 - i) / 2];

      p4est_bal_corner_con_internal (q, &temp, child, balance, &recon);

      if (!recon) {
        add[1 + i] = temp;
      }
    }
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
  p4est_quadrant_t    add[3];
#else
  int                 nextra = 9;
  p8est_quadrant_t    add[9];
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
    p4est_bal_face_con_internal (q, &temp, face, ibalance, &consistent, NULL);
    return consistent;
  }
  else {
    memset (add, -1, nextra * sizeof (p4est_quadrant_t));
    p4est_bal_face_con_internal (q, &temp, face, ibalance, &consistent, add);

    sc_array_resize (seeds, 0);
    if (!consistent) {
      for (i = 0; i < nextra; i++) {
        if (add[i].level != -1) {
          sc_array_resize (seeds, seeds->elem_count + 1);
          s = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);

          *s = add[i];
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

  p4est_bal_corner_con_internal (q, &temp, corner, ibalance, &consistent);
  if (seeds == NULL) {
    return consistent;
  }
  else {
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
  p4est_quadrant_t    add[3];
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
    p8est_bal_edge_con_internal (q, &temp, edge, ibalance, &consistent, NULL);

    return consistent;
  }
  else {
    memset (add, -1, nextra * sizeof (p4est_quadrant_t));

    p8est_bal_edge_con_internal (q, &temp, edge, ibalance, &consistent, add);

    sc_array_resize (seeds, 0);
    if (!consistent) {
      for (i = 0; i < nextra; i++) {
        if (add[i].level != -1) {
          sc_array_resize (seeds, seeds->elem_count + 1);
          s = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);

          *s = add[i];
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
