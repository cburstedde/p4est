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
static inline int
p4est_balance_kernel_1d (p4est_qcoord_t distance, int level)
{
  int                 shift = P4EST_MAXLEVEL - level;
  P4EST_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);
  /* the distance only makes sense if it is an integer number of \a level
   * distances */
  P4EST_ASSERT (distance >= 0);
  P4EST_ASSERT (!(distance & (~(((p4est_qcoord_t) - 1) << shift))));
  distance >>= shift;
  /* The theory says we should use ((distance + 1)&(~1) + 1), but
   * using distance + 1 is equivalent for all distance >= 0 */
  distance++;

  return SC_MAX (0, level - SC_LOG2_32 (distance));
}

/* This is the kernel for 2D balance with face-only balancing */
static inline int
p4est_balance_kernel_2d (p4est_qcoord_t dx, p4est_qcoord_t dy, int level)
{
  int                 shift = P4EST_MAXLEVEL - level;
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

  return SC_MAX (0, level - SC_LOG2_32 (distance));
}

#ifdef P4_TO_P8
/* This is the kernel for 3d balance with face and edge balancing */
static inline int
p8est_balance_kernel_3d_edge (p4est_qcoord_t dx, p4est_qcoord_t dy,
                              p4est_qcoord_t dz, int level)
{
  int                 shift = P4EST_MAXLEVEL - level;
  int                 xbit, ybit, zbit;
  int                 maxbit;
  p4est_qcoord_t      bitwor;
  int                 ret;

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

  /* we want to carry a 1 when there are three 1 bits, so we find places
   * where there is at least one 1 bit and subtract one 1 bit: then if we
   * sum, the binary carry rule will give the correct result */
  bitwor = (dx | dy | dz);
  ret = SC_LOG2_32 (dx + dy + dz - bitwor);
  /* we have to guard against the case when the leading position has one 1
   * bit. */
  ret = SC_MAX (maxbit, ret);
  return SC_MAX (0, level - ret);
}

/* This is the kernel for 3d balance with face balancing only */
static inline int
p8est_balance_kernel_3d_face (p4est_qcoord_t dx, p4est_qcoord_t dy,
                              p4est_qcoord_t dz, int level)
{
  int                 shift = P4EST_MAXLEVEL - level;
  int                 maxbit;
  int                 yzbit, zxbit, xybit;
  p4est_qcoord_t      dyz, dzx, dxy, bitwor;
  int                 ret;

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

  /* this problem is dual to dual to kernel 3d edge */
  dyz = dy + dz;
  dzx = dz + dx;
  dxy = dx + dy;

  yzbit = SC_LOG2_32 (dyz);
  maxbit = yzbit;
  zxbit = SC_LOG2_32 (dzx);
  maxbit = SC_MAX (maxbit, zxbit);
  xybit = SC_LOG2_32 (dxy);
  maxbit = SC_MAX (maxbit, xybit);
  P4EST_ASSERT (maxbit >= 1);

  /* we want to carry a 1 when there are three 1 bits, so we find places
   * where there is at least one 1 bit and subtract one 1 bit: then if we
   * sum, the binary carry rule will give the correct result */
  bitwor = (dyz | dzx | dxy);
  ret = SC_LOG2_32 (dyz + dzx + dxy - bitwor);
  /* we have to guard against the case when the leading position has one 1
   * bit. */
  ret = SC_MAX (maxbit, ret);
  return SC_MAX (0, level - ret);
}
#endif

static void         p4est_bal_corner_con_internal (p4est_quadrant_t const *q,
                                                   p4est_quadrant_t * p,
                                                   int corner,
                                                   int balance,
                                                   int *consisent);

static void         p4est_bal_face_con_internal (p4est_quadrant_t const *q,
                                                 p4est_quadrant_t * p,
                                                 int face, int balance,
                                                 int *consisent,
                                                 p4est_quadrant_t * add);

#ifdef P4_TO_P8
static void         p8est_bal_edge_con_internal (p4est_quadrant_t const *q,
                                                 p4est_quadrant_t * p,
                                                 int edge, int balance,
                                                 int *consistent,
                                                 p4est_quadrant_t * add);
#endif

/* \a corner is the corner of \a p closest to \a q: the corner of \a q closest
 * to \a p is the dual of \a corner */
static void
p4est_bal_corner_con_internal (p4est_quadrant_t const *q,
                               p4est_quadrant_t * p,
                               int corner, int balance, int *consistent)
{
  int                 qlevel = q->level;
  int                 plevel = p->level;
  int                 blevel;
  p4est_qcoord_t      qlen, plen, mask;
  p4est_qcoord_t      dx, dy, dist;
#ifdef P4_TO_P8
  p4est_qcoord_t      dz;
#endif

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (p));

  if (qlevel <= plevel) {
    if (consistent != NULL) {
      *consistent = 1;
    }
    return;
  }

  qlen = P4EST_QUADRANT_LEN (qlevel);
  plen = P4EST_QUADRANT_LEN (plevel);

  dx = (corner & 1) ? ((q->x + qlen) - (p->x + plen)) : p->x - q->x;
  P4EST_ASSERT (dx >= 0);
  dy = (corner & 2) ? ((q->y + qlen) - (p->y + plen)) : p->y - q->y;
  P4EST_ASSERT (dy >= 0);
#ifdef P4_TO_P8
  dz = (corner & 4) ? ((q->z + qlen) - (p->z + plen)) : p->z - q->z;
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

  mask = -1 << (P4EST_MAXLEVEL - blevel);
  p->x = q->x + ((corner & 1) ? -dx : dx);
  p->x &= mask;
  p->y = q->y + ((corner & 2) ? -dy : dy);
  p->y &= mask;
#ifdef P4_TO_P8
  p->z = q->z + ((corner & 4) ? -dz : dz);
  p->z &= mask;
#endif
  p->level = blevel;
  P4EST_ASSERT (p4est_quadrant_is_extended (p));
}

static void
p4est_bal_face_con_internal (p4est_quadrant_t const *q,
                             p4est_quadrant_t * p, int face,
                             int balance, int *consistent,
                             p4est_quadrant_t * add)
{
  int                 qlevel = q->level;
  int                 plevel = p->level;
  int                 blevel;

  int                 child;
#ifdef P4_TO_P8
  int                 edge;
  int                 dual;
  int                 achild = -1;
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
  p4est_qcoord_t      b2mask;
#endif
  p4est_qcoord_t      distance;
  p4est_qcoord_t      qlen, plen, mask, pmask;
  p4est_qcoord_t      b1len;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (p));

  if (qlevel <= plevel) {
    if (consistent != NULL) {
      *consistent = 1;
    }
    return;
  }

  qlen = P4EST_QUADRANT_LEN (qlevel);
  plen = P4EST_QUADRANT_LEN (plevel);

  switch (face) {
  case 0:
    distance = p->x - q->x;
    break;
  case 1:
    distance = (q->x + qlen) - (p->x + plen);
    break;
  case 2:
    distance = p->y - q->y;
    break;
  case 3:
    distance = (q->y + qlen) - (p->y + plen);
    break;
#ifdef P4_TO_P8
  case 4:
    distance = p->z - q->z;
    break;
  case 5:
    distance = (q->z + qlen) - (p->z + plen);
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

  mask = -1 << (P4EST_MAXLEVEL - blevel);
  p->x &= mask;
  p->y &= mask;
#ifdef P4_TO_P8
  p->z &= mask;
#endif
  p->level = blevel;
  P4EST_ASSERT (p4est_quadrant_is_extended (p));

  if (add != NULL) {

    add[nconextra / 2] = *p;

    /* this is the only quad needed if it is only one level smaller than the
     * original quadrant */
    if (blevel == plevel - 1) {
      return;
    }

    mask = -1 << (P4EST_MAXLEVEL - (blevel - 1));
    pmask = -1 << (P4EST_MAXLEVEL - (plevel));
    a = *p;
    a.x &= mask;
    a.y &= mask;
#ifdef P4_TO_P8
    a.z &= mask;
#endif
    a.level = blevel - 1;

    b1len = P4EST_QUADRANT_LEN (blevel - 1);
#ifndef P4_TO_P8
    for (i = -1; i <= 1; i += 2) {
      temp = a;
      /* temp is in a family group one family group over from temp */
      if (face / 2 == 0) {
        temp.y += i * b1len;
      }
      else {
        temp.x += i * b1len;
      }

      if ((temp.x & pmask) != porig.x || (temp.y & pmask) != porig.y) {
        /* only test other descendants of p */
        continue;
      }

      child = p4est_face_corners[face][(1 - i) / 2];

      p4est_bal_corner_con_internal (q, &temp, child, balance, &recon);

      if (!recon) {
        add[1 + i] = temp;
      }
    }
#else
    b2mask = -1 << (P4EST_MAXLEVEL - (blevel - 2));
    if (!balance) {
      achild = p8est_quadrant_child_id (&a);
    }
    for (j = -1; j <= 1; j++) {
      for (i = -1; i <= 1; i++) {
        if (!i & !j) {
          continue;
        }
        temp = a;
        switch (face / 2) {
        case 0:
          temp.y += i * b1len;
          temp.z += j * b1len;
          break;
        case 1:
          temp.x += i * b1len;
          temp.z += j * b1len;
          break;
        case 2:
          temp.x += i * b1len;
          temp.y += j * b1len;
          break;
        default:
          SC_ABORT_NOT_REACHED ();
        }

        if ((temp.x & pmask) != porig.x || (temp.y & pmask) != porig.y ||
            (temp.z & pmask) != porig.z) {
          /* only test other descendants of p */
          continue;
        }

        if (i && j) {

          child = p8est_face_corners[face][(1 - j) + (1 - i) / 2];

          /* for face only balance, we need to check a larger neighbor in one
           * instance */
          if (!balance) {
            dual = p8est_face_corners[face][(1 + j) + (1 + i) / 2];

            if (achild == dual) {
              temp.x &= b2mask;
              temp.y &= b2mask;
              temp.z &= b2mask;
              temp.level = blevel - 2;
            }
          }

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

    if (!balance) {
      for (j = -1; j <= 1; j += 2) {
        for (i = -1; i <= 1; i += 2) {
          if (add[4 + 3 * j + i].level != -1 &&
              add[4 + 3 * j + i].level < blevel) {
            if (add[4 + 3 * j].level != -1 || add[4 + i].level != -1) {
              memset (&(add[4 + 3 * j + i]), -1, sizeof (p4est_quadrant_t));
            }
          }
        }
      }
    }
#endif
  }
}

#ifdef P4_TO_P8
static void
p8est_bal_edge_con_internal (p4est_quadrant_t const *q,
                             p4est_quadrant_t * p, int edge,
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
  p4est_qcoord_t      qlen, plen, mask;
  p4est_qcoord_t      b1len, pmask;
  int                 i;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (p));

  if (qlevel <= plevel) {
    if (consistent != NULL) {
      *consistent = 1;
    }
    return;
  }

  qlen = P4EST_QUADRANT_LEN (qlevel);
  plen = P4EST_QUADRANT_LEN (plevel);

  switch (edge / 4) {
  case 0:
    dx = (edge & 1) ? (q->y + qlen) - (p->y + plen) : p->y - q->y;
    dy = (edge & 2) ? (q->z + qlen) - (p->z + plen) : p->z - q->z;
    break;
  case 1:
    dx = (edge & 1) ? (q->x + qlen) - (p->x + plen) : p->x - q->x;
    dy = (edge & 2) ? (q->z + qlen) - (p->z + plen) : p->z - q->z;
    break;
  case 2:
    dx = (edge & 1) ? (q->x + qlen) - (p->x + plen) : p->x - q->x;
    dy = (edge & 2) ? (q->y + qlen) - (p->y + plen) : p->y - q->y;
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
  mask = -1 << (P4EST_MAXLEVEL - blevel);
  p->x &= mask;
  p->y &= mask;
  p->z &= mask;
  p->level = blevel;
  P4EST_ASSERT (p4est_quadrant_is_extended (p));

  if (add != NULL) {
    add[1] = *p;

    /* this is the only quad needed if it is only one level smaller than the
     * original quadrant */
    if (blevel == plevel - 1) {
      return;
    }

    mask = -1 << (P4EST_MAXLEVEL - (blevel - 1));
    pmask = -1 << (P4EST_MAXLEVEL - (plevel));
    a = *p;
    a.x &= mask;
    a.y &= mask;
    a.z &= mask;
    a.level = blevel - 1;

    b1len = P4EST_QUADRANT_LEN (blevel - 1);
    for (i = -1; i <= 1; i += 2) {
      temp = a;
      /* temp is in a family group one family group over from temp */
      switch (edge / 4) {
      case 0:
        temp.x += i * b1len;
        break;
      case 1:
        temp.y += i * b1len;
        break;
      case 2:
        temp.z += i * b1len;
        break;
      default:
        SC_ABORT_NOT_REACHED ();
      }

      if ((temp.x & pmask) != porig.x || (temp.y & pmask) != porig.y ||
          (temp.z & pmask) != porig.z) {
        /* only test other descendants of p */
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
p4est_balance_seeds_face (p4est_quadrant_t * q,
                          p4est_quadrant_t * p,
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
    return !consistent;
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

    return !consistent;
  }
}

int
p4est_balance_seeds_corner (p4est_quadrant_t * q,
                            p4est_quadrant_t * p,
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
    return !consistent;
  }
  else {
    sc_array_resize (seeds, 0);
    if (!consistent) {
      sc_array_resize (seeds, seeds->elem_count + 1);
      s = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);

      *s = temp;
    }

    return !consistent;
  }
}

#ifdef P4_TO_P8
int
p8est_balance_seeds_edge (p4est_quadrant_t * q,
                          p4est_quadrant_t * p,
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

    return !consistent;
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

    return !consistent;
  }
}
#endif

int
p4est_balance_seeds (p4est_quadrant_t * q, p4est_quadrant_t * p,
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

  /* basic level comparison */
  if (q->level <= p->level + 1) {
    return 0;
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
      /* insulation layer comparison */
      if (diff > pdist) {
        return 0;
      }
      outside[i] = -1;
    }
    else {
      diff = (qc + qdist) - (pc + pdist);
      /* insulation layer comparison */
      if (diff > pdist) {
        return 0;
      }
      if (diff > 0) {
        outside[i] = 1;
      }
    }
    type += (outside[i] ? 1 : 0);
  }

  switch (type) {
  case 0:
    /* q is inside p, so it is its own seed */
    sc_array_resize (seeds, seeds->elem_count + 1);
    s = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);
    *s = *q;
    return 1;
  case 1:
    for (i = 0; i < P4EST_DIM; i++) {
      if (outside[i]) {
        f = 2 * i + (outside[i] > 0 ? 1 : 0);
        return p4est_balance_seeds_face (q, p, f, balance, seeds);
      }
    }
    SC_ABORT_NOT_REACHED ();
    return -1;
  case P4EST_DIM:
    c = 0;
    for (i = 0; i < P4EST_DIM; i++) {
      c += (outside[i] > 0 ? (1 << i) : 0);
    }
    return p4est_balance_seeds_corner (q, p, c, balance, seeds);
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
    return p8est_balance_seeds_edge (q, p, e, balance, seeds);
#endif
  default:
    SC_ABORT_NOT_REACHED ();
    return -1;
  }
}
