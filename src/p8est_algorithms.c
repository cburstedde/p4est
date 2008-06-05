/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2008 Carsten Burstedde, Lucas Wilcox.

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

#include <p4est_to_p8est.h>
#include "p4est_algorithms.c"

bool
p8est_quadrant_is_family (const p4est_quadrant_t * q0,
                          const p4est_quadrant_t * q1,
                          const p4est_quadrant_t * q2,
                          const p4est_quadrant_t * q3,
                          const p4est_quadrant_t * q4,
                          const p4est_quadrant_t * q5,
                          const p4est_quadrant_t * q6,
                          const p4est_quadrant_t * q7)
{
  const int8_t        level = q0->level;
  p4est_qcoord_t      inc;

  P4EST_ASSERT (p4est_quadrant_is_extended (q0));
  P4EST_ASSERT (p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_extended (q2));
  P4EST_ASSERT (p4est_quadrant_is_extended (q3));
  P4EST_ASSERT (p4est_quadrant_is_extended (q4));
  P4EST_ASSERT (p4est_quadrant_is_extended (q5));
  P4EST_ASSERT (p4est_quadrant_is_extended (q6));
  P4EST_ASSERT (p4est_quadrant_is_extended (q7));

  if (level == 0 || level != q1->level ||
      level != q2->level || level != q3->level ||
      level != q4->level || level != q5->level ||
      level != q6->level || level != q7->level) {
    return false;
  }

  inc = P4EST_QUADRANT_LEN (level);
  return ((q0->x + inc == q1->x && q0->y == q1->y && q0->z == q1->z) &&
          (q0->x == q2->x && q0->y + inc == q2->y && q0->z == q2->z) &&
          (q1->x == q3->x && q2->y == q3->y && q0->z == q3->z) &&
          (q0->x == q4->x && q0->y == q4->y && q0->z + inc == q4->z) &&
          (q1->x == q5->x && q1->y == q5->y && q4->z == q5->z) &&
          (q2->x == q6->x && q2->y == q6->y && q4->z == q6->z) &&
          (q3->x == q7->x && q3->y == q7->y && q4->z == q7->z));
}

void
p8est_quadrant_children (const p4est_quadrant_t * q,
                         p4est_quadrant_t * c0, p4est_quadrant_t * c1,
                         p4est_quadrant_t * c2, p4est_quadrant_t * c3,
                         p4est_quadrant_t * c4, p4est_quadrant_t * c5,
                         p4est_quadrant_t * c6, p4est_quadrant_t * c7)
{
  const int8_t        level = q->level + 1;     /* same type */
  const p4est_qcoord_t inc = P4EST_QUADRANT_LEN (level);

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (q->level < P4EST_MAXLEVEL);

  c0->x = q->x;
  c0->y = q->y;
  c0->z = q->z;
  c0->level = level;

  c1->x = c0->x | inc;
  c1->y = c0->y;
  c1->z = c0->z;
  c1->level = level;

  c2->x = c0->x;
  c2->y = c0->y | inc;
  c2->z = c0->z;
  c2->level = level;

  c3->x = c1->x;
  c3->y = c2->y;
  c3->z = c0->z;
  c3->level = level;

  c4->x = c0->x;
  c4->y = c0->y;
  c4->z = c0->z | inc;
  c4->level = level;

  c5->x = c1->x;
  c5->y = c1->y;
  c5->z = c4->z;
  c5->level = level;

  c6->x = c2->x;
  c6->y = c2->y;
  c6->z = c4->z;
  c6->level = level;

  c7->x = c3->x;
  c7->y = c3->y;
  c7->z = c4->z;
  c7->level = level;

  P4EST_ASSERT (p4est_quadrant_is_extended (c0));
  P4EST_ASSERT (p4est_quadrant_is_extended (c1));
  P4EST_ASSERT (p4est_quadrant_is_extended (c2));
  P4EST_ASSERT (p4est_quadrant_is_extended (c3));
  P4EST_ASSERT (p4est_quadrant_is_extended (c4));
  P4EST_ASSERT (p4est_quadrant_is_extended (c5));
  P4EST_ASSERT (p4est_quadrant_is_extended (c6));
  P4EST_ASSERT (p4est_quadrant_is_extended (c7));
  P4EST_ASSERT (p8est_quadrant_is_family (c0, c1, c2, c3, c4, c5, c6, c7));
}

/* EOF p8est_algorithms.c */
