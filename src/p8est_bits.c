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

#include <p4est_to_p8est.h>
#include "p4est_bits.c"

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

bool
p8est_quadrant_is_familyv (const p4est_quadrant_t q[])
{
  const int8_t        level = q[0].level;
  p4est_qcoord_t      inc;

  P4EST_ASSERT (p4est_quadrant_is_extended (&q[0]));
  P4EST_ASSERT (p4est_quadrant_is_extended (&q[1]));
  P4EST_ASSERT (p4est_quadrant_is_extended (&q[2]));
  P4EST_ASSERT (p4est_quadrant_is_extended (&q[3]));
  P4EST_ASSERT (p4est_quadrant_is_extended (&q[4]));
  P4EST_ASSERT (p4est_quadrant_is_extended (&q[5]));
  P4EST_ASSERT (p4est_quadrant_is_extended (&q[6]));
  P4EST_ASSERT (p4est_quadrant_is_extended (&q[7]));

  if (level == 0 || level != q[1].level ||
      level != q[2].level || level != q[3].level ||
      level != q[4].level || level != q[5].level ||
      level != q[6].level || level != q[7].level) {
    return false;
  }

  inc = P4EST_QUADRANT_LEN (level);
  return ((q[0].x + inc == q[1].x && q[0].y == q[1].y && q[0].z == q[1].z) &&
          (q[0].x == q[2].x && q[0].y + inc == q[2].y && q[0].z == q[2].z) &&
          (q[1].x == q[3].x && q[2].y == q[3].y && q[0].z == q[3].z) &&
          (q[0].x == q[4].x && q[0].y == q[4].y && q[0].z + inc == q[4].z) &&
          (q[1].x == q[5].x && q[1].y == q[5].y && q[4].z == q[5].z) &&
          (q[2].x == q[6].x && q[2].y == q[6].y && q[4].z == q[6].z) &&
          (q[3].x == q[7].x && q[3].y == q[7].y && q[4].z == q[7].z));
}

bool
p8est_quadrant_is_familypv (p4est_quadrant_t * q[])
{
  const int8_t        level = q[0]->level;
  p4est_qcoord_t      inc;

  P4EST_ASSERT (p4est_quadrant_is_extended (q[0]));
  P4EST_ASSERT (p4est_quadrant_is_extended (q[1]));
  P4EST_ASSERT (p4est_quadrant_is_extended (q[2]));
  P4EST_ASSERT (p4est_quadrant_is_extended (q[3]));
  P4EST_ASSERT (p4est_quadrant_is_extended (q[4]));
  P4EST_ASSERT (p4est_quadrant_is_extended (q[5]));
  P4EST_ASSERT (p4est_quadrant_is_extended (q[6]));
  P4EST_ASSERT (p4est_quadrant_is_extended (q[7]));

  if (level == 0 || level != q[1]->level ||
      level != q[2]->level || level != q[3]->level ||
      level != q[4]->level || level != q[5]->level ||
      level != q[6]->level || level != q[7]->level) {
    return false;
  }

  inc = P4EST_QUADRANT_LEN (level);
  return ((q[0]->x + inc == q[1]->x && q[0]->y == q[1]->y
           && q[0]->z == q[1]->z) && (q[0]->x == q[2]->x
                                      && q[0]->y + inc == q[2]->y
                                      && q[0]->z == q[2]->z)
          && (q[1]->x == q[3]->x && q[2]->y == q[3]->y && q[0]->z == q[3]->z)
          && (q[0]->x == q[4]->x && q[0]->y == q[4]->y
              && q[0]->z + inc == q[4]->z) && (q[1]->x == q[5]->x
                                               && q[1]->y == q[5]->y
                                               && q[4]->z == q[5]->z)
          && (q[2]->x == q[6]->x && q[2]->y == q[6]->y && q[4]->z == q[6]->z)
          && (q[3]->x == q[7]->x && q[3]->y == q[7]->y
              && q[4]->z == q[7]->z));
}

void
p8est_quadrant_children (const p4est_quadrant_t * q,
                         p4est_quadrant_t * c0, p4est_quadrant_t * c1,
                         p4est_quadrant_t * c2, p4est_quadrant_t * c3,
                         p4est_quadrant_t * c4, p4est_quadrant_t * c5,
                         p4est_quadrant_t * c6, p4est_quadrant_t * c7)
{
  const int8_t        level = (int8_t) (q->level + 1);
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

  /* this also verifies p4est_quadrant_is_extended (c[i]) */
  P4EST_ASSERT (p8est_quadrant_is_family (c0, c1, c2, c3, c4, c5, c6, c7));
}

void
p8est_quadrant_childrenv (const p4est_quadrant_t * q, p4est_quadrant_t c[])
{
  const int8_t        level = (int8_t) (q->level + 1);
  const p4est_qcoord_t inc = P4EST_QUADRANT_LEN (level);

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (q->level < P4EST_MAXLEVEL);

  c[0].x = q->x;
  c[0].y = q->y;
  c[0].z = q->z;
  c[0].level = level;

  c[1].x = c[0].x | inc;
  c[1].y = c[0].y;
  c[1].z = c[0].z;
  c[1].level = level;

  c[2].x = c[0].x;
  c[2].y = c[0].y | inc;
  c[2].z = c[0].z;
  c[2].level = level;

  c[3].x = c[1].x;
  c[3].y = c[2].y;
  c[3].z = c[0].z;
  c[3].level = level;

  c[4].x = c[0].x;
  c[4].y = c[0].y;
  c[4].z = c[0].z | inc;
  c[4].level = level;

  c[5].x = c[1].x;
  c[5].y = c[1].y;
  c[5].z = c[4].z;
  c[5].level = level;

  c[6].x = c[2].x;
  c[6].y = c[2].y;
  c[6].z = c[4].z;
  c[6].level = level;

  c[7].x = c[3].x;
  c[7].y = c[3].y;
  c[7].z = c[4].z;
  c[7].level = level;

  /* this also verifies p4est_quadrant_is_extended (c[i]) */
  P4EST_ASSERT (p8est_quadrant_is_familyv (c));
}

void
p8est_quadrant_transform_face (const p8est_quadrant_t * q,
                               p8est_quadrant_t * r, const int ftransform[])
{
  p4est_qcoord_t      mh, tRmh, Rmh;
  p4est_qcoord_t     *target_xyz[3];
  const p4est_qcoord_t *my_xyz[3];
  const int          *my_axis = &ftransform[0];
  const int          *target_axis = &ftransform[3];
  const int          *edge_reverse = &ftransform[6];

#ifdef P4EST_DEBUG
  int                 i;

  for (i = 0; i < 3; ++i) {
    P4EST_ASSERT (0 <= my_axis[i] && my_axis[i] < 3);
    P4EST_ASSERT (0 <= target_axis[i] && target_axis[i] < 3);
  }
#endif
  P4EST_ASSERT (my_axis[0] != my_axis[1] &&
                my_axis[0] != my_axis[2] && my_axis[1] != my_axis[2]);
  P4EST_ASSERT (target_axis[0] != target_axis[1] &&
                target_axis[0] != target_axis[2] &&
                target_axis[1] != target_axis[2]);
  P4EST_ASSERT (0 <= edge_reverse[0] && edge_reverse[0] < 2);
  P4EST_ASSERT (0 <= edge_reverse[1] && edge_reverse[1] < 2);
  P4EST_ASSERT (0 <= edge_reverse[2] && edge_reverse[2] < 4);
  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (q != r);

  mh = -P4EST_QUADRANT_LEN (q->level);
  Rmh = P4EST_ROOT_LEN + mh;
  tRmh = P4EST_ROOT_LEN + Rmh;

  my_xyz[0] = &q->x;
  my_xyz[1] = &q->y;
  my_xyz[2] = &q->z;

  target_xyz[0] = &r->x;
  target_xyz[1] = &r->y;
  target_xyz[2] = &r->z;
#ifdef P4EST_DEBUG
  r->x = r->y = r->z = (p4est_qcoord_t) P4EST_QCOORD_MIN;
#endif

  *target_xyz[target_axis[0]] =
    !edge_reverse[0] ? *my_xyz[my_axis[0]] : Rmh - *my_xyz[my_axis[0]];
  *target_xyz[target_axis[1]] =
    !edge_reverse[1] ? *my_xyz[my_axis[1]] : Rmh - *my_xyz[my_axis[1]];
  switch (edge_reverse[2]) {
  case 0:
    *target_xyz[target_axis[2]] = mh - *my_xyz[my_axis[2]];
    break;
  case 1:
    *target_xyz[target_axis[2]] = *my_xyz[my_axis[2]] + P4EST_ROOT_LEN;
    break;
  case 2:
    *target_xyz[target_axis[2]] = *my_xyz[my_axis[2]] - P4EST_ROOT_LEN;
    break;
  case 3:
    *target_xyz[target_axis[2]] = tRmh - *my_xyz[my_axis[2]];
    break;
  default:
    SC_CHECK_NOT_REACHED ();
  }

  r->level = q->level;
  P4EST_ASSERT (p4est_quadrant_is_extended (r));
  P4EST_ASSERT ((p4est_quadrant_is_inside_root (q) &&
                 !p4est_quadrant_is_inside_root (r)) ||
                (!p4est_quadrant_is_inside_root (q) &&
                 p4est_quadrant_is_inside_root (r)));
}

bool
p8est_quadrant_touches_edge (const p8est_quadrant_t * q, int edge)
{
  bool                touch;
  bool                inside[2 * P4EST_DIM];
  bool                outside[2 * P4EST_DIM];
  int                 axis, side;
#ifdef P4EST_DEBUG
  int                 incount;
#endif
  p4est_qcoord_t      mh, Rmh;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (0 <= edge && edge < 12);

  axis = edge / 4;
  mh = -P4EST_QUADRANT_LEN (q->level);
  Rmh = P4EST_ROOT_LEN + mh;

  inside[0] = (q->x == 0);
  inside[1] = (q->x == Rmh);
  inside[2] = (q->y == 0);
  inside[3] = (q->y == Rmh);
  inside[4] = (q->z == 0);
  inside[5] = (q->z == Rmh);
  outside[0] = (q->x == mh);
  outside[1] = (q->x == P4EST_ROOT_LEN);
  outside[2] = (q->y == mh);
  outside[3] = (q->y == P4EST_ROOT_LEN);
  outside[4] = (q->z == mh);
  outside[5] = (q->z == P4EST_ROOT_LEN);

  touch = true;
#ifdef P4EST_DEBUG
  incount = 0;
#endif
  if (axis != 0) {
    side = edge % 2;
    touch &= inside[side] || outside[side];
#ifdef P4EST_DEBUG
    incount += inside[side];
#endif
  }
  if (axis != 1) {
    side = (axis == 0) ? (edge % 2) : ((edge / 2) % 2);
    touch &= inside[2 + side] || outside[2 + side];
#ifdef P4EST_DEBUG
    incount += inside[2 + side];
#endif
  }
  if (axis != 2) {
    side = (edge / 2) % 2;
    touch &= inside[4 + side] || outside[4 + side];
#ifdef P4EST_DEBUG
    incount += inside[4 + side];
#endif
  }
  P4EST_ASSERT (!touch || (incount == 0 || incount == 2));
  P4EST_ASSERT (axis != 0 || (q->x >= 0 && q->x < P4EST_ROOT_LEN));
  P4EST_ASSERT (axis != 1 || (q->y >= 0 && q->y < P4EST_ROOT_LEN));
  P4EST_ASSERT (axis != 2 || (q->z >= 0 && q->z < P4EST_ROOT_LEN));

  return touch;
}

void
p8est_quadrant_transform_edge (const p8est_quadrant_t * q,
                               p8est_quadrant_t * r,
                               const p8est_edge_info_t * ei,
                               const p8est_edge_transform_t * et, bool inside)
{
  int                 iaxis;
  p4est_qcoord_t      mh, Rmh;
  p4est_qcoord_t      lshift, rshift;
  p4est_qcoord_t      my_xyz, *target_xyz[3];

  iaxis = (int) ei->iedge / 4;
  P4EST_ASSERT (0 <= ei->iflip && ei->iflip < 2);
  P4EST_ASSERT (0 <= et->naxis[0] && et->naxis[0] < 3);
  P4EST_ASSERT (0 <= et->naxis[1] && et->naxis[1] < 3);
  P4EST_ASSERT (0 <= et->naxis[2] && et->naxis[2] < 3);
  P4EST_ASSERT (et->naxis[0] != et->naxis[1] &&
                et->naxis[0] != et->naxis[2] && et->naxis[1] != et->naxis[2]);
  P4EST_ASSERT (0 <= et->nflip && et->nflip < 2);
  P4EST_ASSERT (p8est_quadrant_touches_edge (q, (int) ei->iedge));
  P4EST_ASSERT (q != r);

  mh = -P4EST_QUADRANT_LEN (q->level);
  Rmh = P4EST_ROOT_LEN + mh;
  lshift = (inside ? 0 : mh);
  rshift = (inside ? Rmh : P4EST_ROOT_LEN);
  target_xyz[0] = &r->x;
  target_xyz[1] = &r->y;
  target_xyz[2] = &r->z;

  /* transform coordinate axis parallel to edge */
  switch (iaxis) {
  case 0:
    my_xyz = q->x;
    break;
  case 1:
    my_xyz = q->y;
    break;
  case 2:
    my_xyz = q->z;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
  }
  if (et->nflip == ei->iflip) {
    *target_xyz[et->naxis[0]] = my_xyz;
  }
  else {
    *target_xyz[et->naxis[0]] = Rmh - my_xyz;
  }

  /* create the other two coordinates */
  switch (et->corners) {
  case 0:
    *target_xyz[et->naxis[1]] = lshift;
    *target_xyz[et->naxis[2]] = lshift;
    break;
  case 1:
    *target_xyz[et->naxis[1]] = rshift;
    *target_xyz[et->naxis[2]] = lshift;
    break;
  case 2:
    *target_xyz[et->naxis[1]] = lshift;
    *target_xyz[et->naxis[2]] = rshift;
    break;
  case 3:
    *target_xyz[et->naxis[1]] = rshift;
    *target_xyz[et->naxis[2]] = rshift;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
  }

  r->level = q->level;
  P4EST_ASSERT (p8est_quadrant_touches_edge (r, (int) et->nedge));
  P4EST_ASSERT (!inside || p4est_quadrant_is_inside_root (r));
}

void
p8est_quadrant_shift_edge (const p8est_quadrant_t * q,
                           p8est_quadrant_t * r, int edge)
{
  int                 outface;
  int                 cid, sid, step[P4EST_DIM];
  p4est_qcoord_t      th;
  p4est_quadrant_t    quad;
  /* *INDENT-OFF* */
  const int           contact[12] = {
    0x14, 0x18, 0x24, 0x28,
    0x11, 0x12, 0x21, 0x22,
    0x05, 0x06, 0x09, 0x0a
  };
  /* *INDENT-ON* */

  P4EST_ASSERT (q != r);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (edge >= 0 && edge < 12);

  P4EST_QUADRANT_INIT (&quad);

  quad = *q;
  for (;;) {
    th = P4EST_LAST_OFFSET (quad.level);
    cid = p4est_quadrant_child_id (&quad);
    switch (edge / 4) {
    case 0:
      sid = 2 * edge + (cid & 0x01);
      step[0] = 0;
      step[1] = 2 * (edge & 0x01) - 1;
      step[2] = (edge & 0x02) - 1;
      break;
    case 1:
      sid = 2 * (edge & 0x02) + (edge & 0x01) + (cid & 0x02);
      step[0] = 2 * (edge & 0x01) - 1;
      step[1] = 0;
      step[2] = (edge & 0x02) - 1;
      break;
    case 2:
      sid = edge - 8 + (cid & 0x04);
      step[0] = 2 * (edge & 0x01) - 1;
      step[1] = (edge & 0x02) - 1;
      step[2] = 0;
      break;
    default:
      SC_CHECK_NOT_REACHED ();
    }
    p4est_quadrant_sibling (&quad, r, sid);
    P4EST_ASSERT (-1 <= step[0] && step[0] <= 1);
    P4EST_ASSERT (-1 <= step[1] && step[1] <= 1);
    P4EST_ASSERT (-1 <= step[2] && step[2] <= 1);

    outface = 0;
    if (step[0] != 0) {
      outface |= ((r->x <= 0) ? 0x01 : 0);
      outface |= ((r->x >= th) ? 0x02 : 0);
    }
    if (step[1] != 0) {
      outface |= ((r->y <= 0) ? 0x04 : 0);
      outface |= ((r->y >= th) ? 0x08 : 0);
    }
    if (step[2] != 0) {
      outface |= ((r->z <= 0) ? 0x10 : 0);
      outface |= ((r->z >= th) ? 0x20 : 0);
    }
    if (outface == contact[edge]) {
      break;
    }
    p4est_quadrant_parent (&quad, &quad);
    quad.x += (p4est_qcoord_t) step[0] * P4EST_QUADRANT_LEN (quad.level);
    quad.y += (p4est_qcoord_t) step[1] * P4EST_QUADRANT_LEN (quad.level);
    quad.z += (p4est_qcoord_t) step[2] * P4EST_QUADRANT_LEN (quad.level);
    P4EST_ASSERT (p4est_quadrant_is_extended (&quad));
  }

  if (step[0] != 0) {
    if (r->x < 0)
      r->x = 0;
    if (r->x >= P4EST_ROOT_LEN)
      r->x = th;
  }
  if (step[1] != 0) {
    if (r->y < 0)
      r->y = 0;
    if (r->y >= P4EST_ROOT_LEN)
      r->y = th;
  }
  if (step[2] != 0) {
    if (r->z < 0)
      r->z = 0;
    if (r->z >= P4EST_ROOT_LEN)
      r->z = th;
  }
  P4EST_ASSERT (p4est_quadrant_is_valid (r));
}

/* EOF p8est_bits.h */
