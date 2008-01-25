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

#include <p4est_mesh.h>
#include <p4est_base.h>
#include <p4est_algorithms.h>

void
p4est_order_local_vertices (p4est_t * p4est,
                            int32_t *
                            num_uniq_local_vertices,
                            int32_t * quadrant_to_local_vertex)
{
  int32_t             Ntotal = p4est->local_num_quadrants * 4;
  int32_t             Ncells = p4est->local_num_quadrants;
  int32_t             i;

  for (i = 0; i < Ncells; ++i) {
    quadrant_to_local_vertex[4 * i + 0] = 4 * i + 0;
    quadrant_to_local_vertex[4 * i + 1] = 4 * i + 1;
    quadrant_to_local_vertex[4 * i + 2] = 4 * i + 2;
    quadrant_to_local_vertex[4 * i + 3] = 4 * i + 3;
  }

  *num_uniq_local_vertices = Ntotal;
}

void
p4est_possible_node_neigbors (p4est_quadrant_t * q,
                              int32_t node,
                              int32_t nnum,
                              int8_t neighbor_rlev,
                              p4est_quadrant_t * neighbor,
                              int32_t * neighbor_node)
{
  p4est_quadrant_t    n;
  int32_t             nnode;
  int32_t             qh = (1 << (P4EST_MAXLEVEL - q->level));
  int32_t             nh =
    (1 << (P4EST_MAXLEVEL - (q->level + neighbor_rlev)));
  int32_t             qx = q->x;
  int32_t             qy = q->y;
  int32_t             cornerx, cornery;

  switch (node) {
  case 0:
    cornerx = qx;
    cornery = qy;
    break;
  case 1:
    cornerx = qx + qh;
    cornery = qy;
    break;
  case 2:
    cornerx = qx;
    cornery = qy + qh;
    break;
  case 3:
    cornerx = qx + qh;
    cornery = qy + qh;
    break;
  default:
    P4EST_ASSERT_NOT_REACHED ();
  }

#ifdef P4EST_HAVE_DEBUG
  /* Check to see if it is possible to construct the neighbor */
  int8_t              qcid = p4est_quadrant_child_id (q);
  P4EST_ASSERT (neighbor_rlev >= 0 || qcid == node);
#endif

  nnode = 3 - nnum;
  n.level = (int8_t) (q->level + neighbor_rlev);
  switch (nnum) {
  case 0:
    n.x = cornerx - nh;
    n.y = cornery - nh;
    break;
  case 1:
    n.x = cornerx;
    n.y = cornery - nh;
    break;
  case 2:
    n.x = cornerx - nh;
    n.y = cornery;
    break;
  case 3:
    n.x = cornerx;
    n.y = cornery;
    break;
  default:
    P4EST_ASSERT_NOT_REACHED ();
  }

  *neighbor = n;
  *neighbor_node = nnode;

  P4EST_ASSERT (p4est_quadrant_is_extended (neighbor));
}

/* EOF p4est_mesh.c */
