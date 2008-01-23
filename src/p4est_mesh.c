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

/* EOF p4est_mesh.c */
