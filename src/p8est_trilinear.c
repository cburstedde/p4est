/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007-2009 Carsten Burstedde, Lucas Wilcox.

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

#include <p8est_trilinear.h>

void
p8est_trilinear_mesh_destroy (trilinear_mesh_t * mesh)
{
  P4EST_ASSERT (mesh->destructor == p8est_trilinear_mesh_destroy);
  P4EST_ASSERT (mesh->extra_info == NULL);

  P4EST_FREE (mesh->elem_table);
  P4EST_FREE (mesh->node_table);
  P4EST_FREE (mesh->fvnid_count_table);
  P4EST_FREE (mesh->fvnid_interval_table);
  sc_mempool_destroy (mesh->sharer_pool);

  P4EST_FREE (mesh->elem_pids);
  P4EST_FREE (mesh->node_pids);

  P4EST_FREE (mesh);
}
