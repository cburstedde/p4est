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

#include <p8est_mesh.h>
#include <p8est_algorithms.h>

#include <p4est_to_p8est.h>

bool
p4est_is_balanced (p4est_t * p4est)
{
  return true;
}

p4est_neighborhood_t *
p4est_neighborhood_new (p4est_t * p4est)
{
  p4est_topidx_t      local_num_trees, flt, nt;
  p4est_locidx_t      local_num_quadrants, lsum;
  p4est_tree_t       *tree;
  p4est_neighborhood_t *nhood;

  P4EST_ASSERT (p4est_is_valid (p4est));

  if (p4est->first_local_tree < 0) {
    flt = 0;
    local_num_trees = 0;
  }
  else {
    flt = p4est->first_local_tree;
    local_num_trees = p4est->last_local_tree - flt + 1; /* type ok */
  }
  local_num_quadrants = p4est->local_num_quadrants;

  nhood = P4EST_ALLOC (p4est_neighborhood_t, 1);
  nhood->cumulative_count = P4EST_ALLOC (p4est_locidx_t, local_num_trees + 1);
  nhood->element_offsets =
    P4EST_ALLOC (p4est_locidx_t, local_num_quadrants + 1);
  nhood->local_neighbors = sc_array_new (sizeof (p4est_locidx_t));

  lsum = 0;
  for (nt = 0; nt < local_num_trees; ++nt) {
    nhood->cumulative_count[nt] = lsum;
    tree = p4est_array_index_topidx (p4est->trees, flt + nt);   /* type ok */
    lsum += (p4est_locidx_t) tree->quadrants.elem_count;        /* type ok */
  }
  P4EST_ASSERT (lsum == local_num_quadrants);
  nhood->cumulative_count[nt] = lsum;

  return nhood;
}

void
p4est_neighborhood_destroy (p4est_neighborhood_t * nhood)
{
  P4EST_FREE (nhood->cumulative_count);
  P4EST_FREE (nhood->element_offsets);
  sc_array_destroy (nhood->local_neighbors);

  P4EST_FREE (nhood);
}

/* EOF p8est_mesh.c */
